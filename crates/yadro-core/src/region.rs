//! Область на плоскости: несколько контуров, внешние — против часовой стрелки,
//! отверстия — по часовой.
//!
//! Инвариант корректной области (его проверяет [`Region::check`]):
//! - контуры не пересекаются и не накладываются; касаться можно только в точке;
//! - слева от каждого ребра обмотка равна 1, справа — 0.
//!
//! Отсюда: внутри области обмотка 1, снаружи 0, площадь — сумма ориентированных
//! площадей контуров. Результаты булевых операций всегда корректны.

use crate::boolean::{self, Edge, Op};
use crate::grid::{GridError, GridPoint};
use crate::offset;
use crate::ring::{Location, Ring};
use crate::seg::{self, Crossing};
use crate::tol::CHORD_TOL_MM;
use std::collections::BTreeSet;
use std::fmt;

/// Что не так с областью.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum RegionIssue {
    /// Два ребра пересекаются или накладываются.
    Crossing {
        /// Первое ребро.
        a: (GridPoint, GridPoint),
        /// Второе ребро.
        b: (GridPoint, GridPoint),
    },
    /// Обмотка возле ребра не 1 слева и 0 справа (неверная ориентация или вложенность).
    Winding {
        /// Ребро.
        edge: (GridPoint, GridPoint),
        /// Обмотка слева.
        left: i64,
        /// Обмотка справа.
        right: i64,
    },
}

impl fmt::Display for RegionIssue {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Crossing { a, b } => {
                write!(f, "рёбра {a:?} и {b:?} пересекаются или накладываются")
            }
            Self::Winding { edge, left, right } => {
                write!(
                    f,
                    "у ребра {edge:?} обмотка слева {left}, справа {right} (нужно 1 и 0)"
                )
            }
        }
    }
}

impl std::error::Error for RegionIssue {}

/// Область: набор контуров с соблюдённым инвариантом.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct Region {
    rings: Vec<Ring>,
}

impl Region {
    /// Пустая область.
    pub fn empty() -> Self {
        Self::default()
    }

    /// Область из контуров с проверкой инварианта.
    pub fn new(rings: Vec<Ring>) -> Result<Self, RegionIssue> {
        let r = Self { rings };
        r.check()?;
        Ok(r)
    }

    /// Область по правилу ненулевой обмотки из любых контуров — самопересекающихся,
    /// наложенных, с любой ориентацией. Внутри — там, где сумма обмоток не ноль.
    pub fn from_rings_nonzero(rings: &[Ring]) -> Self {
        boolean::normalize(rings)
    }

    /// Область из одного контура (по правилу ненулевой обмотки).
    pub fn from_ring(ring: Ring) -> Self {
        Self::from_rings_nonzero(&[ring])
    }

    /// Без проверки: для результатов, корректных по построению.
    pub(crate) fn from_trusted(rings: Vec<Ring>) -> Self {
        Self { rings }
    }

    /// Контуры.
    pub fn rings(&self) -> &[Ring] {
        &self.rings
    }

    /// Пустая ли область.
    pub fn is_empty(&self) -> bool {
        self.rings.is_empty()
    }

    /// Удвоенная площадь, точно (в квадратных шагах сетки).
    pub fn twice_area(&self) -> i128 {
        self.rings.iter().map(Ring::twice_area).sum()
    }

    /// Площадь, мм².
    pub fn area_mm2(&self) -> f64 {
        self.twice_area() as f64 / 2.0 * 1e-6
    }

    /// Положение точки.
    pub fn locate(&self, p: GridPoint) -> Location {
        if self.rings.iter().any(|r| r.on_boundary(p)) {
            Location::Boundary
        } else if self.rings.iter().map(|r| r.winding(p)).sum::<i64>() != 0 {
            Location::Inside
        } else {
            Location::Outside
        }
    }

    /// Объединение.
    pub fn union(&self, other: &Self) -> Self {
        boolean::boolean(self, other, Op::Union)
    }

    /// Пересечение.
    pub fn intersection(&self, other: &Self) -> Self {
        boolean::boolean(self, other, Op::Intersection)
    }

    /// Разность `self ∖ other`.
    pub fn difference(&self, other: &Self) -> Self {
        boolean::boolean(self, other, Op::Difference)
    }

    /// Симметрическая разность.
    pub fn xor(&self, other: &Self) -> Self {
        boolean::boolean(self, other, Op::Xor)
    }

    /// Эквидистанта на `d_mm`: > 0 — наружу, < 0 — внутрь; дуги с допуском хорды
    /// [`CHORD_TOL_MM`]. Подробности и точность — в [`crate::offset`].
    pub fn offset(&self, d_mm: f64) -> Result<Self, GridError> {
        offset::offset(self, d_mm, CHORD_TOL_MM)
    }

    /// Эквидистанта с заданным допуском хорды, мм.
    pub fn offset_tol(&self, d_mm: f64, tol_mm: f64) -> Result<Self, GridError> {
        offset::offset(self, d_mm, tol_mm)
    }

    /// Открытие (opening): сначала внутрь на `r_mm`, потом наружу. Скругляет выпуклые
    /// углы радиусом r и убирает всё, что уже 2r, — перемычки, усы, мелкие выступы.
    pub fn opening(&self, r_mm: f64) -> Result<Self, GridError> {
        self.offset(-r_mm.abs())?.offset(r_mm.abs())
    }

    /// Закрытие (closing): сначала наружу на `r_mm`, потом внутрь. Скругляет вогнутые
    /// углы радиусом r и закрывает щели и отверстия уже 2r. Такой деталь выйдет из-под
    /// фрезы радиуса r: во внутренние углы и узкие пазы фреза не залезет.
    pub fn closing(&self, r_mm: f64) -> Result<Self, GridError> {
        self.offset(r_mm.abs())?.offset(-r_mm.abs())
    }

    /// Убрать острова площадью меньше `min_island_mm2` и заполнить дыры меньше
    /// `min_hole_mm2` (площадь контура — вместе со всем, что внутри). Вложенность
    /// учитывается: убранный остров уносит свои дыры, заполненная дыра — острова в
    /// ней. Так чистят шум топооптимизации.
    pub fn drop_small(&self, min_island_mm2: f64, min_hole_mm2: f64) -> Self {
        let lim = |mm2: f64| mm2.max(0.0) * 2e6; // в удвоенных квадратных шагах сетки
        let (isl, hole) = (lim(min_island_mm2), lim(min_hole_mm2));
        let area: Vec<i128> = self.rings.iter().map(Ring::twice_area).collect();
        let small: Vec<bool> = area
            .iter()
            .map(|&a| {
                if a > 0 {
                    (a as f64) < isl
                } else {
                    ((-a) as f64) < hole
                }
            })
            .collect();
        if !small.contains(&true) {
            return self.clone();
        }
        // Внутри ли контур k контура m: контуры не пересекаются, поэтому хватает одной
        // вершины k не на границе m.
        let inside = |k: usize, m: usize| {
            let ring = &self.rings[m];
            self.rings[k]
                .points()
                .iter()
                .find(|&&q| !ring.on_boundary(q))
                .is_some_and(|&q| ring.winding(q) != 0)
        };
        let n = self.rings.len();
        let parent: Vec<Option<usize>> = (0..n)
            .map(|k| {
                (0..n)
                    .filter(|&m| area[m].abs() > area[k].abs() && inside(k, m))
                    .min_by_key(|&m| area[m].abs())
            })
            .collect();
        let keep = |k: usize| {
            let mut c = Some(k);
            while let Some(i) = c {
                if small[i] {
                    return false;
                }
                c = parent[i];
            }
            true
        };
        let rings: Vec<Ring> = (0..n)
            .filter(|&k| keep(k))
            .map(|k| self.rings[k].clone())
            .collect();
        Self::from_rings_nonzero(&rings)
    }

    /// Все рёбра всех контуров.
    pub fn edges(&self) -> impl Iterator<Item = (GridPoint, GridPoint)> + '_ {
        self.rings.iter().flat_map(Ring::edges)
    }

    /// Проверить инвариант (см. описание модуля): для тестов и для проверки чужих данных.
    pub fn check(&self) -> Result<(), RegionIssue> {
        let edges: Vec<_> = self.edges().collect();
        let mut bad = None;
        boolean::candidate_pairs(&edges, |i, j| {
            if bad.is_some() {
                return;
            }
            let (a, b) = (edges[i], edges[j]);
            if matches!(
                seg::classify(a.0, a.1, b.0, b.1),
                Crossing::Proper | Crossing::Overlap
            ) {
                bad = Some(RegionIssue::Crossing { a, b });
            }
        });
        if let Some(issue) = bad {
            return Err(issue);
        }
        // Делим рёбра в вершинах, лежащих посреди других рёбер (касание вершиной стороны),
        // чтобы середина каждого кусочка была внутренней точкой грани.
        let verts: BTreeSet<GridPoint> = self
            .rings
            .iter()
            .flat_map(|r| r.points().iter().copied())
            .collect();
        let hot = boolean::Hot::new(&verts);
        let mut parts: Vec<Edge> = Vec::new();
        for &(p, q) in &edges {
            let mut mids: Vec<GridPoint> = hot
                .near(p, q)
                .filter(|&c| seg::in_segment_interior(c, p, q))
                .collect();
            mids.sort_by_key(|c| {
                let (dx, dy) = (i128::from(c.x() - p.x()), i128::from(c.y() - p.y()));
                dx * dx + dy * dy
            });
            let mut prev = p;
            for m in mids.into_iter().chain([q]) {
                let (u, v, s) = if prev < m {
                    (prev, m, 1)
                } else {
                    (m, prev, -1)
                };
                parts.push(Edge { u, v, w: [s, 0] });
                prev = m;
            }
        }
        let windings = boolean::edge_windings(&parts);
        for (part, w) in parts.iter().zip(&windings) {
            let left = w[0];
            let right = left - part.w[0];
            // Обмотка «слева» — относительно направления кусочка u → v; приводим к направлению контура.
            let (l, r, e) = if part.w[0] > 0 {
                (left, right, (part.u, part.v))
            } else {
                (right, left, (part.v, part.u))
            };
            if l != 1 || r != 0 {
                return Err(RegionIssue::Winding {
                    edge: e,
                    left: l,
                    right: r,
                });
            }
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn p(x: i64, y: i64) -> GridPoint {
        GridPoint::new(x, y).unwrap()
    }

    #[test]
    fn square_with_hole_is_valid() {
        let outer = Ring::rect(0, 0, 10, 10).unwrap();
        let hole = Ring::rect(3, 3, 6, 6).unwrap().reversed();
        let r = Region::new(vec![outer.clone(), hole]).unwrap();
        assert_eq!(r.twice_area(), 200 - 18);
        assert_eq!(r.locate(p(4, 4)), Location::Outside);
        assert_eq!(r.locate(p(1, 1)), Location::Inside);
        // Отверстие с неправильной ориентацией — ошибка инварианта.
        let bad = Region::new(vec![outer, Ring::rect(3, 3, 6, 6).unwrap()]);
        assert!(matches!(bad, Err(RegionIssue::Winding { .. })));
    }

    #[test]
    fn crossing_rings_are_rejected() {
        let a = Ring::rect(0, 0, 4, 4).unwrap();
        let b = Ring::rect(2, 2, 6, 6).unwrap();
        assert!(matches!(
            Region::new(vec![a, b]),
            Err(RegionIssue::Crossing { .. })
        ));
    }

    #[test]
    fn touching_hole_vertex_on_outer_edge_is_valid() {
        // Отверстие-треугольник касается вершиной середины нижней стороны.
        let outer = Ring::rect(0, 0, 10, 10).unwrap();
        let hole = Ring::new(vec![p(5, 0), p(3, 4), p(7, 4)]).unwrap();
        let hole = if hole.is_ccw() { hole.reversed() } else { hole };
        Region::new(vec![outer, hole]).unwrap();
    }
}
