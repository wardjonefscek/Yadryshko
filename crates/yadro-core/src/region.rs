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
use crate::grid::GridPoint;
use crate::ring::{Location, Ring};
use crate::seg::{self, Crossing};
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

    /// Все рёбра всех контуров.
    pub fn edges(&self) -> impl Iterator<Item = (GridPoint, GridPoint)> + '_ {
        self.rings.iter().flat_map(Ring::edges)
    }

    /// Проверить инвариант (см. описание модуля). Сложность квадратичная — для тестов
    /// и для проверки чужих данных.
    pub fn check(&self) -> Result<(), RegionIssue> {
        let edges: Vec<_> = self.edges().collect();
        for i in 0..edges.len() {
            for j in i + 1..edges.len() {
                let (a, b) = (edges[i], edges[j]);
                match seg::classify(a.0, a.1, b.0, b.1) {
                    Crossing::Proper | Crossing::Overlap => {
                        return Err(RegionIssue::Crossing { a, b });
                    }
                    Crossing::Touch | Crossing::None => {}
                }
            }
        }
        // Делим рёбра в вершинах, лежащих посреди других рёбер (касание вершиной стороны),
        // чтобы середина каждого кусочка была внутренней точкой грани.
        let verts: Vec<GridPoint> = self
            .rings
            .iter()
            .flat_map(|r| r.points().iter().copied())
            .collect();
        let mut parts: Vec<Edge> = Vec::new();
        for &(p, q) in &edges {
            let mut mids: Vec<GridPoint> = verts
                .iter()
                .copied()
                .filter(|&c| seg::in_segment_interior(c, p, q))
                .collect();
            mids.sort_by_key(|c| {
                let (dx, dy) = (i128::from(c.x() - p.x()), i128::from(c.y() - p.y()));
                dx * dx + dy * dy
            });
            mids.dedup();
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
        for k in 0..parts.len() {
            let left = boolean::left_winding(&parts, k)[0];
            let right = left - parts[k].w[0];
            // Обмотка «слева» — относительно направления кусочка u → v; приводим к направлению контура.
            let (l, r) = if parts[k].w[0] > 0 {
                (left, right)
            } else {
                (right, left)
            };
            if l != 1 || r != 0 {
                let e = if parts[k].w[0] > 0 {
                    (parts[k].u, parts[k].v)
                } else {
                    (parts[k].v, parts[k].u)
                };
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
