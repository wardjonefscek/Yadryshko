//! Замкнутый контур на 2D-сетке: площадь, ориентация, положение точки.
//!
//! Контур хранит вершины без повтора первой в конце. При создании убираются
//! повторы подряд, вершины на прямой между соседями и «усы» (заход и возврат
//! по той же прямой) — всё это не меняет ограниченную площадь.

use crate::grid::{GridError, GridPoint};
use crate::predicates::{Sign, orient2d_grid, orient2d_grid_value};
use crate::seg::{Crossing, classify, on_segment};
use crate::vec::Vec2;
use std::fmt;

/// Где точка относительно контура или области.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Location {
    /// Внутри.
    Inside,
    /// На границе.
    Boundary,
    /// Снаружи.
    Outside,
}

/// Почему не удалось построить контур.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum RingError {
    /// После упрощения осталось меньше трёх вершин (все точки на одной прямой).
    TooFewPoints,
    /// Точка за пределом сетки.
    OutOfGrid(GridError),
}

impl fmt::Display for RingError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::TooFewPoints => write!(f, "в контуре меньше трёх вершин не на одной прямой"),
            Self::OutOfGrid(e) => write!(f, "{e}"),
        }
    }
}

impl std::error::Error for RingError {}

/// Замкнутый контур: не меньше трёх вершин, никакие три соседние не на одной прямой.
/// Самопересечения не запрещены — см. [`Ring::is_simple`]. У самопересекающегося
/// контура ориентированная площадь может быть нулевой («бантик»), хотя по правилу
/// ненулевой обмотки он закрашен; поэтому нулевая площадь здесь не ошибка.
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub struct Ring {
    pts: Vec<GridPoint>,
}

impl Ring {
    /// Контур из вершин на сетке (с упрощением, см. описание модуля).
    pub fn new(points: Vec<GridPoint>) -> Result<Self, RingError> {
        let pts = simplify(points);
        if pts.len() < 3 {
            return Err(RingError::TooFewPoints);
        }
        // После упрощения у любых трёх соседних вершин ориентация не нулевая.
        Ok(Self { pts })
    }

    /// Контур из точек в миллиметрах (с округлением на сетку).
    pub fn from_mm(points: &[Vec2]) -> Result<Self, RingError> {
        let pts = points
            .iter()
            .map(|p| GridPoint::from_mm(p.x, p.y))
            .collect::<Result<Vec<_>, _>>()
            .map_err(RingError::OutOfGrid)?;
        Self::new(pts)
    }

    /// Прямоугольник по двум углам (в шагах сетки), обход против часовой.
    pub fn rect(x0: i64, y0: i64, x1: i64, y1: i64) -> Result<Self, RingError> {
        let (xa, xb) = (x0.min(x1), x0.max(x1));
        let (ya, yb) = (y0.min(y1), y0.max(y1));
        let p = |x, y| GridPoint::new(x, y).map_err(RingError::OutOfGrid);
        Self::new(vec![p(xa, ya)?, p(xb, ya)?, p(xb, yb)?, p(xa, yb)?])
    }

    /// Вершины.
    pub fn points(&self) -> &[GridPoint] {
        &self.pts
    }

    /// Число вершин.
    pub fn len(&self) -> usize {
        self.pts.len()
    }

    /// Контур не бывает пустым; метод для порядка (clippy).
    pub fn is_empty(&self) -> bool {
        self.pts.is_empty()
    }

    /// Рёбра: пары соседних вершин, последнее ребро замыкает контур.
    pub fn edges(&self) -> impl Iterator<Item = (GridPoint, GridPoint)> + '_ {
        let n = self.pts.len();
        (0..n).map(move |i| (self.pts[i], self.pts[(i + 1) % n]))
    }

    /// Удвоенная ориентированная площадь, точно (в квадратных шагах сетки).
    /// Больше нуля — обход против часовой стрелки.
    pub fn twice_area(&self) -> i128 {
        let o = self.pts[0];
        self.edges()
            .map(|(a, b)| orient2d_grid_value(o, a, b))
            .sum()
    }

    /// Площадь в мм² (по модулю).
    pub fn area_mm2(&self) -> f64 {
        self.twice_area().unsigned_abs() as f64 / 2.0 * 1e-6
    }

    /// Обход против часовой стрелки.
    pub fn is_ccw(&self) -> bool {
        self.twice_area() > 0
    }

    /// Тот же контур с обратным обходом.
    pub fn reversed(&self) -> Self {
        let mut pts = self.pts.clone();
        pts.reverse();
        Self { pts }
    }

    /// Число обмотки контура вокруг точки `p`, не лежащей на контуре
    /// (алгоритм Сандей: луч вправо, полуоткрытое правило по y).
    pub fn winding(&self, p: GridPoint) -> i64 {
        let mut w = 0;
        for (a, b) in self.edges() {
            if a.y() <= p.y() {
                if b.y() > p.y() && orient2d_grid(a, b, p) == Sign::Positive {
                    w += 1;
                }
            } else if b.y() <= p.y() && orient2d_grid(a, b, p) == Sign::Negative {
                w -= 1;
            }
        }
        w
    }

    /// На контуре ли точка.
    pub fn on_boundary(&self, p: GridPoint) -> bool {
        self.edges().any(|(a, b)| on_segment(p, a, b))
    }

    /// Положение точки по правилу ненулевой обмотки.
    pub fn locate(&self, p: GridPoint) -> Location {
        if self.on_boundary(p) {
            Location::Boundary
        } else if self.winding(p) != 0 {
            Location::Inside
        } else {
            Location::Outside
        }
    }

    /// Нет ли самопересечений: несоседние рёбра не имеют общих точек,
    /// соседние встречаются только в общей вершине.
    pub fn is_simple(&self) -> bool {
        let e: Vec<_> = self.edges().collect();
        let n = e.len();
        for i in 0..n {
            for j in i + 1..n {
                let adjacent = j == i + 1 || (i == 0 && j == n - 1);
                let c = classify(e[i].0, e[i].1, e[j].0, e[j].1);
                match c {
                    Crossing::None => {}
                    Crossing::Touch if adjacent => {}
                    _ => return false,
                }
            }
        }
        true
    }

    /// Габарит: (нижний левый, верхний правый).
    pub fn bbox(&self) -> (GridPoint, GridPoint) {
        let (mut lx, mut ly, mut hx, mut hy) = (i64::MAX, i64::MAX, i64::MIN, i64::MIN);
        for p in &self.pts {
            lx = lx.min(p.x());
            ly = ly.min(p.y());
            hx = hx.max(p.x());
            hy = hy.max(p.y());
        }
        let g = |x, y| GridPoint::new(x, y).expect("вершины на сетке");
        (g(lx, ly), g(hx, hy))
    }

    /// Вершины в миллиметрах.
    pub fn to_mm(&self) -> Vec<Vec2> {
        self.pts.iter().map(|p| p.to_mm()).collect()
    }
}

/// Убрать повторы, вершины на прямой и «усы», пока что-то убирается.
fn simplify(mut pts: Vec<GridPoint>) -> Vec<GridPoint> {
    loop {
        let before = pts.len();
        pts.dedup();
        while pts.len() > 1 && pts.first() == pts.last() {
            pts.pop();
        }
        let n = pts.len();
        if n < 3 {
            return pts;
        }
        // Вершина лишняя, если она на прямой между соседями (прямо или «ус» назад).
        let mut keep = Vec::with_capacity(n);
        let mut removed = false;
        for i in 0..n {
            let prev = if keep.is_empty() {
                pts[(i + n - 1) % n]
            } else {
                *keep.last().expect("не пусто")
            };
            let next = pts[(i + 1) % n];
            if !removed && orient2d_grid(prev, pts[i], next) == Sign::Zero {
                removed = true; // по одной за проход: соседи меняются
                continue;
            }
            keep.push(pts[i]);
        }
        pts = keep;
        if pts.len() == before && !removed {
            return pts;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn p(x: i64, y: i64) -> GridPoint {
        GridPoint::new(x, y).unwrap()
    }

    #[test]
    fn area_and_orientation() {
        let r = Ring::rect(0, 0, 4, 3).unwrap();
        assert_eq!(r.twice_area(), 24);
        assert!(r.is_ccw());
        assert_eq!(r.reversed().twice_area(), -24);
    }

    #[test]
    fn simplification_removes_redundant_vertices() {
        // Повтор, вершина на стороне, «ус» и замыкающая копия первой точки.
        let r = Ring::new(vec![
            p(0, 0),
            p(0, 0),
            p(2, 0),
            p(4, 0),
            p(4, 4),
            p(4, 6),
            p(4, 4),
            p(0, 4),
            p(0, 0),
        ])
        .unwrap();
        assert_eq!(r.points(), &[p(0, 0), p(4, 0), p(4, 4), p(0, 4)]);
        assert_eq!(
            Ring::new(vec![p(0, 0), p(1, 1), p(2, 2)]),
            Err(RingError::TooFewPoints)
        );
    }

    #[test]
    fn locate_points() {
        let r = Ring::rect(0, 0, 4, 4).unwrap();
        assert_eq!(r.locate(p(2, 2)), Location::Inside);
        assert_eq!(r.locate(p(4, 2)), Location::Boundary);
        assert_eq!(r.locate(p(0, 0)), Location::Boundary);
        assert_eq!(r.locate(p(5, 2)), Location::Outside);
        // По часовой — внутри то же самое, обмотка −1.
        assert_eq!(r.reversed().winding(p(2, 2)), -1);
    }

    #[test]
    fn simplicity() {
        assert!(Ring::rect(0, 0, 4, 4).unwrap().is_simple());
        // Несимметричный «бантик»: у симметричного площадь нулевая.
        let bowtie = Ring::new(vec![p(0, 0), p(4, 4), p(4, 0), p(0, 6)]).unwrap();
        assert!(!bowtie.is_simple());
        // «Восьмёрка» с касанием в вершине (2, 2): несоседние рёбра встречаются.
        let pinch = Ring::new(vec![p(0, 0), p(4, 0), p(2, 2), p(4, 4), p(0, 4), p(2, 2)]).unwrap();
        assert!(!pinch.is_simple());
    }
}
