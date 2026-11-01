//! 2D-сетка: целые координаты в микрометрах.
//!
//! Все 2D-построения живут на этой сетке. Так предикаты считаются точно
//! в целых числах, а погрешность округления видна заранее: не больше
//! половины шага (0,5 мкм) по каждой оси.

use crate::tol::{GRID_LIMIT, GRID_PER_MM};
use crate::vec::Vec2;
use std::fmt;

/// Точка 2D-сетки; единица — 1 мкм.
///
/// Инвариант: |x| и |y| не больше [`GRID_LIMIT`]. Поля закрыты,
/// чтобы точку нельзя было создать в обход проверки.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct GridPoint {
    x: i64,
    y: i64,
}

/// Ошибка перевода на сетку.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum GridError {
    /// На входе не число или бесконечность.
    NotFinite,
    /// Значение за пределом сетки.
    OutOfRange,
}

impl fmt::Display for GridError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::NotFinite => write!(f, "на входе не число или бесконечность"),
            Self::OutOfRange => write!(f, "значение за пределом 2D-сетки (±2^40 мкм)"),
        }
    }
}

impl std::error::Error for GridError {}

fn in_range(v: i64) -> bool {
    (-GRID_LIMIT..=GRID_LIMIT).contains(&v)
}

/// Перевести миллиметры в шаги сетки, округляя к ближайшему шагу
/// (ровно посередине — от нуля).
///
/// РИСК: десятичная «середина» вроде 0,0015 мм в двоичном виде может
/// оказаться чуть меньше середины и округлиться вниз. Для геометрии это
/// допустимо, но тесты не должны на это опираться.
pub fn mm_to_grid(mm: f64) -> Result<i64, GridError> {
    if !mm.is_finite() {
        return Err(GridError::NotFinite);
    }
    let g = (mm * GRID_PER_MM).round();
    if g.abs() > GRID_LIMIT as f64 {
        return Err(GridError::OutOfRange);
    }
    Ok(g as i64)
}

/// Перевести шаги сетки в миллиметры.
pub fn grid_to_mm(g: i64) -> f64 {
    g as f64 / GRID_PER_MM
}

impl GridPoint {
    /// Начало координат.
    pub const ORIGIN: Self = Self { x: 0, y: 0 };

    /// Точка из шагов сетки; ошибка, если она за пределом сетки.
    pub fn new(x: i64, y: i64) -> Result<Self, GridError> {
        if in_range(x) && in_range(y) {
            Ok(Self { x, y })
        } else {
            Err(GridError::OutOfRange)
        }
    }

    /// Точка из миллиметров с округлением на сетку.
    pub fn from_mm(x: f64, y: f64) -> Result<Self, GridError> {
        Ok(Self {
            x: mm_to_grid(x)?,
            y: mm_to_grid(y)?,
        })
    }

    /// Координата X в шагах сетки.
    pub fn x(self) -> i64 {
        self.x
    }

    /// Координата Y в шагах сетки.
    pub fn y(self) -> i64 {
        self.y
    }

    /// Точка в миллиметрах.
    pub fn to_mm(self) -> Vec2 {
        Vec2::new(grid_to_mm(self.x), grid_to_mm(self.y))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::testkit::for_cases;

    #[test]
    fn simple_values() {
        assert_eq!(mm_to_grid(1.0), Ok(1000));
        assert_eq!(mm_to_grid(-2.5), Ok(-2500));
        assert_eq!(mm_to_grid(0.0004), Ok(0));
        assert_eq!(mm_to_grid(0.0006), Ok(1));
        assert_eq!(grid_to_mm(1500), 1.5);
        assert_eq!(
            GridPoint::from_mm(12.345, -0.001).map(GridPoint::to_mm),
            Ok(Vec2::new(12.345, -0.001))
        );
    }

    #[test]
    fn rejects_bad_input() {
        assert_eq!(mm_to_grid(f64::NAN), Err(GridError::NotFinite));
        assert_eq!(mm_to_grid(f64::NEG_INFINITY), Err(GridError::NotFinite));
        assert_eq!(mm_to_grid(2.0e9), Err(GridError::OutOfRange));
        assert!(GridPoint::new(GRID_LIMIT, -GRID_LIMIT).is_ok());
        assert_eq!(
            GridPoint::new(GRID_LIMIT + 1, 0),
            Err(GridError::OutOfRange)
        );
        assert_eq!(GridPoint::new(0, i64::MIN), Err(GridError::OutOfRange));
    }

    #[test]
    fn round_trip_through_millimetres() {
        for_cases("round_trip_through_millimetres", 20_000, |rng| {
            let g = if rng.chance(0.2) {
                *rng.pick(&[0, 1, -1, GRID_LIMIT, -GRID_LIMIT])
            } else {
                rng.range_i64(-GRID_LIMIT, GRID_LIMIT)
            };
            assert_eq!(mm_to_grid(grid_to_mm(g)), Ok(g), "g = {g}");
        });
    }
}
