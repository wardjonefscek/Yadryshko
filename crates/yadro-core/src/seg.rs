//! Точные операции с отрезками на 2D-сетке: лежит ли точка на отрезке,
//! как пересекаются два отрезка, где точка пересечения, через какие
//! пиксели проходит отрезок.
//!
//! Всё считается в целых числах без допусков. Координаты сетки не больше
//! 2^40 (см. [`crate::tol::GRID_LIMIT`]), поэтому хватает `i128`.

use crate::grid::GridPoint;
use crate::predicates::{Sign, orient2d_grid};
use std::cmp::Ordering;

/// Лежит ли точка `p` на замкнутом отрезке `ab` (включая концы).
pub fn on_segment(p: GridPoint, a: GridPoint, b: GridPoint) -> bool {
    orient2d_grid(a, b, p) == Sign::Zero && in_box(p, a, b)
}

/// Лежит ли точка `p` строго внутри отрезка `ab` (не в концах).
pub fn in_segment_interior(p: GridPoint, a: GridPoint, b: GridPoint) -> bool {
    p != a && p != b && on_segment(p, a, b)
}

fn in_box(p: GridPoint, a: GridPoint, b: GridPoint) -> bool {
    a.x().min(b.x()) <= p.x()
        && p.x() <= a.x().max(b.x())
        && a.y().min(b.y()) <= p.y()
        && p.y() <= a.y().max(b.y())
}

/// Как пересекаются два отрезка.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Crossing {
    /// Общих точек нет.
    None,
    /// Одна общая точка внутри обоих отрезков (настоящее пересечение).
    Proper,
    /// Одна общая точка, и это конец хотя бы одного отрезка.
    Touch,
    /// Отрезки на одной прямой и общая часть имеет ненулевую длину.
    Overlap,
}

/// Классифицировать пересечение отрезков `ab` и `cd`.
pub fn classify(a: GridPoint, b: GridPoint, c: GridPoint, d: GridPoint) -> Crossing {
    let o1 = orient2d_grid(a, b, c);
    let o2 = orient2d_grid(a, b, d);
    let o3 = orient2d_grid(c, d, a);
    let o4 = orient2d_grid(c, d, b);
    if o1 == Sign::Zero && o2 == Sign::Zero {
        // На одной прямой (или один из отрезков вырожден в точку).
        return collinear_overlap(a, b, c, d);
    }
    if o1 != Sign::Zero && o2 != Sign::Zero && o3 != Sign::Zero && o4 != Sign::Zero {
        return if o1 != o2 && o3 != o4 {
            Crossing::Proper
        } else {
            Crossing::None
        };
    }
    if on_segment(c, a, b) || on_segment(d, a, b) || on_segment(a, c, d) || on_segment(b, c, d) {
        Crossing::Touch
    } else {
        Crossing::None
    }
}

fn collinear_overlap(a: GridPoint, b: GridPoint, c: GridPoint, d: GridPoint) -> Crossing {
    // Проекция на ось с наибольшим разбросом.
    let use_x = (a.x() - b.x()).abs().max((c.x() - d.x()).abs())
        >= (a.y() - b.y()).abs().max((c.y() - d.y()).abs());
    let k = |p: GridPoint| if use_x { p.x() } else { p.y() };
    let (a0, a1) = (k(a).min(k(b)), k(a).max(k(b)));
    let (c0, c1) = (k(c).min(k(d)), k(c).max(k(d)));
    let lo = a0.max(c0);
    let hi = a1.min(c1);
    match lo.cmp(&hi) {
        Ordering::Less => Crossing::Overlap,
        Ordering::Equal => {
            // Одна общая точка проекции; проверим, что точки действительно на одной прямой.
            if on_segment(c, a, b)
                || on_segment(d, a, b)
                || on_segment(a, c, d)
                || on_segment(b, c, d)
            {
                Crossing::Touch
            } else {
                Crossing::None
            }
        }
        Ordering::Greater => Crossing::None,
    }
}

/// Точка с рациональными координатами `(x / den, y / den)`, `den > 0`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct RatPoint {
    /// Числитель x.
    pub x: i128,
    /// Числитель y.
    pub y: i128,
    /// Общий знаменатель, больше нуля.
    pub den: i128,
}

impl RatPoint {
    /// Пиксель сетки, в который попадает точка: полуоткрытый квадрат
    /// `[g − ½, g + ½)` по каждой оси, то есть `g = ⌊v + ½⌋`.
    pub fn pixel(self) -> GridPoint {
        let r = |num: i128| (2 * num + self.den).div_euclid(2 * self.den) as i64;
        GridPoint::new(r(self.x), r(self.y))
            .expect("точка пересечения лежит внутри габарита отрезков")
    }
}

/// Точка настоящего пересечения отрезков `ab` и `cd` (None — если они параллельны).
///
/// РИСК: оценки разрядности. Разности координат до 2^41, знаменатель
/// (векторное произведение) до 2^83, числители x и y до 2^125 — `i128` хватает,
/// и в [`RatPoint::pixel`] удвоенный числитель (до 2^126) тоже помещается.
pub fn intersection(a: GridPoint, b: GridPoint, c: GridPoint, d: GridPoint) -> Option<RatPoint> {
    let (ax, ay) = (i128::from(a.x()), i128::from(a.y()));
    let (rx, ry) = (i128::from(b.x()) - ax, i128::from(b.y()) - ay);
    let (sx, sy) = (i128::from(d.x() - c.x()), i128::from(d.y() - c.y()));
    let (wx, wy) = (i128::from(c.x()) - ax, i128::from(c.y()) - ay);
    let mut den = rx * sy - ry * sx;
    if den == 0 {
        return None;
    }
    let mut t = wx * sy - wy * sx;
    if den < 0 {
        den = -den;
        t = -t;
    }
    Some(RatPoint {
        x: ax * den + t * rx,
        y: ay * den + t * ry,
        den,
    })
}

/// Граница параметра t на отрезке: дробь `num / den` (`den > 0`), строгая или нет.
#[derive(Clone, Copy, Debug)]
struct Bound {
    num: i128,
    den: i128,
    strict: bool,
}

fn frac_cmp(a: Bound, b: Bound) -> Ordering {
    (a.num * b.den).cmp(&(b.num * a.den))
}

/// Проходит ли отрезок `a → b` через пиксель с центром `c` — полуоткрытый квадрат
/// `[cx − ½, cx + ½) × [cy − ½, cy + ½)`, как в snap rounding Хобби.
///
/// Если проходит, возвращает параметр входа в пиксель как дробь `(num, den)` и признак
/// «вход сразу после этого значения»: по ним пиксели упорядочиваются вдоль отрезка.
pub fn hits_pixel(a: GridPoint, b: GridPoint, c: GridPoint) -> Option<(i128, i128, bool)> {
    // Удвоенные координаты: границы пикселя становятся целыми (нечётными).
    let (ax, ay) = (2 * i128::from(a.x()), 2 * i128::from(a.y()));
    let (dx, dy) = (2 * i128::from(b.x() - a.x()), 2 * i128::from(b.y() - a.y()));
    let (cx, cy) = (2 * i128::from(c.x()), 2 * i128::from(c.y()));
    let mut lo = Bound {
        num: 0,
        den: 1,
        strict: false,
    };
    let mut hi = Bound {
        num: 1,
        den: 1,
        strict: false,
    };
    for (p0, d, l, r) in [(ax, dx, cx - 1, cx + 1), (ay, dy, cy - 1, cy + 1)] {
        // Условие l ≤ p0 + t·d < r.
        if d == 0 {
            if !(l <= p0 && p0 < r) {
                return None;
            }
            continue;
        }
        let (low, high) = if d > 0 {
            (
                Bound {
                    num: l - p0,
                    den: d,
                    strict: false,
                },
                Bound {
                    num: r - p0,
                    den: d,
                    strict: true,
                },
            )
        } else {
            (
                Bound {
                    num: p0 - r,
                    den: -d,
                    strict: true,
                },
                Bound {
                    num: p0 - l,
                    den: -d,
                    strict: false,
                },
            )
        };
        lo = match frac_cmp(low, lo) {
            Ordering::Greater => low,
            Ordering::Equal => Bound {
                strict: lo.strict || low.strict,
                ..lo
            },
            Ordering::Less => lo,
        };
        hi = match frac_cmp(high, hi) {
            Ordering::Less => high,
            Ordering::Equal => Bound {
                strict: hi.strict || high.strict,
                ..hi
            },
            Ordering::Greater => hi,
        };
    }
    match frac_cmp(lo, hi) {
        Ordering::Less => Some((lo.num, lo.den, lo.strict)),
        Ordering::Equal if !lo.strict && !hi.strict => Some((lo.num, lo.den, false)),
        _ => None,
    }
}

/// Сравнить два параметра входа из [`hits_pixel`].
pub fn entry_cmp(a: (i128, i128, bool), b: (i128, i128, bool)) -> Ordering {
    (a.0 * b.1).cmp(&(b.0 * a.1)).then(a.2.cmp(&b.2))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::testkit::for_cases;

    fn p(x: i64, y: i64) -> GridPoint {
        GridPoint::new(x, y).unwrap()
    }

    #[test]
    fn classify_basic_cases() {
        assert_eq!(
            classify(p(0, 0), p(4, 4), p(0, 4), p(4, 0)),
            Crossing::Proper
        );
        assert_eq!(
            classify(p(0, 0), p(4, 0), p(2, 0), p(2, 3)),
            Crossing::Touch
        );
        assert_eq!(
            classify(p(0, 0), p(4, 0), p(4, 0), p(6, 2)),
            Crossing::Touch
        );
        assert_eq!(
            classify(p(0, 0), p(4, 0), p(2, 0), p(6, 0)),
            Crossing::Overlap
        );
        assert_eq!(
            classify(p(0, 0), p(4, 0), p(4, 0), p(6, 0)),
            Crossing::Touch
        );
        assert_eq!(classify(p(0, 0), p(4, 0), p(5, 0), p(6, 0)), Crossing::None);
        assert_eq!(classify(p(0, 0), p(4, 0), p(0, 1), p(4, 1)), Crossing::None);
        assert_eq!(
            classify(p(0, 0), p(0, 4), p(0, 1), p(0, 3)),
            Crossing::Overlap
        );
    }

    #[test]
    fn intersection_rounds_to_pixel() {
        // (0,0)–(3,1) и (0,1)–(3,0) пересекаются в (1,5; 0,5) → пиксель (2, 1): ровно на середине — вверх.
        let q = intersection(p(0, 0), p(3, 1), p(0, 1), p(3, 0)).unwrap();
        assert_eq!((q.x * 2, q.y * 2), (3 * q.den, q.den));
        assert_eq!(q.pixel(), p(2, 1));
    }

    #[test]
    fn pixel_hits_are_half_open() {
        // Горизонтальный отрезок y = 0 проходит через пиксели (k, 0).
        assert!(hits_pixel(p(0, 0), p(4, 0), p(2, 0)).is_some());
        assert!(hits_pixel(p(0, 0), p(4, 0), p(2, 1)).is_none());
        // Диагональ (0,0)–(1,1) проходит ровно через общий угол четырёх пикселей (½, ½).
        // Пиксели полуоткрыты, поэтому угол принадлежит только пикселю (1, 1):
        // соседние по диагонали пиксели (1, 0) и (0, 1) отрезок не задевает.
        assert!(hits_pixel(p(0, 0), p(1, 1), p(0, 0)).is_some());
        assert!(hits_pixel(p(0, 0), p(1, 1), p(1, 1)).is_some());
        assert!(hits_pixel(p(0, 0), p(1, 1), p(1, 0)).is_none());
        assert!(hits_pixel(p(0, 0), p(1, 1), p(0, 1)).is_none());
        // Порядок входа вдоль отрезка.
        let e0 = hits_pixel(p(0, 0), p(4, 0), p(1, 0)).unwrap();
        let e1 = hits_pixel(p(0, 0), p(4, 0), p(3, 0)).unwrap();
        assert_eq!(entry_cmp(e0, e1), Ordering::Less);
    }

    #[test]
    fn proper_crossing_point_lies_on_both_segments() {
        for_cases(
            "proper_crossing_point_lies_on_both_segments",
            5_000,
            |rng| {
                let lim = if rng.chance(0.5) {
                    20
                } else {
                    crate::tol::GRID_LIMIT
                };
                let mut r = || p(rng.range_i64(-lim, lim), rng.range_i64(-lim, lim));
                let (a, b, c, d) = (r(), r(), r(), r());
                if classify(a, b, c, d) != Crossing::Proper {
                    return;
                }
                let q = intersection(a, b, c, d).unwrap();
                // Точка на прямой ab: (b − a) × (q − a) = 0 в рациональных числах.
                let on = |u: GridPoint, v: GridPoint| {
                    let (ux, uy) = (i128::from(u.x()) * q.den, i128::from(u.y()) * q.den);
                    let (vx, vy) = (i128::from(v.x() - u.x()), i128::from(v.y() - u.y()));
                    // Делим на den заранее нельзя; проверяем знак без переполнения через деление с остатком.
                    let (qx, qy) = (q.x - ux, q.y - uy);
                    // qx/den, qy/den — смещение от u; коллинеарность: vx·qy − vy·qx = 0.
                    vx.checked_mul(qy)
                        .zip(vy.checked_mul(qx))
                        .map(|(m1, m2)| m1 == m2)
                };
                if let (Some(x), Some(y)) = (on(a, b), on(c, d)) {
                    assert!(
                        x && y,
                        "точка пересечения не на отрезках: {a:?} {b:?} {c:?} {d:?}"
                    );
                }
                // Пиксель пересечения задевают оба отрезка.
                let px = q.pixel();
                assert!(hits_pixel(a, b, px).is_some() && hits_pixel(c, d, px).is_some());
            },
        );
    }
}
