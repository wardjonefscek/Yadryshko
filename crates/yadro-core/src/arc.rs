//! Дуги окружностей: по выпуклости (bulge, как в DXF), по трём точкам, по центру;
//! разбиение на хорды с заданным допуском.
//!
//! Дуги — строительная геометрия в миллиметрах (`f64`). В булевы операции и offset
//! они попадают ломаной: вершины на окружности, прогиб каждой хорды (стрелка) не
//! больше допуска. Вместе с округлением на сетку ломаная отходит от настоящей дуги
//! не больше чем на допуск + 0,71 мкм и почти целиком лежит внутри круга.

use crate::tol::{MAX_ARC_CHORDS, MIN_CHORD_TOL_MM};
use crate::vec::Vec2;
use std::f64::consts::{FRAC_PI_2, TAU};

/// Дуга: центр, радиус, начальный угол и угол раствора (радианы; > 0 — против
/// часовой стрелки). |раствор| не больше 2π.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Arc {
    /// Центр, мм.
    pub center: Vec2,
    /// Радиус, мм.
    pub radius: f64,
    /// Угол начальной точки, рад.
    pub start: f64,
    /// Раствор, рад: > 0 — против часовой стрелки.
    pub sweep: f64,
}

/// Сколько хорд нужно дуге радиуса `radius` с раствором `sweep`, чтобы прогиб каждой
/// был не больше `tol` (единицы у `radius` и `tol` одни и те же). Хорда не длиннее
/// четверти окружности; не больше [`MAX_ARC_CHORDS`].
pub fn chord_count(radius: f64, sweep: f64, tol: f64) -> usize {
    let step = if tol >= radius {
        FRAC_PI_2
    } else {
        (2.0 * (1.0 - tol / radius).acos()).min(FRAC_PI_2)
    };
    let n = (sweep.abs() / step).ceil();
    if n.is_finite() {
        (n as usize).clamp(1, MAX_ARC_CHORDS)
    } else {
        MAX_ARC_CHORDS
    }
}

impl Arc {
    /// Дуга с проверкой: радиус больше нуля, раствор не ноль и не больше 2π по модулю.
    pub fn new(center: Vec2, radius: f64, start: f64, sweep: f64) -> Option<Self> {
        let ok = center.is_finite()
            && radius.is_finite()
            && radius > 0.0
            && start.is_finite()
            && sweep.is_finite()
            && sweep != 0.0
            && sweep.abs() <= TAU;
        ok.then_some(Self {
            center,
            radius,
            start,
            sweep,
        })
    }

    /// Полная окружность против часовой стрелки, начало — справа от центра.
    pub fn circle(center: Vec2, radius: f64) -> Option<Self> {
        Self::new(center, radius, 0.0, TAU)
    }

    /// Дуга из `a` в `b` с выпуклостью `bulge` = tan(раствор / 4), как в DXF:
    /// больше нуля — против часовой стрелки, 1 — полуокружность. `None` для bulge = 0
    /// (это отрезок) и для совпадающих точек.
    pub fn from_bulge(a: Vec2, b: Vec2, bulge: f64) -> Option<Self> {
        let c = a.distance(b);
        if !bulge.is_finite() || bulge == 0.0 || c == 0.0 || !c.is_finite() {
            return None;
        }
        let sweep = 4.0 * bulge.atan();
        let radius = c / (2.0 * (sweep / 2.0).sin().abs());
        let u = (b - a) / c;
        // Центр — на перпендикуляре к хорде через её середину.
        let h = c * (1.0 - bulge * bulge) / (4.0 * bulge);
        let center = (a + b) * 0.5 + u.perp() * h;
        let d = a - center;
        Self::new(center, radius, d.y.atan2(d.x), sweep)
    }

    /// Дуга из `a` через `m` в `b`. `None`, если точки на одной прямой или `a` = `b`.
    pub fn through(a: Vec2, m: Vec2, b: Vec2) -> Option<Self> {
        let (p, q) = (m - a, b - a);
        let d = 2.0 * p.cross(q);
        let scale = p.length_sq().max(q.length_sq());
        if a == b || d.abs() <= 1e-12 * scale {
            return None;
        }
        let (pl, ql) = (p.length_sq(), q.length_sq());
        let center = a + Vec2::new((q.y * pl - p.y * ql) / d, (p.x * ql - q.x * pl) / d);
        let ang = |v: Vec2| (v - center).y.atan2((v - center).x);
        let (sa, sb) = (ang(a), ang(b));
        let sweep = if d > 0.0 {
            (sb - sa).rem_euclid(TAU)
        } else {
            -(sa - sb).rem_euclid(TAU)
        };
        Self::new(center, a.distance(center), sa, sweep)
    }

    /// Точка на угле `start + t · sweep`, `t` от 0 до 1.
    pub fn point(&self, t: f64) -> Vec2 {
        let ang = self.start + t * self.sweep;
        self.center + Vec2::new(ang.cos(), ang.sin()) * self.radius
    }

    /// Начальная точка.
    pub fn start_point(&self) -> Vec2 {
        self.point(0.0)
    }

    /// Конечная точка.
    pub fn end_point(&self) -> Vec2 {
        self.point(1.0)
    }

    /// Длина дуги, мм.
    pub fn length(&self) -> f64 {
        self.radius * self.sweep.abs()
    }

    /// Выпуклость для DXF: tan(раствор / 4).
    pub fn bulge(&self) -> f64 {
        (self.sweep / 4.0).tan()
    }

    /// Против часовой ли стрелки.
    pub fn is_ccw(&self) -> bool {
        self.sweep > 0.0
    }

    /// Та же дуга в обратную сторону.
    pub fn reversed(&self) -> Self {
        Self {
            start: self.start + self.sweep,
            sweep: -self.sweep,
            ..*self
        }
    }

    /// Ориентированная площадь между хордой и дугой: > 0, если дуга против часовой
    /// стрелки (тогда она выпирает вправо от хорды — наружу у контура против часовой).
    pub fn segment_area(&self) -> f64 {
        let t = self.sweep.abs();
        0.5 * self.radius * self.radius * (t - t.sin()) * self.sweep.signum()
    }

    /// Число хорд для допуска `tol_mm` (не меньше [`MIN_CHORD_TOL_MM`]).
    pub fn chord_count(&self, tol_mm: f64) -> usize {
        chord_count(self.radius, self.sweep, tol_mm.max(MIN_CHORD_TOL_MM))
    }

    /// Ломаная по дуге: `n + 1` точек на окружности с обоими концами, прогиб хорд не
    /// больше допуска.
    pub fn flatten(&self, tol_mm: f64) -> Vec<Vec2> {
        let n = self.chord_count(tol_mm);
        (0..=n).map(|i| self.point(i as f64 / n as f64)).collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::testkit::for_cases;
    use std::f64::consts::PI;

    fn close(a: Vec2, b: Vec2) -> bool {
        a.distance(b) < 1e-9
    }

    #[test]
    fn semicircle_from_bulge() {
        let a = Arc::from_bulge(Vec2::new(1.0, 0.0), Vec2::new(-1.0, 0.0), 1.0).unwrap();
        assert!(close(a.center, Vec2::ZERO));
        assert!((a.radius - 1.0).abs() < 1e-12 && (a.sweep - PI).abs() < 1e-12);
        assert!(
            close(a.point(0.5), Vec2::new(0.0, 1.0)),
            "дуга проходит сверху"
        );
        assert!((a.segment_area() - PI / 2.0).abs() < 1e-12);
        assert!(Arc::from_bulge(Vec2::ZERO, Vec2::ZERO, 1.0).is_none());
        assert!(Arc::from_bulge(Vec2::ZERO, Vec2::new(1.0, 0.0), 0.0).is_none());
    }

    #[test]
    fn bulge_round_trip_and_endpoints() {
        for_cases("bulge_round_trip_and_endpoints", 2000, |rng| {
            // Размеры деталей: точки в пределах метра, хорды от микрона до 100 мм.
            let a = Vec2::new(
                rng.range_f64(-1000.0, 1000.0),
                rng.range_f64(-1000.0, 1000.0),
            );
            let len = 10f64.powf(rng.range_f64(-3.0, 2.0));
            let t = rng.range_f64(0.0, TAU);
            let b = a + Vec2::new(t.cos(), t.sin()) * len;
            let bulge =
                10f64.powf(rng.range_f64(-3.0, 1.3)) * if rng.chance(0.5) { 1.0 } else { -1.0 };
            let Some(arc) = Arc::from_bulge(a, b, bulge) else {
                return;
            };
            let scale = 1e-9 * (1.0 + a.length() + arc.radius);
            assert!(arc.start_point().distance(a) < scale);
            assert!(arc.end_point().distance(b) < scale);
            assert!((arc.bulge() - bulge).abs() < 1e-9 * (1.0 + bulge.abs()).powi(2));
            let back = Arc::through(a, arc.point(0.5), b).unwrap();
            assert!((back.sweep - arc.sweep).abs() < 1e-6, "{back:?} / {arc:?}");
        });
    }

    #[test]
    fn chords_within_tolerance() {
        for_cases("chords_within_tolerance", 2000, |rng| {
            let r = 10f64.powf(rng.range_f64(-2.0, 4.0));
            let sweep = rng.range_f64(-TAU, TAU);
            let Some(arc) = Arc::new(Vec2::new(3.0, -7.0), r, rng.range_f64(-PI, PI), sweep) else {
                return;
            };
            let tol = 10f64.powf(rng.range_f64(-3.5, -1.0));
            let pts = arc.flatten(tol);
            let tol = tol.max(MIN_CHORD_TOL_MM);
            for w in pts.windows(2) {
                let mid = (w[0] + w[1]) * 0.5;
                let sag = r - mid.distance(arc.center);
                assert!(sag <= tol * (1.0 + 1e-9) + 1e-12, "прогиб {sag} > {tol}");
                assert!(w[0].distance(w[1]) <= r * 2f64.sqrt() * (1.0 + 1e-9));
            }
        });
    }

    #[test]
    fn three_points_orientation() {
        let ccw = Arc::through(
            Vec2::new(1.0, 0.0),
            Vec2::new(0.0, 1.0),
            Vec2::new(-1.0, 0.0),
        )
        .unwrap();
        assert!(ccw.is_ccw() && (ccw.sweep - PI).abs() < 1e-12);
        let cw = Arc::through(
            Vec2::new(1.0, 0.0),
            Vec2::new(0.0, -1.0),
            Vec2::new(-1.0, 0.0),
        )
        .unwrap();
        assert!(!cw.is_ccw() && (cw.sweep + PI).abs() < 1e-12);
        let big = Arc::through(
            Vec2::new(1.0, 0.0),
            Vec2::new(-1.0, 0.0),
            Vec2::new(0.0, -1.0),
        )
        .unwrap();
        assert!(
            (big.sweep - 1.5 * PI).abs() < 1e-12,
            "через дальнюю сторону — три четверти"
        );
        assert!(Arc::through(Vec2::ZERO, Vec2::new(1.0, 1.0), Vec2::new(2.0, 2.0)).is_none());
    }
}
