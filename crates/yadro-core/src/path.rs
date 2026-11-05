//! Замкнутый контур из отрезков и дуг — как LWPOLYLINE в DXF: вершины и у каждой
//! выпуклость (bulge) участка до следующей вершины. bulge = 0 — отрезок,
//! bulge = tan(раствор дуги / 4): > 0 — дуга против часовой стрелки, < 0 — по часовой.
//!
//! Площадь и длина считаются по настоящим дугам. Для булевых операций контур
//! превращается в [`Ring`] с допуском хорды (см. [`crate::arc`]).

use crate::arc::Arc;
use crate::region::Region;
use crate::ring::{Ring, RingError};
use crate::vec::Vec2;
use std::f64::consts::FRAC_PI_8;

/// Вершина контура и выпуклость участка до следующей вершины.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PathVertex {
    /// Точка, мм.
    pub at: Vec2,
    /// Выпуклость участка к следующей вершине: 0 — отрезок.
    pub bulge: f64,
}

/// Участок контура.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Segment {
    /// Отрезок из первой точки во вторую.
    Line(Vec2, Vec2),
    /// Дуга.
    Arc(Arc),
}

/// Замкнутый контур из отрезков и дуг.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Path {
    verts: Vec<PathVertex>,
}

impl Path {
    /// Пустой контур.
    pub fn new() -> Self {
        Self::default()
    }

    /// Многоугольник: все участки — отрезки.
    pub fn polygon(points: &[Vec2]) -> Self {
        Self {
            verts: points
                .iter()
                .map(|&at| PathVertex { at, bulge: 0.0 })
                .collect(),
        }
    }

    /// Добавить вершину; `bulge` — выпуклость участка от неё к следующей.
    pub fn push(&mut self, at: Vec2, bulge: f64) -> &mut Self {
        self.verts.push(PathVertex { at, bulge });
        self
    }

    /// Окружность против часовой стрелки (две полуокружности).
    pub fn circle(center: Vec2, r: f64) -> Self {
        let mut p = Self::new();
        p.push(center + Vec2::new(r, 0.0), 1.0)
            .push(center - Vec2::new(r, 0.0), 1.0);
        p
    }

    /// Прямоугольник со скруглёнными углами радиуса `r` (не больше половины
    /// меньшей стороны), против часовой стрелки.
    pub fn rounded_rect(x0: f64, y0: f64, x1: f64, y1: f64, r: f64) -> Self {
        let (xa, xb, ya, yb) = (x0.min(x1), x0.max(x1), y0.min(y1), y0.max(y1));
        let r = r.clamp(0.0, 0.5 * (xb - xa).min(yb - ya));
        if r == 0.0 {
            return Self::polygon(&[
                Vec2::new(xa, ya),
                Vec2::new(xb, ya),
                Vec2::new(xb, yb),
                Vec2::new(xa, yb),
            ]);
        }
        let q = FRAC_PI_8.tan(); // четверть окружности
        let mut p = Self::new();
        p.push(Vec2::new(xa + r, ya), 0.0)
            .push(Vec2::new(xb - r, ya), q)
            .push(Vec2::new(xb, ya + r), 0.0)
            .push(Vec2::new(xb, yb - r), q)
            .push(Vec2::new(xb - r, yb), 0.0)
            .push(Vec2::new(xa + r, yb), q)
            .push(Vec2::new(xa, yb - r), 0.0)
            .push(Vec2::new(xa, ya + r), q);
        p
    }

    /// Паз (стадион) вокруг отрезка `a`–`b` радиуса `r`, против часовой стрелки.
    pub fn slot(a: Vec2, b: Vec2, r: f64) -> Self {
        let Some(u) = (b - a).normalized() else {
            return Self::circle(a, r);
        };
        let n = u.perp() * r;
        let mut p = Self::new();
        p.push(a - n, 0.0)
            .push(b - n, 1.0)
            .push(b + n, 0.0)
            .push(a + n, 1.0);
        p
    }

    /// Вершины.
    pub fn vertices(&self) -> &[PathVertex] {
        &self.verts
    }

    /// Число вершин.
    pub fn len(&self) -> usize {
        self.verts.len()
    }

    /// Пустой ли контур.
    pub fn is_empty(&self) -> bool {
        self.verts.is_empty()
    }

    /// Участки по порядку; последний замыкает контур.
    pub fn segments(&self) -> impl Iterator<Item = Segment> + '_ {
        let n = self.verts.len();
        (0..n).map(move |i| {
            let (v, w) = (self.verts[i], self.verts[(i + 1) % n]);
            match Arc::from_bulge(v.at, w.at, v.bulge) {
                Some(arc) => Segment::Arc(arc),
                None => Segment::Line(v.at, w.at),
            }
        })
    }

    /// Ориентированная площадь по настоящим дугам, мм²: > 0 — против часовой.
    pub fn area(&self) -> f64 {
        let o = self.verts.first().map_or(Vec2::ZERO, |v| v.at);
        self.segments()
            .map(|s| match s {
                Segment::Line(a, b) => 0.5 * (a - o).cross(b - o),
                Segment::Arc(arc) => {
                    let (a, b) = (arc.start_point(), arc.end_point());
                    0.5 * (a - o).cross(b - o) + arc.segment_area()
                }
            })
            .sum()
    }

    /// Длина, мм.
    pub fn length(&self) -> f64 {
        self.segments()
            .map(|s| match s {
                Segment::Line(a, b) => a.distance(b),
                Segment::Arc(arc) => arc.length(),
            })
            .sum()
    }

    /// Тот же контур с обратным обходом.
    pub fn reversed(&self) -> Self {
        let n = self.verts.len();
        Self {
            verts: (0..n)
                .map(|j| PathVertex {
                    at: self.verts[n - 1 - j].at,
                    bulge: -self.verts[(2 * n - 2 - j) % n].bulge,
                })
                .collect(),
        }
    }

    /// Ломаная: вершины контура как есть, на дугах — промежуточные точки с прогибом
    /// хорд не больше `tol_mm`. Без повтора первой точки в конце.
    pub fn flatten(&self, tol_mm: f64) -> Vec<Vec2> {
        let mut out = Vec::new();
        for (v, s) in self.verts.iter().zip(self.segments()) {
            out.push(v.at);
            if let Segment::Arc(arc) = s {
                let pts = arc.flatten(tol_mm);
                out.extend_from_slice(&pts[1..pts.len() - 1]);
            }
        }
        out
    }

    /// Контур на сетке с допуском хорды `tol_mm`.
    pub fn to_ring(&self, tol_mm: f64) -> Result<Ring, RingError> {
        Ring::from_mm(&self.flatten(tol_mm))
    }

    /// Область, ограниченная контуром (по правилу ненулевой обмотки).
    pub fn to_region(&self, tol_mm: f64) -> Result<Region, RingError> {
        Ok(Region::from_ring(self.to_ring(tol_mm)?))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::PI;

    #[test]
    fn exact_areas_and_lengths() {
        let c = Path::circle(Vec2::new(5.0, 5.0), 2.0);
        assert!((c.area() - PI * 4.0).abs() < 1e-12);
        assert!((c.length() - 4.0 * PI).abs() < 1e-12);
        assert!((c.reversed().area() + PI * 4.0).abs() < 1e-12);
        let rr = Path::rounded_rect(0.0, 0.0, 10.0, 6.0, 1.5);
        assert!((rr.area() - (60.0 - (4.0 - PI) * 2.25)).abs() < 1e-12);
        let s = Path::slot(Vec2::ZERO, Vec2::new(8.0, 0.0), 1.0);
        assert!((s.area() - (16.0 + PI)).abs() < 1e-12);
        assert!((s.length() - (16.0 + 2.0 * PI)).abs() < 1e-12);
    }

    #[test]
    fn reversed_keeps_shape() {
        let rr = Path::rounded_rect(0.0, 0.0, 10.0, 6.0, 1.5);
        let back = rr.reversed();
        assert!((back.area() + rr.area()).abs() < 1e-12);
        assert!((back.length() - rr.length()).abs() < 1e-12);
        assert_eq!(back.reversed(), rr);
    }

    #[test]
    fn flattened_area_close_to_exact() {
        let rr = Path::rounded_rect(0.0, 0.0, 10.0, 6.0, 1.5);
        for tol in [0.0005, 0.002, 0.05] {
            let ring = rr.to_ring(tol).unwrap();
            let flat = ring.area_mm2();
            // Хорды внутри дуг: площадь меньше, но не больше чем на длину × допуск
            // (плюс округление на сетку).
            let diff = rr.area() - flat;
            assert!(diff > -rr.length() * 0.000_71, "tol {tol}: {diff}");
            assert!(diff < rr.length() * (tol + 0.000_71), "tol {tol}: {diff}");
        }
    }
}
