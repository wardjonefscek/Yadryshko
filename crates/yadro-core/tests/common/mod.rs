//! Общее для тестов ядра: короткие конструкторы, случайные области, расстояния.
#![allow(dead_code)]

use yadro_core::testkit::Rng;
use yadro_core::tol::GRID_LIMIT;
use yadro_core::{GridPoint, Region, Ring};

pub fn p(x: i64, y: i64) -> GridPoint {
    GridPoint::new(x, y).unwrap()
}

pub fn rect(x0: i64, y0: i64, x1: i64, y1: i64) -> Region {
    Region::from_ring(Ring::rect(x0, y0, x1, y1).unwrap())
}

pub fn poly(pts: &[(i64, i64)]) -> Region {
    Region::from_ring(Ring::new(pts.iter().map(|&(x, y)| p(x, y)).collect()).unwrap())
}

pub fn rings_text(r: &Region) -> String {
    r.rings()
        .iter()
        .map(|ring| {
            format!(
                "{:?}",
                ring.points()
                    .iter()
                    .map(|q| (q.x(), q.y()))
                    .collect::<Vec<_>>()
            )
        })
        .collect::<Vec<_>>()
        .join("; ")
}

/// Звёздчатый многоугольник вокруг центра.
pub fn star(rng: &mut Rng, cx: i64, cy: i64, rmin: f64, rmax: f64, k: usize) -> Ring {
    let mut ang: Vec<f64> = (0..k)
        .map(|_| rng.range_f64(0.0, std::f64::consts::TAU))
        .collect();
    ang.sort_by(f64::total_cmp);
    let pts: Vec<GridPoint> = ang
        .iter()
        .map(|&t| {
            let r = rng.range_f64(rmin, rmax);
            let clamp = |v: f64| (v.round() as i64).clamp(-GRID_LIMIT, GRID_LIMIT);
            p(
                clamp(cx as f64 + r * t.cos()),
                clamp(cy as f64 + r * t.sin()),
            )
        })
        .collect();
    Ring::new(pts).unwrap_or_else(|_| Ring::rect(cx, cy, cx + 1, cy + 1).unwrap())
}

/// Многоугольник с вершинами на мелкой решётке — много совпадений и наложений.
pub fn lattice_ring(rng: &mut Rng, span: i64, step: i64, k: usize) -> Ring {
    loop {
        let pts: Vec<GridPoint> = (0..k)
            .map(|_| {
                p(
                    rng.range_i64(-span, span) * step,
                    rng.range_i64(-span, span) * step,
                )
            })
            .collect();
        if let Ok(r) = Ring::new(pts) {
            return r;
        }
    }
}

pub fn random_region(rng: &mut Rng, scale: i64) -> Region {
    let c = || 0;
    match rng.below(5) {
        0 => {
            let k = rng.range_i64(3, 12) as usize;
            Region::from_ring(star(rng, c(), c(), scale as f64 * 0.1, scale as f64, k))
        }
        1 => {
            let k = rng.range_i64(3, 8) as usize;
            Region::from_ring(lattice_ring(rng, 4, (scale / 4).max(1), k))
        }
        2 => {
            // Объединение нескольких прямоугольников на решётке.
            let step = (scale / 5).max(1);
            let rings: Vec<Ring> = (0..rng.range_i64(1, 4))
                .map(|_| {
                    let (x0, y0) = (rng.range_i64(-5, 4), rng.range_i64(-5, 4));
                    let (x1, y1) = (x0 + rng.range_i64(1, 4), y0 + rng.range_i64(1, 4));
                    Ring::rect(x0 * step, y0 * step, x1 * step, y1 * step).unwrap()
                })
                .collect();
            Region::from_rings_nonzero(&rings)
        }
        3 => {
            // С отверстием.
            let k = rng.range_i64(4, 10) as usize;
            let outer = Region::from_ring(star(rng, 0, 0, scale as f64 * 0.6, scale as f64, k));
            let inner =
                Region::from_ring(star(rng, 0, 0, scale as f64 * 0.1, scale as f64 * 0.4, k));
            outer.difference(&inner)
        }
        _ => {
            // Произвольный, возможно самопересекающийся контур.
            let k = rng.range_i64(3, 9) as usize;
            Region::from_ring(lattice_ring(
                rng,
                (scale).clamp(1, 8),
                (scale / 8).max(1),
                k,
            ))
        }
    }
}

/// Расстояние от точки до отрезка с точными промежуточными величинами.
pub fn dist(q: GridPoint, a: GridPoint, b: GridPoint) -> f64 {
    let (ax, ay, bx, by, qx, qy) = (
        a.x() as i128,
        a.y() as i128,
        b.x() as i128,
        b.y() as i128,
        q.x() as i128,
        q.y() as i128,
    );
    let (dx, dy) = (bx - ax, by - ay);
    let t = (qx - ax) * dx + (qy - ay) * dy;
    let l2 = dx * dx + dy * dy;
    if l2 == 0 || t <= 0 {
        return ((qx - ax) as f64).hypot((qy - ay) as f64);
    }
    if t >= l2 {
        return ((qx - bx) as f64).hypot((qy - by) as f64);
    }
    let c = (qx - ax) * dy - (qy - ay) * dx;
    c.unsigned_abs() as f64 / (l2 as f64).sqrt()
}

/// Расстояние от точки до границы области (до ближайшего ребра).
pub fn dist_to_boundary(q: GridPoint, r: &Region) -> f64 {
    r.edges()
        .map(|(a, b)| dist(q, a, b))
        .fold(f64::INFINITY, f64::min)
}
