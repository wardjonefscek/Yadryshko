//! Восстановление дуг: ломаная контура → отрезки и дуги с допуском.
//!
//! Булевы операции, эквидистанта и контур из плотности дают ломаные: дуга там —
//! сотни коротких хорд. Для фрезы и для DXF лучше настоящие дуги: G-кода меньше,
//! станок идёт ровнее. Подгонка жадная: от каждой вершины берётся самый длинный
//! отрезок или дуга, которые укладываются в допуск, — и так до конца контура.
//!
//! Гарантия: каждая точка ломаной не дальше `tol` от подогнанного контура, и
//! наоборот (расстояние Хаусдорфа ≤ `tol`). Почему хватает проверять вершины и рёбра:
//! - отрезок: все вершины ломаной в полосе ширины 2·tol вокруг него, а полоса
//!   выпуклая, значит, в ней и рёбра; ломаная идёт от одного конца отрезка к
//!   другому внутри полосы и пересекает каждый перпендикуляр к нему;
//! - дуга: все вершины в кольце R ± tol, их углы идут по порядку в сторону дуги,
//!   а каждое ребро не заходит внутрь кольца глубже R − tol (проверяется ближайшая
//!   к центру точка ребра); значит, ломаная пересекает каждый луч из центра в
//!   пределах дуги внутри кольца.
//!
//! Вершины подогнанного контура — вершины исходной ломаной (кроме целой
//! окружности: она становится двумя полуокружностями с точным центром, а в DXF —
//! объектом CIRCLE). Вершины, общие для
//! нескольких контуров области или лежащие на чужом контуре, сохраняются, чтобы
//! касания остались на месте. Контуры, которые ближе друг к другу чем 2·tol, после
//! подгонки могут задеть друг друга на долю допуска — это ниже точности станка.
//! Деталь тоньше допуска подгонка не трогает: если контур схлопнулся бы или
//! потерял больше половины площади, остаётся исходная ломаная.

use crate::arc::Arc;
use crate::grid::GridPoint;
use crate::path::Path;
use crate::region::Region;
use crate::ring::Ring;
use crate::vec::Vec2;
use std::collections::BTreeMap;
use std::f64::consts::{PI, TAU};

/// Наибольший раствор одной дуги: чуть больше полуокружности. Полная окружность —
/// две дуги, так её понимают все CAM.
const MAX_SWEEP: f64 = PI + 0.1;

/// Расстояние от точки до отрезка.
fn seg_dist(p: Vec2, a: Vec2, b: Vec2) -> f64 {
    let ab = b - a;
    let l2 = ab.length_sq();
    if l2 == 0.0 {
        return p.distance(a);
    }
    let t = ((p - a).dot(ab) / l2).clamp(0.0, 1.0);
    p.distance(a + ab * t)
}

/// Все внутренние точки в полосе tol вокруг отрезка первая–последняя.
fn line_ok(s: &[Vec2], tol: f64) -> bool {
    let (a, b) = (s[0], s[s.len() - 1]);
    a != b && s[1..s.len() - 1].iter().all(|&p| seg_dist(p, a, b) <= tol)
}

/// Дуга через первую, среднюю и последнюю точки, если в неё укладываются все.
fn arc_fit(s: &[Vec2], tol: f64) -> Option<Arc> {
    if s.len() < 3 {
        return None;
    }
    let arc = Arc::through(s[0], s[s.len() / 2], s[s.len() - 1])?;
    if arc.sweep.abs() > MAX_SWEEP || !arc.radius.is_finite() || arc.radius <= tol {
        return None;
    }
    let (c, r) = (arc.center, arc.radius);
    let dir = arc.sweep.signum();
    let mut prev = -1e-12;
    for (k, &p) in s.iter().enumerate() {
        if (p.distance(c) - r).abs() > tol {
            return None;
        }
        let d = p - c;
        let mut phi = (dir * (d.y.atan2(d.x) - arc.start)).rem_euclid(TAU);
        if k == s.len() - 1 {
            phi = arc.sweep.abs();
        } else if k > 0 && phi > arc.sweep.abs() + 1e-12 {
            return None; // вышла за конец дуги
        }
        if k > 0 && phi <= prev {
            return None; // углы должны идти по порядку
        }
        prev = phi;
    }
    // Рёбра не заходят внутрь кольца глубже допуска.
    let inward_ok = s.windows(2).all(|w| r - seg_dist(c, w[0], w[1]) <= tol);
    inward_ok.then_some(arc)
}

/// Угол поворота в вершине, рад (по модулю).
fn turn(a: Vec2, b: Vec2, c: Vec2) -> f64 {
    let (u, v) = (b - a, c - b);
    u.cross(v).atan2(u.dot(v)).abs()
}

/// Вся ломаная — окружность в пределах допуска? Тогда (центр, радиус).
/// Центр — по трём далеко разнесённым вершинам, проверка — как у дуги: вершины в
/// кольце R ± tol, рёбра не глубже R − tol, обход ровно один оборот в одну сторону.
fn circle_fit(pts: &[Vec2], tol: f64) -> Option<(Vec2, f64)> {
    let n = pts.len();
    if n < 8 {
        return None;
    }
    let arc = Arc::through(pts[0], pts[n / 3], pts[2 * n / 3])?;
    let (c, r) = (arc.center, arc.radius);
    if r <= tol {
        return None;
    }
    let mut total = 0.0;
    let mut sign = 0.0;
    for k in 0..n {
        let (a, b) = (pts[k], pts[(k + 1) % n]);
        if (a.distance(c) - r).abs() > tol || r - seg_dist(c, a, b) > tol {
            return None;
        }
        let step = (a - c).cross(b - c).atan2((a - c).dot(b - c));
        if step == 0.0 || (sign != 0.0 && step.signum() != sign) {
            return None;
        }
        sign = step.signum();
        total += step;
    }
    ((total.abs() - TAU).abs() < 1e-6).then_some((c, r))
}

/// Подогнать замкнутую ломаную (мм); `pinned[k]` — вершину k нельзя убирать.
fn fit_points(pts: &[Vec2], pinned: &[bool], tol: f64) -> Path {
    let n = pts.len();
    if n <= 3 {
        return Path::polygon(pts);
    }
    // Целая окружность — две полуокружности (в DXF это CIRCLE). Вершины такого
    // контура уже не из ломаной, поэтому только если ни одна не закреплена.
    if !pinned.contains(&true)
        && let Some((c, r)) = circle_fit(pts, tol)
    {
        let p = Path::circle(c, r);
        let ccw = pts
            .iter()
            .enumerate()
            .map(|(k, a)| a.cross(pts[(k + 1) % n]))
            .sum::<f64>()
            > 0.0;
        return if ccw { p } else { p.reversed() };
    }
    // Начало: закреплённая вершина, иначе самый острый угол.
    let start = (0..n).find(|&k| pinned[k]).unwrap_or_else(|| {
        (0..n)
            .max_by(|&i, &j| {
                let t = |k: usize| turn(pts[(k + n - 1) % n], pts[k], pts[(k + 1) % n]);
                t(i).total_cmp(&t(j))
            })
            .expect("вершины есть")
    });
    let q: Vec<Vec2> = (0..=n).map(|k| pts[(start + k) % n]).collect();
    let pin: Vec<bool> = (0..=n).map(|k| pinned[(start + k) % n]).collect();
    // Куски: (начало, конец, дуга или отрезок).
    let mut pieces: Vec<(usize, usize, Option<Arc>)> = Vec::new();
    let mut i = 0;
    while i < n {
        let mut line_end = i + 1;
        while line_end < n && !pin[line_end] && line_ok(&q[i..=line_end + 1], tol) {
            line_end += 1;
        }
        let mut arc = None;
        let mut j = i + 2;
        while j <= n {
            match arc_fit(&q[i..=j], tol) {
                Some(a) => arc = Some((j, a)),
                None => break,
            }
            if pin[j] {
                break;
            }
            j += 1;
        }
        // Закреплённая вершина внутри дуги недопустима: дуга кончается не дальше неё.
        let arc = arc.filter(|&(j, _)| !(i + 1..j).any(|k| pin[k]));
        match arc {
            Some((j, a)) if j > line_end => {
                pieces.push((i, j, Some(a)));
                i = j;
            }
            _ => {
                pieces.push((i, line_end, None));
                i = line_end;
            }
        }
    }
    // Шов: если начало не закреплено и первый и последний кусок сливаются — слить.
    if pieces.len() >= 2 && !pin[0] {
        let (ls, _, la) = pieces[pieces.len() - 1];
        let (_, fe, fa) = pieces[0];
        let seam: Vec<Vec2> = q[ls..n].iter().chain(&q[..=fe]).copied().collect();
        let merged = match (la, fa) {
            (None, None) => line_ok(&seam, tol).then_some(None),
            (Some(_), Some(_)) => arc_fit(&seam, tol).map(Some),
            _ => None,
        };
        if let Some(m) = merged {
            // Новый первый кусок начинается с ls.
            pieces.pop();
            pieces[0] = (ls, fe, m);
        }
    }
    let mut path = Path::new();
    for &(s, _, a) in &pieces {
        path.push(q[s % n], a.map_or(0.0, |a| a.bulge()));
    }
    // Деталь тоньше допуска могла схлопнуться (два отрезка туда-обратно) или
    // вывернуться. Тогда — ломаная как есть: в файле не должно быть пустых контуров.
    let area = |p: &[Vec2]| -> f64 {
        (0..p.len())
            .map(|k| p[k].cross(p[(k + 1) % p.len()]))
            .sum::<f64>()
            * 0.5
    };
    let (a0, a1) = (area(pts), path.area());
    if path.len() < 3 && pieces.iter().all(|p| p.2.is_none())
        || a1 * a0 <= 0.0
        || a1.abs() < 0.5 * a0.abs()
    {
        return Path::polygon(pts);
    }
    path
}

/// Подогнать один контур (все вершины можно убирать).
pub fn fit_ring(ring: &Ring, tol_mm: f64) -> Path {
    let pts = ring.to_mm();
    fit_points(&pts, &vec![false; pts.len()], tol_mm)
}

/// Подогнать все контуры области. Вершины, общие для нескольких контуров или
/// лежащие на чужом контуре, остаются на месте.
pub fn fit_region(region: &Region, tol_mm: f64) -> Vec<Path> {
    let mut count: BTreeMap<GridPoint, usize> = BTreeMap::new();
    for r in region.rings() {
        for &p in r.points() {
            *count.entry(p).or_default() += 1;
        }
    }
    let boxes: Vec<_> = region.rings().iter().map(Ring::bbox).collect();
    region
        .rings()
        .iter()
        .enumerate()
        .map(|(k, ring)| {
            let pinned: Vec<bool> = ring
                .points()
                .iter()
                .map(|&p| {
                    count[&p] > 1
                        || region.rings().iter().enumerate().any(|(m, other)| {
                            let (lo, hi) = boxes[m];
                            m != k
                                && (lo.x()..=hi.x()).contains(&p.x())
                                && (lo.y()..=hi.y()).contains(&p.y())
                                && other.on_boundary(p)
                        })
                })
                .collect();
            fit_points(&ring.to_mm(), &pinned, tol_mm)
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn square_stays_four_lines() {
        let r = Ring::rect(0, 0, 10_000, 10_000).unwrap();
        let p = fit_ring(&r, 0.005);
        assert_eq!(p.len(), 4);
        assert!(p.vertices().iter().all(|v| v.bulge == 0.0));
    }

    #[test]
    fn circle_comes_back_as_arcs() {
        let c = Path::circle(Vec2::new(3.0, -2.0), 7.5);
        let ring = c.to_ring(0.002).unwrap();
        assert!(ring.len() > 100);
        let p = fit_ring(&ring, 0.005);
        assert!(p.len() <= 3, "окружность — две-три дуги, а не {}", p.len());
        assert!((p.area() - c.area()).abs() < 7.5 * TAU * 0.005);
        for s in p.segments() {
            if let crate::path::Segment::Arc(a) = s {
                assert!((a.radius - 7.5).abs() < 0.005, "радиус {}", a.radius);
                assert!(a.center.distance(Vec2::new(3.0, -2.0)) < 0.005);
            } else {
                panic!("в окружности не должно быть отрезков");
            }
        }
    }

    #[test]
    fn rounded_rect_recovers_corners() {
        let rr = Path::rounded_rect(0.0, 0.0, 40.0, 20.0, 3.0);
        let p = fit_ring(&rr.to_ring(0.002).unwrap(), 0.005);
        let arcs = p.vertices().iter().filter(|v| v.bulge != 0.0).count();
        assert_eq!((p.len(), arcs), (8, 4), "{:?}", p.vertices());
    }
}
