//! Эквидистанта (offset) области на расстояние d: d > 0 — наружу, d < 0 — внутрь.
//!
//! Наружу результат — сумма Минковского A ⊕ круг радиуса d: сама A, полосы ширины d
//! снаружи (справа) от каждого ребра и сектора радиуса d у выпуклых вершин. Точка не
//! дальше d от A либо лежит в A, либо её ближайшая точка на границе — внутри ребра
//! (тогда точка в полосе), либо в выпуклой вершине (тогда в секторе). Внутрь — то же
//! для дополнения: из A вычитаются полосы слева от рёбер и сектора у вогнутых вершин.
//!
//! Все кусочки собираются одной булевой операцией, поэтому результат корректен, как
//! любой результат булевых операций, и само собой разбирается со схлопнувшимися
//! отверстиями, слившимися деталями и перемычками, порвавшимися при сужении.
//!
//! Точность. Дуги секторов — ломаные с вершинами на окружности и прогибом не больше
//! допуска хорды `tol`; углы полос — точки на расстоянии d, округлённые на сетку.
//! Булева операция сдвигает рёбра ещё не больше чем на 0,71 мкм. Поэтому граница
//! результата отстоит от границы A на расстояние от |d| − tol − 1,42 мкм до
//! |d| + 1,42 мкм (√2 шага сетки): на дугах хорды чуть срезают круг, так что наружу
//! результат может быть чуть меньше точного, а внутрь — чуть больше.

use crate::arc::chord_count;
use crate::boolean;
use crate::grid::{GridError, GridPoint};
use crate::predicates::{Sign, orient2d_grid};
use crate::region::Region;
use crate::tol::{GRID_LIMIT, GRID_PER_MM, MIN_CHORD_TOL_MM};

/// Единичная нормаль к ребру p → q: справа (`right`) или слева.
fn normal(p: GridPoint, q: GridPoint, right: bool) -> (f64, f64) {
    let (dx, dy) = ((q.x() - p.x()) as f64, (q.y() - p.y()) as f64);
    let len = dx.hypot(dy);
    if right {
        (dy / len, -dx / len)
    } else {
        (-dy / len, dx / len)
    }
}

/// Точка сетки, ближайшая к p + e·n.
fn shift(p: GridPoint, n: (f64, f64), e: f64) -> Result<GridPoint, GridError> {
    let round = |v: f64| {
        let r = v.round();
        if r.abs() > GRID_LIMIT as f64 {
            Err(GridError::OutOfRange)
        } else {
            Ok(r as i64)
        }
    };
    GridPoint::new(
        round(p.x() as f64 + e * n.0)?,
        round(p.y() as f64 + e * n.1)?,
    )
}

/// Сектор у вершины `v` радиуса `e` от направления `a` до `b` против часовой стрелки
/// (раствор меньше π). Крайние точки считаются той же функцией, что и углы полос,
/// поэтому совпадают с ними точно.
///
/// РИСК: если считать крайние точки иначе (например, через угол), они разойдутся с
/// углами полос на доли микрона, и между кусочками останутся щели.
fn sector(
    v: GridPoint,
    a: (f64, f64),
    b: (f64, f64),
    e: f64,
    tol: f64,
) -> Result<Vec<GridPoint>, GridError> {
    let sweep = (a.0 * b.1 - a.1 * b.0)
        .atan2(a.0 * b.0 + a.1 * b.1)
        .max(0.0);
    let n = chord_count(e, sweep, tol);
    let start = a.1.atan2(a.0);
    let mut pts = vec![v, shift(v, a, e)?];
    for k in 1..n {
        let t = start + sweep * k as f64 / n as f64;
        pts.push(shift(v, (t.cos(), t.sin()), e)?);
    }
    pts.push(shift(v, b, e)?);
    Ok(pts)
}

/// Эквидистанта области на `d_mm` (> 0 — наружу, < 0 — внутрь) с допуском хорды
/// `tol_mm` для дуг. Ошибка — только если результат выходит за предел сетки или
/// расстояние не число.
pub fn offset(region: &Region, d_mm: f64, tol_mm: f64) -> Result<Region, GridError> {
    if !d_mm.is_finite() || !tol_mm.is_finite() {
        return Err(GridError::NotFinite);
    }
    let d = d_mm * GRID_PER_MM;
    if d == 0.0 || region.is_empty() {
        return Ok(region.clone());
    }
    let outward = d > 0.0;
    let e = d.abs();
    let tol = tol_mm.max(MIN_CHORD_TOL_MM) * GRID_PER_MM;
    let mut loops: Vec<Vec<GridPoint>> = Vec::new();
    for ring in region.rings() {
        let pts = ring.points();
        let n = pts.len();
        // Полосы вдоль рёбер; все кусочки против часовой стрелки.
        for i in 0..n {
            let (p, q) = (pts[i], pts[(i + 1) % n]);
            let s = normal(p, q, outward);
            let (p1, q1) = (shift(p, s, e)?, shift(q, s, e)?);
            loops.push(if outward {
                vec![p, p1, q1, q]
            } else {
                vec![p, q, q1, p1]
            });
        }
        // Сектора: наружу — у выпуклых вершин (поворот влево), внутрь — у вогнутых.
        for i in 0..n {
            let (prev, v, next) = (pts[(i + n - 1) % n], pts[i], pts[(i + 1) % n]);
            match (outward, orient2d_grid(prev, v, next)) {
                (true, Sign::Positive) => {
                    loops.push(sector(
                        v,
                        normal(prev, v, true),
                        normal(v, next, true),
                        e,
                        tol,
                    )?);
                }
                (false, Sign::Negative) => {
                    loops.push(sector(
                        v,
                        normal(v, next, false),
                        normal(prev, v, false),
                        e,
                        tol,
                    )?);
                }
                _ => {}
            }
        }
    }
    Ok(if outward {
        boolean::combine_loops(region, &loops, |a, b| a || b)
    } else {
        boolean::combine_loops(region, &loops, |a, b| a && !b)
    })
}
