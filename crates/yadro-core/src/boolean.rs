//! Булевы операции над областями: объединение, пересечение, разность, XOR.
//!
//! Метод — snap rounding Хобби (Hobby, 1999) на сетке 1 мкм:
//!
//! 1. «Горячие» пиксели — пиксели всех вершин и всех точек пересечения рёбер.
//! 2. Каждое ребро заменяется ломаной через центры горячих пикселей, которые оно
//!    задевает, в порядке вдоль ребра. По теореме Хобби такие ломаные не пересекаются,
//!    они могут только сходиться в вершинах или совпадать.
//! 3. Рёбра делятся в вершинах, совпадающие рёбра сливаются, у каждого остаётся вклад
//!    в обмотку каждого операнда.
//! 4. Для каждого ребра считается обмотка обоих операндов слева от него — лучом из точки,
//!    символически сдвинутой с середины ребра на бесконечно малое расстояние влево.
//!    Справа обмотка меньше на вклад ребра.
//! 5. Ребро — граница результата, если правило операции даёт разный ответ слева и справа.
//!    Из таких рёбер собираются контуры: внутренность всегда слева.
//!
//! Внутренность операнда — по правилу ненулевой обмотки, поэтому на вход годятся любые
//! замкнутые контуры, даже самопересекающиеся. Результат — корректная [`Region`].
//!
//! Точность: вершины результата — вершины входа или точки пересечения, округлённые
//! на сетку; рёбра смещаются не больше чем на √2/2 шага (0,71 мкм). Предикаты
//! считаются точно, допусков нет.
//!
//! Сложность пока квадратичная по числу рёбер: для деталей каркаса (сотни рёбер)
//! это миллисекунды. Заметка на будущее: заменить перебор пар заметающей прямой.

use crate::grid::GridPoint;
use crate::region::Region;
use crate::ring::Ring;
use crate::seg::{self, Crossing};
use std::cmp::Ordering;
use std::collections::{BTreeMap, BTreeSet};

/// Булева операция.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum Op {
    /// A ∪ B.
    Union,
    /// A ∩ B.
    Intersection,
    /// A ∖ B.
    Difference,
    /// A △ B (симметрическая разность).
    Xor,
}

impl Op {
    /// Правило: внутри ли результат, если точка внутри A — `a`, внутри B — `b`.
    pub fn apply(self, a: bool, b: bool) -> bool {
        match self {
            Op::Union => a || b,
            Op::Intersection => a && b,
            Op::Difference => a && !b,
            Op::Xor => a != b,
        }
    }
}

/// Направленный отрезок входа и номер операнда (0 — A, 1 — B).
#[derive(Clone, Copy, Debug)]
pub(crate) struct Input {
    pub a: GridPoint,
    pub b: GridPoint,
    pub owner: usize,
}

/// Ребро расположения: `u < v`; `w[k]` — вклад операнда k в обмотку при обходе u → v.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) struct Edge {
    pub u: GridPoint,
    pub v: GridPoint,
    pub w: [i64; 2],
}

/// Объединить по произвольному правилу: результат внутри там, где `rule(в A, в B)`.
pub fn combine(a: &Region, b: &Region, rule: impl Fn(bool, bool) -> bool) -> Region {
    let mut inputs = Vec::new();
    for (owner, reg) in [(0, a), (1, b)] {
        for r in reg.rings() {
            inputs.extend(r.edges().map(|(p, q)| Input { a: p, b: q, owner }));
        }
    }
    from_inputs(&inputs, |w: [i64; 2]| rule(w[0] != 0, w[1] != 0))
}

/// Булева операция.
pub fn boolean(a: &Region, b: &Region, op: Op) -> Region {
    combine(a, b, |x, y| op.apply(x, y))
}

/// Область по правилу ненулевой обмотки из произвольных контуров.
pub(crate) fn normalize(rings: &[Ring]) -> Region {
    let inputs: Vec<Input> = rings
        .iter()
        .flat_map(|r| {
            r.edges().map(|(p, q)| Input {
                a: p,
                b: q,
                owner: 0,
            })
        })
        .collect();
    from_inputs(&inputs, |w: [i64; 2]| w[0] != 0)
}

fn from_inputs(inputs: &[Input], inside: impl Fn([i64; 2]) -> bool) -> Region {
    if inputs.is_empty() {
        return Region::empty();
    }
    let edges = arrangement(inputs);
    let mut directed = Vec::new();
    for k in 0..edges.len() {
        let left = left_winding(&edges, k);
        let e = edges[k];
        let right = [left[0] - e.w[0], left[1] - e.w[1]];
        match (inside(left), inside(right)) {
            (true, false) => directed.push((e.u, e.v)),
            (false, true) => directed.push((e.v, e.u)),
            _ => {}
        }
    }
    Region::from_trusted(assemble(&directed))
}

// ───────────────────────────── Snap rounding ─────────────────────────────

/// Сколько раз можно добавить горячие пиксели, если после округления нашлись
/// пересечения. По теореме Хобби хватает одного раунда; остальные — страховка.
/// РИСК: если раундов не хватило, операция паникует — такой случай сразу в каталог.
const MAX_ROUNDS: usize = 8;

fn bbox(a: GridPoint, b: GridPoint) -> (i64, i64, i64, i64) {
    (
        a.x().min(b.x()),
        a.y().min(b.y()),
        a.x().max(b.x()),
        a.y().max(b.y()),
    )
}

fn boxes_meet(p: (i64, i64, i64, i64), q: (i64, i64, i64, i64)) -> bool {
    p.0 <= q.2 && q.0 <= p.2 && p.1 <= q.3 && q.1 <= p.3
}

/// Пары отрезков с пересекающимися габаритами (сортировка по левому краю).
fn candidate_pairs(segs: &[(GridPoint, GridPoint)], mut f: impl FnMut(usize, usize)) {
    let boxes: Vec<_> = segs.iter().map(|&(a, b)| bbox(a, b)).collect();
    let mut order: Vec<usize> = (0..segs.len()).collect();
    order.sort_by_key(|&i| boxes[i].0);
    for (k, &i) in order.iter().enumerate() {
        for &j in &order[k + 1..] {
            if boxes[j].0 > boxes[i].2 {
                break;
            }
            if boxes_meet(boxes[i], boxes[j]) {
                f(i.min(j), i.max(j));
            }
        }
    }
}

/// Пиксели настоящих пересечений.
fn crossing_pixels(segs: &[(GridPoint, GridPoint)]) -> Vec<GridPoint> {
    let mut out = Vec::new();
    candidate_pairs(segs, |i, j| {
        let ((a, b), (c, d)) = (segs[i], segs[j]);
        if seg::classify(a, b, c, d) == Crossing::Proper {
            out.push(
                seg::intersection(a, b, c, d)
                    .expect("настоящее пересечение не параллельно")
                    .pixel(),
            );
        }
    });
    out
}

pub(crate) fn arrangement(inputs: &[Input]) -> Vec<Edge> {
    let segs: Vec<_> = inputs.iter().map(|e| (e.a, e.b)).collect();
    let mut hot: BTreeSet<GridPoint> = inputs.iter().flat_map(|e| [e.a, e.b]).collect();
    hot.extend(crossing_pixels(&segs));
    for _ in 0..MAX_ROUNDS {
        let edges = snap(inputs, &hot);
        let pairs: Vec<_> = edges.iter().map(|e| (e.u, e.v)).collect();
        let extra = crossing_pixels(&pairs);
        if extra.is_empty() {
            return edges;
        }
        hot.extend(extra);
    }
    panic!("ядро: snap rounding не сошёлся за {MAX_ROUNDS} раундов — это ошибка ядра, нужен тест");
}

/// Горячие пиксели, отсортированные по x, для выборки по габариту.
struct Hot {
    pts: Vec<GridPoint>,
}

impl Hot {
    fn new(set: &BTreeSet<GridPoint>) -> Self {
        // GridPoint упорядочены по (x, y) — BTreeSet уже отсортирован.
        Self {
            pts: set.iter().copied().collect(),
        }
    }

    /// Пиксели, центры которых в габарите отрезка, расширенном на 1 шаг.
    fn near(&self, a: GridPoint, b: GridPoint) -> impl Iterator<Item = GridPoint> + '_ {
        let (lx, ly, hx, hy) = bbox(a, b);
        let start = self.pts.partition_point(|p| p.x() < lx - 1);
        self.pts[start..]
            .iter()
            .take_while(move |p| p.x() <= hx + 1)
            .filter(move |p| p.y() >= ly - 1 && p.y() <= hy + 1)
            .copied()
    }
}

fn snap(inputs: &[Input], hot_set: &BTreeSet<GridPoint>) -> Vec<Edge> {
    let hot = Hot::new(hot_set);
    let mut merged: BTreeMap<(GridPoint, GridPoint), [i64; 2]> = BTreeMap::new();
    let mut add = |p: GridPoint, q: GridPoint, owner: usize| {
        if p == q {
            return;
        }
        let (key, s) = if p < q { ((p, q), 1) } else { ((q, p), -1) };
        merged.entry(key).or_insert([0, 0])[owner] += s;
    };
    for e in inputs {
        if e.a == e.b {
            continue;
        }
        // Горячие пиксели вдоль ребра в порядке входа.
        let mut hits: Vec<((i128, i128, bool), GridPoint)> = hot
            .near(e.a, e.b)
            .filter_map(|c| seg::hits_pixel(e.a, e.b, c).map(|t| (t, c)))
            .collect();
        hits.sort_by(|x, y| seg::entry_cmp(x.0, y.0));
        let mut path: Vec<GridPoint> = hits.into_iter().map(|h| h.1).collect();
        path.dedup();
        debug_assert_eq!(path.first(), Some(&e.a));
        debug_assert_eq!(path.last(), Some(&e.b));
        for w in path.windows(2) {
            // Делим звено во всех горячих центрах, лежащих внутри него: геометрия не меняется,
            // зато в расположении нет вершин посреди рёбер.
            let (p, q) = (w[0], w[1]);
            let mut mids: Vec<GridPoint> = hot
                .near(p, q)
                .filter(|&c| seg::in_segment_interior(c, p, q))
                .collect();
            mids.sort_by_key(|c| {
                let (dx, dy) = (i128::from(c.x() - p.x()), i128::from(c.y() - p.y()));
                dx * dx + dy * dy
            });
            let mut prev = p;
            for m in mids {
                add(prev, m, e.owner);
                prev = m;
            }
            add(prev, q, e.owner);
        }
    }
    merged
        .into_iter()
        .filter(|(_, w)| *w != [0, 0])
        .map(|((u, v), w)| Edge { u, v, w })
        .collect()
}

// ───────────────────────────── Обмотка ─────────────────────────────

/// Обмотка обоих операндов в точке, бесконечно близкой к середине ребра `k` слева от него.
///
/// Точка P = M + ε·n, где M — середина ребра, n — нормаль влево, ε → +0. Сравнения с P
/// раскрываются символически: сначала по M, при равенстве — по знаку члена при ε.
/// Это корректно, потому что ни одна вершина не лежит внутри ребра и рёбра не пересекаются:
/// при малом ε точка P строго внутри грани слева от ребра.
///
/// РИСК: символический сдвиг держится на том, что `u < v` (тогда n.y ≥ 0) и что в
/// расположении нет вершин внутри рёбер. Сложность O(n) на ребро, O(n²) на операцию —
/// для контуров из тысяч рёбер (TopOpt) понадобится заметающая прямая.
pub(crate) fn left_winding(edges: &[Edge], k: usize) -> [i64; 2] {
    let e = edges[k];
    // Удвоенные координаты: середина целая.
    let mx = i128::from(e.u.x()) + i128::from(e.v.x());
    let my = i128::from(e.u.y()) + i128::from(e.v.y());
    let (nx, ny) = (
        -(i128::from(e.v.y() - e.u.y())),
        i128::from(e.v.x() - e.u.x()),
    );
    // y ≤ P.y ?
    let le = |y: i128| y < my || (y == my && ny >= 0);
    let mut w = [0i64; 2];
    for f in edges {
        let (ax, ay) = (2 * i128::from(f.u.x()), 2 * i128::from(f.u.y()));
        let (bx, by) = (2 * i128::from(f.v.x()), 2 * i128::from(f.v.y()));
        let (lo, hi) = (le(ay), le(by));
        if lo == hi {
            continue;
        }
        // Знак ориентации (A, B, P) с символическим ε.
        let o = (bx - ax) * (my - ay) - (by - ay) * (mx - ax);
        let s = if o != 0 {
            o.signum()
        } else {
            ((bx - ax) * ny - (by - ay) * nx).signum()
        };
        if lo && !hi && s > 0 {
            // Ребро идёт вверх, P слева от него.
            w[0] += f.w[0];
            w[1] += f.w[1];
        } else if !lo && hi && s < 0 {
            // Ребро идёт вниз, P справа от него.
            w[0] -= f.w[0];
            w[1] -= f.w[1];
        }
    }
    w
}

// ───────────────────────────── Сборка контуров ─────────────────────────────

fn sub(a: GridPoint, b: GridPoint) -> (i128, i128) {
    (
        i128::from(a.x()) - i128::from(b.x()),
        i128::from(a.y()) - i128::from(b.y()),
    )
}

fn cross(a: (i128, i128), b: (i128, i128)) -> i128 {
    a.0 * b.1 - a.1 * b.0
}

/// Порядок направлений по углу, отсчитанному по часовой стрелке от `r`, в (0, 2π].
fn cw_from(r: (i128, i128), a: (i128, i128), b: (i128, i128)) -> Ordering {
    let half = |v: (i128, i128)| {
        let c = cross(r, v);
        let d = r.0 * v.0 + r.1 * v.1;
        if c < 0 {
            0
        } else if c == 0 && d < 0 {
            1
        } else if c > 0 {
            2
        } else {
            3
        }
    };
    // В одной половине раньше то, от которого второе повёрнуто по часовой (векторное произведение < 0).
    half(a).cmp(&half(b)).then_with(|| cross(a, b).cmp(&0))
}

/// Собрать контуры из направленных рёбер границы (внутренность слева).
///
/// После ребра u → v берётся то выходящее из v ребро, которое первым встречается
/// при повороте по часовой стрелке от направления назад (v → u). Для границы области
/// входящие и выходящие рёбра вокруг вершины чередуются, поэтому правило задаёт
/// взаимно однозначное «следующее ребро», и касающиеся в вершине контуры остаются раздельными.
fn assemble(directed: &[(GridPoint, GridPoint)]) -> Vec<Ring> {
    let mut out: BTreeMap<GridPoint, Vec<usize>> = BTreeMap::new();
    for (i, &(a, _)) in directed.iter().enumerate() {
        out.entry(a).or_default().push(i);
    }
    let next: Vec<usize> = directed
        .iter()
        .map(|&(a, b)| {
            let back = sub(a, b);
            *out[&b]
                .iter()
                .min_by(|&&x, &&y| cw_from(back, sub(directed[x].1, b), sub(directed[y].1, b)))
                .expect("у вершины границы есть выходящее ребро")
        })
        .collect();
    let mut used = vec![false; directed.len()];
    let mut rings = Vec::new();
    for s in 0..directed.len() {
        if used[s] {
            continue;
        }
        let mut pts = Vec::new();
        let mut cur = s;
        while !used[cur] {
            used[cur] = true;
            pts.push(directed[cur].0);
            cur = next[cur];
        }
        debug_assert_eq!(cur, s, "следующее ребро должно замыкать контур");
        if let Ok(r) = Ring::new(pts) {
            rings.push(r);
        }
    }
    rings
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn op_rules() {
        assert!(Op::Union.apply(true, false));
        assert!(!Op::Intersection.apply(true, false));
        assert!(Op::Difference.apply(true, false) && !Op::Difference.apply(true, true));
        assert!(Op::Xor.apply(false, true) && !Op::Xor.apply(true, true));
    }

    #[test]
    fn clockwise_order() {
        // От направления «вниз» по часовой: запад, север, восток.
        let down = (0, -1);
        let (w, n, e) = ((-1, 0), (0, 1), (1, 0));
        assert_eq!(cw_from(down, w, n), Ordering::Less);
        assert_eq!(cw_from(down, n, e), Ordering::Less);
        assert_eq!(cw_from(down, w, e), Ordering::Less);
    }
}
