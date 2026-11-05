//! Ворота ступени 1, шаг 2: булевы операции.
//!
//! 1. Каталог вырожденных случаев — каждый с известным ответом.
//! 2. Случайные пары областей (полный прогон — 100 000): результат корректен
//!    ([`Region::check`]), выполняются точные тождества площадей, а случайные точки
//!    вдали от границ классифицируются так же, как по исходным областям.
//!
//! Тождества точные: все четыре операции и «A, округлённая в паре с B» строятся
//! на одном и том же расположении рёбер, поэтому площади складываются без погрешности.

mod common;

use common::*;
use yadro_core::boolean::combine;
use yadro_core::svg::{Style, Svg};
use yadro_core::testkit::{Rng, for_cases};
use yadro_core::tol::GRID_LIMIT;
use yadro_core::{Location, Op, Region, Ring};

/// Карточка ошибки: картинка A, B и результата в папке временных файлов тестов.
fn card(name: &str, a: &Region, b: &Region, r: &Region) -> String {
    let mut s = Svg::new();
    s.region(a, Style::filled("#1f5fbf", 0.25), Style::line("#1f5fbf"));
    s.region(b, Style::filled("#c0392b", 0.25), Style::line("#c0392b"));
    s.region(
        r,
        Style::line("#202020").width(2.0),
        Style::line("#27ae60").width(2.0),
    );
    let path = format!("{}/{name}.svg", env!("CARGO_TARGET_TMPDIR"));
    let _ = std::fs::write(&path, s.render(900.0));
    path
}

/// Все проверки пары: корректность, тождества. Возвращает [U, I, D, X].
fn check_pair(name: &str, a: &Region, b: &Region) -> [Region; 4] {
    let ops = [Op::Union, Op::Intersection, Op::Difference, Op::Xor];
    let res = ops.map(|op| yadro_core::boolean::boolean(a, b, op));
    let fail = |what: String, r: &Region| -> ! {
        let path = card(name, a, b, r);
        panic!(
            "{what}\nA = {}\nB = {}\nкарточка: {path}",
            rings_text(a),
            rings_text(b)
        )
    };
    for (op, r) in ops.iter().zip(&res) {
        if let Err(e) = r.check() {
            fail(format!("{op:?}: результат некорректен: {e}"), r);
        }
    }
    let a1 = combine(a, b, |x, _| x);
    let b1 = combine(a, b, |_, y| y);
    let ba = combine(a, b, |x, y| y && !x);
    let [u, i, d, x] = res.each_ref().map(Region::twice_area);
    let (sa, sb, sba) = (a1.twice_area(), b1.twice_area(), ba.twice_area());
    if u + i != sa + sb || d + i != sa || x != u - i || u != d + sba + i {
        fail(
            format!(
                "тождества площадей нарушены: U={u} I={i} D={d} X={x} A'={sa} B'={sb} B∖A={sba}"
            ),
            &res[0],
        );
    }
    // Округление не уводит площадь операнда дальше, чем на периметр × 0,71 шага.
    for (orig, snapped) in [(a, sa), (b, sb)] {
        let per: f64 = orig
            .edges()
            .map(|(q, r)| ((q.x() - r.x()) as f64).hypot((q.y() - r.y()) as f64))
            .sum();
        let diff = (orig.twice_area() - snapped).unsigned_abs() as f64 / 2.0;
        if diff > per * 0.71 + 1.0 {
            fail(
                format!("округление сдвинуло площадь на {diff} при периметре {per}"),
                &res[0],
            );
        }
    }
    res
}

// ───────────────────────────── Каталог вырожденных случаев ─────────────────────────────

#[test]
fn catalog_identical() {
    let a = rect(0, 0, 10, 10);
    let [u, i, d, x] = check_pair("catalog_identical", &a, &a);
    assert_eq!((u.twice_area(), i.twice_area()), (200, 200));
    assert!(d.is_empty() && x.is_empty());
}

#[test]
fn catalog_shared_edge() {
    let [u, i, d, x] = check_pair(
        "catalog_shared_edge",
        &rect(0, 0, 10, 10),
        &rect(10, 0, 20, 10),
    );
    assert_eq!(u.twice_area(), 400);
    assert_eq!(u.rings().len(), 1);
    assert_eq!(
        u.rings()[0].len(),
        4,
        "общая сторона исчезает, остаётся прямоугольник"
    );
    assert!(i.is_empty());
    assert_eq!((d.twice_area(), x.twice_area()), (200, 400));
}

#[test]
fn catalog_shared_vertex() {
    let [u, i, ..] = check_pair(
        "catalog_shared_vertex",
        &rect(0, 0, 10, 10),
        &rect(10, 10, 20, 20),
    );
    assert_eq!(u.twice_area(), 400);
    assert_eq!(u.rings().len(), 2, "касание в вершине — два контура");
    assert!(i.is_empty());
}

#[test]
fn catalog_t_junction() {
    let tri = poly(&[(5, 10), (7, 14), (3, 14)]);
    let [u, i, d, _] = check_pair("catalog_t_junction", &rect(0, 0, 10, 10), &tri);
    assert_eq!(u.twice_area(), 200 + 16);
    assert!(i.is_empty());
    assert_eq!(d.twice_area(), 200);
}

#[test]
fn catalog_collinear_overlap() {
    let [u, i, d, x] = check_pair(
        "catalog_collinear_overlap",
        &rect(0, 0, 10, 10),
        &rect(5, 0, 15, 10),
    );
    assert_eq!(
        (
            u.twice_area(),
            i.twice_area(),
            d.twice_area(),
            x.twice_area()
        ),
        (300, 100, 100, 200)
    );
    assert_eq!(u.rings()[0].len(), 4);
    // Частичное наложение сторон без пересечения внутренностей.
    let [u2, i2, ..] = check_pair(
        "catalog_partial_edge",
        &rect(0, 0, 10, 10),
        &rect(5, 10, 15, 20),
    );
    assert_eq!(u2.twice_area(), 400);
    assert_eq!(u2.rings().len(), 1);
    assert!(i2.is_empty());
}

#[test]
fn catalog_nested() {
    let [u, i, d, x] = check_pair("catalog_nested", &rect(0, 0, 10, 10), &rect(3, 3, 6, 6));
    assert_eq!(
        (
            u.twice_area(),
            i.twice_area(),
            d.twice_area(),
            x.twice_area()
        ),
        (200, 18, 182, 182)
    );
    assert_eq!(d.rings().len(), 2);
    assert_eq!(
        d.rings().iter().filter(|r| !r.is_ccw()).count(),
        1,
        "одно отверстие по часовой"
    );
    assert_eq!(d.locate(p(4, 4)), Location::Outside);
}

#[test]
fn catalog_nested_touching_edge() {
    let [_, i, d, _] = check_pair(
        "catalog_nested_touching_edge",
        &rect(0, 0, 10, 10),
        &rect(0, 3, 4, 6),
    );
    assert_eq!(i.twice_area(), 24);
    assert_eq!(d.twice_area(), 176);
    assert_eq!(d.rings().len(), 1, "вырез у края — один контур с выемкой");
}

#[test]
fn catalog_hole_touching_outer_at_vertex() {
    let tri = poly(&[(5, 0), (7, 4), (3, 4)]);
    let [_, _, d, _] = check_pair(
        "catalog_hole_touching_outer_at_vertex",
        &rect(0, 0, 10, 10),
        &tri,
    );
    assert_eq!(d.twice_area(), 200 - 16);
    d.check().unwrap();
}

#[test]
fn catalog_plus() {
    let [u, i, ..] = check_pair("catalog_plus", &rect(0, 4, 12, 8), &rect(4, 0, 8, 12));
    assert_eq!((u.twice_area(), i.twice_area()), (160, 32));
    assert_eq!(u.rings().len(), 1);
    assert_eq!(u.rings()[0].len(), 12);
}

#[test]
fn catalog_shallow_crossings_need_rounding() {
    // Почти параллельные рёбра: точки пересечения не на сетке и рядом друг с другом.
    let a = poly(&[(0, 0), (1000, 1), (1000, 3), (0, 2)]);
    let b = poly(&[(0, 1), (1000, 0), (1000, 2), (0, 3)]);
    check_pair("catalog_shallow_crossings_need_rounding", &a, &b);
    let fan: Vec<Region> = (0..6)
        .map(|k| poly(&[(0, k), (997, 7 - k), (997, 8 - k), (0, k + 1)]))
        .collect();
    let mut acc = Region::empty();
    for f in &fan {
        check_pair("catalog_fan", &acc, f);
        acc = acc.union(f);
    }
}

#[test]
fn catalog_thin_sliver() {
    let sliver = poly(&[(0, 0), (1000, 0), (1000, 1)]);
    let [_, i, ..] = check_pair("catalog_thin_sliver", &sliver, &rect(500, -5, 501, 5));
    assert!(i.twice_area() <= 2);
}

#[test]
fn catalog_huge_coordinates() {
    let h = GRID_LIMIT / 2;
    let square = rect(-h, -h, h, h);
    let diamond = poly(&[
        (0, -GRID_LIMIT),
        (GRID_LIMIT, 0),
        (0, GRID_LIMIT),
        (-GRID_LIMIT, 0),
    ]);
    let [u, i, ..] = check_pair("catalog_huge_coordinates", &square, &diamond);
    assert_eq!(i.twice_area(), square.twice_area());
    assert_eq!(u.twice_area(), diamond.twice_area());
    let skew = poly(&[
        (-GRID_LIMIT, -GRID_LIMIT + 7),
        (GRID_LIMIT, -GRID_LIMIT),
        (GRID_LIMIT - 3, GRID_LIMIT),
        (-GRID_LIMIT + 11, GRID_LIMIT - 5),
    ]);
    check_pair("catalog_huge_skew", &skew, &diamond);
}

#[test]
fn catalog_bowtie_nonzero() {
    // Симметричный «бантик»: ориентированная площадь 0, но по ненулевой обмотке — два треугольника.
    let r = poly(&[(0, 0), (4, 4), (4, 0), (0, 4)]);
    r.check().unwrap();
    assert_eq!(r.rings().len(), 2);
    assert_eq!(r.twice_area(), 16);
}

#[test]
fn catalog_pentagram_nonzero() {
    let star = poly(&[(0, 100), (59, -81), (-95, 31), (95, 31), (-59, -81)]);
    star.check().unwrap();
    assert_eq!(
        star.locate(p(0, 0)),
        Location::Inside,
        "в центре обмотка 2 — внутри"
    );
    assert_eq!(star.rings().len(), 1);
    assert_eq!(star.rings()[0].len(), 10);
}

#[test]
fn catalog_empty_operands() {
    let a = rect(0, 0, 10, 10);
    let e = Region::empty();
    let [u, i, d, x] = check_pair("catalog_empty_operands", &a, &e);
    assert_eq!(
        (
            u.twice_area(),
            i.twice_area(),
            d.twice_area(),
            x.twice_area()
        ),
        (200, 0, 200, 200)
    );
    assert!(e.difference(&a).is_empty());
    assert!(e.union(&e).is_empty());
}

#[test]
fn catalog_opposite_duplicates_cancel() {
    let r = Ring::rect(0, 0, 10, 10).unwrap();
    assert!(Region::from_rings_nonzero(&[r.clone(), r.reversed()]).is_empty());
    // А одинаково ориентированные — обмотка 2, это та же область.
    assert_eq!(
        Region::from_rings_nonzero(&[r.clone(), r]).twice_area(),
        200
    );
}

#[test]
fn catalog_unit_shift() {
    let [u, i, ..] = check_pair(
        "catalog_unit_shift",
        &rect(0, 0, 1000, 1000),
        &rect(1, 1, 1001, 1001),
    );
    assert_eq!(i.twice_area(), 2 * 999 * 999);
    assert_eq!(u.twice_area(), 2 * (2 * 1_000_000 - 999 * 999));
}

// ───────────────────────────── Случайные пары ─────────────────────────────

#[test]
fn random_pairs_agree() {
    for_cases("random_pairs_agree", 2_000, |rng| {
        let scale = *rng.pick(&[8, 40, 1_000, 1_000_000, GRID_LIMIT / 4]);
        let a = random_region(rng, scale);
        let mut b = random_region(rng, scale);
        // Иногда B — сдвинутая на шаг копия A: почти совпадающие рёбра.
        if rng.chance(0.15) {
            let (sx, sy) = (rng.range_i64(-1, 1), rng.range_i64(-1, 1));
            let rings: Vec<Ring> = a
                .rings()
                .iter()
                .map(|r| {
                    Ring::new(
                        r.points()
                            .iter()
                            .map(|q| p(q.x() + sx, q.y() + sy))
                            .collect(),
                    )
                    .unwrap()
                })
                .collect();
            b = Region::from_rings_nonzero(&rings);
        }
        let res = check_pair("random_pairs_agree", &a, &b);
        // Независимая сверка случайными точками вдали от границ.
        let edges: Vec<_> = a.edges().chain(b.edges()).collect();
        if edges.is_empty() {
            return;
        }
        let lim = (scale + 2).min(GRID_LIMIT);
        for _ in 0..24 {
            let q = if rng.chance(0.3) {
                // Рядом с вершиной.
                let v = edges[rng.below(edges.len() as u64) as usize].0;
                p(
                    (v.x() + rng.range_i64(-3, 3)).clamp(-GRID_LIMIT, GRID_LIMIT),
                    (v.y() + rng.range_i64(-3, 3)).clamp(-GRID_LIMIT, GRID_LIMIT),
                )
            } else {
                p(rng.range_i64(-lim, lim), rng.range_i64(-lim, lim))
            };
            if edges.iter().any(|&(u, v)| dist(q, u, v) < 3.0) {
                continue;
            }
            let ina = a.locate(q) == Location::Inside;
            let inb = b.locate(q) == Location::Inside;
            for (op, r) in [Op::Union, Op::Intersection, Op::Difference, Op::Xor]
                .iter()
                .zip(&res)
            {
                let got = r.locate(q);
                let want = if op.apply(ina, inb) {
                    Location::Inside
                } else {
                    Location::Outside
                };
                if got != want {
                    let path = card("random_pairs_agree", &a, &b, r);
                    panic!(
                        "{op:?}: точка {q:?} — {got:?}, ожидалось {want:?}\nA = {}\nB = {}\nкарточка: {path}",
                        rings_text(&a),
                        rings_text(&b)
                    );
                }
            }
        }
    });
}

#[test]
fn associativity_of_union_on_lattice() {
    // (A ∪ B) ∪ C и A ∪ (B ∪ C) на решётке без округлений совпадают по площади точно.
    for_cases("associativity_of_union_on_lattice", 500, |rng| {
        let r = |rng: &mut Rng| {
            let (x0, y0) = (rng.range_i64(-6, 5), rng.range_i64(-6, 5));
            rect(
                x0 * 10,
                y0 * 10,
                (x0 + rng.range_i64(1, 5)) * 10,
                (y0 + rng.range_i64(1, 5)) * 10,
            )
        };
        let (a, b, c) = (r(rng), r(rng), r(rng));
        let l = a.union(&b).union(&c);
        let rr = a.union(&b.union(&c));
        l.check().unwrap();
        rr.check().unwrap();
        assert_eq!(l.twice_area(), rr.twice_area());
    });
}
