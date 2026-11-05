//! Ворота ступени 1, шаг 3: дуги и эквидистанта (offset).
//!
//! Независимая проверка — расстояния. Знаковое расстояние s(q) до исходной области
//! (снаружи > 0, внутри < 0) считается прямо по рёбрам. У эквидистанты на d:
//! - каждая вершина результата лежит там, где s ≈ d (с допуском хорды и округлением);
//! - точки с s заметно меньше d — внутри результата, заметно больше — снаружи;
//! - наружу результат содержит исходную область, внутрь — содержится в ней;
//! - при большем d результат не меньше (вложенность).
//!
//! Плюс каталог случаев с известным ответом: точные квадраты, схлопнувшиеся
//! отверстия, слившиеся детали, порванные перемычки, скругления фрезой.

mod common;

use common::*;
use std::f64::consts::PI;
use yadro_core::card::shape_of;
use yadro_core::testkit::{Rng, for_cases};
use yadro_core::tol::{GRID_LIMIT, MIN_CHORD_TOL_MM};
use yadro_core::{GridPoint, Location, Path, Region, Vec2};

/// Запас на округление, мкм: углы кусочков округляются на сетку (до 0,71 мкм), и
/// snap rounding сдвигает границу ещё до 0,71 мкм — всего √2 ≈ 1,42, плюс погрешность f64.
const SLACK_VERTEX: f64 = 1.5;
/// Запас для точек-проб: граница результата и проба округляются независимо, мкм.
const SLACK_PROBE: f64 = 2.0;

/// Знаковое расстояние до области: снаружи > 0, внутри < 0, на границе 0.
fn signed_dist(q: GridPoint, r: &Region) -> f64 {
    match r.locate(q) {
        Location::Boundary => 0.0,
        Location::Inside => -dist_to_boundary(q, r),
        Location::Outside => dist_to_boundary(q, r),
    }
}

/// Все проверки эквидистанты на `d` мкм с допуском хорды `tol` мкм. Возвращает результат.
fn check_offset(name: &str, r: &Region, d: f64, tol: f64, rng: Option<&mut Rng>) -> Region {
    let o = r
        .offset_tol(d / 1000.0, tol / 1000.0)
        .expect("в пределах сетки");
    let tol = tol.max(MIN_CHORD_TOL_MM * 1000.0);
    let fail = |what: String| -> ! {
        panic!(
            "{name}: {what}\nd = {d} мкм, допуск {tol} мкм\nA = {}",
            rings_text(r)
        )
    };
    if let Err(e) = o.check() {
        fail(format!("результат некорректен: {e}"));
    }
    if d == 0.0 {
        if o != *r {
            fail("при d = 0 область должна остаться той же".into());
        }
        return o;
    }
    // Где может лежать граница результата по знаковому расстоянию.
    let (lo, hi) = if d > 0.0 { (d - tol, d) } else { (d, d + tol) };
    for q in o
        .rings()
        .iter()
        .flat_map(|ring| ring.points().iter().copied())
    {
        let s = signed_dist(q, r);
        if s < lo - SLACK_VERTEX || s > hi + SLACK_VERTEX {
            fail(format!(
                "вершина {q:?}: расстояние {s}, нужно от {lo} до {hi}"
            ));
        }
    }
    // Вложенность в исходную область.
    if d > tol + 4.0 && !r.difference(&o).is_empty() {
        fail("наружу: исходная область должна целиком войти в результат".into());
    }
    if d < -4.0 && !o.difference(r).is_empty() {
        fail("внутрь: результат должен целиком лежать в исходной области".into());
    }
    if let Some(rng) = rng {
        // Пробы: случайные точки в габарите с запасом и точки возле вершин.
        let verts: Vec<GridPoint> = r
            .rings()
            .iter()
            .chain(o.rings())
            .flat_map(|ring| ring.points().iter().copied())
            .collect();
        if verts.is_empty() {
            return o;
        }
        let (mut lx, mut ly, mut hx, mut hy) = (i64::MAX, i64::MAX, i64::MIN, i64::MIN);
        for v in &verts {
            lx = lx.min(v.x());
            ly = ly.min(v.y());
            hx = hx.max(v.x());
            hy = hy.max(v.y());
        }
        let m = d.abs() as i64 + 10;
        let clamp = |v: i64| v.clamp(-GRID_LIMIT, GRID_LIMIT);
        for _ in 0..32 {
            let q = if rng.chance(0.4) {
                let v = verts[rng.below(verts.len() as u64) as usize];
                let j = (d.abs() as i64 / 2).max(4);
                p(
                    clamp(v.x() + rng.range_i64(-j, j)),
                    clamp(v.y() + rng.range_i64(-j, j)),
                )
            } else {
                p(
                    clamp(rng.range_i64(lx - m, hx + m)),
                    clamp(rng.range_i64(ly - m, hy + m)),
                )
            };
            let s = signed_dist(q, r);
            let want = if s < lo - SLACK_PROBE {
                Location::Inside
            } else if s > hi + SLACK_PROBE {
                Location::Outside
            } else {
                continue;
            };
            let got = o.locate(q);
            if got != want {
                fail(format!(
                    "проба {q:?}: расстояние {s} — {got:?}, ожидалось {want:?}"
                ));
            }
        }
    }
    o
}

fn mm(r: &Region) -> f64 {
    r.area_mm2()
}

// ───────────────────────────── Каталог ─────────────────────────────

#[test]
fn zero_offset_is_identity() {
    let r = poly(&[
        (0, 0),
        (5000, 0),
        (5000, 3000),
        (2000, 3000),
        (2000, 6000),
        (0, 6000),
    ]);
    assert_eq!(r.offset(0.0).unwrap(), r);
    check_offset("zero_offset_is_identity", &r, 0.0, 2.0, None);
}

#[test]
fn square_inward_is_exact() {
    let r = rect(0, 0, 10_000, 10_000);
    let o = check_offset("square_inward_is_exact", &r, -1000.0, 2.0, None);
    assert!(o.xor(&rect(1000, 1000, 9000, 9000)).is_empty());
    assert_eq!(o.rings().len(), 1);
    assert_eq!(o.rings()[0].len(), 4, "острые углы остаются острыми");
    // Почти до конца — квадрат 2 × 2 мкм, до конца и дальше — пусто.
    let tiny = r.offset(-4.999).unwrap();
    assert!(tiny.xor(&rect(4999, 4999, 5001, 5001)).is_empty());
    assert!(r.offset(-5.0).unwrap().is_empty());
    assert!(r.offset(-7.0).unwrap().is_empty());
}

#[test]
fn square_outward_has_round_corners() {
    let r = rect(0, 0, 10_000, 10_000);
    let (d, tol) = (1000.0, 2.0);
    let o = check_offset("square_outward_has_round_corners", &r, d, tol, None);
    assert_eq!(o.rings().len(), 1);
    let exact = 100.0 + 40.0 * 1.0 + PI * 1.0;
    let got = mm(&o);
    // Хорды внутри дуг: площадь меньше, но не больше чем на длину дуг × допуск.
    let deficit = exact - got;
    assert!(deficit > -0.01, "{got} > {exact}");
    assert!(
        deficit < 2.0 * PI * 1.0 * tol / 1000.0 + 0.01,
        "{got} ≪ {exact}"
    );
    // Стороны остались прямыми и ровно на расстоянии d.
    assert_eq!(o.locate(p(5000, 11_000)), Location::Boundary);
    assert_eq!(o.locate(p(-1000, 5000)), Location::Boundary);
}

#[test]
fn hole_shrinks_then_closes() {
    let r = rect(0, 0, 10_000, 10_000).difference(&rect(4500, 4500, 5500, 5500));
    // Углы отверстия — вогнутые для детали: дуг нет, отверстие остаётся квадратом.
    let o = check_offset("hole_shrinks", &r, 400.0, 2.0, None);
    assert_eq!(o.rings().len(), 2);
    let hole = rect(4900, 4900, 5100, 5100);
    assert!(o.difference(&hole).xor(&o).is_empty() || o.intersection(&hole).is_empty());
    assert!(
        o.intersection(&hole).is_empty(),
        "отверстие 0,2 мм осталось"
    );
    assert!(!o.locate(p(5000, 5000)).eq(&Location::Inside));
    let closed = check_offset("hole_closes", &r, 600.0, 2.0, None);
    assert_eq!(closed.rings().len(), 1, "отверстие 1 мм закрылось");
    assert_eq!(closed.locate(p(5000, 5000)), Location::Inside);
}

#[test]
fn gap_closes_and_neck_breaks() {
    let two = rect(0, 0, 10_000, 10_000).union(&rect(11_000, 0, 21_000, 10_000));
    assert_eq!(
        check_offset("gap_open", &two, 400.0, 2.0, None)
            .rings()
            .len(),
        2
    );
    assert_eq!(
        check_offset("gap_closes", &two, 600.0, 2.0, None)
            .rings()
            .len(),
        1
    );
    let bell = two.union(&rect(9000, 4500, 12_000, 5500));
    assert_eq!(
        check_offset("neck_holds", &bell, -400.0, 2.0, None)
            .rings()
            .len(),
        1
    );
    assert_eq!(
        check_offset("neck_breaks", &bell, -600.0, 2.0, None)
            .rings()
            .len(),
        2
    );
}

#[test]
fn l_shape_both_ways() {
    let l = poly(&[
        (0, 0),
        (8000, 0),
        (8000, 3000),
        (3000, 3000),
        (3000, 8000),
        (0, 8000),
    ]);
    let mut rng = Rng::new(7);
    for d in [-1400.0, -700.0, -10.0, 10.0, 700.0, 2500.0] {
        check_offset("l_shape_both_ways", &l, d, 2.0, Some(&mut rng));
    }
}

#[test]
fn opening_rounds_convex_corners() {
    let r = rect(0, 0, 10_000, 10_000);
    let o = r.opening(2.0).unwrap();
    o.check().unwrap();
    let exact = 100.0 - (4.0 - PI) * 4.0;
    assert!((mm(&o) - exact).abs() < 0.05, "{} vs {exact}", mm(&o));
    assert_eq!(o.locate(p(5000, 0)), Location::Boundary, "стороны на месте");
    assert_eq!(o.locate(p(100, 100)), Location::Outside, "угол срезан");
}

#[test]
fn closing_is_what_a_cutter_leaves() {
    // Внутренний угол L-детали фреза радиуса 2 мм скругляет: + (1 − π/4) r².
    let l = poly(&[
        (0, 0),
        (8000, 0),
        (8000, 3000),
        (3000, 3000),
        (3000, 8000),
        (0, 8000),
    ]);
    let c = l.closing(2.0).unwrap();
    c.check().unwrap();
    let exact = mm(&l) + (1.0 - PI / 4.0) * 4.0;
    assert!((mm(&c) - exact).abs() < 0.05, "{} vs {exact}", mm(&c));
    assert_eq!(c.locate(p(3100, 3100)), Location::Inside, "угол заполнен");
    // Паз 1 мм фреза радиуса 0,6 мм не прорежет, радиуса 0,4 мм — прорежет.
    let slotted = rect(0, 0, 10_000, 10_000).difference(&rect(4500, 3000, 5500, 10_000));
    // Закрытый карман 1 мм заполняется целиком (с точностью до хорд в углах).
    let pocket = rect(0, 0, 10_000, 10_000).difference(&rect(4500, 3000, 5500, 7000));
    let filled = pocket.closing(0.6).unwrap();
    assert!(filled.xor(&rect(0, 0, 10_000, 10_000)).twice_area() < 2 * 20_000);
    // У открытого паза фреза заходит в устье и оставляет полукруглую выемку.
    let mouth = slotted.closing(0.6).unwrap();
    assert_eq!(mouth.locate(p(5000, 9990)), Location::Outside);
    assert_eq!(mouth.locate(p(5000, 9000)), Location::Inside);
    let cut = slotted.closing(0.4).unwrap();
    assert_eq!(cut.locate(p(5000, 6000)), Location::Outside);
    assert_eq!(
        cut.locate(p(4510, 3010)),
        Location::Inside,
        "в углу паза остаётся радиус"
    );
}

#[test]
fn sharp_tip_and_tiny_offsets() {
    // Остриё в 1°: снаружи — почти полукруг на вершине, внутри остриё укорачивается.
    let tip = poly(&[(0, 0), (100_000, -873), (100_000, 873)]);
    let mut rng = Rng::new(11);
    for d in [-300.0, -50.0, 0.5, 1.0, 3.0, 50.0, 3000.0] {
        check_offset("sharp_tip_and_tiny_offsets", &tip, d, 2.0, Some(&mut rng));
    }
}

#[test]
fn huge_offset_of_tiny_shape_is_a_circle() {
    let t = poly(&[(0, 0), (10, 0), (0, 10)]);
    let o = check_offset("huge_offset_of_tiny_shape", &t, 100_000.0, 5.0, None);
    let exact = PI * 100.0 * 100.0;
    assert!(
        (mm(&o) - exact).abs() / exact < 1e-3,
        "{} vs {exact}",
        mm(&o)
    );
}

#[test]
fn circles_from_paths() {
    let c = Path::circle(Vec2::new(0.0, 0.0), 5.0)
        .to_region(0.002)
        .unwrap();
    let grow = c.offset(1.0).unwrap();
    let shrink = c.offset(-1.0).unwrap();
    assert!((mm(&grow) - PI * 36.0).abs() < 0.05, "{}", mm(&grow));
    assert!((mm(&shrink) - PI * 16.0).abs() < 0.05, "{}", mm(&shrink));
    let mut rng = Rng::new(3);
    check_offset("circles_from_paths", &c, 1000.0, 2.0, Some(&mut rng));
    check_offset("circles_from_paths", &c, -1000.0, 2.0, Some(&mut rng));
}

#[test]
fn touching_squares() {
    let r = rect(0, 0, 10_000, 10_000).union(&rect(10_000, 10_000, 20_000, 20_000));
    assert_eq!(r.rings().len(), 2);
    let mut rng = Rng::new(5);
    let out = check_offset("touching_squares", &r, 500.0, 2.0, Some(&mut rng));
    assert_eq!(out.rings().len(), 1);
    let inn = check_offset("touching_squares", &r, -500.0, 2.0, Some(&mut rng));
    assert_eq!(inn.rings().len(), 2);
}

#[test]
fn near_grid_limit() {
    let h = GRID_LIMIT - 5_000_000;
    let r = rect(h - 10_000_000, h - 10_000_000, h, h);
    let mut rng = Rng::new(9);
    check_offset("near_grid_limit", &r, 2_000_000.0, 50.0, Some(&mut rng));
    assert!(
        r.offset(10_000.0).is_err(),
        "за предел сетки — ошибка, а не паника"
    );
    assert!(r.offset(f64::NAN).is_err());
}

// ───────────────────────────── Случайные области ─────────────────────────────

/// Один случайный случай: область, расстояние, допуск; все проверки и вложенность.
fn random_offset_case(rng: &mut Rng) {
    let scale = *rng.pick(&[40, 1_000, 100_000, 10_000_000]);
    let r = random_region(rng, scale);
    let mag = if rng.chance(0.15) {
        rng.range_f64(0.5, 6.0)
    } else {
        scale as f64 * 10f64.powf(rng.range_f64(-3.0, 0.0))
    };
    let d = if rng.chance(0.5) { mag } else { -mag };
    let tol = (mag * 10f64.powf(rng.range_f64(-5.0, -2.0))).max(0.5);
    let seed = rng.next_u64();
    check_shape_or_card("random_offsets_agree_with_distances", &shape_of(&r), |r| {
        let mut probes = Rng::new(seed);
        let o = check_offset(
            "random_offsets_agree_with_distances",
            r,
            d,
            tol,
            Some(&mut probes),
        );
        // Вложенность: при большем d результат не меньше.
        let step = tol + 5.0;
        let o2 = r
            .offset_tol((d + step) / 1000.0, tol / 1000.0)
            .expect("в пределах сетки");
        assert!(
            o.difference(&o2).is_empty(),
            "вложенность нарушена: d = {d}, d₂ = {}",
            d + step
        );
    });
}

#[test]
fn random_offsets_agree_with_distances() {
    for_cases(
        "random_offsets_agree_with_distances",
        300,
        random_offset_case,
    );
}

/// Случаи, на которых когда-то падала проверка (запас на округление был меньше √2).
#[test]
fn regression_offset_seeds() {
    for seed in [
        11032491998542159587,
        3788776528944386978,
        9689675955521663328,
    ] {
        random_offset_case(&mut Rng::new(seed));
    }
}

#[test]
fn random_paths_flatten_within_tolerance() {
    for_cases("random_paths_flatten_within_tolerance", 500, |rng| {
        let c = Vec2::new(rng.range_f64(-500.0, 500.0), rng.range_f64(-500.0, 500.0));
        let size = 10f64.powf(rng.range_f64(-1.0, 2.5));
        let path = match rng.below(3) {
            0 => Path::circle(c, size),
            1 => Path::rounded_rect(
                c.x,
                c.y,
                c.x + size * rng.range_f64(0.2, 3.0),
                c.y + size,
                size * rng.range_f64(0.0, 0.6),
            ),
            _ => {
                let t = rng.range_f64(0.0, 2.0 * PI);
                Path::slot(
                    c,
                    c + Vec2::new(t.cos(), t.sin()) * size,
                    size * rng.range_f64(0.05, 1.0),
                )
            }
        };
        let tol = 10f64.powf(rng.range_f64(-3.3, -1.0));
        let ring = path.to_ring(tol).unwrap();
        let exact = path.area();
        let flat = ring.twice_area() as f64 / 2.0 * 1e-6;
        let deficit = exact - flat;
        let len = path.length();
        let t = tol.max(MIN_CHORD_TOL_MM);
        assert!(
            deficit > -len * 0.000_71 - 1e-9 && deficit < len * (t + 0.000_71) + 1e-9,
            "площадь {flat} против {exact}, длина {len}, допуск {t}"
        );
        assert!(ring.is_ccw());
    });
}
