//! Ворота ступени 1, шаг 4б: контур из поля плотности TopOpt.
//!
//! Независимая проверка — сами отсчёты: каждый центр элемента с плотностью
//! заметно выше порога должен оказаться внутри области, заметно ниже — снаружи.
//! Плюс точные случаи: поле из единиц — ровно расчётная область, линейное поле —
//! ровно полуплоскость, один элемент — ромб известной площади; седла, острова,
//! дыры, гладкость границы с уплотнением.

mod common;

use common::*;
use yadro_core::testkit::{Rng, for_cases};
use yadro_core::{Airfoil, ContourOptions, DensityField, Location, Region, Vec2};

fn opts(threshold: f64, refine: usize) -> ContourOptions {
    ContourOptions {
        threshold,
        refine,
        ..ContourOptions::default()
    }
}

fn grid_point(v: Vec2) -> yadro_core::GridPoint {
    yadro_core::GridPoint::from_mm(v.x, v.y).unwrap()
}

/// Каждый центр элемента, достаточно далёкий от порога, классифицирован верно.
fn check_samples(f: &DensityField, r: &Region, o: &ContourOptions) {
    let (nx, ny) = f.size();
    let (mut lo, mut hi) = (f64::INFINITY, f64::NEG_INFINITY);
    for j in 0..ny {
        for i in 0..nx {
            lo = lo.min(f.get(i, j));
            hi = hi.max(f.get(i, j));
        }
    }
    // Пересечение с ребром соседней клетки не ближе 3 мкм к отсчёту: |v − t| / разброс
    // (с запасом 1,25 на выбросы Catmull-Rom) × шаг ≥ 0,003 мм.
    let step = f.cell_mm() / o.refine.max(1) as f64;
    let spread =
        1.25 * (hi - lo).max(1e-12) + (o.threshold - lo).abs().max((hi - o.threshold).abs());
    let margin = 0.003 / step * spread;
    for j in 0..ny {
        for i in 0..nx {
            let v = f.get(i, j);
            if (v - o.threshold).abs() < margin {
                continue;
            }
            let want = if v >= o.threshold {
                Location::Inside
            } else {
                Location::Outside
            };
            let got = r.locate(grid_point(f.center(i, j)));
            assert_eq!(got, want, "элемент ({i}, {j}) с плотностью {v}");
        }
    }
}

fn domain_region(f: &DensityField) -> Region {
    let (a, b) = f.domain();
    let (a, b) = (grid_point(a), grid_point(b));
    rect(a.x(), a.y(), b.x(), b.y())
}

#[test]
fn full_field_is_exactly_the_domain() {
    let f = DensityField::from_fn(7, 5, 1.25, Vec2::new(-3.5, 2.0), |_, _| 1.0).unwrap();
    for refine in [1, 2, 3, 5] {
        let r = f.contour(&opts(0.5, refine)).unwrap();
        r.check().unwrap();
        assert!(r.xor(&domain_region(&f)).is_empty(), "refine {refine}");
        assert_eq!(r.rings()[0].len(), 4, "углы области острые");
    }
}

#[test]
fn empty_field_is_empty() {
    let f = DensityField::from_fn(6, 4, 1.0, Vec2::ZERO, |_, _| 0.1).unwrap();
    assert!(f.contour(&opts(0.5, 1)).unwrap().is_empty());
    assert!(f.contour(&opts(0.5, 4)).unwrap().is_empty());
}

#[test]
fn single_element_is_a_diamond() {
    let f = DensityField::from_fn(5, 5, 1.0, Vec2::ZERO, |i, j| {
        if (i, j) == (2, 2) { 1.0 } else { 0.0 }
    })
    .unwrap();
    let r = f.contour(&opts(0.5, 1)).unwrap();
    assert_eq!(r.rings().len(), 1);
    assert_eq!(r.rings()[0].len(), 4);
    assert_eq!(
        r.twice_area(),
        1_000_000,
        "ромб с полудиагоналями 0,5 мм — 0,5 мм²"
    );
    // Отсев островов меньше 1 мм².
    let mut o = opts(0.5, 1);
    o.min_island_mm2 = 1.0;
    assert!(f.contour(&o).unwrap().is_empty());
}

#[test]
fn linear_fields_give_exact_half_planes() {
    // Плотность растёт по x: порог 0,5 ровно на x = 5 мм.
    let f = DensityField::from_fn(10, 4, 1.0, Vec2::ZERO, |i, _| (i as f64 + 0.5) / 10.0).unwrap();
    for refine in [1, 2, 4] {
        let r = f.contour(&opts(0.5, refine)).unwrap();
        assert!(
            r.xor(&rect(5000, 0, 10_000, 4000)).is_empty(),
            "refine {refine}"
        );
    }
    // По диагонали: материал там, где x + y ≥ 6 мм. За краем поле продолжается
    // постоянным по нормали к краю, поэтому в полуклетке у края линия загибается;
    // между центрами элементов — ровно прямая.
    let g = DensityField::from_fn(8, 8, 1.0, Vec2::ZERO, |i, j| {
        (i as f64 + j as f64 + 1.0) / 12.0
    })
    .unwrap();
    let r = g.contour(&opts(0.5, 1)).unwrap();
    let inner = rect(500, 500, 7500, 7500);
    let want = inner.intersection(&poly(&[
        (6000, 0),
        (8000, 0),
        (8000, 8000),
        (0, 8000),
        (0, 6000),
    ]));
    assert!(r.intersection(&inner).xor(&want).is_empty());
    // У края — линия под прямым углом к краю: на нижней кромке x = 5,5 мм.
    assert_eq!(
        r.locate(grid_point(Vec2::new(5.5, 0.0))),
        Location::Boundary
    );
}

#[test]
fn saddles_follow_the_bilinear_field() {
    let f = DensityField::from_fn(6, 6, 1.0, Vec2::ZERO, |i, j| ((i + j) % 2) as f64).unwrap();
    // Порог выше седловой точки (0,5): клетки-единицы — отдельные острова.
    let apart = f.contour(&opts(0.6, 1)).unwrap();
    apart.check().unwrap();
    assert_eq!(apart.rings().len(), 18);
    assert!(apart.rings().iter().all(|r| r.is_ccw()));
    // Порог на седле: единицы соединены по диагонали в одну деталь.
    let joined = f.contour(&opts(0.5, 1)).unwrap();
    joined.check().unwrap();
    assert_eq!(joined.rings().iter().filter(|r| r.is_ccw()).count(), 1);
    check_samples(&f, &joined, &opts(0.5, 1));
}

#[test]
fn small_islands_and_holes_are_dropped() {
    // Квадрат 10 × 10 элементов с дырой в один элемент и отдельной соринкой.
    let f = DensityField::from_fn(20, 14, 1.0, Vec2::ZERO, |i, j| {
        let block = (1..11).contains(&i) && (1..11).contains(&j) && (i, j) != (5, 5);
        let speck = (i, j) == (16, 7);
        if block || speck { 1.0 } else { 0.0 }
    })
    .unwrap();
    let raw = f.contour(&opts(0.5, 1)).unwrap();
    assert_eq!(raw.rings().len(), 3, "деталь, дыра и соринка");
    let mut o = opts(0.5, 1);
    o.min_island_mm2 = 1.0;
    o.min_hole_mm2 = 1.0;
    let clean = f.contour(&o).unwrap();
    clean.check().unwrap();
    assert_eq!(clean.rings().len(), 1);
    assert_eq!(
        clean.locate(grid_point(f.center(5, 5))),
        Location::Inside,
        "дыра заполнена"
    );
    assert_eq!(
        clean.locate(grid_point(f.center(16, 7))),
        Location::Outside,
        "соринка убрана"
    );
}

#[test]
fn island_inside_a_hole_survives_hole_rules() {
    // Кольцо с островом внутри: большой остров остаётся, маленькая дыра вокруг
    // него не «съедает» его при заполнении (дыра большая — не заполняется).
    let f = DensityField::from_fn(21, 21, 1.0, Vec2::ZERO, |i, j| {
        let d = ((i as f64 - 10.0).powi(2) + (j as f64 - 10.0).powi(2)).sqrt();
        if !(3.0..7.0).contains(&d) { 1.0 } else { 0.0 }
    })
    .unwrap();
    let r = f.contour(&opts(0.5, 2)).unwrap();
    r.check().unwrap();
    assert_eq!(
        r.rings().iter().filter(|r| r.is_ccw()).count(),
        2,
        "внешняя деталь и остров"
    );
    let filled = r.drop_small(0.0, 1000.0);
    filled.check().unwrap();
    assert_eq!(filled.rings().len(), 1, "дыра заполнена вместе с островом");
}

#[test]
fn nested_small_things_go_together() {
    // Большой блок; в нём дыра 3 × 3 с островом в один элемент. Отдельно — малый
    // остров 4 × 4 со своей дырой в один элемент.
    let f = DensityField::from_fn(30, 16, 1.0, Vec2::ZERO, |i, j| {
        let big = (1..14).contains(&i) && (1..14).contains(&j);
        let pocket = (5..8).contains(&i) && (5..8).contains(&j);
        let isle = (i, j) == (6, 6);
        let small = (18..22).contains(&i) && (5..9).contains(&j) && (i, j) != (19, 6);
        if (big && !pocket) || isle || small {
            1.0
        } else {
            0.0
        }
    })
    .unwrap();
    let r = f.contour(&opts(0.5, 1)).unwrap();
    assert_eq!(
        r.rings().len(),
        5,
        "блок, дыра, островок, малый остров, его дыра"
    );
    // Малый остров уходит вместе со своей дырой — дыра не превращается в соринку.
    let a = r.drop_small(20.0, 0.0);
    a.check().unwrap();
    assert_eq!(
        a.rings().len(),
        2,
        "блок и его дыра; островок в дыре тоже мал"
    );
    // Дыра в блоке заполняется вместе с островком в ней.
    let b = r.drop_small(0.0, 20.0);
    b.check().unwrap();
    assert_eq!(b.rings().len(), 2, "блок сплошной, малый остров тоже");
    assert_eq!(b.locate(grid_point(f.center(5, 5))), Location::Inside);
    // Регрессия: так когда-то ложилась проверка — дыра внутри убранного острова.
    let c = r.drop_small(20.0, 20.0);
    assert_eq!(c.rings().len(), 1);
}

#[test]
fn refine_makes_circles_rounder() {
    // Конус: плотность падает с расстоянием от центра, порог — окружность радиуса 8 мм.
    let (c, rad) = (Vec2::new(12.0, 12.0), 8.0);
    let f = DensityField::from_fn(24, 24, 1.0, Vec2::ZERO, |i, j| {
        let p = Vec2::new(i as f64 + 0.5, j as f64 + 0.5);
        0.5 + (rad - p.distance(c)) / 6.0
    })
    .unwrap();
    let err = |refine: usize| {
        let r = f.contour(&opts(0.5, refine)).unwrap();
        r.check().unwrap();
        r.rings()[0]
            .points()
            .iter()
            .map(|q| (q.to_mm().distance(c) - rad).abs())
            .fold(0.0, f64::max)
    };
    let (e1, e4) = (err(1), err(4));
    assert!(e1 < 0.05, "без уплотнения отклонение {e1} мм");
    assert!(
        e4 < e1 / 2.0,
        "уплотнение должно сглаживать: {e4} против {e1}"
    );
}

#[test]
fn rib_from_density_inside_an_airfoil() {
    // Поле «как из TopOpt» в прямоугольнике под нервюру, обрезанное профилем.
    let a = Airfoil::naca4_code("2412", 61, true).unwrap();
    let rib = a
        .to_region(&Airfoil::placement(200.0, 0.0, 0.0, Vec2::ZERO))
        .unwrap();
    let f = DensityField::from_fn(100, 20, 2.0, Vec2::new(0.0, -10.0), |i, j| {
        let (x, y) = (i as f64 * 2.0 + 1.0, j as f64 * 2.0 - 9.0);
        let truss = ((x / 25.0 * std::f64::consts::PI).sin() * 8.0 - y).abs() < 2.0;
        let skin = !(-6.0..10.0).contains(&y);
        if truss || skin { 1.0 } else { 0.05 }
    })
    .unwrap();
    let mut o = opts(0.5, 3);
    o.min_island_mm2 = 2.0;
    o.min_hole_mm2 = 2.0;
    let part = f.contour(&o).unwrap().intersection(&rib);
    part.check().unwrap();
    assert!(part.area_mm2() > 0.2 * rib.area_mm2() && part.area_mm2() < rib.area_mm2());
}

/// Случаи, на которых когда-то ложилась проверка: дыра внутри убранного острова
/// превращалась в соринку; порог площади сравнивался с округлением вниз.
#[test]
fn regression_density_seeds() {
    for seed in [
        4445012358684040281,
        12581959153181244188,
        15812848193996606630,
    ] {
        random_field_case(&mut Rng::new(seed));
    }
}

#[test]
fn random_fields_classify_samples() {
    for_cases("random_fields_classify_samples", 400, random_field_case);
}

fn random_field_case(rng: &mut Rng) {
    let (nx, ny) = (rng.range_i64(1, 24) as usize, rng.range_i64(1, 24) as usize);
    let cell = *rng.pick(&[0.25, 0.5, 1.0, 2.5]);
    let origin = Vec2::new(
        rng.range_i64(-20_000, 20_000) as f64 * 0.001,
        rng.range_i64(-20_000, 20_000) as f64 * 0.001,
    );
    let kind = rng.below(4);
    let bumps: Vec<(f64, f64, f64, f64)> = (0..rng.range_i64(1, 5))
        .map(|_| {
            (
                rng.range_f64(0.0, nx as f64),
                rng.range_f64(0.0, ny as f64),
                rng.range_f64(0.5, 6.0),
                rng.range_f64(-1.0, 1.5),
            )
        })
        .collect();
    let values: Vec<f64> = (0..nx * ny)
        .map(|k| {
            let (i, j) = ((k % nx) as f64 + 0.5, (k / nx) as f64 + 0.5);
            match kind {
                0 => rng.unit_f64(),
                1 => f64::from(u8::from(rng.chance(0.5))),
                2 => *rng.pick(&[0.0, 0.5, 1.0]),
                _ => bumps
                    .iter()
                    .map(|&(x, y, s, a)| a * (-((i - x).powi(2) + (j - y).powi(2)) / (s * s)).exp())
                    .sum(),
            }
        })
        .collect();
    let f = DensityField::new(nx, ny, cell, origin, values).unwrap();
    let o = opts(
        *rng.pick(&[0.5, 0.3, 0.7, 0.01]),
        rng.range_i64(1, 3) as usize,
    );
    let r = f.contour(&o).unwrap();
    r.check().unwrap();
    assert!(
        r.difference(&domain_region(&f)).is_empty(),
        "за пределами области"
    );
    check_samples(&f, &r, &o);
    // Отсев по площади: не остаётся ни островов, ни дыр меньше порогов.
    let (mi, mh) = (
        rng.range_f64(0.0, 4.0) * cell * cell,
        rng.range_f64(0.0, 4.0) * cell * cell,
    );
    let clean = r.drop_small(mi, mh);
    clean.check().unwrap();
    let areas = |g: &Region| {
        g.rings()
            .iter()
            .map(|x| format!("{:.6}", x.twice_area() as f64 / 2e6))
            .collect::<Vec<_>>()
            .join(", ")
    };
    let info = format!("до: [{}]\nпосле: [{}]", areas(&r), areas(&clean));
    for ring in clean.rings() {
        let a = ring.area_mm2();
        if ring.is_ccw() {
            assert!(a >= mi * (1.0 - 1e-9), "остался остров {a} < {mi}\n{info}");
        } else {
            assert!(a >= mh * (1.0 - 1e-9), "осталась дыра {a} < {mh}\n{info}");
        }
    }
}
