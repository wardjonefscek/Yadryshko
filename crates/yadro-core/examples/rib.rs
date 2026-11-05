//! Пример ступени 1: заготовка нервюры булевыми операциями.
//!
//! Профиль NACA 2412 (хорда 200 мм) минус облегчающие отверстия и паз под
//! лонжерон, плюс язычок крепления. Картинка — в `target/rib.svg`.
//!
//! Запуск: `cargo run -p yadro-core --example rib`

use yadro_core::svg::{Style, Svg};
use yadro_core::{Region, Ring, Vec2};

/// Профиль NACA 4-й серии с замкнутой задней кромкой, точки против часовой.
fn naca4(m: f64, p: f64, t: f64, chord: f64, n: usize) -> Vec<Vec2> {
    let camber = |x: f64| {
        if x < p {
            (
                m / (p * p) * (2.0 * p * x - x * x),
                2.0 * m / (p * p) * (p - x),
            )
        } else {
            (
                m / ((1.0 - p) * (1.0 - p)) * (1.0 - 2.0 * p + 2.0 * p * x - x * x),
                2.0 * m / ((1.0 - p) * (1.0 - p)) * (p - x),
            )
        }
    };
    let thick = |x: f64| {
        5.0 * t
            * (0.2969 * x.sqrt() - 0.1260 * x - 0.3516 * x * x + 0.2843 * x.powi(3)
                - 0.1036 * x.powi(4))
    };
    let side = |x: f64, upper: bool| {
        let (yc, dy) = camber(x);
        let (yt, th) = (thick(x), dy.atan());
        let s = if upper { 1.0 } else { -1.0 };
        Vec2::new(
            (x - s * yt * th.sin()) * chord,
            (yc + s * yt * th.cos()) * chord,
        )
    };
    // Косинусная разбивка: гуще у кромок.
    let xs: Vec<f64> = (0..=n)
        .map(|i| 0.5 * (1.0 - (std::f64::consts::PI * i as f64 / n as f64).cos()))
        .collect();
    let mut pts: Vec<Vec2> = xs.iter().rev().map(|&x| side(x, false)).collect();
    pts.extend(xs.iter().skip(1).map(|&x| side(x, true)));
    pts.reverse(); // нижняя поверхность от носка к хвосту, затем верхняя обратно
    pts
}

fn circle(cx: f64, cy: f64, r: f64, n: usize) -> Ring {
    let pts: Vec<Vec2> = (0..n)
        .map(|i| {
            let a = std::f64::consts::TAU * i as f64 / n as f64;
            Vec2::new(cx + r * a.cos(), cy + r * a.sin())
        })
        .collect();
    Ring::from_mm(&pts).expect("окружность на сетке")
}

fn rect(x0: f64, y0: f64, x1: f64, y1: f64) -> Region {
    Region::from_ring(
        Ring::from_mm(&[
            Vec2::new(x0, y0),
            Vec2::new(x1, y0),
            Vec2::new(x1, y1),
            Vec2::new(x0, y1),
        ])
        .expect("на сетке"),
    )
}

fn main() {
    let profile = Region::from_ring(
        Ring::from_mm(&naca4(0.02, 0.4, 0.12, 200.0, 80)).expect("профиль на сетке"),
    );
    let holes = Region::from_rings_nonzero(&[
        circle(28.0, 2.5, 6.0, 64),
        circle(90.0, 2.0, 7.5, 64),
        circle(130.0, 1.2, 5.5, 64),
    ]);
    let spar = rect(56.0, -4.0, 62.0, 30.0);
    let tab = rect(195.0, -3.0, 210.0, 3.0);

    let rib = profile.difference(&holes).difference(&spar).union(&tab);
    rib.check().expect("результат корректен");

    println!("Профиль:            {:9.3} мм²", profile.area_mm2());
    println!("Нервюра:            {:9.3} мм²", rib.area_mm2());
    println!(
        "Контуров:           {} (отверстий: {})",
        rib.rings().len(),
        rib.rings().iter().filter(|r| !r.is_ccw()).count()
    );
    println!(
        "Вершин:             {}",
        rib.rings().iter().map(Ring::len).sum::<usize>()
    );

    let mut svg = Svg::new();
    svg.region(&profile, Style::line("#9aa5b1"), Style::line("#9aa5b1"));
    svg.region(
        &rib,
        Style::filled("#1f5fbf", 0.35).width(1.5),
        Style::line("#c0392b").width(1.5),
    );
    let out = concat!(env!("CARGO_MANIFEST_DIR"), "/../../target/rib.svg");
    svg.save(out, 1200.0).expect("запись SVG");
    println!("Картинка:           target/rib.svg");
}
