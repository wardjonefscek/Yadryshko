//! Пример ступени 1: нервюра булевыми операциями, дугами и эквидистантой.
//!
//! Профиль NACA 2412 (хорда 200 мм) минус облегчающие отверстия (настоящие
//! окружности) и паз под лонжерон со скруглённым дном, плюс язычок крепления.
//! Для фрезы ⌀3 мм: траектория центра фрезы (эквидистанта наружу на 1,5 мм) и то,
//! что фреза оставит (closing: внутренние углы получают радиус 1,5 мм).
//! Картинка — в `target/rib.svg`.
//!
//! Запуск: `cargo run -p yadro-core --example rib`

use yadro_core::svg::{Style, Svg};
use yadro_core::{Path, Region, Ring, Vec2};

/// Допуск хорды для дуг, мм.
const TOL: f64 = 0.002;
/// Радиус фрезы, мм.
const CUTTER: f64 = 1.5;

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

fn region(path: &Path) -> Region {
    path.to_region(TOL).expect("контур на сетке")
}

fn main() {
    let profile = Region::from_ring(
        Ring::from_mm(&naca4(0.02, 0.4, 0.12, 200.0, 80)).expect("профиль на сетке"),
    );
    let holes = region(&Path::circle(Vec2::new(28.0, 2.5), 6.0))
        .union(&region(&Path::circle(Vec2::new(90.0, 2.0), 7.5)))
        .union(&region(&Path::circle(Vec2::new(130.0, 1.2), 5.5)));
    // Паз под лонжерон 6 мм сверху, дно скруглено радиусом 1 мм.
    let spar = region(&Path::rounded_rect(56.0, -4.0, 62.0, 30.0, 1.0));
    let tab = region(&Path::rounded_rect(195.0, -3.0, 210.0, 3.0, 1.0));

    let rib = profile.difference(&holes).difference(&spar).union(&tab);
    rib.check().expect("результат корректен");
    let tool_path = rib.offset(CUTTER).expect("эквидистанта");
    let machined = rib.closing(CUTTER).expect("closing");

    println!("Профиль:              {:9.3} мм²", profile.area_mm2());
    println!("Нервюра:              {:9.3} мм²", rib.area_mm2());
    println!(
        "Контуров:             {} (отверстий: {})",
        rib.rings().len(),
        rib.rings().iter().filter(|r| !r.is_ccw()).count()
    );
    println!(
        "Вершин:               {}",
        rib.rings().iter().map(Ring::len).sum::<usize>()
    );
    println!(
        "После фрезы ⌀{:.0} мм:   {:9.3} мм² (+{:.3} мм² во внутренних углах)",
        2.0 * CUTTER,
        machined.area_mm2(),
        machined.area_mm2() - rib.area_mm2()
    );

    let mut svg = Svg::new();
    svg.region(&profile, Style::line("#9aa5b1"), Style::line("#9aa5b1"));
    svg.region(
        &rib,
        Style::filled("#1f5fbf", 0.35).width(1.5),
        Style::line("#c0392b").width(1.5),
    );
    svg.region(
        &tool_path,
        Style::line("#27ae60").width(1.0),
        Style::line("#27ae60").width(1.0),
    );
    let out = concat!(env!("CARGO_MANIFEST_DIR"), "/../../target/rib.svg");
    svg.save(out, 1200.0).expect("запись SVG");
    println!("Картинка:             target/rib.svg");
}
