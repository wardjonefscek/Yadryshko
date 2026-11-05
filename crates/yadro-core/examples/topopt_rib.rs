//! Пример шага 4: нервюра Clark Y из файла базы UIUC с облегчением по полю
//! плотности «как из TopOpt».
//!
//! 1. Профиль читается из clarky.dat, хорда 200 мм, задняя кромка не тоньше 0,6 мм.
//! 2. Поле плотности на сетке 1 мм: обшивка вдоль контура, стенка у лонжерона и
//!    раскосы фермы (так
//!    выглядит типичный результат топооптимизации нервюры; настоящий даст Cad+).
//! 3. Контур по порогу 0,5 с уплотнением ×3, без островов и дыр меньше 4 мм²,
//!    обрезка профилем, паз под лонжерон, проверка под фрезу ⌀3 мм.
//!
//! Картинка — в `target/topopt_rib.svg`. Запуск:
//! `cargo run -p yadro-core --release --example topopt_rib`

use yadro_core::svg::{Style, Svg};
use yadro_core::{Airfoil, ContourOptions, DensityField, Path, Vec2};

const CLARK_Y: &str = include_str!("../tests/data/clarky.dat");
const CHORD: f64 = 200.0;
const CELL: f64 = 1.0;

fn main() {
    let foil = Airfoil::parse_dat(CLARK_Y)
        .expect("файл профиля")
        .with_min_te_gap(0.6 / CHORD);
    let xf = Airfoil::placement(CHORD, 0.0, 0.0, Vec2::ZERO);
    let rib = foil.to_region(&xf).expect("профиль на сетке");

    // Поле плотности над прямоугольником вокруг профиля. Плотность меняется плавно
    // (как после фильтра TopOpt): ramp(s) — от 0,02 до 1 на полосе шириной 1,5 мм.
    let ramp = |s: f64| (0.5 + s / 1.5).clamp(0.02, 1.0);
    let surf = |x: f64| {
        let u = (x / CHORD).clamp(0.0, 1.0);
        (foil.upper_y(u) * CHORD, foil.lower_y(u) * CHORD)
    };
    // Узлы фермы: попеременно у низа и у верха, шаг 32 мм.
    let nodes: Vec<Vec2> = (0..5)
        .map(|k| {
            let x = 22.0 + 32.0 * k as f64;
            let (up, lo) = surf(x);
            Vec2::new(x, if k % 2 == 0 { lo + 1.5 } else { up - 1.5 })
        })
        .collect();
    let seg_dist = |p: Vec2, a: Vec2, b: Vec2| {
        let ab = b - a;
        let t = ((p - a).dot(ab) / ab.length_sq()).clamp(0.0, 1.0);
        p.distance(a + ab * t)
    };
    let (y0, y1) = (-6.0, 26.0);
    let (nx, ny) = ((CHORD / CELL) as usize, ((y1 - y0) / CELL) as usize);
    let field = DensityField::from_fn(nx, ny, CELL, Vec2::new(0.0, y0), |i, j| {
        let p = Vec2::new((i as f64 + 0.5) * CELL, y0 + (j as f64 + 0.5) * CELL);
        let (up, lo) = surf(p.x);
        let skin = (up - p.y).min(p.y - lo); // глубина под обшивкой
        let truss = nodes
            .windows(2)
            .map(|w| seg_dist(p, w[0], w[1]))
            .fold(f64::INFINITY, f64::min);
        // Стенка вокруг лонжерона: нагрузка с него идёт в нервюру здесь.
        let web = (p.x - 51.0).abs();
        ramp(2.5 - skin)
            .max(ramp(1.8 - truss))
            .max(ramp(18.0 - p.x))
            .max(ramp(7.0 - web))
    })
    .expect("поле");
    let opts = ContourOptions {
        threshold: 0.5,
        refine: 3,
        min_island_mm2: 4.0,
        min_hole_mm2: 4.0,
    };
    let lightened = field.contour(&opts).expect("контур").intersection(&rib);
    // Паз под лонжерон 6 мм сверху до средней линии, дно скруглено.
    let mid = foil.camber(0.255) * CHORD;
    let spar = Path::rounded_rect(48.0, mid, 54.0, 30.0, 1.0)
        .to_region(0.002)
        .expect("паз");
    let part = lightened.difference(&spar);
    let machined = part.closing(1.5).expect("фреза");
    part.check().expect("деталь корректна");

    println!("Профиль:          {}", foil.name());
    let (xt, t) = foil.max_thickness();
    println!(
        "Толщина:          {:.1}% на {:.0}% хорды",
        t * 100.0,
        xt * 100.0
    );
    println!("Сплошная:         {:9.1} мм²", rib.area_mm2());
    println!(
        "Облегчённая:      {:9.1} мм² ({:.0}% массы)",
        part.area_mm2(),
        100.0 * part.area_mm2() / rib.area_mm2()
    );
    println!(
        "Контуров:         {} (отверстий: {})",
        part.rings().len(),
        part.rings().iter().filter(|r| !r.is_ccw()).count()
    );
    println!(
        "После фрезы ⌀3:   {:9.1} мм² (+{:.2} мм² во внутренних углах)",
        machined.area_mm2(),
        machined.area_mm2() - part.area_mm2()
    );

    let mut svg = Svg::new();
    svg.region(&rib, Style::line("#9aa5b1"), Style::line("#9aa5b1"));
    svg.region(
        &part,
        Style::filled("#1f5fbf", 0.4).width(1.2),
        Style::line("#c0392b").width(1.2),
    );
    let out = concat!(env!("CARGO_MANIFEST_DIR"), "/../../target/topopt_rib.svg");
    svg.save(out, 1400.0).expect("запись SVG");
    println!("Картинка:         target/topopt_rib.svg");
}
