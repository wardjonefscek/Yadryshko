//! Ворота ступени 1, шаг 4а: профили крыла.
//!
//! - Настоящий файл базы UIUC (Clark Y) читается, толщина совпадает с паспортной.
//! - Selig и Lednicer одного профиля дают одно и то же.
//! - Повёрнутый, растянутый и сдвинутый файл приводится обратно к хорде 1.
//! - Мусор в оформлении (CRLF, табуляция, запятые, комментарии, повторы точек)
//!   не мешает, а настоящие ошибки называются с номером строки.
//! - Площадь контура NACA совпадает с интегралом формулы толщины.

use yadro_core::airfoil::AirfoilError;
use yadro_core::testkit::{Rng, for_cases};
use yadro_core::{Airfoil, Vec2, Xform2};

const CLARK_Y: &str = include_str!("data/clarky.dat");

fn same(a: &Airfoil, b: &Airfoil, tol: f64) {
    for x in (0..=200).map(|i| i as f64 / 200.0) {
        let (du, dl) = (
            (a.upper_y(x) - b.upper_y(x)).abs(),
            (a.lower_y(x) - b.lower_y(x)).abs(),
        );
        assert!(
            du <= tol && dl <= tol,
            "x = {x}: верх {du}, низ {dl}, допуск {tol}"
        );
    }
}

/// Контуры совпадают по точкам (одинаковое число, каждая не дальше `tol`).
fn same_outline(a: &Airfoil, b: &Airfoil, tol: f64) {
    let (pa, pb) = (a.outline(&Xform2::IDENTITY), b.outline(&Xform2::IDENTITY));
    assert_eq!(pa.len(), pb.len(), "разное число точек");
    for (p, q) in pa.iter().zip(&pb) {
        assert!(p.distance(*q) <= tol, "{p:?} / {q:?}, допуск {tol}");
    }
}

fn lednicer(a: &Airfoil, decimals: usize) -> String {
    let mut s = format!(
        "{}\n{}. {}.\n\n",
        a.name(),
        a.upper().len(),
        a.lower().len()
    );
    for p in a.upper() {
        s.push_str(&format!("{:.*} {:.*}\n", decimals, p.x, decimals, p.y));
    }
    s.push('\n');
    for p in a.lower() {
        s.push_str(&format!("{:.*} {:.*}\n", decimals, p.x, decimals, p.y));
    }
    s
}

/// Simpson по [0, 1] на `n` (чётном) отрезках.
fn simpson(f: impl Fn(f64) -> f64, n: usize) -> f64 {
    let h = 1.0 / n as f64;
    let s: f64 = (0..=n)
        .map(|i| {
            let w = if i == 0 || i == n {
                1.0
            } else if i % 2 == 1 {
                4.0
            } else {
                2.0
            };
            w * f(i as f64 * h)
        })
        .sum();
    s * h / 3.0
}

#[test]
fn clark_y_from_uiuc() {
    let a = Airfoil::parse_dat(CLARK_Y).unwrap();
    assert_eq!(a.name(), "CLARK Y AIRFOIL");
    let (x, t) = a.max_thickness();
    assert!((t - 0.117).abs() < 0.001, "толщина {t}");
    assert!((0.25..0.32).contains(&x), "место толщины {x}");
    assert!(
        (a.te_gap() - 0.0011986).abs() < 1e-6,
        "зазор кромки {}",
        a.te_gap()
    );
    assert!(
        a.upper().iter().all(|p| p.y >= -1e-12),
        "верх Clark Y не ниже хорды"
    );
    // Нервюра 200 мм, носок вверх на 3° вокруг четверти хорды.
    let xf = Airfoil::placement(200.0, 3.0, 0.25, Vec2::new(0.0, 0.0));
    let r = a.to_region(&xf).unwrap();
    r.check().unwrap();
    assert_eq!(r.rings().len(), 1);
    let area = r.area_mm2();
    assert!(
        (area - 0.0815 * 200.0 * 200.0).abs() < 0.01 * area,
        "площадь {area}"
    );
}

#[test]
fn selig_and_lednicer_agree() {
    let a = Airfoil::naca4_code("2412", 61, false).unwrap();
    let s = Airfoil::parse_dat(&a.to_selig(12)).unwrap();
    let l = Airfoil::parse_dat(&lednicer(&a, 12)).unwrap();
    same(&a, &s, 1e-11);
    same(&a, &l, 1e-11);
    assert_eq!(s.name(), "NACA 2412");
    assert_eq!(l.name(), "NACA 2412");
}

fn transformed_text(a: &Airfoil, xf: &Xform2) -> String {
    let mut text = format!("{} в мм\n", a.name());
    for p in a.outline(xf) {
        text.push_str(&format!("{:.10} {:.10}\n", p.x, p.y));
    }
    text
}

#[test]
fn transformed_file_normalizes_back() {
    // Симметричный профиль: носок однозначен, приведение возвращает исходник.
    let sym = Airfoil::naca4_code("0012", 81, true).unwrap();
    let xf = Xform2::scaling(250.0, 250.0)
        .then(&Xform2::rotation(7f64.to_radians()))
        .then(&Xform2::translation(Vec2::new(-30.0, 12.5)));
    same(
        &sym,
        &Airfoil::parse_dat(&transformed_text(&sym, &xf)).unwrap(),
        1e-9,
    );
    // С вогнутостью носок файла не помечен; приведение не зависит от того, как файл
    // повёрнут, растянут и сдвинут.
    let a = Airfoil::naca4_code("4415", 81, true).unwrap();
    let xf2 = Xform2::scaling(80.0, 80.0)
        .then(&Xform2::rotation(-12f64.to_radians()))
        .then(&Xform2::translation(Vec2::new(500.0, -3.0)));
    let b1 = Airfoil::parse_dat(&transformed_text(&a, &xf)).unwrap();
    let b2 = Airfoil::parse_dat(&transformed_text(&a, &xf2)).unwrap();
    for (p, q) in b1
        .outline(&Xform2::IDENTITY)
        .iter()
        .zip(b2.outline(&Xform2::IDENTITY))
    {
        assert!(p.distance(q) < 1e-9, "{p:?} / {q:?}");
    }
    // А от исходника отличается лишь чуть: носок по формуле NACA и самая далёкая
    // от кромки точка расходятся на доли процента хорды.
    same(&a, &b1, 0.01);
}

#[test]
fn formatting_noise_is_harmless() {
    let a = Airfoil::naca4_code("0012", 41, true).unwrap();
    let pts = a.outline(&Xform2::IDENTITY);
    let mut text = String::from("# комментарий\r\n  NACA 0012 с мусором  \r\n\r\n");
    for (k, p) in pts.iter().enumerate() {
        let line = match k % 4 {
            0 => format!("{:.10}\t{:.10}  ", p.x, p.y),
            1 => format!("{:.10}, {:.10}", p.x, p.y),
            2 => format!("  {:.10e} {:.10e}", p.x, p.y).replace('e', "D"),
            _ => format!("{:.10} {:.10}", p.x, p.y),
        };
        text.push_str(&line);
        text.push_str("\r\n");
        if k == pts.len() / 2 {
            // Повтор носка и комментарий посреди данных.
            text.push_str(&format!("{:.10} {:.10}\r\n# носок\r\n\r\n", p.x, p.y));
        }
    }
    // Замыкающая копия первой точки.
    text.push_str(&format!("{:.10} {:.10}\r\n", pts[0].x, pts[0].y));
    let b = Airfoil::parse_dat(&text).unwrap();
    assert_eq!(b.name(), "NACA 0012 с мусором");
    same(&a, &b, 1e-9);
}

#[test]
fn reversed_direction_is_detected() {
    let a = Airfoil::naca4_code("2412", 51, true).unwrap();
    let mut pts = a.outline(&Xform2::IDENTITY);
    pts.reverse(); // сначала низ
    let text: String = pts.iter().map(|p| format!("{} {}\n", p.x, p.y)).collect();
    let b = Airfoil::parse_dat(&text).unwrap();
    same(&a, &b, 1e-12);
    assert!(b.upper_y(0.3) > b.lower_y(0.3));
}

#[test]
fn errors_name_the_line() {
    assert_eq!(Airfoil::parse_dat(""), Err(AirfoilError::Empty));
    assert_eq!(
        Airfoil::parse_dat("Только название\n"),
        Err(AirfoilError::Empty)
    );
    assert_eq!(
        Airfoil::parse_dat("X\n1 0\n0.5 0.1\nоops\n0 0\n"),
        Err(AirfoilError::BadLine {
            line: 4,
            text: "оops".into()
        })
    );
    assert_eq!(
        Airfoil::parse_dat("X\n1 0\n0.5 NaN\n"),
        Err(AirfoilError::NotFinite { line: 3 })
    );
    assert_eq!(
        Airfoil::parse_dat("X\n1 0\n0.5 0.1 0.2\n"),
        Err(AirfoilError::BadLine {
            line: 3,
            text: "0.5 0.1 0.2".into()
        })
    );
    assert!(matches!(
        Airfoil::parse_dat("X\n1 0\n0 0\n1 0\n"),
        Err(AirfoilError::TooFewPoints(_))
    ));
    assert_eq!(
        Airfoil::parse_dat("X\n3. 3.\n0 0\n0.5 0.1\n1 0\n0 0\n1 0\n"),
        Err(AirfoilError::LednicerCounts {
            declared: 6,
            found: 5
        })
    );
    let v = |x, y| Vec2::new(x, y);
    assert_eq!(
        Airfoil::from_surfaces(
            "X",
            vec![v(0.0, 0.0), v(0.5, 0.1), v(1.0, 0.0)],
            vec![v(1.0, 0.0), v(0.5, -0.1), v(0.0, 0.0)]
        ),
        Err(AirfoilError::ZeroChord)
    );
}

#[test]
fn naca_area_matches_formula() {
    for (code, t) in [("0012", 0.12), ("0006", 0.06), ("0021", 0.21)] {
        let a = Airfoil::naca4_code(code, 401, true).unwrap();
        let yt = |x: f64| {
            5.0 * t
                * (0.2969 * x.sqrt() - 0.1260 * x - 0.3516 * x * x + 0.2843 * x.powi(3)
                    - 0.1036 * x.powi(4))
        };
        // √x у носка: интеграл по замене x = u², dx = 2u du.
        let exact = simpson(|u| 2.0 * yt(u * u) * 2.0 * u, 20_000);
        let chord = 1000.0;
        let r = a
            .to_region(&Airfoil::placement(chord, 0.0, 0.0, Vec2::ZERO))
            .unwrap();
        let got = r.area_mm2() / (chord * chord);
        assert!(
            (got - exact).abs() < 2e-5,
            "{code}: площадь {got}, интеграл {exact}"
        );
    }
}

#[test]
fn placement_keeps_chord_and_pivot() {
    let a = Airfoil::naca4_code("2412", 61, true).unwrap();
    let at = Vec2::new(100.0, -20.0);
    let xf = Airfoil::placement(180.0, 5.0, 0.25, at);
    let le = xf.apply_point(Vec2::ZERO);
    let te = xf.apply_point(Vec2::new(1.0, 0.0));
    assert!((le.distance(te) - 180.0).abs() < 1e-9);
    assert!(xf.apply_point(Vec2::new(0.25, 0.0)).distance(at) < 1e-12);
    assert!(le.y > te.y, "положительная крутка поднимает носок");
    // Носок профиля в контуре стоит ровно там, куда переносится точка (0, 0).
    assert!(a.outline(&xf).iter().any(|p| p.distance(le) < 1e-9));
}

#[test]
fn thick_trailing_edge_for_printing() {
    let a = Airfoil::naca4_code("2412", 81, true).unwrap();
    let chord = 150.0;
    let b = a.with_min_te_gap(0.6 / chord); // кромка 0,6 мм — полторы ширины сопла
    let r = b
        .to_region(&Airfoil::placement(chord, 0.0, 0.0, Vec2::ZERO))
        .unwrap();
    r.check().unwrap();
    let pts = b.outline(&Airfoil::placement(chord, 0.0, 0.0, Vec2::ZERO));
    let gap = pts[0].distance(pts[pts.len() - 1]);
    assert!((gap - 0.6).abs() < 1e-9, "кромка {gap} мм");
}

#[test]
fn random_naca_round_trips() {
    for_cases("random_naca_round_trips", 1000, |rng: &mut Rng| {
        let m = rng.range_i64(0, 9) as f64 / 100.0;
        let p = rng.range_i64(1, 9) as f64 / 10.0;
        let t = rng.range_i64(6, 24) as f64 / 100.0;
        let n = rng.range_i64(15, 160) as usize;
        let a = Airfoil::naca4(m, p, t, n, rng.chance(0.5));
        let decimals = rng.range_i64(6, 12) as usize;
        let text = if rng.chance(0.5) {
            a.to_selig(decimals)
        } else {
            lednicer(&a, decimals)
        };
        let b = Airfoil::parse_dat(&text).unwrap();
        // Файл в долях хорды берётся как есть: точки совпадают с точностью до
        // округления на decimals знаков (по y у носка сравнивать нельзя — там
        // поверхность почти вертикальна, и любое округление x меняет y сильно).
        same_outline(&a, &b, 10f64.powi(-(decimals as i32)));
        assert_eq!(b.name(), a.name());
        // Деталь: случайная хорда, крутка, место.
        let chord = 10f64.powf(rng.range_f64(1.0, 3.0));
        let at = Vec2::new(rng.range_f64(-500.0, 500.0), rng.range_f64(-500.0, 500.0));
        let xf = Airfoil::placement(
            chord,
            rng.range_f64(-10.0, 10.0),
            rng.range_f64(0.0, 1.0),
            at,
        );
        let r = b.to_region(&xf).unwrap();
        r.check().unwrap();
        let unit = b
            .to_region(&Airfoil::placement(1000.0, 0.0, 0.0, Vec2::ZERO))
            .unwrap();
        let ratio = r.area_mm2() / chord.powi(2) / (unit.area_mm2() / 1e6);
        // Сетка 1 мкм: при хорде 10 мм относительная погрешность площади до ~1e-3.
        assert!(
            (ratio - 1.0).abs() < 2e-3,
            "площадь не масштабируется: {ratio}"
        );
    });
}
