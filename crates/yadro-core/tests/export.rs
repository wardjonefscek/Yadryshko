//! Ворота ступени 1, шаг 5: восстановление дуг, DXF и SVG.
//!
//! - Подогнанный контур не дальше допуска от ломаной и обратно (Хаусдорф),
//!   проверяется выборкой точек в обе стороны — независимо от устройства подгонки.
//! - Окружности, скругления, пазы и углы эквидистанты возвращаются дугами
//!   с верными радиусами; общие вершины контуров остаются на месте.
//! - DXF и SVG читаются обратно в то же самое (свой разбор SVG здесь, в тесте).
//! - Чтение DXF чужих программ: LWPOLYLINE, POLYLINE, LINE+ARC вразброс, CIRCLE,
//!   комментарии и лишние объекты; ошибки — с номером строки.

mod common;

use common::*;
use std::f64::consts::{PI, TAU};
use yadro_core::dxf::{self, DxfError};
use yadro_core::fit::{fit_region, fit_ring};
use yadro_core::testkit::{Rng, for_cases};
use yadro_core::{Arc, ContourOptions, DensityField, Path, Region, Ring, Segment, Sheet, Vec2};

// ───────────────────────────── Расстояния ─────────────────────────────

fn seg_dist(p: Vec2, a: Vec2, b: Vec2) -> f64 {
    let ab = b - a;
    let l2 = ab.length_sq();
    if l2 == 0.0 {
        return p.distance(a);
    }
    let t = ((p - a).dot(ab) / l2).clamp(0.0, 1.0);
    p.distance(a + ab * t)
}

fn arc_dist(p: Vec2, a: &Arc) -> f64 {
    let d = p - a.center;
    let rel = if a.sweep > 0.0 {
        (d.y.atan2(d.x) - a.start).rem_euclid(TAU)
    } else {
        (a.start - d.y.atan2(d.x)).rem_euclid(TAU)
    };
    if rel <= a.sweep.abs() {
        (d.length() - a.radius).abs()
    } else {
        p.distance(a.start_point()).min(p.distance(a.end_point()))
    }
}

/// Точки по контуру: на каждом участке `k` штук.
fn samples(path: &Path, k: usize) -> Vec<Vec2> {
    let mut out = Vec::new();
    for s in path.segments() {
        for i in 0..k {
            let t = i as f64 / k as f64;
            out.push(match s {
                Segment::Line(a, b) => a + (b - a) * t,
                Segment::Arc(a) => a.point(t),
            });
        }
    }
    out
}

/// Участки, разложенные по клеткам сетки: быстрый ответ «есть ли участок ближе
/// `tol`» для точки. Без этого проверка на контурах в тысячи вершин квадратична.
/// Участок кладётся во все клетки, до которых от него не дальше `tol`.
struct Near {
    cell: f64,
    segs: Vec<Segment>,
    grid: std::collections::HashMap<(i64, i64), Vec<usize>>,
}

impl Near {
    fn new(segs: Vec<Segment>, tol: f64) -> Self {
        let len = |s: &Segment| match s {
            Segment::Line(a, b) => a.distance(*b),
            Segment::Arc(a) => a.length(),
        };
        let mean = segs.iter().map(len).sum::<f64>() / segs.len().max(1) as f64;
        let cell = mean.max(tol * 4.0).max(1e-6);
        let mut grid: std::collections::HashMap<(i64, i64), Vec<usize>> = Default::default();
        let step = cell / 2.0;
        for (k, s) in segs.iter().enumerate() {
            // Точки вдоль участка через полклетки; вокруг каждой — запас tol + полшага.
            let n = (len(s) / step).ceil().max(1.0) as usize;
            let m = tol + step;
            let c = |v: f64| (v / cell).floor() as i64;
            let mut cells = std::collections::BTreeSet::new();
            for i in 0..=n {
                let t = i as f64 / n as f64;
                let q = match s {
                    Segment::Line(a, b) => *a + (*b - *a) * t,
                    Segment::Arc(a) => a.point(t),
                };
                for x in c(q.x - m)..=c(q.x + m) {
                    for y in c(q.y - m)..=c(q.y + m) {
                        cells.insert((x, y));
                    }
                }
            }
            for key in cells {
                grid.entry(key).or_default().push(k);
            }
        }
        Self { cell, segs, grid }
    }

    /// Расстояние до ближайшего участка, если он ближе допуска (иначе может быть
    /// бесконечность — значит, ближе допуска участков нет).
    fn dist(&self, p: Vec2) -> f64 {
        let key = (
            (p.x / self.cell).floor() as i64,
            (p.y / self.cell).floor() as i64,
        );
        self.grid.get(&key).map_or(f64::INFINITY, |ks| {
            ks.iter()
                .map(|&k| match self.segs[k] {
                    Segment::Line(a, b) => seg_dist(p, a, b),
                    Segment::Arc(a) => arc_dist(p, &a),
                })
                .fold(f64::INFINITY, f64::min)
        })
    }
}

/// Расстояние Хаусдорфа между ломаной и контуром (выборкой); больше допуска — «бесконечно».
fn hausdorff(pts: &[Vec2], path: &Path, tol: f64) -> f64 {
    let n = pts.len();
    let to_path = Near::new(path.segments().collect(), tol);
    let poly: Vec<Segment> = (0..n)
        .map(|i| Segment::Line(pts[i], pts[(i + 1) % n]))
        .collect();
    let to_poly = Near::new(poly, tol);
    let there = (0..n)
        .flat_map(|i| {
            let (a, b) = (pts[i], pts[(i + 1) % n]);
            (0..4).map(move |k| a + (b - a) * (k as f64 / 4.0))
        })
        .map(|q| to_path.dist(q))
        .fold(0.0, f64::max);
    let back = samples(path, 24)
        .into_iter()
        .map(|q| to_poly.dist(q))
        .fold(0.0, f64::max);
    there.max(back)
}

fn check_fit(ring: &Ring, path: &Path, tol: f64) {
    let pts = ring.to_mm();
    let h = hausdorff(&pts, path, tol);
    assert!(h <= tol * (1.0 + 1e-6) + 1e-9, "Хаусдорф {h} > {tol}");
    assert!(path.len() <= ring.len(), "вершин стало больше");
    // Вершины — из ломаной; исключение — целая окружность (две полуокружности).
    if yadro_core::sheet::is_circle(path).is_none() {
        let key = |p: Vec2| ((p.x * 1e6).round() as i64, (p.y * 1e6).round() as i64);
        let set: std::collections::HashSet<_> = pts.iter().map(|&p| key(p)).collect();
        for v in path.vertices() {
            assert!(set.contains(&key(v.at)), "вершина {:?} не из ломаной", v.at);
        }
    }
}

fn arcs(p: &Path) -> Vec<Arc> {
    p.segments()
        .filter_map(|s| match s {
            Segment::Arc(a) => Some(a),
            Segment::Line(..) => None,
        })
        .collect()
}

// ───────────────────────────── Восстановление дуг ─────────────────────────────

#[test]
fn fit_recovers_known_shapes() {
    let tol = 0.005;
    let cases: Vec<(Path, usize, usize, f64)> = vec![
        // (контур, участков, дуг, радиус)
        (Path::circle(Vec2::new(12.5, -4.0), 9.0), 2, 2, 9.0),
        (Path::circle(Vec2::new(0.0, 0.0), 0.3).reversed(), 2, 2, 0.3),
        (Path::rounded_rect(0.0, 0.0, 60.0, 25.0, 4.0), 8, 4, 4.0),
        (
            Path::slot(Vec2::ZERO, Vec2::new(30.0, 10.0), 3.0),
            4,
            2,
            3.0,
        ),
    ];
    for (path, len, n_arcs, r) in cases {
        let ring = path.to_ring(0.002).unwrap();
        let fitted = fit_ring(&ring, tol);
        check_fit(&ring, &fitted, tol);
        assert_eq!(fitted.len(), len, "{:?}", fitted.vertices());
        let a = arcs(&fitted);
        assert_eq!(a.len(), n_arcs);
        for arc in a {
            assert!(
                (arc.radius - r).abs() <= tol,
                "радиус {} вместо {r}",
                arc.radius
            );
        }
        assert!((fitted.area() - path.area()).abs() < path.length() * tol);
        assert_eq!(fitted.area() > 0.0, path.area() > 0.0);
    }
    // Короткая прямая между скруглениями (0,2 мм): подгонка может отдать её дуге
    // большого радиуса — это в допуске; участков всё равно не больше восьми, а
    // крупные дуги — радиуса скругления.
    let small = Path::rounded_rect(-5.0, 2.0, 5.0, 3.0, 0.4);
    let ring = small.to_ring(0.002).unwrap();
    let fitted = fit_ring(&ring, tol);
    check_fit(&ring, &fitted, tol);
    assert!(fitted.len() <= 8);
    for arc in arcs(&fitted).iter().filter(|a| a.sweep.abs() > 0.5) {
        assert!((arc.radius - 0.4).abs() <= tol, "радиус {}", arc.radius);
    }
    // Многоугольник без дуг остаётся многоугольником.
    let poly = poly(&[
        (0, 0),
        (7000, 1000),
        (9000, 6000),
        (3000, 9000),
        (-2000, 4000),
    ]);
    let f = fit_ring(&poly.rings()[0], tol);
    assert_eq!(f.len(), 5);
    assert!(arcs(&f).is_empty());
}

#[test]
fn fit_leaves_slivers_alone() {
    // Щепка толщиной 2 мкм при допуске 50 мкм: подгонка могла бы схлопнуть её в два
    // отрезка туда-обратно; остаётся ломаная с той же площадью.
    let sliver = poly(&[
        (0, 0),
        (5000, 1),
        (10_000, 0),
        (10_000, 2),
        (5000, 3),
        (0, 2),
    ]);
    let ring = &sliver.rings()[0];
    let p = fit_ring(ring, 0.05);
    assert_eq!(p.len(), ring.len());
    assert!((p.area() - ring.area_mm2()).abs() < 1e-12);
}

#[test]
fn fit_offset_corners_become_arcs() {
    let plate = rect(0, 0, 40_000, 20_000).difference(&rect(10_000, 5000, 30_000, 15_000));
    let grown = plate.offset(2.0).unwrap();
    let paths = fit_region(&grown, 0.005);
    assert_eq!(paths.len(), 2);
    let outer = paths.iter().find(|p| p.area() > 0.0).unwrap();
    let hole = paths.iter().find(|p| p.area() < 0.0).unwrap();
    // Снаружи: 4 стороны и 4 дуги радиуса 2; отверстие — прямоугольник с острыми углами.
    assert_eq!(
        (outer.len(), arcs(outer).len()),
        (8, 4),
        "{:?}",
        outer.vertices()
    );
    assert!(arcs(outer).iter().all(|a| (a.radius - 2.0).abs() < 0.005));
    assert_eq!((hole.len(), arcs(hole).len()), (4, 0));
    for (ring, path) in grown.rings().iter().zip(&paths) {
        check_fit(ring, path, 0.005);
    }
}

#[test]
fn fit_keeps_touching_vertices() {
    // Два квадрата касаются в вершине; отверстие касается внешнего контура вершиной.
    let touch = rect(0, 0, 10_000, 10_000).union(&rect(10_000, 10_000, 20_000, 20_000));
    let hole = poly(&[(5000, 0), (7000, 4000), (3000, 4000)]);
    let plate = rect(-30_000, -5000, -10_000, 15_000).difference(&poly(&[
        (-20_000, -5000),
        (-17_000, 0),
        (-23_000, 0),
    ]));
    for r in [touch, rect(0, 0, 10_000, 10_000).difference(&hole), plate] {
        let paths = fit_region(&r, 0.05);
        for (ring, path) in r.rings().iter().zip(&paths) {
            check_fit(ring, path, 0.05);
        }
        // Вершина касания (общая для контуров или дважды в одном) осталась столько
        // же раз.
        let all: Vec<Vec2> = r.rings().iter().flat_map(|g| g.to_mm()).collect();
        let fitted: Vec<Vec2> = paths
            .iter()
            .flat_map(|p| p.vertices().iter().map(|v| v.at))
            .collect();
        for p in &all {
            let times = all.iter().filter(|q| q.distance(*p) < 1e-9).count();
            if times > 1 {
                let kept = fitted.iter().filter(|q| q.distance(*p) < 1e-9).count();
                assert_eq!(kept, times, "вершина касания {p:?}");
            }
        }
    }
}

#[test]
fn fit_shrinks_density_contours() {
    let (c, rad) = (Vec2::new(15.0, 15.0), 10.0);
    let f = DensityField::from_fn(30, 30, 1.0, Vec2::ZERO, |i, j| {
        let p = Vec2::new(i as f64 + 0.5, j as f64 + 0.5);
        0.5 + (rad - p.distance(c)) / 6.0
    })
    .unwrap();
    let opts = ContourOptions {
        refine: 4,
        ..ContourOptions::default()
    };
    let region = f.contour(&opts).unwrap();
    let ring = &region.rings()[0];
    let tol = 0.02;
    let fitted = fit_ring(ring, tol);
    check_fit(ring, &fitted, tol);
    assert!(
        fitted.len() * 5 < ring.len(),
        "{} вершин из {}",
        fitted.len(),
        ring.len()
    );
}

// ───────────────────────────── DXF ─────────────────────────────

fn sample_sheet() -> Sheet {
    let mut s = Sheet::new();
    s.add_path("CONTOUR", Path::rounded_rect(0.0, 0.0, 120.0, 40.0, 5.0))
        .add_path("HOLES", Path::circle(Vec2::new(20.0, 20.0), 6.5))
        .add_path(
            "HOLES",
            Path::slot(Vec2::new(50.0, 15.0), Vec2::new(90.0, 25.0), 4.0).reversed(),
        )
        .add_path(
            "POCKETS",
            Path::polygon(&[
                Vec2::new(100.0, 5.0),
                Vec2::new(112.0, 5.0),
                Vec2::new(106.0, 30.0),
            ]),
        );
    s
}

/// Контуры совпадают: окружность (CIRCLE) — по центру и радиусу, остальное — по вершинам.
/// Контуры сравниваются по слоям: SVG пишет слой группой, так что порядок между
/// слоями может поменяться, а внутри слоя — нет.
fn same_sheet(a: &Sheet, b: &Sheet, tol: f64) {
    assert_eq!(a.items().len(), b.items().len());
    let mut sa: Vec<_> = a.items().iter().collect();
    let mut sb: Vec<_> = b.items().iter().collect();
    let order = |l: &str| a.layers().iter().position(|x| *x == l);
    sa.sort_by_key(|it| order(&it.layer));
    sb.sort_by_key(|it| order(&it.layer));
    for (x, y) in sa.into_iter().zip(sb) {
        assert_eq!(x.layer, y.layer);
        let (nx, ny) = (
            Near::new(x.path.segments().collect(), tol),
            Near::new(y.path.segments().collect(), tol),
        );
        let h = samples(&x.path, 32)
            .into_iter()
            .map(|q| ny.dist(q))
            .chain(samples(&y.path, 32).into_iter().map(|q| nx.dist(q)))
            .fold(0.0, f64::max);
        assert!(h <= tol, "слой {}: расхождение {h}", x.layer);
        // Площадь сдвигается не больше чем на расхождение × длину контура.
        let (ax, ay) = (x.path.area().abs(), y.path.area().abs());
        let allowed = 2.0 * tol * x.path.length() + 1e-9 * ax.max(1.0);
        assert!(
            (ax - ay).abs() <= allowed,
            "слой {}: площадь {ax} против {ay}, допуск {allowed}",
            x.layer
        );
    }
}

#[test]
fn dxf_round_trip() {
    let s = sample_sheet();
    let text = s.to_dxf();
    assert!(text.contains("AC1009") && text.contains("CIRCLE") && text.contains("POLYLINE"));
    let back = Sheet::from_dxf(&text).unwrap();
    same_sheet(&s, &back, 2e-6);
    assert_eq!(back.layers(), vec!["CONTOUR", "HOLES", "POCKETS"]);
    // CRLF и пробелы вокруг кодов не мешают.
    let crlf = text
        .replace('\n', "\r\n")
        .replace("\r\n10\r\n", "\r\n 10 \r\n");
    same_sheet(&s, &Sheet::from_dxf(&crlf).unwrap(), 2e-6);
}

#[test]
fn dxf_from_other_programs() {
    // LWPOLYLINE (R2000) с выпуклостью: прямоугольник 10 × 4 со скруглённым правым торцом.
    let lw = "0\nSECTION\n2\nENTITIES\n0\nLWPOLYLINE\n5\n2F\n100\nAcDbEntity\n8\nPLATE\n100\nAcDbPolyline\n90\n4\n70\n1\n\
              10\n0\n20\n0\n10\n10\n20\n0\n42\n1\n10\n10\n20\n4\n10\n0\n20\n4\n0\nENDSEC\n0\nEOF\n";
    let s = Sheet::from_dxf(lw).unwrap();
    assert_eq!(s.items().len(), 1);
    assert_eq!(s.items()[0].layer, "PLATE");
    assert!((s.items()[0].path.area() - (40.0 + PI * 2.0)).abs() < 1e-9);
    // LINE и ARC вразброс и в разные стороны — паз 20 × 6 мм; плюс TEXT и комментарий.
    let parts = "999\nэкспорт из другой программы\n0\nSECTION\n2\nENTITIES\n\
                 0\nARC\n8\nS\n10\n20\n20\n3\n40\n3\n50\n270\n51\n90\n\
                 0\nTEXT\n8\nS\n10\n0\n20\n0\n40\n2\n1\nподпись\n\
                 0\nLINE\n8\nS\n10\n0\n20\n6\n11\n20\n21\n6\n\
                 0\nLINE\n8\nS\n10\n20\n20\n0\n11\n0\n21\n0\n\
                 0\nARC\n8\nS\n10\n0\n20\n3\n40\n3\n50\n90\n51\n270\n0\nENDSEC\n0\nEOF\n";
    let s = Sheet::from_dxf(parts).unwrap();
    assert_eq!(s.items().len(), 1, "4 куска — один контур");
    assert!((s.items()[0].path.area().abs() - (120.0 + PI * 9.0)).abs() < 1e-9);
    // POLYLINE без флага замкнутости, но концы совпадают; CIRCLE; файл без секций.
    let r12 = "0\nPOLYLINE\n8\nA\n66\n1\n70\n0\n0\nVERTEX\n10\n0\n20\n0\n0\nVERTEX\n10\n5\n20\n0\n\
               0\nVERTEX\n10\n0\n20\n5\n0\nVERTEX\n10\n0\n20\n0\n0\nSEQEND\n0\nCIRCLE\n8\nB\n10\n1\n20\n1\n40\n0.5\n0\nEOF\n";
    let s = Sheet::from_dxf(r12).unwrap();
    assert_eq!(s.items().len(), 2);
    assert!((s.items()[0].path.area() - 12.5).abs() < 1e-12);
    assert!((s.items()[1].path.area() - PI * 0.25).abs() < 1e-12);
    // Область по чёт-нечет: окружность внутри треугольника — дыра, а не двойная заливка.
    let r = s.to_region(0.002);
    r.check().unwrap();
    assert!((r.area_mm2() - (12.5 - PI * 0.25)).abs() < 0.01);
}

#[test]
fn dxf_errors_are_named() {
    assert_eq!(
        Sheet::from_dxf("0\nLINE\nоops\n1\n"),
        Err(DxfError::BadCode { line: 3 })
    );
    assert_eq!(
        Sheet::from_dxf("0\nCIRCLE\n10\nабв\n20\n0\n40\n1\n"),
        Err(DxfError::BadNumber { line: 4 })
    );
    assert_eq!(Sheet::from_dxf("0\nLINE\n10\n"), Err(DxfError::Truncated));
    let open =
        "0\nLINE\n8\nL\n10\n0\n20\n0\n11\n5\n21\n0\n0\nLINE\n8\nL\n10\n5\n20\n0\n11\n5\n21\n5\n";
    assert!(matches!(
        Sheet::from_dxf(open),
        Err(DxfError::Unclosed { ref layer, .. }) if layer == "L"
    ));
    // Зазор меньше допуска сцепки — замкнуто; больше — нет.
    let gap = |g: f64| {
        format!(
            "0\nLINE\n10\n0\n20\n0\n11\n5\n21\n0\n0\nLINE\n10\n5\n20\n0\n11\n0\n21\n5\n\
             0\nLINE\n10\n0\n20\n5\n11\n{g}\n21\n0\n"
        )
    };
    assert!(Sheet::from_dxf(&gap(0.0005)).is_ok());
    assert!(dxf::read(&gap(0.0005), 0.0001).is_err());
}

// ───────────────────────────── SVG ─────────────────────────────

/// Разбор нашего SVG: размеры, слои (id групп) и пути `M L A Z` обратно в контуры.
fn parse_svg(text: &str) -> (f64, f64, [f64; 4], Sheet) {
    let attr = |s: &str, name: &str| -> String {
        let k = s.find(&format!("{name}=\"")).expect(name) + name.len() + 2;
        s[k..k + s[k..].find('"').unwrap()].to_string()
    };
    let svg_tag = &text[text.find("<svg").unwrap()..];
    let w: f64 = attr(svg_tag, "width")
        .trim_end_matches("mm")
        .parse()
        .unwrap();
    let h: f64 = attr(svg_tag, "height")
        .trim_end_matches("mm")
        .parse()
        .unwrap();
    let vb: Vec<f64> = attr(svg_tag, "viewBox")
        .split_whitespace()
        .map(|v| v.parse().unwrap())
        .collect();
    let mut sheet = Sheet::new();
    let mut layer = String::new();
    for line in text.lines() {
        if line.starts_with("<g ") {
            layer = attr(line, "id");
        }
        if !line.starts_with("<path ") {
            continue;
        }
        let d = attr(line, "d");
        let toks: Vec<String> = d.split_whitespace().map(String::from).collect();
        let mut k = 0;
        let mut verts: Vec<(Vec2, f64)> = Vec::new();
        let num = |s: &str| -> f64 { s.trim_start_matches(['M', 'L', 'A']).parse().unwrap() };
        let mut cur = Vec2::ZERO;
        while k < toks.len() {
            let t = &toks[k];
            if t.starts_with('M') || t.starts_with('L') {
                cur = Vec2::new(num(t), -num(&toks[k + 1]));
                verts.push((cur, 0.0));
                k += 2;
            } else if t.starts_with('A') {
                let r = num(t);
                let large = toks[k + 3] == "1";
                let sweep_cw = toks[k + 4] == "1";
                let end = Vec2::new(num(&toks[k + 5]), -num(&toks[k + 6]));
                let half = (cur.distance(end) / (2.0 * r)).min(1.0).asin();
                let theta = if large { TAU - 2.0 * half } else { 2.0 * half };
                let theta = if sweep_cw { -theta } else { theta };
                verts.last_mut().unwrap().1 = (theta / 4.0).tan();
                cur = end;
                verts.push((cur, 0.0));
                k += 7;
            } else {
                assert_eq!(t, "Z");
                k += 1;
            }
        }
        // Последний участок вернулся в начало: точка-повтор не нужна.
        let (last, _) = verts.pop().unwrap();
        assert!(last.distance(verts[0].0) < 1e-9, "путь не замкнут");
        let mut p = Path::new();
        for (at, bulge) in verts {
            p.push(at, bulge);
        }
        sheet.add_path(&layer, p);
    }
    (w, h, [vb[0], vb[1], vb[2], vb[3]], sheet)
}

#[test]
fn svg_round_trip_in_mm() {
    let s = sample_sheet();
    let text = s.to_svg();
    let (w, h, vb, back) = parse_svg(&text);
    let (lo, hi) = s.bounds().unwrap();
    assert!((w - (hi.x - lo.x + 2.0)).abs() < 1e-6 && (h - (hi.y - lo.y + 2.0)).abs() < 1e-6);
    assert_eq!(vb[2], w, "1 единица рисунка = 1 мм");
    assert!((vb[0] - (lo.x - 1.0)).abs() < 1e-6 && (vb[1] - (-hi.y - 1.0)).abs() < 1e-6);
    same_sheet(&s, &back, 5e-4);
}

// ───────────────────────────── От области до файла ─────────────────────────────

#[test]
fn region_to_dxf_and_back() {
    // Пластина с отверстиями-окружностями и пазом — через булевы операции.
    let tol = 0.002;
    let plate = Path::rounded_rect(0.0, 0.0, 80.0, 30.0, 3.0)
        .to_region(tol)
        .unwrap();
    let holes = Path::circle(Vec2::new(15.0, 15.0), 5.0)
        .to_region(tol)
        .unwrap()
        .union(
            &Path::slot(Vec2::new(35.0, 15.0), Vec2::new(65.0, 15.0), 4.0)
                .to_region(tol)
                .unwrap(),
        );
    let part = plate.difference(&holes);
    let mut sheet = Sheet::new();
    sheet.add_region("CONTOUR", "HOLES", &part, 0.005);
    assert_eq!(sheet.layers(), vec!["CONTOUR", "HOLES"]);
    // Окружность вернулась окружностью — в DXF это CIRCLE.
    let text = sheet.to_dxf();
    assert_eq!(text.matches("\nCIRCLE\n").count(), 1, "{text}");
    let back = Sheet::from_dxf(&text).unwrap();
    let r = back.to_region(tol);
    r.check().unwrap();
    let diff = r.xor(&part).area_mm2();
    let per: f64 = part
        .edges()
        .map(|(a, b)| a.to_mm().distance(b.to_mm()))
        .sum();
    assert!(
        diff < per * 0.012,
        "расхождение площади {diff} мм² при периметре {per}"
    );
}

#[test]
fn random_regions_export_within_tolerance() {
    for_cases(
        "random_regions_export_within_tolerance",
        150,
        |rng: &mut Rng| {
            let scale = *rng.pick(&[1_000, 50_000, 2_000_000]);
            let base = random_region(rng, scale);
            // С дугами: часть случаев — эквидистанта со скруглёнными углами.
            let r = if rng.chance(0.5) {
                base.offset(scale as f64 * 1e-3 * rng.range_f64(0.02, 0.3))
                    .unwrap()
            } else {
                base
            };
            // Допуск подгонки не меньше допуска хорды (2 мкм) с округлением — иначе дуги
            // не узнаются, и контур остаётся ломаной (это тоже проверяется).
            let tol = *rng.pick(&[0.002, 0.005, 0.02, 0.05]);
            let paths = fit_region(&r, tol);
            assert_eq!(paths.len(), r.rings().len());
            for (ring, path) in r.rings().iter().zip(&paths) {
                check_fit(ring, path, tol);
                assert_eq!(path.area() > 0.0, ring.is_ccw(), "ориентация сохранилась");
            }
            let mut sheet = Sheet::new();
            for p in &paths {
                sheet.add_path(if p.area() > 0.0 { "OUT" } else { "IN" }, p.clone());
            }
            // DXF: 6 знаков (1 нм) у координат, центра и радиуса окружности — до 2 нм.
            same_sheet(&sheet, &Sheet::from_dxf(&sheet.to_dxf()).unwrap(), 2e-6);
            same_sheet(&sheet, &parse_svg(&sheet.to_svg()).3, 5e-4);
        },
    );
}

#[test]
fn shapes_survive_regions() {
    // Region → данные → Region: то же самое (данные карточек и каталога).
    let r = rect(0, 0, 10_000, 10_000).difference(&rect(2000, 2000, 4000, 4000));
    let s = yadro_core::card::shape_of(&r);
    assert!(yadro_core::card::region_of(&s).xor(&r).is_empty());
    let _: Region = yadro_core::card::region_of(&vec![vec![(0, 0), (1, 1), (2, 2)]]);
}
