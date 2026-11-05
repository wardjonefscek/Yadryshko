//! Лист для изготовления: контуры по слоям → DXF (R12) и SVG в мм.
//!
//! Контуры хранятся как [`Path`] — отрезки и настоящие дуги. Область ядра
//! (ломаная) попадает на лист через восстановление дуг ([`crate::fit`]).
//!
//! DXF — версия R12 (AC1009): её читают все CAM, лазеры и станки. Контур —
//! замкнутая POLYLINE с выпуклостью (bulge) у вершин, полная окружность — CIRCLE.
//! SVG — в мм (ширина и высота с единицами, 1 единица рисунка = 1 мм), ось Y
//! перевёрнута, как принято в SVG, дуги — командой `A` кусками не больше четверти
//! окружности (так их форма не зависит от округления чисел в файле). Такой файл без
//! масштабирования открывают программы для лазера и плоттера.

use crate::dxf;
use crate::fit;
use crate::path::{Path, Segment};
use crate::region::Region;
use crate::vec::Vec2;

/// Контур на листе и его слой.
#[derive(Clone, Debug, PartialEq)]
pub struct SheetItem {
    /// Слой (в DXF — LAYER, в SVG — группа с этим id).
    pub layer: String,
    /// Контур.
    pub path: Path,
}

/// Подпись на листе (в CAM слой подписей обычно не режут).
#[derive(Clone, Debug, PartialEq)]
pub struct SheetText {
    /// Слой.
    pub layer: String,
    /// Левый нижний угол первой буквы, мм.
    pub at: Vec2,
    /// Высота букв, мм.
    pub height: f64,
    /// Текст.
    pub text: String,
}

/// Лист: контуры и подписи по слоям.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Sheet {
    items: Vec<SheetItem>,
    texts: Vec<SheetText>,
}

/// Число в файл: 6 знаков после запятой (нанометры), без хвостовых нулей.
pub(crate) fn num(x: f64) -> String {
    let s = format!("{x:.6}");
    let s = s.trim_end_matches('0').trim_end_matches('.');
    if s == "-0" || s.is_empty() {
        "0".to_string()
    } else {
        s.to_string()
    }
}

/// Цвета слоёв в SVG: тот же порядок, что номера цветов в DXF.
const SVG_COLORS: [&str; 6] = [
    "#000000", "#d62728", "#1f77b4", "#2ca02c", "#9467bd", "#ff7f0e",
];

/// Полная ли это окружность: две вершины, обе дуги — полуокружности в одну сторону.
pub fn is_circle(p: &Path) -> Option<(Vec2, f64)> {
    let v = p.vertices();
    let semi = |b: f64| (b.abs() - 1.0).abs() < 1e-12;
    if v.len() == 2 && semi(v[0].bulge) && semi(v[1].bulge) && v[0].bulge * v[1].bulge > 0.0 {
        Some(((v[0].at + v[1].at) * 0.5, v[0].at.distance(v[1].at) * 0.5))
    } else {
        None
    }
}

impl Sheet {
    /// Пустой лист.
    pub fn new() -> Self {
        Self::default()
    }

    /// Добавить контур в слой.
    pub fn add_path(&mut self, layer: &str, path: Path) -> &mut Self {
        if !path.is_empty() {
            self.items.push(SheetItem {
                layer: layer.to_string(),
                path,
            });
        }
        self
    }

    /// Добавить область: дуги восстанавливаются с допуском `fit_tol_mm`; внешние
    /// контуры — в слой `outer`, отверстия — в `holes` (можно в тот же).
    pub fn add_region(
        &mut self,
        outer: &str,
        holes: &str,
        region: &Region,
        fit_tol_mm: f64,
    ) -> &mut Self {
        for path in fit::fit_region(region, fit_tol_mm) {
            let layer = if path.area() >= 0.0 { outer } else { holes };
            self.add_path(layer, path);
        }
        self
    }

    /// Добавить подпись.
    pub fn add_text(&mut self, layer: &str, at: Vec2, height: f64, text: &str) -> &mut Self {
        if !text.is_empty() && height > 0.0 {
            self.texts.push(SheetText {
                layer: layer.to_string(),
                at,
                height,
                text: text.to_string(),
            });
        }
        self
    }

    /// Контуры.
    pub fn items(&self) -> &[SheetItem] {
        &self.items
    }

    /// Подписи.
    pub fn texts(&self) -> &[SheetText] {
        &self.texts
    }

    /// Слои в порядке появления: сначала слои контуров, потом слои подписей.
    pub fn layers(&self) -> Vec<&str> {
        let mut out: Vec<&str> = Vec::new();
        let names = self
            .items
            .iter()
            .map(|it| it.layer.as_str())
            .chain(self.texts.iter().map(|t| t.layer.as_str()));
        for name in names {
            if !out.contains(&name) {
                out.push(name);
            }
        }
        out
    }

    /// Габарит по настоящим дугам: (левый нижний, правый верхний), мм.
    pub fn bounds(&self) -> Option<(Vec2, Vec2)> {
        let mut lo = Vec2::new(f64::INFINITY, f64::INFINITY);
        let mut hi = Vec2::new(f64::NEG_INFINITY, f64::NEG_INFINITY);
        let mut add = |p: Vec2| {
            lo = Vec2::new(lo.x.min(p.x), lo.y.min(p.y));
            hi = Vec2::new(hi.x.max(p.x), hi.y.max(p.y));
        };
        for it in &self.items {
            for s in it.path.segments() {
                match s {
                    Segment::Line(a, _) => add(a),
                    Segment::Arc(arc) => {
                        add(arc.start_point());
                        add(arc.end_point());
                        // Крайние точки окружности, попавшие в дугу.
                        for k in 0..4 {
                            let ang = k as f64 * std::f64::consts::FRAC_PI_2;
                            let rel = if arc.sweep > 0.0 {
                                (ang - arc.start).rem_euclid(std::f64::consts::TAU)
                            } else {
                                (arc.start - ang).rem_euclid(std::f64::consts::TAU)
                            };
                            if rel <= arc.sweep.abs() {
                                add(arc.center + Vec2::new(ang.cos(), ang.sin()) * arc.radius);
                            }
                        }
                    }
                }
            }
        }
        // Подпись — примерно: ширина буквы 0,6 высоты.
        for tx in &self.texts {
            add(tx.at);
            let w = 0.6 * tx.height * tx.text.chars().count() as f64;
            add(tx.at + Vec2::new(w, tx.height));
        }
        lo.x.is_finite().then_some((lo, hi))
    }

    /// Текст DXF R12.
    pub fn to_dxf(&self) -> String {
        dxf::write(self)
    }

    /// Прочитать DXF (см. [`crate::dxf::read`]).
    pub fn from_dxf(text: &str) -> Result<Self, dxf::DxfError> {
        dxf::read(text, dxf::JOIN_TOL_MM)
    }

    /// Область, ограниченная контурами листа, по правилу чёт-нечет (как понимают
    /// вложенные контуры CAM): дуги — ломаной с допуском `tol_mm`.
    pub fn to_region(&self, tol_mm: f64) -> Region {
        let mut out = Region::empty();
        for it in &self.items {
            if let Ok(r) = it.path.to_region(tol_mm) {
                out = out.xor(&r);
            }
        }
        out
    }

    /// Текст SVG в мм.
    pub fn to_svg(&self) -> String {
        let (lo, hi) = self.bounds().unwrap_or((Vec2::ZERO, Vec2::new(1.0, 1.0)));
        let m = 1.0; // поле 1 мм вокруг
        let (w, h) = (hi.x - lo.x + 2.0 * m, hi.y - lo.y + 2.0 * m);
        let mut s = String::from("<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n");
        s.push_str(&format!(
            "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"{}mm\" height=\"{}mm\" viewBox=\"{} {} {} {}\">\n",
            num(w),
            num(h),
            num(lo.x - m),
            num(-hi.y - m),
            num(w),
            num(h)
        ));
        for (k, layer) in self.layers().iter().enumerate() {
            s.push_str(&format!(
                "<g id=\"{}\" fill=\"none\" stroke=\"{}\" stroke-width=\"0.1\">\n",
                escape(layer),
                SVG_COLORS[k % SVG_COLORS.len()]
            ));
            for it in self.items.iter().filter(|it| it.layer == *layer) {
                s.push_str(&format!("<path d=\"{}\"/>\n", svg_d(&it.path)));
            }
            for tx in self.texts.iter().filter(|tx| tx.layer == *layer) {
                s.push_str(&format!(
                    "<text x=\"{}\" y=\"{}\" font-size=\"{}\" font-family=\"sans-serif\" fill=\"{}\" stroke=\"none\">{}</text>\n",
                    num(tx.at.x),
                    num(-tx.at.y),
                    num(tx.height),
                    SVG_COLORS[k % SVG_COLORS.len()],
                    escape(&tx.text)
                ));
            }
            s.push_str("</g>\n");
        }
        s.push_str("</svg>\n");
        s
    }
}

/// Данные пути SVG: Y перевёрнута; дуга против часовой стрелки в мм — по часовой
/// на экране, то есть флаг sweep = 0.
fn svg_d(p: &Path) -> String {
    let v = p.vertices();
    let pt = |q: Vec2| format!("{} {}", num(q.x), num(-q.y));
    let mut d = format!("M{}", pt(v[0].at));
    for (k, s) in p.segments().enumerate() {
        let end = v[(k + 1) % v.len()].at; // конец — ровно следующая вершина
        match s {
            Segment::Line(..) => d.push_str(&format!(" L{}", pt(end))),
            Segment::Arc(arc) => {
                // В SVG дуга задаётся концами и радиусом; около полуокружности середина
                // очень чувствительна к округлению (сдвиг ≈ √(2·R·δ)). Поэтому дуга
                // пишется кусками не больше четверти окружности, а радиус — с 9 знаками.
                let n = (arc.sweep.abs() / std::f64::consts::FRAC_PI_2 - 1e-9)
                    .ceil()
                    .max(1.0) as usize;
                let sweep = u8::from(arc.sweep < 0.0);
                let r = format!("{:.9}", arc.radius);
                for i in 1..=n {
                    let q = if i == n {
                        end
                    } else {
                        arc.point(i as f64 / n as f64)
                    };
                    d.push_str(&format!(" A{r} {r} 0 0 {sweep} {}", pt(q)));
                }
            }
        }
    }
    d.push_str(" Z");
    d
}

/// Экранировать текст для атрибута XML.
fn escape(s: &str) -> String {
    s.replace('&', "&amp;")
        .replace('<', "&lt;")
        .replace('>', "&gt;")
        .replace('"', "&quot;")
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn numbers_are_short() {
        assert_eq!(num(1.0), "1");
        assert_eq!(num(-0.0000001), "0");
        assert_eq!(num(12.345), "12.345");
        assert_eq!(num(-3.1234567), "-3.123457");
    }

    #[test]
    fn circle_is_detected() {
        let c = Path::circle(Vec2::new(1.0, 2.0), 3.0);
        assert_eq!(is_circle(&c), Some((Vec2::new(1.0, 2.0), 3.0)));
        assert_eq!(
            is_circle(&Path::rounded_rect(0.0, 0.0, 1.0, 1.0, 0.5)),
            None
        );
    }

    #[test]
    fn texts_go_to_both_files() {
        let mut s = Sheet::new();
        s.add_path("CONTOUR", Path::circle(Vec2::ZERO, 5.0))
            .add_text("LABELS", Vec2::new(-5.0, -8.0), 2.0, "Нервюра 1 <z = 0>");
        assert_eq!(s.layers(), vec!["CONTOUR", "LABELS"]);
        let dxf = s.to_dxf();
        assert!(dxf.contains("\nTEXT\n8\nLABELS\n"));
        assert!(dxf.contains("\\U+041D\\U+0435"), "кириллица — через \\U+");
        let svg = s.to_svg();
        assert!(svg.contains(">Нервюра 1 &lt;z = 0&gt;</text>"));
        // Чтение DXF подписи пропускает.
        let back = Sheet::from_dxf(&dxf).unwrap();
        assert_eq!(back.items().len(), 1);
        let (lo, _) = s.bounds().unwrap();
        assert!((lo.y + 8.0).abs() < 1e-12, "подпись входит в габарит");
    }

    #[test]
    fn bounds_follow_arcs() {
        let mut s = Sheet::new();
        s.add_path("A", Path::circle(Vec2::new(10.0, 0.0), 2.0));
        let (lo, hi) = s.bounds().unwrap();
        assert!(lo.distance(Vec2::new(8.0, -2.0)) < 1e-12);
        assert!(hi.distance(Vec2::new(12.0, 2.0)) < 1e-12);
    }
}
