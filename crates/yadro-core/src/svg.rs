//! Отладочная картинка 2D в SVG: контуры, точки и подписи.
//!
//! Ось Y направлена вверх, как в CAD; в SVG она смотрит вниз, поэтому при
//! выводе координата Y меняет знак. Толщина линий, размер точек и шрифта
//! заданы в пикселях картинки и не зависят от размеров модели: и пластина
//! в метр, и щель в микрометр видны одинаково чётко.

use crate::grid::GridPoint;
use crate::vec::Vec2;
use std::fmt::Write as _;
use std::io;
use std::path::Path;

/// Цвета для типичных ролей в отладке.
pub mod palette {
    /// Первое тело (A).
    pub const A: &str = "#2563eb";
    /// Второе тело (B).
    pub const B: &str = "#16a34a";
    /// Результат операции.
    pub const RESULT: &str = "#111827";
    /// Ошибка или расхождение.
    pub const ERROR: &str = "#dc2626";
    /// Вспомогательные пометки.
    pub const NOTE: &str = "#6b7280";
}

/// Стиль линии и заливки.
#[derive(Clone, Debug, PartialEq)]
pub struct Style {
    /// Цвет линии.
    pub stroke: String,
    /// Толщина линии, пиксели картинки.
    pub width_px: f64,
    /// Цвет заливки; `None` — без заливки.
    pub fill: Option<String>,
    /// Непрозрачность заливки, от 0 до 1.
    pub fill_opacity: f64,
}

impl Style {
    /// Линия заданного цвета толщиной 1,5 пикселя.
    pub fn line(color: &str) -> Self {
        Self {
            stroke: color.into(),
            width_px: 1.5,
            fill: None,
            fill_opacity: 0.0,
        }
    }

    /// Линия с полупрозрачной заливкой того же цвета.
    pub fn filled(color: &str, opacity: f64) -> Self {
        Self {
            fill: Some(color.into()),
            fill_opacity: opacity,
            ..Self::line(color)
        }
    }

    /// Та же линия другой толщины.
    pub fn width(mut self, px: f64) -> Self {
        self.width_px = px;
        self
    }
}

enum Item {
    Path {
        points: Vec<Vec2>,
        closed: bool,
        style: Style,
    },
    Dot {
        at: Vec2,
        color: String,
    },
    Label {
        at: Vec2,
        text: String,
        color: String,
    },
}

/// Картинка для отладки: элементы набираются по одному, потом выводятся в SVG.
#[derive(Default)]
pub struct Svg {
    items: Vec<Item>,
}

fn escape(text: &str) -> String {
    text.replace('&', "&amp;")
        .replace('<', "&lt;")
        .replace('>', "&gt;")
        .replace('"', "&quot;")
}

impl Svg {
    /// Пустая картинка.
    pub fn new() -> Self {
        Self::default()
    }

    /// Ломаная.
    pub fn polyline(&mut self, points: &[Vec2], style: Style) -> &mut Self {
        self.items.push(Item::Path {
            points: points.to_vec(),
            closed: false,
            style,
        });
        self
    }

    /// Замкнутый контур.
    pub fn polygon(&mut self, points: &[Vec2], style: Style) -> &mut Self {
        self.items.push(Item::Path {
            points: points.to_vec(),
            closed: true,
            style,
        });
        self
    }

    /// Замкнутый контур на 2D-сетке.
    pub fn grid_polygon(&mut self, points: &[GridPoint], style: Style) -> &mut Self {
        let points: Vec<Vec2> = points.iter().map(|p| p.to_mm()).collect();
        self.polygon(&points, style)
    }

    /// Контур на 2D-сетке.
    pub fn ring(&mut self, ring: &crate::ring::Ring, style: Style) -> &mut Self {
        self.grid_polygon(ring.points(), style)
    }

    /// Область: внешние контуры стилем `outer`, отверстия стилем `hole`.
    pub fn region(
        &mut self,
        region: &crate::region::Region,
        outer: Style,
        hole: Style,
    ) -> &mut Self {
        for r in region.rings() {
            let s = if r.is_ccw() {
                outer.clone()
            } else {
                hole.clone()
            };
            self.ring(r, s);
        }
        self
    }

    /// Точка (кружок радиусом 3 пикселя).
    pub fn dot(&mut self, at: Vec2, color: &str) -> &mut Self {
        self.items.push(Item::Dot {
            at,
            color: color.into(),
        });
        self
    }

    /// Подпись рядом с точкой.
    pub fn label(&mut self, at: Vec2, text: &str, color: &str) -> &mut Self {
        self.items.push(Item::Label {
            at,
            text: text.into(),
            color: color.into(),
        });
        self
    }

    /// Габарит всех элементов; `None`, если картинка пуста.
    pub fn bounds(&self) -> Option<(Vec2, Vec2)> {
        self.items
            .iter()
            .flat_map(|item| match item {
                Item::Path { points, .. } => points.clone(),
                Item::Dot { at, .. } | Item::Label { at, .. } => vec![*at],
            })
            .fold(None, |acc, p| match acc {
                None => Some((p, p)),
                Some((lo, hi)) => Some((lo.min(p), hi.max(p))),
            })
    }

    /// SVG-текст шириной `width_px` пикселей.
    pub fn render(&self, width_px: f64) -> String {
        let (lo, hi) = self.bounds().unwrap_or((Vec2::ZERO, Vec2::new(1.0, 1.0)));
        let size = (hi.x - lo.x).max(hi.y - lo.y).max(1e-9);
        let margin = size * 0.05;
        let (left, top) = (lo.x - margin, -hi.y - margin);
        let (w, h) = (hi.x - lo.x + 2.0 * margin, hi.y - lo.y + 2.0 * margin);
        let unit = w / width_px; // миллиметров в одном пикселе
        let height_px = width_px * h / w;
        let mut s = String::new();
        let _ = writeln!(
            s,
            r#"<svg xmlns="http://www.w3.org/2000/svg" width="{width_px}" height="{height_px:.1}" viewBox="{left} {top} {w} {h}">"#
        );
        let _ = writeln!(
            s,
            r#"<rect x="{left}" y="{top}" width="{w}" height="{h}" fill="white"/>"#
        );
        for item in &self.items {
            match item {
                Item::Path {
                    points,
                    closed,
                    style,
                } => {
                    let pts: Vec<String> =
                        points.iter().map(|p| format!("{},{}", p.x, -p.y)).collect();
                    let tag = if *closed { "polygon" } else { "polyline" };
                    let fill = match &style.fill {
                        Some(c) => format!(r#"fill="{c}" fill-opacity="{}""#, style.fill_opacity),
                        None => r#"fill="none""#.to_string(),
                    };
                    let _ = writeln!(
                        s,
                        r#"<{tag} points="{}" stroke="{}" stroke-width="{}" vector-effect="non-scaling-stroke" stroke-linejoin="round" {fill}/>"#,
                        pts.join(" "),
                        style.stroke,
                        style.width_px
                    );
                }
                Item::Dot { at, color } => {
                    let _ = writeln!(
                        s,
                        r#"<circle cx="{}" cy="{}" r="{}" fill="{color}"/>"#,
                        at.x,
                        -at.y,
                        3.0 * unit
                    );
                }
                Item::Label { at, text, color } => {
                    let _ = writeln!(
                        s,
                        r#"<text x="{}" y="{}" font-size="{}" font-family="sans-serif" fill="{color}">{}</text>"#,
                        at.x + 4.0 * unit,
                        -at.y - 4.0 * unit,
                        12.0 * unit,
                        escape(text)
                    );
                }
            }
        }
        s.push_str("</svg>\n");
        s
    }

    /// Сохранить картинку в файл.
    pub fn save(&self, path: impl AsRef<Path>, width_px: f64) -> io::Result<()> {
        std::fs::write(path, self.render(width_px))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn y_axis_points_up() {
        let mut svg = Svg::new();
        svg.dot(Vec2::new(0.0, 10.0), palette::ERROR);
        assert!(svg.render(400.0).contains(r#"cy="-10""#));
    }

    #[test]
    fn view_box_wraps_content_with_margin() {
        let square = [
            Vec2::new(0.0, 0.0),
            Vec2::new(10.0, 0.0),
            Vec2::new(10.0, 10.0),
            Vec2::new(0.0, 10.0),
        ];
        let mut svg = Svg::new();
        svg.polygon(&square, Style::filled(palette::A, 0.2));
        let text = svg.render(400.0);
        assert!(text.contains(r#"viewBox="-0.5 -10.5 11 11""#), "{text}");
        assert!(
            text.contains(r#"<polygon points="0,-0 10,-0 10,-10 0,-10""#),
            "{text}"
        );
    }

    #[test]
    fn labels_are_escaped_and_grid_points_converted() {
        let mut svg = Svg::new();
        svg.label(Vec2::ZERO, r#"A<B & "C""#, palette::NOTE);
        let g = |x, y| GridPoint::new(x, y).unwrap();
        svg.grid_polygon(
            &[g(0, 0), g(1500, 0), g(0, 2000)],
            Style::line(palette::B).width(2.0),
        );
        let text = svg.render(300.0);
        assert!(text.contains("A&lt;B &amp; &quot;C&quot;"));
        assert!(text.contains("1.5,-0"));
        assert!(text.contains(r#"stroke-width="2""#));
        assert_eq!(svg.bounds(), Some((Vec2::ZERO, Vec2::new(1.5, 2.0))));
    }

    #[test]
    fn empty_picture_is_still_valid() {
        let text = Svg::new().render(100.0);
        assert!(text.starts_with("<svg"));
        assert!(text.ends_with("</svg>\n"));
    }
}
