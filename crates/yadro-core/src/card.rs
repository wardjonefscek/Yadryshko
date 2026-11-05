//! Карточка ошибки: что упало, как повторить, картинка, данные для каталога.
//!
//! Когда случайная проверка падает, пример сначала ужимается ([`shrink_shapes`]:
//! убираются контуры и вершины, координаты становятся короче, пока проверка
//! продолжает падать), а потом сохраняется карточка — страница HTML:
//! - текст падения и команда, которая повторяет ровно этот случай;
//! - картинка исходных данных и результата (ужатого и, для сравнения, исходного);
//! - данные ужатого примера кодом на Rust — их можно сразу вставить в каталог
//!   вырожденных случаев (правило журнала: найденная ошибка сначала становится тестом).
//!
//! Карточки кладутся в `target/cards` (или в папку из `YADRO_CARDS`).

use crate::grid::GridPoint;
use crate::region::Region;
use crate::ring::Ring;
use crate::svg::{Style, Svg};
use crate::vec::Vec2;
use std::path::PathBuf;

/// Набор контуров на сетке: списки вершин (x, y) в микрометрах.
pub type Shape = Vec<Vec<(i64, i64)>>;

/// Контуры области как данные.
pub fn shape_of(region: &Region) -> Shape {
    region
        .rings()
        .iter()
        .map(|r| r.points().iter().map(|p| (p.x(), p.y())).collect())
        .collect()
}

/// Область из данных по правилу ненулевой обмотки; негодные контуры пропускаются.
pub fn region_of(shape: &Shape) -> Region {
    let rings: Vec<Ring> = shape
        .iter()
        .filter_map(|r| {
            let pts: Option<Vec<GridPoint>> =
                r.iter().map(|&(x, y)| GridPoint::new(x, y).ok()).collect();
            Ring::new(pts?).ok()
        })
        .collect();
    Region::from_rings_nonzero(&rings)
}

/// Сколько всего вершин.
pub fn vertex_count(shape: &Shape) -> usize {
    shape.iter().map(Vec::len).sum()
}

/// Варианты поменьше: без контура, без вершины, с координатами вдвое меньше,
/// сдвинутые к нулю, с округлёнными координатами.
pub fn smaller_shapes(s: &Shape) -> Vec<Shape> {
    let mut out = Vec::new();
    for k in 0..s.len() {
        let mut c = s.clone();
        c.remove(k);
        out.push(c);
    }
    for k in 0..s.len() {
        if s[k].len() > 3 {
            for v in 0..s[k].len() {
                let mut c = s.clone();
                c[k].remove(v);
                out.push(c);
            }
        }
    }
    let max = s
        .iter()
        .flatten()
        .map(|&(x, y)| x.abs().max(y.abs()))
        .max()
        .unwrap_or(0);
    if max > 1 {
        out.push(map(s, |x| x.div_euclid(2)));
    }
    let (mx, my) = (
        s.iter().flatten().map(|p| p.0).min().unwrap_or(0),
        s.iter().flatten().map(|p| p.1).min().unwrap_or(0),
    );
    if mx != 0 || my != 0 {
        out.push(
            s.iter()
                .map(|r| r.iter().map(|&(x, y)| (x - mx, y - my)).collect())
                .collect(),
        );
    }
    // Короче числа: каждую координату — к ближайшему круглому значению.
    for k in 0..s.len() {
        for v in 0..s[k].len() {
            for axis in 0..2 {
                let val = if axis == 0 { s[k][v].0 } else { s[k][v].1 };
                if let Some(r) = rounder(val) {
                    let mut c = s.clone();
                    if axis == 0 {
                        c[k][v].0 = r;
                    } else {
                        c[k][v].1 = r;
                    }
                    out.push(c);
                }
            }
        }
    }
    out
}

/// Ближайшее к `v` число с меньшим числом значащих цифр, если оно другое.
fn rounder(v: i64) -> Option<i64> {
    if v == 0 {
        return None;
    }
    let mut p = 10i64;
    while p <= v.abs() {
        let r = (v as f64 / p as f64).round() as i64 * p;
        if r != v {
            return Some(r);
        }
        p = p.saturating_mul(10);
    }
    Some(0)
}

fn map(s: &Shape, f: impl Fn(i64) -> i64) -> Shape {
    s.iter()
        .map(|r| r.iter().map(|&(x, y)| (f(x), f(y))).collect())
        .collect()
}

/// Варианты поменьше для пары: каждую по отдельности и обе сразу (вдвое меньше,
/// сдвиг к нулю), чтобы взаимное расположение сохранялось.
pub fn smaller_pairs(p: &(Shape, Shape)) -> Vec<(Shape, Shape)> {
    let (a, b) = p;
    let mut out = Vec::new();
    let max = a
        .iter()
        .chain(b)
        .flatten()
        .map(|&(x, y)| x.abs().max(y.abs()))
        .max()
        .unwrap_or(0);
    if max > 1 {
        out.push((map(a, |x| x.div_euclid(2)), map(b, |x| x.div_euclid(2))));
    }
    let all = a.iter().chain(b).flatten();
    let (mx, my) = all.fold((i64::MAX, i64::MAX), |m, p| (m.0.min(p.0), m.1.min(p.1)));
    if mx != i64::MAX && (mx != 0 || my != 0) {
        let t = |s: &Shape| -> Shape {
            s.iter()
                .map(|r| r.iter().map(|&(x, y)| (x - mx, y - my)).collect())
                .collect()
        };
        out.push((t(a), t(b)));
    }
    out.extend(smaller_shapes(a).into_iter().map(|c| (c, b.clone())));
    out.extend(smaller_shapes(b).into_iter().map(|c| (a.clone(), c)));
    out
}

/// Ужать пару фигур, пока `fails` остаётся верным (см. [`crate::testkit::shrink`]).
pub fn shrink_pair(
    start: (Shape, Shape),
    fails: impl FnMut(&(Shape, Shape)) -> bool,
) -> (Shape, Shape) {
    crate::testkit::shrink(start, smaller_pairs, fails, 4000)
}

/// Ужать одну фигуру.
pub fn shrink_shapes(start: Shape, fails: impl FnMut(&Shape) -> bool) -> Shape {
    crate::testkit::shrink(start, smaller_shapes, fails, 4000)
}

/// Фигура кодом на Rust для каталога.
pub fn shape_code(name: &str, s: &Shape) -> String {
    let rings: Vec<String> = s
        .iter()
        .map(|r| {
            let pts: Vec<String> = r.iter().map(|(x, y)| format!("({x}, {y})")).collect();
            format!("        &[{}],", pts.join(", "))
        })
        .collect();
    format!(
        "    let {name}: &[&[(i64, i64)]] = &[\n{}\n    ];",
        rings.join("\n")
    )
}

/// Карточка ошибки.
#[derive(Clone, Debug, Default)]
pub struct Card {
    title: String,
    message: String,
    command: Option<String>,
    layers: Vec<(String, String, Shape)>,
    code: Vec<(String, Shape)>,
    notes: Vec<String>,
}

fn escape(s: &str) -> String {
    s.replace('&', "&amp;")
        .replace('<', "&lt;")
        .replace('>', "&gt;")
}

impl Card {
    /// Карточка с заголовком (обычно — имя теста). Команда повтора берётся из
    /// текущего случая `for_cases`, если он идёт.
    pub fn new(title: &str) -> Self {
        Self {
            title: title.to_string(),
            command: crate::testkit::replay_command(),
            ..Self::default()
        }
    }

    /// Текст падения.
    pub fn message(mut self, text: &str) -> Self {
        self.message = text.to_string();
        self
    }

    /// Слой картинки.
    pub fn shape(mut self, label: &str, color: &str, s: &Shape) -> Self {
        self.layers
            .push((label.to_string(), color.to_string(), s.clone()));
        self
    }

    /// Данные для каталога (кодом на Rust).
    pub fn code(mut self, name: &str, s: &Shape) -> Self {
        self.code.push((name.to_string(), s.clone()));
        self
    }

    /// Пометка (пустая не добавляется).
    pub fn note(mut self, text: &str) -> Self {
        if !text.is_empty() {
            self.notes.push(text.to_string());
        }
        self
    }

    /// Картинка слоёв.
    fn picture(&self) -> String {
        let mut svg = Svg::new();
        let small = self
            .layers
            .iter()
            .map(|l| vertex_count(&l.2))
            .sum::<usize>()
            <= 24;
        for (_, color, s) in &self.layers {
            for r in s {
                let pts: Vec<Vec2> = r
                    .iter()
                    .map(|&(x, y)| Vec2::new(x as f64 * 1e-3, y as f64 * 1e-3))
                    .collect();
                svg.polygon(&pts, Style::filled(color, 0.18).width(1.6));
                if small {
                    for (k, p) in pts.iter().enumerate() {
                        svg.dot(*p, color);
                        svg.label(*p, &format!("{:?}", r[k]), color);
                    }
                }
            }
        }
        svg.render(760.0)
    }

    /// Страница HTML.
    pub fn html(&self) -> String {
        let legend: Vec<String> = self
            .layers
            .iter()
            .map(|(label, color, s)| {
                format!(
                    "<span><i style=\"background:{color}\"></i>{} — контуров {}, вершин {}</span>",
                    escape(label),
                    s.len(),
                    vertex_count(s)
                )
            })
            .collect();
        let code: Vec<String> = self.code.iter().map(|(n, s)| shape_code(n, s)).collect();
        let mut h = String::from(
            "<!doctype html>\n<html lang=\"ru\"><head><meta charset=\"utf-8\">\
             <meta name=\"viewport\" content=\"width=device-width, initial-scale=1\">\n",
        );
        h.push_str(&format!(
            "<title>Карточка: {}</title>\n",
            escape(&self.title)
        ));
        h.push_str(
            "<style>body{font:15px/1.45 system-ui,sans-serif;margin:24px auto;max-width:820px;\
             padding:0 16px;color:#1b1b1b;background:#fff}h1{font-size:22px}h2{font-size:16px;\
             margin-top:22px}pre{background:#f4f4f4;padding:10px 12px;overflow-x:auto;\
             border-radius:6px}.legend span{margin-right:18px;display:inline-block}\
             .legend i{display:inline-block;width:12px;height:12px;margin-right:6px;\
             border-radius:2px}svg{max-width:100%;height:auto;border:1px solid #ddd}</style>\n\
             </head><body>\n",
        );
        h.push_str(&format!(
            "<h1>Карточка ошибки: {}</h1>\n",
            escape(&self.title)
        ));
        h.push_str(&format!(
            "<h2>Что упало</h2>\n<pre>{}</pre>\n",
            escape(&self.message)
        ));
        if let Some(c) = &self.command {
            h.push_str(&format!(
                "<h2>Как повторить</h2>\n<pre>{}</pre>\n",
                escape(c)
            ));
        }
        for n in &self.notes {
            h.push_str(&format!("<p>{}</p>\n", escape(n)));
        }
        h.push_str("<h2>Картинка</h2>\n");
        h.push_str(&format!(
            "<div class=\"legend\">{}</div>\n",
            legend.join("")
        ));
        h.push_str(&self.picture());
        if !code.is_empty() {
            h.push_str(&format!(
                "<h2>Для каталога</h2>\n<pre>{}</pre>\n",
                escape(&code.join("\n"))
            ));
        }
        h.push_str("</body></html>\n");
        h
    }

    /// Сохранить в `target/cards` (или в `YADRO_CARDS`); вернуть путь.
    pub fn save(&self) -> std::io::Result<PathBuf> {
        let dir = std::env::var_os("YADRO_CARDS").map_or_else(
            || PathBuf::from(concat!(env!("CARGO_MANIFEST_DIR"), "/../../target/cards")),
            PathBuf::from,
        );
        std::fs::create_dir_all(&dir)?;
        let slug: String = self
            .title
            .chars()
            .map(|c| {
                if c.is_ascii_alphanumeric() || c == '_' {
                    c
                } else {
                    '-'
                }
            })
            .collect();
        let tag = self.command.as_deref().map_or_else(
            || format!("{:x}", self.message.len() * 31 + vertex_count_all(self)),
            |c| {
                c.split_whitespace()
                    .next()
                    .and_then(|s| s.strip_prefix("YADRO_SEED="))
                    .unwrap_or("0")
                    .to_string()
            },
        );
        let path = dir.join(format!("{slug}-{tag}.html"));
        std::fs::write(&path, self.html())?;
        Ok(path)
    }
}

fn vertex_count_all(c: &Card) -> usize {
    c.layers.iter().map(|l| vertex_count(&l.2)).sum()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn shrink_finds_a_small_example() {
        // «Ошибка»: есть вершина с x > 1000. Минимум — один треугольник.
        let start: Shape = vec![
            (0..12)
                .map(|k| (k * 371 - 900, (k * k * 37) % 1500))
                .collect(),
            vec![(5000, 5000), (6000, 5000), (6000, 7000), (5000, 7000)],
        ];
        let fails = |s: &Shape| s.iter().flatten().any(|p| p.0 > 1000);
        let small = shrink_shapes(start, fails);
        assert!(fails(&small));
        assert_eq!(small.len(), 1);
        assert_eq!(small[0].len(), 3);
    }

    #[test]
    fn rounder_shortens_numbers() {
        assert_eq!(rounder(1234), Some(1230));
        assert_eq!(rounder(1230), Some(1200));
        assert_eq!(rounder(1000), Some(0));
        assert_eq!(rounder(7), Some(0));
        assert_eq!(rounder(0), None);
        assert_eq!(rounder(-1250), Some(-1300));
    }

    #[test]
    fn card_has_everything() {
        let s: Shape = vec![vec![(0, 0), (1000, 0), (0, 1000)]];
        let html = Card::new("пример")
            .message("сломалось <здесь>")
            .shape("A", "#1f5fbf", &s)
            .code("a", &s)
            .html();
        assert!(html.contains("сломалось &lt;здесь&gt;"));
        assert!(html.contains("<svg"));
        // В HTML «&» экранируется.
        assert!(html.contains("let a: &amp;[&amp;[(i64, i64)]] = &amp;["));
        assert!(html.contains("(1000, 0)"));
    }
}
