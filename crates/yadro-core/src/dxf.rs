//! DXF: запись листа в R12 и чтение контуров из файлов других программ.
//!
//! Запись: контуры — POLYLINE с выпуклостью или CIRCLE, подписи — TEXT.
//!
//! Чтение понимает POLYLINE (с VERTEX и выпуклостью), LWPOLYLINE, CIRCLE, а также
//! LINE и ARC — их концы сцепляются в замкнутые контуры (с зазором до
//! [`JOIN_TOL_MM`]). Остальные объекты (TEXT, DIMENSION, INSERT…) пропускаются.
//! Пустые строки по краям, CRLF и комментарии (код 999) не мешают.

use crate::path::Path;
use crate::sheet::{Sheet, is_circle, num};
use crate::vec::Vec2;
use std::fmt;

/// Зазор, до которого концы отрезков и дуг считаются совпавшими, мм.
pub const JOIN_TOL_MM: f64 = 0.001;

/// Почему не удалось прочитать DXF.
#[derive(Clone, Debug, PartialEq)]
pub enum DxfError {
    /// Строка кода группы — не целое число.
    BadCode {
        /// Номер строки (с 1).
        line: usize,
    },
    /// Значение числовой группы — не число.
    BadNumber {
        /// Номер строки (с 1).
        line: usize,
    },
    /// После кода группы нет значения.
    Truncated,
    /// Отрезки и дуги не сцепляются в замкнутый контур.
    Unclosed {
        /// Слой.
        layer: String,
        /// Где оборвалась цепочка, мм.
        at: Vec2,
    },
}

impl fmt::Display for DxfError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::BadCode { line } => write!(f, "DXF, строка {line}: код группы — не число"),
            Self::BadNumber { line } => write!(f, "DXF, строка {line}: ожидалось число"),
            Self::Truncated => write!(f, "DXF оборван: после кода группы нет значения"),
            Self::Unclosed { layer, at } => write!(
                f,
                "DXF, слой {layer}: контур не замкнут у точки ({:.3}; {:.3})",
                at.x, at.y
            ),
        }
    }
}

impl std::error::Error for DxfError {}

/// Записать лист в DXF R12.
pub fn write(sheet: &Sheet) -> String {
    let (lo, hi) = sheet.bounds().unwrap_or((Vec2::ZERO, Vec2::ZERO));
    let mut s = String::new();
    let mut g = |code: u32, value: &str| {
        s.push_str(&format!("{code}\n{value}\n"));
    };
    g(0, "SECTION");
    g(2, "HEADER");
    g(9, "$ACADVER");
    g(1, "AC1009");
    g(9, "$INSUNITS");
    g(70, "4");
    g(9, "$EXTMIN");
    g(10, &num(lo.x));
    g(20, &num(lo.y));
    g(9, "$EXTMAX");
    g(10, &num(hi.x));
    g(20, &num(hi.y));
    g(0, "ENDSEC");
    g(0, "SECTION");
    g(2, "TABLES");
    g(0, "TABLE");
    g(2, "LAYER");
    let layers = sheet.layers();
    g(70, &layers.len().to_string());
    for (k, name) in layers.iter().enumerate() {
        g(0, "LAYER");
        g(2, name);
        g(70, "0");
        g(62, &[7, 1, 5, 3, 6, 30][k % 6].to_string());
        g(6, "CONTINUOUS");
    }
    g(0, "ENDTAB");
    g(0, "ENDSEC");
    g(0, "SECTION");
    g(2, "ENTITIES");
    for it in sheet.items() {
        if let Some((c, r)) = is_circle(&it.path) {
            g(0, "CIRCLE");
            g(8, &it.layer);
            g(10, &num(c.x));
            g(20, &num(c.y));
            g(30, "0");
            g(40, &num(r));
            continue;
        }
        g(0, "POLYLINE");
        g(8, &it.layer);
        g(66, "1");
        g(10, "0");
        g(20, "0");
        g(30, "0");
        g(70, "1");
        for v in it.path.vertices() {
            g(0, "VERTEX");
            g(8, &it.layer);
            g(10, &num(v.at.x));
            g(20, &num(v.at.y));
            g(30, "0");
            if v.bulge != 0.0 {
                g(42, &format!("{:.12}", v.bulge));
            }
        }
        g(0, "SEQEND");
        g(8, &it.layer);
    }
    for tx in sheet.texts() {
        g(0, "TEXT");
        g(8, &tx.layer);
        g(10, &num(tx.at.x));
        g(20, &num(tx.at.y));
        g(30, "0");
        g(40, &num(tx.height));
        g(1, &dxf_text(&tx.text));
    }
    g(0, "ENDSEC");
    g(0, "EOF");
    s
}

/// Текст для R12: не-ASCII символы — как `\U+XXXX` (так их понимают AutoCAD и
/// большинство CAM), переводы строк — пробелами.
fn dxf_text(s: &str) -> String {
    s.chars()
        .map(|c| match c {
            '\n' | '\r' => " ".to_string(),
            c if c.is_ascii() => c.to_string(),
            c => format!("\\U+{:04X}", u32::from(c)),
        })
        .collect()
}

/// Пары (код, значение, номер строки значения).
fn pairs(text: &str) -> Result<Vec<(i32, String, usize)>, DxfError> {
    let lines: Vec<&str> = text.lines().collect();
    let mut out = Vec::new();
    let mut k = 0;
    while k < lines.len() {
        let code = lines[k].trim();
        if code.is_empty() && k + 1 >= lines.len() {
            break; // пустой хвост
        }
        let code: i32 = code
            .parse()
            .map_err(|_| DxfError::BadCode { line: k + 1 })?;
        let value = lines.get(k + 1).ok_or(DxfError::Truncated)?;
        out.push((code, value.trim().to_string(), k + 2));
        k += 2;
    }
    Ok(out)
}

/// Кусок, который ещё надо сцепить: из `a` в `b` с выпуклостью.
struct Piece {
    a: Vec2,
    b: Vec2,
    bulge: f64,
    layer: String,
}

/// Объект DXF: тип и его группы.
struct Entity {
    kind: String,
    groups: Vec<(i32, String, usize)>,
}

impl Entity {
    fn num(&self, code: i32) -> Result<Option<f64>, DxfError> {
        match self.groups.iter().find(|g| g.0 == code) {
            None => Ok(None),
            Some((_, v, line)) => v
                .parse::<f64>()
                .ok()
                .filter(|x| x.is_finite())
                .map(Some)
                .ok_or(DxfError::BadNumber { line: *line }),
        }
    }

    fn need(&self, code: i32) -> Result<f64, DxfError> {
        Ok(self.num(code)?.unwrap_or(0.0))
    }

    fn layer(&self) -> String {
        self.groups
            .iter()
            .find(|g| g.0 == 8)
            .map_or_else(|| "0".to_string(), |g| g.1.clone())
    }
}

/// Разбить секцию ENTITIES (или весь файл, если секций нет) на объекты.
fn entities(p: &[(i32, String, usize)]) -> Vec<Entity> {
    let start = p
        .windows(2)
        .position(|w| w[0].0 == 2 && w[0].1 == "ENTITIES")
        .map_or(0, |i| i + 1);
    let mut out: Vec<Entity> = Vec::new();
    for g in &p[start..] {
        if g.0 == 0 {
            if g.1 == "ENDSEC" || g.1 == "EOF" {
                break;
            }
            out.push(Entity {
                kind: g.1.clone(),
                groups: Vec::new(),
            });
        } else if g.0 != 999
            && let Some(e) = out.last_mut()
        {
            e.groups.push(g.clone());
        }
    }
    out
}

/// Прочитать контуры из DXF; концы отрезков и дуг ближе `join_tol_mm` сцепляются.
pub fn read(text: &str, join_tol_mm: f64) -> Result<Sheet, DxfError> {
    let p = pairs(text)?;
    let ents = entities(&p);
    let mut sheet = Sheet::new();
    let mut pieces: Vec<Piece> = Vec::new();
    let mut k = 0;
    while k < ents.len() {
        let e = &ents[k];
        k += 1;
        match e.kind.as_str() {
            "CIRCLE" => {
                let (c, r) = (Vec2::new(e.need(10)?, e.need(20)?), e.need(40)?);
                if r > 0.0 {
                    sheet.add_path(&e.layer(), Path::circle(c, r));
                }
            }
            "LINE" => pieces.push(Piece {
                a: Vec2::new(e.need(10)?, e.need(20)?),
                b: Vec2::new(e.need(11)?, e.need(21)?),
                bulge: 0.0,
                layer: e.layer(),
            }),
            "ARC" => {
                let (c, r) = (Vec2::new(e.need(10)?, e.need(20)?), e.need(40)?);
                let (s, t) = (e.need(50)?.to_radians(), e.need(51)?.to_radians());
                let mut sweep = (t - s).rem_euclid(std::f64::consts::TAU);
                if sweep == 0.0 {
                    sweep = std::f64::consts::TAU;
                }
                let at = |a: f64| c + Vec2::new(a.cos(), a.sin()) * r;
                pieces.push(Piece {
                    a: at(s),
                    b: at(s + sweep),
                    bulge: (sweep / 4.0).tan(),
                    layer: e.layer(),
                });
            }
            "LWPOLYLINE" => {
                let mut verts: Vec<(Vec2, f64)> = Vec::new();
                let mut x = None;
                for (code, v, line) in &e.groups {
                    let val = || {
                        v.parse::<f64>()
                            .ok()
                            .filter(|x| x.is_finite())
                            .ok_or(DxfError::BadNumber { line: *line })
                    };
                    match code {
                        10 => x = Some(val()?),
                        20 => verts.push((Vec2::new(x.take().unwrap_or(0.0), val()?), 0.0)),
                        42 => {
                            if let Some(last) = verts.last_mut() {
                                last.1 = val()?;
                            }
                        }
                        _ => {}
                    }
                }
                let closed = e.need(70)? as i64 & 1 == 1;
                polyline(
                    &mut sheet,
                    &mut pieces,
                    &e.layer(),
                    verts,
                    closed,
                    join_tol_mm,
                );
            }
            "POLYLINE" => {
                let closed = e.need(70)? as i64 & 1 == 1;
                let mut verts = Vec::new();
                while k < ents.len() && ents[k].kind == "VERTEX" {
                    let v = &ents[k];
                    verts.push((Vec2::new(v.need(10)?, v.need(20)?), v.need(42)?));
                    k += 1;
                }
                if k < ents.len() && ents[k].kind == "SEQEND" {
                    k += 1;
                }
                polyline(
                    &mut sheet,
                    &mut pieces,
                    &e.layer(),
                    verts,
                    closed,
                    join_tol_mm,
                );
            }
            _ => {}
        }
    }
    chain(&mut sheet, pieces, join_tol_mm)?;
    Ok(sheet)
}

/// Полилиния: замкнутая (или с совпавшими концами) — сразу контур, иначе — куски.
fn polyline(
    sheet: &mut Sheet,
    pieces: &mut Vec<Piece>,
    layer: &str,
    mut verts: Vec<(Vec2, f64)>,
    closed: bool,
    tol: f64,
) {
    if verts.len() < 2 {
        return;
    }
    let ends_meet = verts[0].0.distance(verts[verts.len() - 1].0) <= tol;
    if closed || ends_meet {
        if !closed {
            verts.pop(); // последняя точка повторяет первую
        }
        let mut p = Path::new();
        for (at, b) in verts {
            p.push(at, b);
        }
        sheet.add_path(layer, p);
    } else {
        for w in verts.windows(2) {
            pieces.push(Piece {
                a: w[0].0,
                b: w[1].0,
                bulge: w[0].1,
                layer: layer.to_string(),
            });
        }
    }
}

/// Сцепить куски в замкнутые контуры (внутри одного слоя).
fn chain(sheet: &mut Sheet, mut pieces: Vec<Piece>, tol: f64) -> Result<(), DxfError> {
    while let Some(first) = pieces.pop() {
        let start = first.a;
        let layer = first.layer.clone();
        let mut path = Path::new();
        let mut end = first.b;
        path.push(first.a, first.bulge);
        while end.distance(start) > tol {
            // Следующий кусок: начинается в конце цепочки или кончается там (тогда развернуть).
            let found = pieces.iter().position(|q| {
                q.layer == layer && (q.a.distance(end) <= tol || q.b.distance(end) <= tol)
            });
            let Some(i) = found else {
                return Err(DxfError::Unclosed { layer, at: end });
            };
            let q = pieces.swap_remove(i);
            let (a, b, bulge) = if q.a.distance(end) <= tol {
                (q.a, q.b, q.bulge)
            } else {
                (q.b, q.a, -q.bulge)
            };
            path.push(a, bulge);
            end = b;
        }
        sheet.add_path(&layer, path);
    }
    Ok(())
}
