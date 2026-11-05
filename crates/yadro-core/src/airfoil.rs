//! Профили крыла: чтение файлов .dat (база UIUC, форматы Selig и Lednicer),
//! NACA 4-й серии, толщина и средняя линия, контур нервюры в мм.
//!
//! Профиль хранится в долях хорды: носок в (0, 0), середина задней кромки в (1, 0).
//! Файл, который уже так записан (носок и кромка не дальше 0,001 от этих точек), —
//! а так записаны файлы базы UIUC — берётся как есть, без перестроения. Иначе
//! (файл в мм, повёрнутый, сдвинутый) профиль приводится к хорде 1, и носком
//! считается точка, самая далёкая от середины задней кромки: так приведение не
//! зависит от того, как файл был повёрнут и растянут.
//!
//! Форматы:
//! - Selig — точки от задней кромки по верху к носку и по низу обратно;
//! - Lednicer — после названия строка с числами точек верха и низа, затем верх и
//!   низ отдельно, каждый от носка к задней кромке.
//!
//! Строки, начинающиеся с `#`, и пустые строки пропускаются; числа можно разделять
//! пробелами, табуляцией или запятыми; экспонента может быть с `D` (Фортран).

use crate::region::Region;
use crate::ring::{Ring, RingError};
use crate::vec::Vec2;
use crate::xform::Xform2;
use std::fmt;

/// Почему не удалось прочитать профиль.
#[derive(Clone, Debug, PartialEq)]
pub enum AirfoilError {
    /// В файле нет точек.
    Empty,
    /// Строка данных не из двух чисел.
    BadLine {
        /// Номер строки (с 1).
        line: usize,
        /// Сама строка.
        text: String,
    },
    /// Не число или бесконечность.
    NotFinite {
        /// Номер строки (с 1).
        line: usize,
    },
    /// Слишком мало точек: нужно хотя бы по две на каждую поверхность.
    TooFewPoints(usize),
    /// Lednicer: заявлено одно число точек, а в файле другое.
    LednicerCounts {
        /// Сколько заявлено (верх + низ).
        declared: usize,
        /// Сколько есть.
        found: usize,
    },
    /// Хорда нулевая: носок совпадает с задней кромкой.
    ZeroChord,
}

impl fmt::Display for AirfoilError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Empty => write!(f, "в файле профиля нет точек"),
            Self::BadLine { line, text } => {
                write!(f, "строка {line}: ожидались два числа, а здесь «{text}»")
            }
            Self::NotFinite { line } => write!(f, "строка {line}: не число или бесконечность"),
            Self::TooFewPoints(n) => write!(f, "в профиле {n} точек — слишком мало"),
            Self::LednicerCounts { declared, found } => write!(
                f,
                "формат Lednicer: заявлено {declared} точек, а в файле {found}"
            ),
            Self::ZeroChord => write!(f, "хорда профиля нулевая"),
        }
    }
}

impl std::error::Error for AirfoilError {}

/// Профиль, приведённый к хорде 1.
#[derive(Clone, Debug, PartialEq)]
pub struct Airfoil {
    name: String,
    /// Верх от носка к задней кромке.
    upper: Vec<Vec2>,
    /// Низ от носка к задней кромке.
    lower: Vec<Vec2>,
}

/// Числа строки или `None`, если строка — не числа.
fn numbers(line: &str) -> Option<Vec<f64>> {
    let toks: Vec<&str> = line
        .split(|c: char| c.is_whitespace() || c == ',' || c == ';')
        .filter(|t| !t.is_empty())
        .collect();
    if toks.is_empty() {
        return None;
    }
    toks.iter()
        .map(|t| t.replace(['D', 'd'], "e").parse::<f64>().ok())
        .collect()
}

/// Похоже ли на файл в долях хорды: носок у (0, 0), середина задней кромки у (1, 0).
///
/// РИСК: от порога зависит, перестраивается ли файл. 0,001 хорды — с запасом для
/// файлов UIUC (у них ровно 0 и 1) и мало для файла в мм или повёрнутого.
fn is_unit(le: Vec2, te: Vec2) -> bool {
    const NEAR: f64 = 1e-3;
    le.length() <= NEAR && te.distance(Vec2::new(1.0, 0.0)) <= NEAR
}

/// y на ломаной поверхности при данном x (линейно). Поверхность идёт от носка к
/// задней кромке; берётся участок, ближайший к задней кромке, — у носка ломаная
/// может чуть заходить назад по x. За пределами — по крайней точке.
fn surface_y(s: &[Vec2], x: f64) -> f64 {
    for w in s.windows(2).rev() {
        let (a, b) = (w[0], w[1]);
        let (lo, hi) = if a.x <= b.x { (a, b) } else { (b, a) };
        if lo.x <= x && x <= hi.x {
            if hi.x == lo.x {
                return 0.5 * (lo.y + hi.y);
            }
            return lo.y + (hi.y - lo.y) * (x - lo.x) / (hi.x - lo.x);
        }
    }
    if x <= s[0].x {
        s[0].y
    } else {
        s[s.len() - 1].y
    }
}

impl Airfoil {
    /// Профиль из верха и низа (каждый от носка к задней кромке) с приведением к
    /// хорде 1. Если «верх» оказался ниже «низа», они меняются местами.
    pub fn from_surfaces(
        name: &str,
        mut upper: Vec<Vec2>,
        mut lower: Vec<Vec2>,
    ) -> Result<Self, AirfoilError> {
        upper.dedup();
        lower.dedup();
        if upper.len() < 2 || lower.len() < 2 || upper.len() + lower.len() < 5 {
            return Err(AirfoilError::TooFewPoints(upper.len() + lower.len()));
        }
        // Носок общий: если поверхности начинаются в разных точках, носок — середина.
        let le = (upper[0] + lower[0]) * 0.5;
        let te = (upper[upper.len() - 1] + lower[lower.len() - 1]) * 0.5;
        let c = te - le;
        let len = c.length();
        if len == 0.0 || !len.is_finite() {
            return Err(AirfoilError::ZeroChord);
        }
        // Файл уже в долях хорды — координаты берутся как есть, без перестроения.
        let xf = if is_unit(le, te) {
            Xform2::IDENTITY
        } else {
            Xform2::translation(-le)
                .then(&Xform2::rotation(-c.y.atan2(c.x)))
                .then(&Xform2::scaling(1.0 / len, 1.0 / len))
        };
        let map = |s: &[Vec2]| s.iter().map(|&p| xf.apply_point(p)).collect::<Vec<_>>();
        let (mut up, mut lo) = (map(&upper), map(&lower));
        // Ориентация: площадь контура «задняя кромка → верх → носок → низ» > 0.
        let mut ring: Vec<Vec2> = up.iter().rev().copied().collect();
        ring.extend(lo.iter().skip(1).copied());
        let area: f64 = (0..ring.len())
            .map(|i| ring[i].cross(ring[(i + 1) % ring.len()]))
            .sum();
        if area < 0.0 {
            std::mem::swap(&mut up, &mut lo);
        }
        Ok(Self {
            name: name.trim().to_string(),
            upper: up,
            lower: lo,
        })
    }

    /// Прочитать файл .dat (Selig или Lednicer, определяется сам).
    pub fn parse_dat(text: &str) -> Result<Self, AirfoilError> {
        let mut name = String::new();
        let mut pts: Vec<(Vec2, usize)> = Vec::new();
        for (k, raw) in text.lines().enumerate() {
            let line = raw.trim();
            if line.is_empty() || line.starts_with('#') {
                continue;
            }
            match numbers(line) {
                Some(v) if v.len() == 2 => {
                    if !v[0].is_finite() || !v[1].is_finite() {
                        return Err(AirfoilError::NotFinite { line: k + 1 });
                    }
                    pts.push((Vec2::new(v[0], v[1]), k + 1));
                }
                _ if pts.is_empty() => {
                    if name.is_empty() {
                        name = line.to_string();
                    }
                }
                _ => {
                    return Err(AirfoilError::BadLine {
                        line: k + 1,
                        text: line.to_string(),
                    });
                }
            }
        }
        let Some(&(first, _)) = pts.first() else {
            return Err(AirfoilError::Empty);
        };
        let pts: Vec<Vec2> = pts.into_iter().map(|(p, _)| p).collect();
        let is_count = |v: f64| v >= 2.0 && v.fract() == 0.0;
        // Lednicer: первая «точка» — два целых числа точек. Если сумма не сходится,
        // а остальные точки похожи на хорду 1, это ошибка в файле; иначе это Selig
        // в мм с целой первой точкой.
        if is_count(first.x) && is_count(first.y) {
            let (nu, nl) = (first.x as usize, first.y as usize);
            let rest = &pts[1..];
            if rest.len() != nu + nl {
                let unit = rest.iter().all(|p| p.x.abs() <= 1.5 && p.y.abs() <= 1.5);
                if unit {
                    return Err(AirfoilError::LednicerCounts {
                        declared: nu + nl,
                        found: rest.len(),
                    });
                }
                return Self::parse_selig(&name, pts);
            }
            let mut up = rest[..nu].to_vec();
            let mut lo = rest[nu..].to_vec();
            // Обычно от носка; если наоборот — развернуть.
            if up[0].x > up[up.len() - 1].x {
                up.reverse();
            }
            if lo[0].x > lo[lo.len() - 1].x {
                lo.reverse();
            }
            return Self::from_surfaces(&name, up, lo);
        }
        Self::parse_selig(&name, pts)
    }

    /// Selig: носок — точка, самая далёкая от середины задней кромки.
    fn parse_selig(name: &str, pts: Vec<Vec2>) -> Result<Self, AirfoilError> {
        let mut p = pts;
        p.dedup();
        if p.len() < 5 {
            return Err(AirfoilError::TooFewPoints(p.len()));
        }
        let te = (p[0] + p[p.len() - 1]) * 0.5;
        let near_origin = (0..p.len())
            .min_by(|&i, &j| p[i].length().total_cmp(&p[j].length()))
            .expect("точки есть");
        let le = if is_unit(p[near_origin], te) {
            near_origin
        } else {
            (0..p.len())
                .max_by(|&i, &j| p[i].distance(te).total_cmp(&p[j].distance(te)))
                .expect("точки есть")
        };
        let up: Vec<Vec2> = p[..=le].iter().rev().copied().collect();
        let lo: Vec<Vec2> = p[le..].to_vec();
        Self::from_surfaces(name, up, lo)
    }

    /// NACA 4-й серии: `m` — наибольшая вогнутость, `p` — её место, `t` — толщина
    /// (доли хорды); `n` точек на поверхность (косинусная разбивка, гуще у кромок).
    /// `closed_te` — замкнутая задняя кромка (коэффициент −0,1036 вместо −0,1015).
    pub fn naca4(m: f64, p: f64, t: f64, n: usize, closed_te: bool) -> Self {
        let n = n.max(3);
        let k4 = if closed_te { -0.1036 } else { -0.1015 };
        let yt = |x: f64| {
            5.0 * t
                * (0.2969 * x.sqrt() - 0.1260 * x - 0.3516 * x * x
                    + 0.2843 * x.powi(3)
                    + k4 * x.powi(4))
        };
        let camber = |x: f64| {
            if m == 0.0 || p <= 0.0 || p >= 1.0 {
                (0.0, 0.0)
            } else if x < p {
                (
                    m / (p * p) * (2.0 * p * x - x * x),
                    2.0 * m / (p * p) * (p - x),
                )
            } else {
                let q = (1.0 - p) * (1.0 - p);
                (
                    m / q * (1.0 - 2.0 * p + 2.0 * p * x - x * x),
                    2.0 * m / q * (p - x),
                )
            }
        };
        let side = |x: f64, s: f64| {
            let (yc, dy) = camber(x);
            let th = dy.atan();
            let y = yt(x);
            Vec2::new(x - s * y * th.sin(), yc + s * y * th.cos())
        };
        let xs: Vec<f64> = (0..n)
            .map(|i| 0.5 * (1.0 - (std::f64::consts::PI * i as f64 / (n - 1) as f64).cos()))
            .collect();
        let name = format!("NACA {:.0}{:.0}{:02.0}", m * 100.0, p * 10.0, t * 100.0);
        Self {
            name,
            upper: xs.iter().map(|&x| side(x, 1.0)).collect(),
            lower: xs.iter().map(|&x| side(x, -1.0)).collect(),
        }
    }

    /// NACA 4-й серии по коду вроде «2412»; `None`, если код не из четырёх цифр.
    pub fn naca4_code(code: &str, n: usize, closed_te: bool) -> Option<Self> {
        let d: Vec<u32> = code
            .trim()
            .chars()
            .map(|c| c.to_digit(10))
            .collect::<Option<_>>()?;
        if d.len() != 4 {
            return None;
        }
        Some(Self::naca4(
            f64::from(d[0]) / 100.0,
            f64::from(d[1]) / 10.0,
            f64::from(d[2] * 10 + d[3]) / 100.0,
            n,
            closed_te,
        ))
    }

    /// Название.
    pub fn name(&self) -> &str {
        &self.name
    }

    /// Верх от носка к задней кромке (хорда 1).
    pub fn upper(&self) -> &[Vec2] {
        &self.upper
    }

    /// Низ от носка к задней кромке (хорда 1).
    pub fn lower(&self) -> &[Vec2] {
        &self.lower
    }

    /// y верха при данном x (доли хорды).
    pub fn upper_y(&self, x: f64) -> f64 {
        surface_y(&self.upper, x)
    }

    /// y низа при данном x (доли хорды).
    pub fn lower_y(&self, x: f64) -> f64 {
        surface_y(&self.lower, x)
    }

    /// Толщина по вертикали при данном x (доли хорды).
    pub fn thickness(&self, x: f64) -> f64 {
        self.upper_y(x) - self.lower_y(x)
    }

    /// Средняя линия при данном x (доли хорды).
    pub fn camber(&self, x: f64) -> f64 {
        0.5 * (self.upper_y(x) + self.lower_y(x))
    }

    /// Наибольшая толщина и её место: (x, толщина), доли хорды.
    pub fn max_thickness(&self) -> (f64, f64) {
        self.upper
            .iter()
            .chain(&self.lower)
            .map(|p| p.x.clamp(0.0, 1.0))
            .map(|x| (x, self.thickness(x)))
            .fold(
                (0.0, f64::NEG_INFINITY),
                |a, b| if b.1 > a.1 { b } else { a },
            )
    }

    /// Зазор задней кромки (верх минус низ), доли хорды.
    pub fn te_gap(&self) -> f64 {
        self.upper[self.upper.len() - 1].y - self.lower[self.lower.len() - 1].y
    }

    /// Профиль с задней кромкой не тоньше `min_gap` (доли хорды): недостающее
    /// добавляется линейно по x — половина вверх, половина вниз, носок не меняется.
    /// Нужно для печати и фрезы: кромка тоньше сопла или фрезы не получится.
    pub fn with_min_te_gap(&self, min_gap: f64) -> Self {
        let gap = self.te_gap();
        if gap >= min_gap {
            return self.clone();
        }
        let d = 0.5 * (min_gap - gap);
        let shift = |s: &[Vec2], k: f64| {
            s.iter()
                .map(|p| Vec2::new(p.x, p.y + k * d * p.x.clamp(0.0, 1.0)))
                .collect()
        };
        Self {
            name: self.name.clone(),
            upper: shift(&self.upper, 1.0),
            lower: shift(&self.lower, -1.0),
        }
    }

    /// Преобразование, которое ставит профиль с хордой `chord_mm` на место:
    /// поворот на `twist_deg` (> 0 — носок вверх) вокруг точки хорды на доле `pivot`
    /// (обычно 0,25) и перенос этой точки в `at`.
    pub fn placement(chord_mm: f64, twist_deg: f64, pivot: f64, at: Vec2) -> Xform2 {
        Xform2::scaling(chord_mm, chord_mm)
            .then(&Xform2::translation(Vec2::new(-pivot * chord_mm, 0.0)))
            .then(&Xform2::rotation(-twist_deg.to_radians()))
            .then(&Xform2::translation(at))
    }

    /// Контур в мм: задняя кромка → верх → носок → низ → задняя кромка (против
    /// часовой стрелки, если `xf` не зеркалит). `xf` применяется к профилю с хордой 1,
    /// масштаб задаёт оно же (см. [`Airfoil::placement`]).
    pub fn outline(&self, xf: &Xform2) -> Vec<Vec2> {
        self.upper
            .iter()
            .rev()
            .chain(self.lower.iter().skip(1))
            .map(|&p| xf.apply_point(p))
            .collect()
    }

    /// Контур на сетке.
    pub fn to_ring(&self, xf: &Xform2) -> Result<Ring, RingError> {
        Ring::from_mm(&self.outline(xf))
    }

    /// Область профиля (по правилу ненулевой обмотки — перехлёст у кромки не мешает).
    pub fn to_region(&self, xf: &Xform2) -> Result<Region, RingError> {
        Ok(Region::from_ring(self.to_ring(xf)?))
    }

    /// Текст в формате Selig с `decimals` знаками после запятой.
    pub fn to_selig(&self, decimals: usize) -> String {
        let mut s = format!("{}\n", self.name);
        for p in self.outline(&Xform2::IDENTITY) {
            s.push_str(&format!("{:.*} {:.*}\n", decimals, p.x, decimals, p.y));
        }
        s
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn naca_basics() {
        let a = Airfoil::naca4_code("0012", 101, true).unwrap();
        assert_eq!(a.name(), "NACA 0012");
        let (x, t) = a.max_thickness();
        assert!((x - 0.3).abs() < 0.02, "место толщины {x}");
        assert!((t - 0.1205).abs() < 0.001, "толщина {t}");
        assert!(a.te_gap().abs() < 1e-12);
        assert!(a.camber(0.5).abs() < 1e-12);
        let b = Airfoil::naca4_code("2412", 101, false).unwrap();
        assert!(b.te_gap() > 0.002 && b.te_gap() < 0.003);
        assert!((b.camber(0.4) - 0.02).abs() < 1e-3);
        assert!(Airfoil::naca4_code("24x2", 10, true).is_none());
    }

    #[test]
    fn te_gap_thickening() {
        let a = Airfoil::naca4_code("2412", 81, true).unwrap();
        let b = a.with_min_te_gap(0.004);
        assert!((b.te_gap() - 0.004).abs() < 1e-12);
        assert!(
            (b.thickness(0.0) - a.thickness(0.0)).abs() < 1e-12,
            "носок на месте"
        );
        assert_eq!(b.with_min_te_gap(0.001), b, "уже толще — без изменений");
    }

    #[test]
    fn numbers_parsing() {
        assert_eq!(numbers("1.0 2.5"), Some(vec![1.0, 2.5]));
        assert_eq!(numbers(" 1,\t-2.5D-01 "), Some(vec![1.0, -0.25]));
        assert_eq!(numbers("CLARK Y"), None);
        assert_eq!(numbers("NACA 2412"), None);
    }
}
