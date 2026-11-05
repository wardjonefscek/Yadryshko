//! Контур из поля плотности TopOpt: где плотность не меньше порога — материал.
//!
//! Поле — плотности элементов прямоугольной сетки (значения в центрах элементов).
//! Граница строится методом marching squares по линейной интерполяции между
//! центрами, поэтому получается гладкой, а не «лесенкой» по краям элементов:
//! - каждый отсчёт с плотностью ≥ порога оказывается внутри области, ниже — снаружи;
//! - седловые клетки (материал по диагонали) решаются по значению в седловой точке
//!   билинейной интерполяции (asymptotic decider), так что связность честная;
//! - `refine > 1` — поле сначала уплотняется бикубической интерполяцией
//!   Catmull-Rom (она проходит через исходные отсчёты), граница становится плавнее;
//! - края расчётной области точные: поле продолжается за край тем же значением,
//!   а результат обрезается по прямоугольнику области;
//! - мелкие острова и дыры (шум оптимизации) можно отбросить по площади.
//!
//! Контуры собираются в [`Region`] обычной нормализацией, поэтому результат
//! корректен, как любой результат булевых операций.

use crate::grid::{GridError, GridPoint};
use crate::region::Region;
use crate::ring::Ring;
use crate::vec::Vec2;
use std::fmt;

/// Почему поле плотности не годится.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum DensityError {
    /// Нет ни одного элемента.
    Empty,
    /// Размер элемента не положительное конечное число.
    BadCell,
    /// Число значений не равно nx · ny.
    SizeMismatch {
        /// Сколько нужно.
        expected: usize,
        /// Сколько дали.
        got: usize,
    },
    /// Значение элемента — не число или бесконечность.
    NotFinite {
        /// Номер по x.
        i: usize,
        /// Номер по y.
        j: usize,
    },
    /// Контур не помещается на 2D-сетку.
    Grid(GridError),
}

impl fmt::Display for DensityError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Empty => write!(f, "в поле плотности нет элементов"),
            Self::BadCell => write!(f, "размер элемента должен быть положительным числом"),
            Self::SizeMismatch { expected, got } => {
                write!(f, "нужно {expected} значений плотности, а дано {got}")
            }
            Self::NotFinite { i, j } => write!(f, "плотность элемента ({i}, {j}) — не число"),
            Self::Grid(e) => write!(f, "{e}"),
        }
    }
}

impl std::error::Error for DensityError {}

/// Параметры построения контура.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ContourOptions {
    /// Порог плотности: где не меньше — материал. Обычно 0,5.
    pub threshold: f64,
    /// Во сколько раз уплотнить поле перед построением (1 — как есть).
    pub refine: usize,
    /// Острова меньше этой площади отбрасываются, мм².
    pub min_island_mm2: f64,
    /// Дыры меньше этой площади заполняются, мм².
    pub min_hole_mm2: f64,
}

impl Default for ContourOptions {
    fn default() -> Self {
        Self {
            threshold: 0.5,
            refine: 1,
            min_island_mm2: 0.0,
            min_hole_mm2: 0.0,
        }
    }
}

/// Поле плотности: `nx × ny` элементов размера `cell_mm`, левый нижний угол — `origin`.
/// Элемент (i, j) — i-й по x слева, j-й по y снизу.
#[derive(Clone, Debug, PartialEq)]
pub struct DensityField {
    nx: usize,
    ny: usize,
    cell: f64,
    origin: Vec2,
    values: Vec<f64>,
}

/// Сколько отсчётов добавить с каждой стороны: два повторяющих край (чтобы
/// интерполяция внутри области не видела «пустоты» снаружи) и один пустой (чтобы
/// контуры замыкались).
///
/// РИСК: Catmull-Rom берёт по два отсчёта с каждой стороны; с одним повторяющим
/// рядом «пустота» просочится внутрь области на полклетки от края.
const PAD: usize = 3;

/// Веса Catmull-Rom для параметра t ∈ [0, 1] между отсчётами 1 и 2 из четырёх.
fn catmull_rom(t: f64) -> [f64; 4] {
    let (t2, t3) = (t * t, t * t * t);
    [
        0.5 * (-t3 + 2.0 * t2 - t),
        0.5 * (3.0 * t3 - 5.0 * t2 + 2.0),
        0.5 * (-3.0 * t3 + 4.0 * t2 + t),
        0.5 * (t3 - t2),
    ]
}

impl DensityField {
    /// Поле из значений по строкам: сначала все i при j = 0 (нижний ряд), потом j = 1…
    pub fn new(
        nx: usize,
        ny: usize,
        cell_mm: f64,
        origin: Vec2,
        values: Vec<f64>,
    ) -> Result<Self, DensityError> {
        if nx == 0 || ny == 0 {
            return Err(DensityError::Empty);
        }
        if !(cell_mm.is_finite() && cell_mm > 0.0) || !origin.is_finite() {
            return Err(DensityError::BadCell);
        }
        if values.len() != nx * ny {
            return Err(DensityError::SizeMismatch {
                expected: nx * ny,
                got: values.len(),
            });
        }
        if let Some(k) = values.iter().position(|v| !v.is_finite()) {
            return Err(DensityError::NotFinite {
                i: k % nx,
                j: k / nx,
            });
        }
        Ok(Self {
            nx,
            ny,
            cell: cell_mm,
            origin,
            values,
        })
    }

    /// Поле из функции `f(i, j)`.
    pub fn from_fn(
        nx: usize,
        ny: usize,
        cell_mm: f64,
        origin: Vec2,
        f: impl Fn(usize, usize) -> f64,
    ) -> Result<Self, DensityError> {
        let values = (0..ny)
            .flat_map(|j| (0..nx).map(move |i| (i, j)))
            .map(|(i, j)| f(i, j))
            .collect();
        Self::new(nx, ny, cell_mm, origin, values)
    }

    /// Число элементов по x и по y.
    pub fn size(&self) -> (usize, usize) {
        (self.nx, self.ny)
    }

    /// Размер элемента, мм.
    pub fn cell_mm(&self) -> f64 {
        self.cell
    }

    /// Плотность элемента (i, j).
    pub fn get(&self, i: usize, j: usize) -> f64 {
        self.values[i + j * self.nx]
    }

    /// Центр элемента (i, j), мм.
    pub fn center(&self, i: usize, j: usize) -> Vec2 {
        self.origin + Vec2::new((i as f64 + 0.5) * self.cell, (j as f64 + 0.5) * self.cell)
    }

    /// Расчётная область: (левый нижний, правый верхний угол), мм.
    pub fn domain(&self) -> (Vec2, Vec2) {
        let size = Vec2::new(self.nx as f64 * self.cell, self.ny as f64 * self.cell);
        (self.origin, self.origin + size)
    }

    /// Отсчёт расширенной сетки: индексы со сдвигом на PAD; снаружи — «пусто».
    fn padded(&self, a: usize, b: usize, empty: f64) -> f64 {
        let (w, h) = (self.nx + 2 * PAD, self.ny + 2 * PAD);
        if a == 0 || b == 0 || a == w - 1 || b == h - 1 {
            return empty;
        }
        let i = a.saturating_sub(PAD).min(self.nx - 1);
        let j = b.saturating_sub(PAD).min(self.ny - 1);
        self.get(i, j)
    }

    /// Область, где плотность не меньше порога.
    pub fn contour(&self, opts: &ContourOptions) -> Result<Region, DensityError> {
        let t = opts.threshold;
        let r = opts.refine.max(1);
        let empty = t.min(0.0) - 1.0;
        let (cw, ch) = (self.nx + 2 * PAD, self.ny + 2 * PAD);
        let coarse: Vec<f64> = (0..ch)
            .flat_map(|b| (0..cw).map(move |a| (a, b)))
            .map(|(a, b)| self.padded(a, b, empty))
            .collect();
        // Уплотнённая сетка отсчётов: шаг cell / r.
        let (w, h) = ((cw - 1) * r + 1, (ch - 1) * r + 1);
        let at = |a: isize, b: isize| -> f64 {
            let a = a.clamp(0, cw as isize - 1) as usize;
            let b = b.clamp(0, ch as isize - 1) as usize;
            coarse[a + b * cw]
        };
        let mut v = vec![0.0; w * h];
        for fb in 0..h {
            let (bb, tb) = ((fb / r) as isize, (fb % r) as f64 / r as f64);
            let wy = catmull_rom(tb);
            for fa in 0..w {
                let (ba, ta) = ((fa / r) as isize, (fa % r) as f64 / r as f64);
                v[fa + fb * w] = if fa % r == 0 && fb % r == 0 {
                    at(ba, bb)
                } else {
                    let wx = catmull_rom(ta);
                    let mut s = 0.0;
                    for (dy, wyk) in wy.iter().enumerate() {
                        for (dx, wxk) in wx.iter().enumerate() {
                            s += wyk * wxk * at(ba + dx as isize - 1, bb + dy as isize - 1);
                        }
                    }
                    s
                };
            }
        }
        // Положение отсчёта (fa, fb), мм: отсчёт (PAD·r, PAD·r) — центр элемента (0, 0).
        let step = self.cell / r as f64;
        let base = self.center(0, 0) - Vec2::new(PAD as f64 * self.cell, PAD as f64 * self.cell);
        let pos = |fa: usize, fb: usize| base + Vec2::new(fa as f64 * step, fb as f64 * step);
        let inside = |k: usize| v[k] >= t;

        // Рёбра сетки: горизонтальное (a, b)–(a+1, b) — 2k, вертикальное (a, b)–(a, b+1) — 2k+1.
        let crossing = |edge: usize| -> Vec2 {
            let k = edge / 2;
            let (a, b) = (k % w, k / w);
            let (p, q, kq) = if edge.is_multiple_of(2) {
                (pos(a, b), pos(a + 1, b), k + 1)
            } else {
                (pos(a, b), pos(a, b + 1), k + w)
            };
            let (va, vb) = (v[k], v[kq]);
            let s = if va == vb {
                0.5
            } else {
                ((t - va) / (vb - va)).clamp(0.0, 1.0)
            };
            p + (q - p) * s
        };
        let none = usize::MAX;
        let mut next = vec![none; 2 * w * h];
        let mut starts = Vec::new();
        for b in 0..h - 1 {
            for a in 0..w - 1 {
                let k = a + b * w;
                // Углы против часовой: bl, br, tr, tl; рёбра: низ, право, верх, лево.
                let corners = [k, k + 1, k + 1 + w, k + w];
                let edges = [2 * k, 2 * (k + 1) + 1, 2 * (k + w), 2 * k + 1];
                let ins = corners.map(inside);
                // Начало — где по обходу внутри → снаружи, конец — снаружи → внутри.
                let mut st = [false; 4];
                let mut en = [false; 4];
                let mut count = 0;
                for e in 0..4 {
                    let (c0, c1) = (ins[e], ins[(e + 1) % 4]);
                    st[e] = c0 && !c1;
                    en[e] = !c0 && c1;
                    count += usize::from(c0 != c1);
                }
                if count == 0 {
                    continue;
                }
                // Седло: центр внутри — начало соединяется со следующим концом по обходу,
                // снаружи — с предыдущим. Для двух пересечений оба правила совпадают.
                let center_inside = if count == 4 {
                    let [bl, br, tr, tl] = corners.map(|c| v[c]);
                    let den = bl + tr - br - tl;
                    let saddle = (bl * tr - br * tl) / den;
                    saddle >= t
                } else {
                    true
                };
                for e in 0..4 {
                    if !st[e] {
                        continue;
                    }
                    let mut f = e;
                    loop {
                        f = if center_inside {
                            (f + 1) % 4
                        } else {
                            (f + 3) % 4
                        };
                        if en[f] {
                            break;
                        }
                    }
                    next[edges[e]] = edges[f];
                    starts.push(edges[e]);
                }
            }
        }
        // Обход: каждое пересечение — начало ровно одного отрезка и конец ровно одного.
        let mut rings = Vec::new();
        let mut used = vec![false; next.len()];
        for &s in &starts {
            if used[s] {
                continue;
            }
            let mut pts = Vec::new();
            let mut e = s;
            while !used[e] {
                used[e] = true;
                let p = crossing(e);
                pts.push(GridPoint::from_mm(p.x, p.y).map_err(DensityError::Grid)?);
                e = next[e];
                debug_assert!(e != none, "контур marching squares должен замыкаться");
            }
            if let Ok(ring) = Ring::new(pts) {
                rings.push(ring);
            }
        }
        let (lo, hi) = self.domain();
        let g = |p: Vec2| GridPoint::from_mm(p.x, p.y).map_err(DensityError::Grid);
        let (glo, ghi) = (g(lo)?, g(hi)?);
        let domain = match Ring::rect(glo.x(), glo.y(), ghi.x(), ghi.y()) {
            Ok(ring) => Region::from_ring(ring),
            Err(_) => return Ok(Region::empty()),
        };
        let shape = Region::from_rings_nonzero(&rings).intersection(&domain);
        Ok(shape.drop_small(opts.min_island_mm2, opts.min_hole_mm2))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn catmull_rom_weights() {
        for t in [0.0, 0.25, 0.5, 1.0] {
            let w = catmull_rom(t);
            assert!((w.iter().sum::<f64>() - 1.0).abs() < 1e-15, "сумма весов 1");
            // Линейная функция воспроизводится точно.
            let lin: f64 = w.iter().enumerate().map(|(k, w)| w * k as f64).sum();
            assert!((lin - (1.0 + t)).abs() < 1e-15);
        }
        assert_eq!(catmull_rom(0.0), [0.0, 1.0, 0.0, 0.0]);
    }

    #[test]
    fn input_checks() {
        let o = Vec2::ZERO;
        assert_eq!(
            DensityField::new(0, 3, 1.0, o, vec![]),
            Err(DensityError::Empty)
        );
        assert_eq!(
            DensityField::new(2, 2, 0.0, o, vec![0.0; 4]),
            Err(DensityError::BadCell)
        );
        assert_eq!(
            DensityField::new(2, 2, 1.0, o, vec![0.0; 3]),
            Err(DensityError::SizeMismatch {
                expected: 4,
                got: 3
            })
        );
        assert_eq!(
            DensityField::new(2, 2, 1.0, o, vec![0.0, 0.0, f64::NAN, 0.0]),
            Err(DensityError::NotFinite { i: 0, j: 1 })
        );
    }
}
