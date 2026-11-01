//! Точные предикаты: знак ориентации без ошибок округления.
//!
//! Для `f64` сначала идёт быстрый расчёт с оценкой погрешности (фильтр
//! Шевчука). Если знак не гарантирован — точки почти вырождены, — считаем
//! точно в своей длинной арифметике. Ответ всегда точный для тех чисел,
//! что пришли на вход: никаких допусков.
//!
//! На 2D-сетке (целые микрометры) считаем сразу точно в `i128`.

use crate::bigint::BigInt;
use crate::exact::to_common_integers;
use crate::grid::GridPoint;
use crate::trace::{self, Predicate};
use crate::vec::{Vec2, Vec3};
use std::cmp::Ordering;

/// Знак результата предиката.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum Sign {
    /// Меньше нуля.
    Negative,
    /// Ровно ноль: вырожденный случай.
    Zero,
    /// Больше нуля.
    Positive,
}

impl Sign {
    /// Знак целого числа.
    pub fn of_i128(v: i128) -> Self {
        match v.cmp(&0) {
            Ordering::Less => Self::Negative,
            Ordering::Equal => Self::Zero,
            Ordering::Greater => Self::Positive,
        }
    }

    /// Знак длинного целого.
    pub fn of_bigint(v: &BigInt) -> Self {
        Self::of_i128(i128::from(v.signum()))
    }

    /// Знак числа `f64`; ±0 и нечисло дают `Zero`.
    pub fn of_f64(v: f64) -> Self {
        match v.partial_cmp(&0.0) {
            Some(Ordering::Less) => Self::Negative,
            Some(Ordering::Greater) => Self::Positive,
            _ => Self::Zero,
        }
    }

    /// −1, 0 или +1.
    pub fn to_i32(self) -> i32 {
        match self {
            Self::Negative => -1,
            Self::Zero => 0,
            Self::Positive => 1,
        }
    }

    /// Противоположный знак.
    pub fn flip(self) -> Self {
        match self {
            Self::Negative => Self::Positive,
            Self::Zero => Self::Zero,
            Self::Positive => Self::Negative,
        }
    }
}

/// Удвоенная ориентированная площадь треугольника `abc` на 2D-сетке, точно.
///
/// Больше нуля, если обход a → b → c идёт против часовой стрелки
/// (точка c левее направленной прямой a → b); ноль — точки на одной прямой.
pub fn orient2d_grid_value(a: GridPoint, b: GridPoint, c: GridPoint) -> i128 {
    // РИСК: точность держится на пределе сетки GRID_LIMIT = 2^40 (см. tol.rs):
    // разности до 2^41, произведения до 2^82 — до переполнения i128 далеко.
    let (ax, ay) = (i128::from(a.x()), i128::from(a.y()));
    let (bx, by) = (i128::from(b.x()), i128::from(b.y()));
    let (cx, cy) = (i128::from(c.x()), i128::from(c.y()));
    (bx - ax) * (cy - ay) - (by - ay) * (cx - ax)
}

/// Ориентация трёх точек 2D-сетки: `Positive` — против часовой стрелки.
pub fn orient2d_grid(a: GridPoint, b: GridPoint, c: GridPoint) -> Sign {
    Sign::of_i128(orient2d_grid_value(a, b, c))
}

/// Половина машинного эпсилона, 2^−53 (так её определяет Шевчук).
const EPS: f64 = f64::EPSILON / 2.0;
/// Оценка погрешности быстрого orient2d (у Шевчука — ccwerrboundA).
const ORIENT2D_BOUND: f64 = (3.0 + 16.0 * EPS) * EPS;
/// Оценка погрешности быстрого orient3d (у Шевчука — o3derrboundA).
const ORIENT3D_BOUND: f64 = (7.0 + 56.0 * EPS) * EPS;

const fn pow2(e: i32) -> f64 {
    f64::from_bits(((1023 + e) as u64) << 52)
}

/// Разности координат, при которых фильтру можно верить.
///
/// РИСК: оценки погрешности Шевчука верны только без переполнения и без
/// денормализованных чисел. Для 2D (выражение второй степени) берём
/// |d| ∈ [2^−450, 2^450], для 3D (третьей степени) — [2^−300, 2^300]:
/// все произведения остаются в [2^−900, 2^900]. Вне диапазона — сразу
/// точный путь.
const RANGE_2D: (f64, f64) = (pow2(-450), pow2(450));
const RANGE_3D: (f64, f64) = (pow2(-300), pow2(300));

/// Все разности — нули или в безопасном диапазоне (нечисло сюда не проходит).
fn filter_applies(diffs: &[f64], (lo, hi): (f64, f64)) -> bool {
    diffs
        .iter()
        .all(|&d| d == 0.0 || (lo..=hi).contains(&d.abs()))
}

/// Ориентация трёх точек на плоскости (`f64`), точно.
///
/// `Positive` — обход a → b → c против часовой стрелки, `Zero` — точки
/// на одной прямой. Координаты должны быть конечными числами.
pub fn orient2d(a: Vec2, b: Vec2, c: Vec2) -> Sign {
    let (acx, bcx) = (a.x - c.x, b.x - c.x);
    let (acy, bcy) = (a.y - c.y, b.y - c.y);
    if filter_applies(&[acx, bcx, acy, bcy], RANGE_2D) {
        let left = acx * bcy;
        let right = acy * bcx;
        let det = left - right;
        let bound = ORIENT2D_BOUND * (left.abs() + right.abs());
        let fast = if det > bound {
            Some(Sign::Positive)
        } else if -det > bound {
            Some(Sign::Negative)
        } else if bound == 0.0 {
            // Оба произведения ровно нули: в диапазоне фильтра это значит,
            // что ноль и точный ответ.
            Some(Sign::Zero)
        } else {
            None
        };
        if let Some(sign) = fast {
            trace::note_fast(Predicate::Orient2d);
            return sign;
        }
    }
    let sign = orient2d_exact(a, b, c);
    trace::note_exact(Predicate::Orient2d, &[a.x, a.y, b.x, b.y, c.x, c.y], sign);
    sign
}

/// Точный путь [`orient2d`]: длинная арифметика без округлений.
pub fn orient2d_exact(a: Vec2, b: Vec2, c: Vec2) -> Sign {
    Sign::of_bigint(&orient2d_exact_value(a, b, c))
}

/// Точное значение определителя [`orient2d`] в целых числах общего масштаба
/// (см. [`to_common_integers`]): настоящее значение, умноженное на 2^−2s.
/// Знак тот же, что у настоящего; само значение нужно для проверок.
pub fn orient2d_exact_value(a: Vec2, b: Vec2, c: Vec2) -> BigInt {
    let n = to_common_integers(&[a.x, a.y, b.x, b.y, c.x, c.y])
        .expect("orient2d: координата — не число или бесконечность");
    let [ax, ay, bx, by, cx, cy] = n.as_slice() else {
        unreachable!("шесть координат")
    };
    &(ax - cx) * &(by - cy) - &(ay - cy) * &(bx - cx)
}

/// Положение точки `d` относительно плоскости через `a`, `b`, `c` (`f64`), точно.
///
/// Это знак определителя det(a − d, b − d, c − d). `Positive` — точка d
/// лежит ниже плоскости, если смотреть сверху так, чтобы обход a → b → c
/// шёл против часовой стрелки (соглашение Шевчука). `Zero` — все четыре
/// точки в одной плоскости. Координаты должны быть конечными числами.
pub fn orient3d(a: Vec3, b: Vec3, c: Vec3, d: Vec3) -> Sign {
    let (adx, bdx, cdx) = (a.x - d.x, b.x - d.x, c.x - d.x);
    let (ady, bdy, cdy) = (a.y - d.y, b.y - d.y, c.y - d.y);
    let (adz, bdz, cdz) = (a.z - d.z, b.z - d.z, c.z - d.z);
    if filter_applies(&[adx, bdx, cdx, ady, bdy, cdy, adz, bdz, cdz], RANGE_3D) {
        let (bdxcdy, cdxbdy) = (bdx * cdy, cdx * bdy);
        let (cdxady, adxcdy) = (cdx * ady, adx * cdy);
        let (adxbdy, bdxady) = (adx * bdy, bdx * ady);
        let det = adz * (bdxcdy - cdxbdy) + bdz * (cdxady - adxcdy) + cdz * (adxbdy - bdxady);
        let permanent = (bdxcdy.abs() + cdxbdy.abs()) * adz.abs()
            + (cdxady.abs() + adxcdy.abs()) * bdz.abs()
            + (adxbdy.abs() + bdxady.abs()) * cdz.abs();
        let bound = ORIENT3D_BOUND * permanent;
        let fast = if det > bound {
            Some(Sign::Positive)
        } else if -det > bound {
            Some(Sign::Negative)
        } else if bound == 0.0 {
            // Все слагаемые ровно нули: в диапазоне фильтра — точный ноль.
            Some(Sign::Zero)
        } else {
            None
        };
        if let Some(sign) = fast {
            trace::note_fast(Predicate::Orient3d);
            return sign;
        }
    }
    let sign = orient3d_exact(a, b, c, d);
    let inputs = [a.x, a.y, a.z, b.x, b.y, b.z, c.x, c.y, c.z, d.x, d.y, d.z];
    trace::note_exact(Predicate::Orient3d, &inputs, sign);
    sign
}

/// Точный путь [`orient3d`]: длинная арифметика без округлений.
pub fn orient3d_exact(a: Vec3, b: Vec3, c: Vec3, d: Vec3) -> Sign {
    Sign::of_bigint(&orient3d_exact_value(a, b, c, d))
}

/// Точное значение определителя [`orient3d`] в целых числах общего масштаба
/// (см. [`to_common_integers`]): настоящее значение, умноженное на 2^−3s.
/// Знак тот же, что у настоящего; само значение нужно для проверок.
pub fn orient3d_exact_value(a: Vec3, b: Vec3, c: Vec3, d: Vec3) -> BigInt {
    let n = to_common_integers(&[a.x, a.y, a.z, b.x, b.y, b.z, c.x, c.y, c.z, d.x, d.y, d.z])
        .expect("orient3d: координата — не число или бесконечность");
    let [ax, ay, az, bx, by, bz, cx, cy, cz, dx, dy, dz] = n.as_slice() else {
        unreachable!("двенадцать координат")
    };
    let (adx, ady, adz) = (ax - dx, ay - dy, az - dz);
    let (bdx, bdy, bdz) = (bx - dx, by - dy, bz - dz);
    let (cdx, cdy, cdz) = (cx - dx, cy - dy, cz - dz);
    let m1 = &bdx * &cdy - &cdx * &bdy;
    let m2 = &cdx * &ady - &adx * &cdy;
    let m3 = &adx * &bdy - &bdx * &ady;
    &adz * &m1 + &bdz * &m2 + &cdz * &m3
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::exact::decompose;
    use crate::testkit::{Rng, for_cases};
    use crate::tol::GRID_LIMIT;

    // --- Независимая сверка точного пути по модулю простых чисел -------------
    //
    // Точный путь считает определитель длинной арифметикой. Здесь то же
    // значение считается по модулю двух простых чисел — только на u128, без
    // BigInt и без его сдвигов. Совпадение остатков по обоим модулям значит,
    // что длинное значение верно вместе со знаком: случайно совпасть оба
    // остатка могут с вероятностью около 2^−91.

    const PRIMES: [u64; 2] = [(1 << 61) - 1, 1_000_000_007];

    fn mul_mod(a: u64, b: u64, p: u64) -> u64 {
        ((u128::from(a) * u128::from(b)) % u128::from(p)) as u64
    }

    fn sub_mod(a: u64, b: u64, p: u64) -> u64 {
        ((u128::from(a) + u128::from(p) - u128::from(b)) % u128::from(p)) as u64
    }

    fn pow2_mod(mut e: u64, p: u64) -> u64 {
        let (mut base, mut r) = (2 % p, 1 % p);
        while e > 0 {
            if e & 1 == 1 {
                r = mul_mod(r, base, p);
            }
            base = mul_mod(base, base, p);
            e >>= 1;
        }
        r
    }

    /// Координаты по модулю `p` в том же общем масштабе, что у точного пути.
    fn residues(xs: &[f64], p: u64) -> Vec<u64> {
        let parts: Vec<_> = xs.iter().map(|&x| decompose(x).unwrap()).collect();
        let scale = parts
            .iter()
            .filter(|q| q.1 != 0)
            .map(|q| q.2)
            .min()
            .unwrap_or(0);
        parts
            .iter()
            .map(|&(neg, m, e)| {
                let r = mul_mod(m % p, pow2_mod((e - scale) as u64, p), p);
                if neg { sub_mod(0, r, p) } else { r }
            })
            .collect()
    }

    #[test]
    fn orient2d_exact_value_passes_modular_check() {
        for_cases("orient2d_exact_value_passes_modular_check", 20_000, |rng| {
            let [a, b, c] = triangle(rng);
            let value = orient2d_exact_value(a, b, c);
            for p in PRIMES {
                let r = residues(&[a.x, a.y, b.x, b.y, c.x, c.y], p);
                let (sub, mul) = (|x, y| sub_mod(x, y, p), |x, y| mul_mod(x, y, p));
                let want = sub(
                    mul(sub(r[0], r[4]), sub(r[3], r[5])),
                    mul(sub(r[1], r[5]), sub(r[2], r[4])),
                );
                assert_eq!(value.rem_u64(p), want, "p={p} a={a:?} b={b:?} c={c:?}");
            }
        });
    }

    #[test]
    fn orient3d_exact_value_passes_modular_check() {
        for_cases("orient3d_exact_value_passes_modular_check", 20_000, |rng| {
            let [a, b, c, d] = tetrahedron(rng);
            let value = orient3d_exact_value(a, b, c, d);
            let xs = [a.x, a.y, a.z, b.x, b.y, b.z, c.x, c.y, c.z, d.x, d.y, d.z];
            for p in PRIMES {
                let r = residues(&xs, p);
                let (sub, mul) = (|x, y| sub_mod(x, y, p), |x, y| mul_mod(x, y, p));
                let add = |x: u64, y: u64| ((u128::from(x) + u128::from(y)) % u128::from(p)) as u64;
                // Разложение по первой строке — не тот порядок, что в точном пути.
                let row = |i: usize| [sub(r[i], r[9]), sub(r[i + 1], r[10]), sub(r[i + 2], r[11])];
                let ([ax, ay, az], [bx, by, bz], [cx, cy, cz]) = (row(0), row(3), row(6));
                let want = add(
                    sub(
                        mul(ax, sub(mul(by, cz), mul(bz, cy))),
                        mul(ay, sub(mul(bx, cz), mul(bz, cx))),
                    ),
                    mul(az, sub(mul(bx, cy), mul(by, cx))),
                );
                assert_eq!(
                    value.rem_u64(p),
                    want,
                    "p={p} a={a:?} b={b:?} c={c:?} d={d:?}"
                );
            }
        });
    }

    /// Случайный треугольник: треть — общего вида, треть ровно вырождена,
    /// треть почти вырождена (третья точка в паре ulp от прямой).
    fn triangle(rng: &mut Rng) -> [Vec2; 3] {
        loop {
            let a = rng.geom_vec2();
            let b = rng.geom_vec2();
            let c = match rng.below(3) {
                0 => rng.geom_vec2(),
                1 => a + (b - a) * rng.dyadic_t(),
                _ => {
                    let t = rng.unit_f64();
                    rng.nudge2(a + (b - a) * t)
                }
            };
            if c.is_finite() {
                return [a, b, c];
            }
        }
    }

    /// Случайный тетраэдр: общий, ровно плоский или почти плоский.
    fn tetrahedron(rng: &mut Rng) -> [Vec3; 4] {
        loop {
            let a = rng.geom_vec3();
            let b = rng.geom_vec3();
            let c = rng.geom_vec3();
            let d = match rng.below(3) {
                0 => rng.geom_vec3(),
                1 => a + (b - a) * rng.dyadic_t() + (c - a) * rng.dyadic_t(),
                _ => {
                    let (s, t) = (rng.unit_f64(), rng.unit_f64());
                    rng.nudge3(a + (b - a) * s + (c - a) * t)
                }
            };
            if d.is_finite() {
                return [a, b, c, d];
            }
        }
    }

    #[test]
    fn orient2d_fast_path_agrees_with_exact() {
        for_cases("orient2d_fast_path_agrees_with_exact", 20_000, |rng| {
            let [a, b, c] = triangle(rng);
            assert_eq!(
                orient2d(a, b, c),
                orient2d_exact(a, b, c),
                "a={a:?} b={b:?} c={c:?}"
            );
        });
    }

    #[test]
    fn orient2d_respects_symmetries() {
        for_cases("orient2d_respects_symmetries", 20_000, |rng| {
            let [a, b, c] = triangle(rng);
            let s = orient2d(a, b, c);
            assert_eq!(orient2d(b, c, a), s, "a={a:?} b={b:?} c={c:?}");
            assert_eq!(orient2d(c, a, b), s);
            assert_eq!(orient2d(b, a, c), s.flip());
            assert_eq!(orient2d(a, c, b), s.flip());
            assert_eq!(orient2d(c, b, a), s.flip());
        });
    }

    /// Независимая сверка без нашей длинной арифметики: все координаты —
    /// целые с общим множителем 2^k, точный ответ считается прямо в i128.
    #[test]
    fn orient2d_matches_i128_on_scaled_integers() {
        for_cases("orient2d_matches_i128_on_scaled_integers", 20_000, |rng| {
            let k = rng.range_i64(-40, 40) as i32;
            let span = if rng.chance(0.5) { 8 } else { 1 << 30 };
            let mut w: [i64; 6] = std::array::from_fn(|_| rng.range_i64(-span, span));
            if rng.chance(0.3) {
                // c = a + t·(b − a): ровно на прямой.
                let t = rng.range_i64(-3, 3);
                w[4] = w[0] + t * (w[2] - w[0]);
                w[5] = w[1] + t * (w[3] - w[1]);
            }
            let scale = 2f64.powi(k);
            let p = |i: usize| Vec2::new(w[i] as f64 * scale, w[i + 1] as f64 * scale);
            let v = |i: usize| i128::from(w[i]);
            let want = (v(0) - v(4)) * (v(3) - v(5)) - (v(1) - v(5)) * (v(2) - v(4));
            assert_eq!(
                orient2d(p(0), p(2), p(4)),
                Sign::of_i128(want),
                "w={w:?} k={k}"
            );
        });
    }

    #[test]
    fn orient3d_fast_path_agrees_with_exact() {
        for_cases("orient3d_fast_path_agrees_with_exact", 20_000, |rng| {
            let [a, b, c, d] = tetrahedron(rng);
            assert_eq!(
                orient3d(a, b, c, d),
                orient3d_exact(a, b, c, d),
                "a={a:?} b={b:?} c={c:?} d={d:?}"
            );
        });
    }

    #[test]
    fn orient3d_respects_symmetries() {
        for_cases("orient3d_respects_symmetries", 20_000, |rng| {
            let [a, b, c, d] = tetrahedron(rng);
            let s = orient3d(a, b, c, d);
            // Чётные перестановки знак сохраняют, нечётные — меняют.
            assert_eq!(orient3d(b, c, a, d), s, "a={a:?} b={b:?} c={c:?} d={d:?}");
            assert_eq!(orient3d(d, c, b, a), s);
            assert_eq!(orient3d(b, a, c, d), s.flip());
            assert_eq!(orient3d(a, b, d, c), s.flip());
            assert_eq!(orient3d(d, b, c, a), s.flip());
        });
    }

    #[test]
    fn orient3d_matches_i128_on_scaled_integers() {
        for_cases("orient3d_matches_i128_on_scaled_integers", 20_000, |rng| {
            let k = rng.range_i64(-30, 30) as i32;
            let span = if rng.chance(0.5) { 8 } else { 1 << 20 };
            let mut w: [i64; 12] = std::array::from_fn(|_| rng.range_i64(-span, span));
            if rng.chance(0.3) {
                // d = a + s·(b − a) + t·(c − a): ровно в плоскости.
                let (s, t) = (rng.range_i64(-2, 2), rng.range_i64(-2, 2));
                for axis in 0..3 {
                    w[9 + axis] =
                        w[axis] + s * (w[3 + axis] - w[axis]) + t * (w[6 + axis] - w[axis]);
                }
            }
            let scale = 2f64.powi(k);
            let p = |i: usize| {
                Vec3::new(
                    w[i] as f64 * scale,
                    w[i + 1] as f64 * scale,
                    w[i + 2] as f64 * scale,
                )
            };
            let v = |i: usize| i128::from(w[i]);
            let (adx, ady, adz) = (v(0) - v(9), v(1) - v(10), v(2) - v(11));
            let (bdx, bdy, bdz) = (v(3) - v(9), v(4) - v(10), v(5) - v(11));
            let (cdx, cdy, cdz) = (v(6) - v(9), v(7) - v(10), v(8) - v(11));
            let want = adx * (bdy * cdz - bdz * cdy) - ady * (bdx * cdz - bdz * cdx)
                + adz * (bdx * cdy - bdy * cdx);
            assert_eq!(
                orient3d(p(0), p(3), p(6), p(9)),
                Sign::of_i128(want),
                "w={w:?} k={k}"
            );
        });
    }

    fn grid_coord(rng: &mut Rng) -> i64 {
        if rng.chance(0.3) {
            *rng.pick(&[-GRID_LIMIT, GRID_LIMIT, 0, 1, -1])
        } else {
            rng.range_i64(-GRID_LIMIT, GRID_LIMIT)
        }
    }

    #[test]
    fn orient2d_grid_is_exact_even_at_the_limits() {
        for_cases("orient2d_grid_is_exact_even_at_the_limits", 20_000, |rng| {
            let p: [GridPoint; 3] =
                std::array::from_fn(|_| GridPoint::new(grid_coord(rng), grid_coord(rng)).unwrap());
            let big = |q: GridPoint| {
                (
                    BigInt::from_i128(q.x().into()),
                    BigInt::from_i128(q.y().into()),
                )
            };
            let ((ax, ay), (bx, by), (cx, cy)) = (big(p[0]), big(p[1]), big(p[2]));
            let want = &(&bx - &ax) * &(&cy - &ay) - &(&by - &ay) * &(&cx - &ax);
            assert_eq!(
                BigInt::from_i128(orient2d_grid_value(p[0], p[1], p[2])),
                want
            );
        });
    }

    /// Тесты чего-то стоят, только если генераторы часто дают трудные
    /// случаи: проверяем, что точный путь срабатывает заметную долю раз.
    #[test]
    fn generators_often_hit_the_exact_path() {
        let mut rng = Rng::new(42);
        trace::reset_counters();
        for _ in 0..4000 {
            let [a, b, c] = triangle(&mut rng);
            orient2d(a, b, c);
            let [a, b, c, d] = tetrahedron(&mut rng);
            orient3d(a, b, c, d);
        }
        let k = trace::counters();
        let share2 = k.orient2d_exact as f64 / 4000.0;
        let share3 = k.orient3d_exact as f64 / 4000.0;
        assert!(
            (0.2..0.9).contains(&share2),
            "2D: точный путь в {share2:.2} случаев"
        );
        assert!(
            (0.2..0.9).contains(&share3),
            "3D: точный путь в {share3:.2} случаев"
        );
    }

    #[test]
    fn conventions_on_simple_shapes() {
        let o = Vec2::ZERO;
        assert_eq!(
            orient2d(o, Vec2::new(1.0, 0.0), Vec2::new(0.0, 1.0)),
            Sign::Positive
        );
        assert_eq!(
            orient2d(o, Vec2::new(0.0, 1.0), Vec2::new(1.0, 0.0)),
            Sign::Negative
        );
        let (a, b, c) = (Vec3::ZERO, Vec3::X, Vec3::Y);
        // a → b → c против часовой, если смотреть сверху (с +Z).
        assert_eq!(orient3d(a, b, c, Vec3::Z), Sign::Negative, "над плоскостью");
        assert_eq!(
            orient3d(a, b, c, -Vec3::Z),
            Sign::Positive,
            "под плоскостью"
        );
        assert_eq!(orient3d(a, b, c, Vec3::new(5.0, -7.0, 0.0)), Sign::Zero);
        let g = |x, y| GridPoint::new(x, y).unwrap();
        assert_eq!(orient2d_grid(g(0, 0), g(10, 0), g(3, 5)), Sign::Positive);
        assert_eq!(orient2d_grid(g(0, 0), g(10, 10), g(-4, -4)), Sign::Zero);
    }

    #[test]
    fn degenerate_input_takes_exact_path_and_is_logged() {
        trace::reset_counters();
        trace::start_capture();
        let collinear = [Vec2::ZERO, Vec2::new(1.0, 1.0), Vec2::new(3.0, 3.0)];
        assert_eq!(
            orient2d(collinear[0], collinear[1], collinear[2]),
            Sign::Zero
        );
        assert_eq!(
            orient2d(Vec2::ZERO, Vec2::new(1.0, 0.0), Vec2::new(0.0, 1.0)),
            Sign::Positive
        );
        let c = trace::counters();
        assert_eq!((c.orient2d_fast, c.orient2d_exact), (1, 1));
        let log = trace::take_captured();
        assert_eq!(log.len(), 1);
        assert_eq!(log[0].predicate, Predicate::Orient2d);
        assert_eq!(log[0].inputs, vec![0.0, 0.0, 1.0, 1.0, 3.0, 3.0]);
        assert_eq!(log[0].sign, Sign::Zero);
    }

    #[test]
    fn extreme_magnitudes() {
        let (big, tiny) = (1e300, f64::from_bits(1));
        let (a, b) = (Vec2::new(-big, -big), Vec2::new(big, big));
        // Разности переполняют фильтр — ответ всё равно точный.
        assert_eq!(orient2d(a, b, Vec2::ZERO), Sign::Zero);
        assert_eq!(orient2d(a, b, Vec2::new(0.0, tiny)), Sign::Positive);
        assert_eq!(orient2d(a, b, Vec2::new(tiny, 0.0)), Sign::Negative);
        let m = f64::MAX;
        let up = Vec3::new(0.0, 0.0, tiny);
        assert_eq!(
            orient3d(
                Vec3::new(-m, 0.0, 0.0),
                Vec3::new(m, 0.0, 0.0),
                Vec3::new(0.0, m, 0.0),
                up
            ),
            Sign::Negative
        );
    }
}
