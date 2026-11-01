//! Своя длинная арифметика: целые числа любой длины со знаком.
//!
//! Нужна для запасного точного пути предикатов. Он срабатывает редко —
//! только когда быстрый расчёт в `f64` не может гарантировать знак, —
//! поэтому здесь важнее простота и очевидная правильность, чем скорость.

use std::cmp::Ordering;
use std::ops::{Add, Mul, Neg, Sub};

/// Целое число любой длины со знаком.
///
/// Модуль хранится словами по 64 бита, младшее слово первым.
/// Инвариант: старших нулевых слов нет; у нуля пустой модуль и знак «плюс».
#[derive(Clone, Debug, Default, PartialEq, Eq, Hash)]
pub struct BigInt {
    neg: bool,
    mag: Vec<u64>,
}

impl BigInt {
    /// Ноль.
    pub fn zero() -> Self {
        Self::default()
    }

    /// Число из `i128`.
    pub fn from_i128(v: i128) -> Self {
        let m = v.unsigned_abs();
        Self::from_limbs(v < 0, vec![m as u64, (m >> 64) as u64])
    }

    /// Число (−1)^`neg` · `m` · 2^`shift`.
    pub fn from_shifted(neg: bool, m: u64, shift: u32) -> Self {
        let words = (shift / 64) as usize;
        let bits = shift % 64;
        let mut mag = vec![0u64; words + 2];
        mag[words] = m << bits;
        if bits != 0 {
            mag[words + 1] = m >> (64 - bits);
        }
        Self::from_limbs(neg, mag)
    }

    /// Число из знака и слов модуля (младшее слово первым).
    pub fn from_limbs(neg: bool, mut limbs: Vec<u64>) -> Self {
        trim(&mut limbs);
        Self {
            neg: neg && !limbs.is_empty(),
            mag: limbs,
        }
    }

    /// Ноль ли это.
    pub fn is_zero(&self) -> bool {
        self.mag.is_empty()
    }

    /// Знак: −1, 0 или +1.
    pub fn signum(&self) -> i32 {
        match (self.mag.is_empty(), self.neg) {
            (true, _) => 0,
            (false, true) => -1,
            (false, false) => 1,
        }
    }

    /// Число значащих бит модуля (у нуля — 0).
    pub fn bit_len(&self) -> u64 {
        match self.mag.last() {
            None => 0,
            Some(&top) => (self.mag.len() as u64 - 1) * 64 + u64::from(64 - top.leading_zeros()),
        }
    }

    /// Значение как `i128`, если оно туда помещается.
    pub fn to_i128(&self) -> Option<i128> {
        if self.mag.len() > 2 {
            return None;
        }
        let lo = u128::from(self.mag.first().copied().unwrap_or(0));
        let hi = u128::from(self.mag.get(1).copied().unwrap_or(0));
        let m = lo | (hi << 64);
        if self.neg {
            // Модуль 2^127 — это ровно i128::MIN.
            (m <= 1u128 << 127).then(|| (m as i128).wrapping_neg())
        } else {
            i128::try_from(m).ok()
        }
    }

    /// Остаток от деления на `p > 0`, всегда в [0, p) — и для отрицательных чисел.
    ///
    /// Нужен независимой проверке: значение, посчитанное длинной арифметикой,
    /// сверяется с тем же выражением, посчитанным по модулю простого числа.
    pub fn rem_u64(&self, p: u64) -> u64 {
        assert!(p > 0, "rem_u64(0)");
        let p128 = u128::from(p);
        let r = self
            .mag
            .iter()
            .rev()
            .fold(0u128, |r, &limb| ((r << 64) | u128::from(limb)) % p128) as u64;
        if self.neg && r != 0 { p - r } else { r }
    }
}

fn trim(mag: &mut Vec<u64>) {
    while mag.last() == Some(&0) {
        mag.pop();
    }
}

fn cmp_mag(a: &[u64], b: &[u64]) -> Ordering {
    a.len()
        .cmp(&b.len())
        .then_with(|| a.iter().rev().cmp(b.iter().rev()))
}

fn add_mag(a: &[u64], b: &[u64]) -> Vec<u64> {
    let (long, short) = if a.len() >= b.len() { (a, b) } else { (b, a) };
    let mut out = Vec::with_capacity(long.len() + 1);
    let mut carry = false;
    for (i, &x) in long.iter().enumerate() {
        let (s1, c1) = x.overflowing_add(short.get(i).copied().unwrap_or(0));
        let (s2, c2) = s1.overflowing_add(u64::from(carry));
        out.push(s2);
        carry = c1 || c2;
    }
    if carry {
        out.push(1);
    }
    out
}

/// |a| − |b| при |a| ≥ |b|.
fn sub_mag(a: &[u64], b: &[u64]) -> Vec<u64> {
    let mut out = Vec::with_capacity(a.len());
    let mut borrow = false;
    for (i, &x) in a.iter().enumerate() {
        let (d1, b1) = x.overflowing_sub(b.get(i).copied().unwrap_or(0));
        let (d2, b2) = d1.overflowing_sub(u64::from(borrow));
        out.push(d2);
        borrow = b1 || b2;
    }
    debug_assert!(!borrow, "sub_mag: уменьшаемое меньше вычитаемого");
    trim(&mut out);
    out
}

fn mul_mag(a: &[u64], b: &[u64]) -> Vec<u64> {
    if a.is_empty() || b.is_empty() {
        return Vec::new();
    }
    let mut out = vec![0u64; a.len() + b.len()];
    for (i, &x) in a.iter().enumerate() {
        let mut carry = 0u128;
        for (j, &y) in b.iter().enumerate() {
            // Без переполнения: (2^64 − 1)² + 2·(2^64 − 1) = 2^128 − 1.
            let t = u128::from(x) * u128::from(y) + u128::from(out[i + j]) + carry;
            out[i + j] = t as u64;
            carry = t >> 64;
        }
        // Это слово строка i ещё не трогала, а прежние строки до него не доходят.
        out[i + b.len()] = carry as u64;
    }
    trim(&mut out);
    out
}

fn signed_add(an: bool, a: &[u64], bn: bool, b: &[u64]) -> BigInt {
    if an == bn {
        return BigInt::from_limbs(an, add_mag(a, b));
    }
    match cmp_mag(a, b) {
        Ordering::Equal => BigInt::zero(),
        Ordering::Greater => BigInt::from_limbs(an, sub_mag(a, b)),
        Ordering::Less => BigInt::from_limbs(bn, sub_mag(b, a)),
    }
}

impl Add for &BigInt {
    type Output = BigInt;
    fn add(self, o: &BigInt) -> BigInt {
        signed_add(self.neg, &self.mag, o.neg, &o.mag)
    }
}

impl Sub for &BigInt {
    type Output = BigInt;
    fn sub(self, o: &BigInt) -> BigInt {
        signed_add(self.neg, &self.mag, !o.neg, &o.mag)
    }
}

impl Mul for &BigInt {
    type Output = BigInt;
    fn mul(self, o: &BigInt) -> BigInt {
        BigInt::from_limbs(self.neg != o.neg, mul_mag(&self.mag, &o.mag))
    }
}

impl Neg for &BigInt {
    type Output = BigInt;
    fn neg(self) -> BigInt {
        BigInt::from_limbs(!self.neg, self.mag.clone())
    }
}

macro_rules! forward_owned {
    ($tr:ident, $f:ident) => {
        impl $tr for BigInt {
            type Output = BigInt;
            fn $f(self, o: BigInt) -> BigInt {
                (&self).$f(&o)
            }
        }
    };
}

forward_owned!(Add, add);
forward_owned!(Sub, sub);
forward_owned!(Mul, mul);

impl Neg for BigInt {
    type Output = BigInt;
    fn neg(mut self) -> BigInt {
        self.neg = !self.neg && !self.mag.is_empty();
        self
    }
}

impl Ord for BigInt {
    fn cmp(&self, o: &Self) -> Ordering {
        match (self.signum(), o.signum()) {
            (a, b) if a != b => a.cmp(&b),
            (-1, _) => cmp_mag(&o.mag, &self.mag),
            _ => cmp_mag(&self.mag, &o.mag),
        }
    }
}

impl PartialOrd for BigInt {
    fn partial_cmp(&self, o: &Self) -> Option<Ordering> {
        Some(self.cmp(o))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::testkit::{Rng, for_cases};

    fn random_i128(rng: &mut Rng) -> i128 {
        ((u128::from(rng.next_u64()) << 64) | u128::from(rng.next_u64())) as i128
    }

    /// Случайное длинное число: до пяти слов, часто нули и «все единицы» —
    /// на них ломаются переносы.
    fn random_big(rng: &mut Rng) -> BigInt {
        let n = rng.below(6) as usize;
        let limbs = (0..n)
            .map(|_| match rng.below(4) {
                0 => 0,
                1 => u64::MAX,
                _ => rng.next_u64(),
            })
            .collect();
        BigInt::from_limbs(rng.chance(0.5), limbs)
    }

    #[test]
    fn i128_round_trip() {
        for v in [
            0,
            1,
            -1,
            i128::MAX,
            i128::MIN,
            1 << 64,
            -(1 << 64),
            u64::MAX as i128,
        ] {
            assert_eq!(BigInt::from_i128(v).to_i128(), Some(v), "{v}");
        }
        for_cases("i128_round_trip", 5000, |rng| {
            let v = random_i128(rng);
            assert_eq!(BigInt::from_i128(v).to_i128(), Some(v));
        });
        assert_eq!(BigInt::from_shifted(false, 1, 127).to_i128(), None);
        assert_eq!(
            BigInt::from_shifted(true, 1, 127).to_i128(),
            Some(i128::MIN)
        );
    }

    #[test]
    fn operations_match_i128() {
        for_cases("operations_match_i128", 20_000, |rng| {
            let (a, b) = (rng.next_u64() as i64, rng.next_u64() as i64);
            let (x, y) = (BigInt::from_i128(a.into()), BigInt::from_i128(b.into()));
            let (a, b) = (i128::from(a), i128::from(b));
            assert_eq!((&x + &y).to_i128(), Some(a + b));
            assert_eq!((&x - &y).to_i128(), Some(a - b));
            assert_eq!((&x * &y).to_i128(), Some(a * b));
            assert_eq!((-&x).to_i128(), Some(-a));
            assert_eq!(x.cmp(&y), a.cmp(&b));
            assert_eq!(i128::from(x.signum()), a.signum());
        });
    }

    #[test]
    fn algebra_holds_for_long_numbers() {
        for_cases("algebra_holds_for_long_numbers", 5000, |rng| {
            let (a, b, c) = (random_big(rng), random_big(rng), random_big(rng));
            assert_eq!(&a + &b, &b + &a);
            assert_eq!(&a * &b, &b * &a);
            assert_eq!(&(&a + &b) + &c, &a + &(&b + &c));
            assert_eq!(&(&a * &b) * &c, &a * &(&b * &c));
            assert_eq!(&a * &(&b + &c), &a * &b + &a * &c);
            assert_eq!(&(&a + &b) - &b, a);
            assert_eq!(&a - &a, BigInt::zero());
            assert_eq!(-(-a.clone()), a);
            assert!((&a * &a).signum() >= 0);
            assert_eq!(&(&a + &b) * &(&a - &b), &a * &a - &b * &b);
            assert_eq!(a.cmp(&b), (&a - &b).signum().cmp(&0));
        });
    }

    #[test]
    fn shifted_values() {
        for_cases("shifted_values", 20_000, |rng| {
            let m = rng.next_u64() >> 11;
            let s = rng.below(70) as u32;
            let neg = rng.chance(0.5);
            let want = i128::from(m) << s;
            let want = if neg { -want } else { want };
            assert_eq!(
                BigInt::from_shifted(neg, m, s).to_i128(),
                Some(want),
                "m={m} s={s}"
            );
        });
        assert_eq!(BigInt::from_shifted(false, 1, 2000).bit_len(), 2001);
        assert_eq!(BigInt::from_shifted(true, 0, 500), BigInt::zero());
        assert_eq!(BigInt::from_shifted(false, u64::MAX, 64).bit_len(), 128);
    }

    #[test]
    fn remainder_matches_i128() {
        for_cases("remainder_matches_i128", 20_000, |rng| {
            let v = random_i128(rng);
            let p = *rng.pick(&[(1u64 << 61) - 1, 1_000_000_007, 2, 3, u64::MAX]);
            let want = v.rem_euclid(i128::from(p)) as u64;
            assert_eq!(BigInt::from_i128(v).rem_u64(p), want, "v={v} p={p}");
        });
    }

    #[test]
    fn zero_has_one_form() {
        let z = BigInt::from_limbs(true, vec![0, 0, 0]);
        assert_eq!(z, BigInt::zero());
        assert!(z.is_zero());
        assert_eq!(z.signum(), 0);
        assert_eq!(-z.clone(), z);
        assert_eq!(BigInt::from_i128(5) - BigInt::from_i128(5), BigInt::zero());
    }
}
