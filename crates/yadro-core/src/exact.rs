//! Точный вид чисел `f64`: каждое конечное число — целое, умноженное на степень двойки.

use crate::bigint::BigInt;

/// Разложение конечного `f64`: x = (−1)^neg · m · 2^e. Для нечисел и бесконечностей — `None`.
///
/// У нуля `m = 0`. Денормализованные числа (меньше 2^−1022) раскладываются тоже точно.
pub fn decompose(x: f64) -> Option<(bool, u64, i32)> {
    if !x.is_finite() {
        return None;
    }
    let bits = x.to_bits();
    let neg = bits >> 63 == 1;
    let exp = ((bits >> 52) & 0x7ff) as i32;
    let frac = bits & ((1u64 << 52) - 1);
    Some(if exp == 0 {
        (neg, frac, -1074)
    } else {
        (neg, frac | (1u64 << 52), exp - 1075)
    })
}

/// Перевести числа в целые с общим масштабом: x_i = n_i · 2^s.
///
/// Масштаб s — наименьший показатель среди ненулевых чисел. Все числа
/// умножаются на одно и то же положительное 2^−s, поэтому знак
/// однородного выражения от разностей координат (любого предиката)
/// не меняется. `None`, если среди чисел есть нечисло или бесконечность.
pub fn to_common_integers(xs: &[f64]) -> Option<Vec<BigInt>> {
    let parts = xs
        .iter()
        .map(|&x| decompose(x))
        .collect::<Option<Vec<_>>>()?;
    let scale = parts
        .iter()
        .filter(|p| p.1 != 0)
        .map(|p| p.2)
        .min()
        .unwrap_or(0);
    Some(
        parts
            .into_iter()
            .map(|(neg, m, e)| {
                if m == 0 {
                    BigInt::zero()
                } else {
                    // Не больше 971 − (−1074) = 2045 бит сдвига.
                    BigInt::from_shifted(neg, m, (e - scale) as u32)
                }
            })
            .collect(),
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::testkit::for_cases;

    #[test]
    fn known_values() {
        assert_eq!(decompose(1.0), Some((false, 1 << 52, -52)));
        assert_eq!(decompose(0.5), Some((false, 1 << 52, -53)));
        assert_eq!(decompose(-3.0), Some((true, 3 << 51, -51)));
        assert_eq!(decompose(0.0), Some((false, 0, -1074)));
        assert_eq!(decompose(-0.0), Some((true, 0, -1074)));
        assert_eq!(decompose(f64::from_bits(1)), Some((false, 1, -1074)));
        assert_eq!(decompose(f64::MAX), Some((false, (1 << 53) - 1, 971)));
        assert_eq!(decompose(f64::NAN), None);
        assert_eq!(decompose(f64::INFINITY), None);
    }

    #[test]
    fn decomposition_rebuilds_the_number() {
        for_cases("decomposition_rebuilds_the_number", 20_000, |rng| {
            let exp = rng.range_i64(200, 1800) as u64;
            let x = f64::from_bits((rng.next_u64() & 0x800F_FFFF_FFFF_FFFF) | (exp << 52));
            let (neg, m, e) = decompose(x).unwrap();
            let rebuilt = m as f64 * 2f64.powi(e);
            assert_eq!(if neg { -rebuilt } else { rebuilt }, x);
        });
    }

    #[test]
    fn common_integers_keep_ratios() {
        let n = to_common_integers(&[0.5, 3.0, -0.25, 0.0]).unwrap();
        // Масштаб 2^−54: 0,5 → 2^53, 3 → 3·2^54, −0,25 → −2^52.
        assert_eq!(n[0].to_i128(), Some(1 << 53));
        assert_eq!(n[1].to_i128(), Some(3 << 54));
        assert_eq!(n[2].to_i128(), Some(-(1 << 52)));
        assert!(n[3].is_zero());
        assert_eq!(to_common_integers(&[1.0, f64::NAN]), None);
        assert!(
            to_common_integers(&[0.0, -0.0])
                .unwrap()
                .iter()
                .all(BigInt::is_zero)
        );
    }
}
