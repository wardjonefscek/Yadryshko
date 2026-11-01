//! Инструменты для тестов: воспроизводимые случайные случаи.
//!
//! Своя замена proptest: проект собирается без интернета и без чужих
//! библиотек. Каждый случай получает свой сид; если случай падает,
//! в сообщении есть команда, которая повторяет ровно его.
//!
//! Переменные окружения: `YADRO_CASES` — число случаев на тест,
//! `YADRO_SEED` — начальный сид.

use crate::vec::{Vec2, Vec3};
use std::panic::{self, AssertUnwindSafe};

/// Генератор SplitMix64: быстрый и воспроизводимый.
#[derive(Clone, Debug)]
pub struct Rng {
    state: u64,
}

impl Rng {
    /// Генератор с заданным сидом.
    pub fn new(seed: u64) -> Self {
        Self { state: seed }
    }

    /// Следующее 64-битное число.
    pub fn next_u64(&mut self) -> u64 {
        self.state = self.state.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.state;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }

    /// Целое из [0, n); n > 0.
    pub fn below(&mut self, n: u64) -> u64 {
        assert!(n > 0, "below(0)");
        ((u128::from(self.next_u64()) * u128::from(n)) >> 64) as u64
    }

    /// Целое из [lo, hi] включительно.
    pub fn range_i64(&mut self, lo: i64, hi: i64) -> i64 {
        assert!(lo <= hi, "range_i64: lo > hi");
        let span = (i128::from(hi) - i128::from(lo) + 1) as u128;
        if span > u128::from(u64::MAX) {
            return self.next_u64() as i64;
        }
        (i128::from(lo) + i128::from(self.below(span as u64))) as i64
    }

    /// Число из [0, 1).
    pub fn unit_f64(&mut self) -> f64 {
        (self.next_u64() >> 11) as f64 * (f64::EPSILON / 2.0)
    }

    /// Число из [lo, hi).
    pub fn range_f64(&mut self, lo: f64, hi: f64) -> f64 {
        lo + (hi - lo) * self.unit_f64()
    }

    /// Истина с вероятностью `p`.
    pub fn chance(&mut self, p: f64) -> bool {
        self.unit_f64() < p
    }

    /// Случайный элемент непустого среза.
    pub fn pick<'a, T>(&mut self, items: &'a [T]) -> &'a T {
        &items[self.below(items.len() as u64) as usize]
    }

    /// Число для проверки геометрии: смесь «удобных» и «злых» значений —
    /// мелкая решётка (частые вырождения), обычные размеры в мм и любые
    /// порядки, от денормализованных до огромных.
    pub fn geom_f64(&mut self) -> f64 {
        match self.below(6) {
            0 => self.range_i64(-8, 8) as f64,
            1 => self.range_i64(-64, 64) as f64 * 2f64.powi(self.range_i64(-20, 20) as i32),
            2 => self.range_f64(-1000.0, 1000.0),
            3 => self.range_i64(-1_000_000, 1_000_000) as f64 * 0.001,
            4 => {
                // Любой порядок: случайные биты без бесконечностей и нечисел.
                let exp = self.range_i64(0, 2046) as u64;
                f64::from_bits((self.next_u64() & 0x800F_FFFF_FFFF_FFFF) | (exp << 52))
            }
            _ => {
                // Крошечные: денормализованные и около 2^−1022.
                let exp = self.range_i64(0, 3) as u64;
                f64::from_bits((self.next_u64() & 0x800F_FFFF_FFFF_FFFF) | (exp << 52))
            }
        }
    }

    /// Точка на плоскости из [`Self::geom_f64`].
    pub fn geom_vec2(&mut self) -> Vec2 {
        Vec2::new(self.geom_f64(), self.geom_f64())
    }

    /// Точка в пространстве из [`Self::geom_f64`].
    pub fn geom_vec3(&mut self) -> Vec3 {
        Vec3::new(self.geom_f64(), self.geom_f64(), self.geom_f64())
    }

    /// Параметр точки на прямой: «круглая» двоичная дробь, чтобы точка
    /// ложилась на прямую ровно.
    pub fn dyadic_t(&mut self) -> f64 {
        *self.pick(&[0.0, 1.0, 0.5, 2.0, -1.0, 0.25, 0.75, 3.0, -0.5])
    }

    /// Сдвинуть число на 0–3 шага представления (ulp) в любую сторону.
    pub fn nudge(&mut self, x: f64) -> f64 {
        let steps = self.range_i64(-3, 3);
        (0..steps.unsigned_abs()).fold(x, |v, _| {
            if steps > 0 {
                v.next_up()
            } else {
                v.next_down()
            }
        })
    }

    /// [`Self::nudge`] для каждой координаты точки на плоскости.
    pub fn nudge2(&mut self, p: Vec2) -> Vec2 {
        Vec2::new(self.nudge(p.x), self.nudge(p.y))
    }

    /// [`Self::nudge`] для каждой координаты точки в пространстве.
    pub fn nudge3(&mut self, p: Vec3) -> Vec3 {
        Vec3::new(self.nudge(p.x), self.nudge(p.y), self.nudge(p.z))
    }
}

/// Число случаев на тест: `YADRO_CASES` или `default`.
pub fn cases(default: u64) -> u64 {
    std::env::var("YADRO_CASES")
        .ok()
        .and_then(|v| v.parse().ok())
        .unwrap_or(default)
}

/// Начальный сид: `YADRO_SEED` или постоянный, чтобы прогоны повторялись.
pub fn base_seed() -> u64 {
    std::env::var("YADRO_SEED")
        .ok()
        .and_then(|v| v.parse().ok())
        .unwrap_or(0x5EED_2026_1030)
}

/// Прогнать проверку на случайных случаях.
///
/// `name` — имя тестовой функции: по нему строится команда повтора.
/// При падении сообщает номер и сид случая.
pub fn for_cases(name: &str, default_cases: u64, mut check: impl FnMut(&mut Rng)) {
    let n = cases(default_cases);
    let base = base_seed();
    for i in 0..n {
        let seed = base.wrapping_add(i.wrapping_mul(0xD1B5_4A32_D192_ED03));
        let mut rng = Rng::new(seed);
        if let Err(payload) = panic::catch_unwind(AssertUnwindSafe(|| check(&mut rng))) {
            let msg = payload
                .downcast_ref::<String>()
                .map(String::as_str)
                .or_else(|| payload.downcast_ref::<&str>().copied())
                .unwrap_or("(без текста)");
            panic!(
                "{name}: упал случай {i} из {n}.\n\
                 Повторить только его: YADRO_SEED={seed} YADRO_CASES=1 cargo test -p yadro-core {name}\n\
                 {msg}"
            );
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn generator_is_reproducible_and_in_range() {
        let (mut a, mut b) = (Rng::new(7), Rng::new(7));
        for _ in 0..1000 {
            assert_eq!(a.next_u64(), b.next_u64());
            let v = a.range_i64(-3, 5);
            assert!((-3..=5).contains(&v));
            let u = a.unit_f64();
            assert!((0.0..1.0).contains(&u));
            assert!(a.geom_f64().is_finite());
            b.range_i64(-3, 5);
            b.unit_f64();
            b.geom_f64();
        }
        assert_eq!(
            Rng::new(1).range_i64(i64::MIN, i64::MAX),
            Rng::new(1).next_u64() as i64
        );
    }

    #[test]
    fn nudge_moves_by_a_few_ulps() {
        let mut rng = Rng::new(3);
        for _ in 0..1000 {
            let x = 1.0;
            let y = rng.nudge(x);
            assert!((y - x).abs() <= 3.0 * f64::EPSILON);
        }
    }
}
