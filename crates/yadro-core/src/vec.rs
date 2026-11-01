//! Векторы на `f64`: [`Vec2`] и [`Vec3`]. Точка и направление — один тип.

use std::ops::{Add, AddAssign, Div, Mul, Neg, Sub, SubAssign};

/// Вектор или точка на плоскости, мм.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct Vec2 {
    /// Координата X.
    pub x: f64,
    /// Координата Y.
    pub y: f64,
}

/// Вектор или точка в пространстве, мм.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct Vec3 {
    /// Координата X.
    pub x: f64,
    /// Координата Y.
    pub y: f64,
    /// Координата Z.
    pub z: f64,
}

/// Общие операции для обоих типов векторов.
macro_rules! impl_vector {
    ($t:ident { $($f:ident),+ }) => {
        impl Add for $t {
            type Output = Self;
            fn add(self, o: Self) -> Self {
                Self { $($f: self.$f + o.$f),+ }
            }
        }

        impl Sub for $t {
            type Output = Self;
            fn sub(self, o: Self) -> Self {
                Self { $($f: self.$f - o.$f),+ }
            }
        }

        impl Mul<f64> for $t {
            type Output = Self;
            fn mul(self, k: f64) -> Self {
                Self { $($f: self.$f * k),+ }
            }
        }

        impl Mul<$t> for f64 {
            type Output = $t;
            fn mul(self, v: $t) -> $t {
                v * self
            }
        }

        impl Div<f64> for $t {
            type Output = Self;
            fn div(self, k: f64) -> Self {
                Self { $($f: self.$f / k),+ }
            }
        }

        impl Neg for $t {
            type Output = Self;
            fn neg(self) -> Self {
                Self { $($f: -self.$f),+ }
            }
        }

        impl AddAssign for $t {
            fn add_assign(&mut self, o: Self) {
                $(self.$f += o.$f;)+
            }
        }

        impl SubAssign for $t {
            fn sub_assign(&mut self, o: Self) {
                $(self.$f -= o.$f;)+
            }
        }

        impl $t {
            /// Скалярное произведение.
            pub fn dot(self, o: Self) -> f64 {
                0.0 $(+ self.$f * o.$f)+
            }

            /// Квадрат длины.
            pub fn length_sq(self) -> f64 {
                self.dot(self)
            }

            /// Расстояние до точки `o`.
            pub fn distance(self, o: Self) -> f64 {
                (self - o).length()
            }

            /// Единичный вектор того же направления; `None` для нулевого вектора и нечисел.
            pub fn normalized(self) -> Option<Self> {
                let len = self.length();
                (len > 0.0 && len.is_finite()).then(|| self / len)
            }

            /// Все компоненты — конечные числа.
            pub fn is_finite(self) -> bool {
                true $(&& self.$f.is_finite())+
            }

            /// Покомпонентный минимум.
            pub fn min(self, o: Self) -> Self {
                Self { $($f: self.$f.min(o.$f)),+ }
            }

            /// Покомпонентный максимум.
            pub fn max(self, o: Self) -> Self {
                Self { $($f: self.$f.max(o.$f)),+ }
            }
        }
    };
}

impl_vector!(Vec2 { x, y });
impl_vector!(Vec3 { x, y, z });

impl Vec2 {
    /// Нулевой вектор.
    pub const ZERO: Self = Self::new(0.0, 0.0);

    /// Создать вектор.
    pub const fn new(x: f64, y: f64) -> Self {
        Self { x, y }
    }

    /// Косое произведение `self × o`: больше нуля, если `o` повёрнут
    /// от `self` против часовой стрелки.
    pub fn cross(self, o: Self) -> f64 {
        self.x * o.y - self.y * o.x
    }

    /// Длина.
    pub fn length(self) -> f64 {
        self.x.hypot(self.y)
    }

    /// Поворот на 90° против часовой стрелки.
    pub fn perp(self) -> Self {
        Self::new(-self.y, self.x)
    }

    /// Точка в пространстве с заданной высотой `z`.
    pub fn with_z(self, z: f64) -> Vec3 {
        Vec3::new(self.x, self.y, z)
    }
}

impl Vec3 {
    /// Нулевой вектор.
    pub const ZERO: Self = Self::new(0.0, 0.0, 0.0);
    /// Единичный вектор оси X.
    pub const X: Self = Self::new(1.0, 0.0, 0.0);
    /// Единичный вектор оси Y.
    pub const Y: Self = Self::new(0.0, 1.0, 0.0);
    /// Единичный вектор оси Z.
    pub const Z: Self = Self::new(0.0, 0.0, 1.0);

    /// Создать вектор.
    pub const fn new(x: f64, y: f64, z: f64) -> Self {
        Self { x, y, z }
    }

    /// Векторное произведение (правило правой руки).
    pub fn cross(self, o: Self) -> Self {
        Self::new(
            self.y * o.z - self.z * o.y,
            self.z * o.x - self.x * o.z,
            self.x * o.y - self.y * o.x,
        )
    }

    /// Длина.
    pub fn length(self) -> f64 {
        self.x.hypot(self.y).hypot(self.z)
    }

    /// Проекция на плоскость XY.
    pub fn xy(self) -> Vec2 {
        Vec2::new(self.x, self.y)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn cross_products_follow_right_hand_rule() {
        assert_eq!(Vec3::X.cross(Vec3::Y), Vec3::Z);
        assert_eq!(Vec3::Y.cross(Vec3::Z), Vec3::X);
        assert_eq!(Vec2::new(1.0, 0.0).cross(Vec2::new(0.0, 1.0)), 1.0);
        assert_eq!(Vec2::new(1.0, 0.0).perp(), Vec2::new(0.0, 1.0));
    }

    #[test]
    fn lengths_and_normalization() {
        assert_eq!(Vec2::new(3.0, 4.0).length(), 5.0);
        assert_eq!(Vec3::new(2.0, 3.0, 6.0).length(), 7.0);
        assert_eq!(Vec2::ZERO.normalized(), None);
        assert_eq!(Vec2::new(f64::NAN, 1.0).normalized(), None);
        assert_eq!(Vec2::new(0.0, 5.0).normalized(), Some(Vec2::new(0.0, 1.0)));
        assert_eq!(Vec3::new(1.0, 2.0, 2.0).distance(Vec3::ZERO), 3.0);
    }

    #[test]
    fn arithmetic_and_bounds() {
        let a = Vec3::new(1.0, -2.0, 3.0);
        let b = Vec3::new(-1.0, 5.0, 0.5);
        assert_eq!(a + b - b, a);
        assert_eq!(2.0 * a, a * 2.0);
        assert_eq!(-a / 2.0, Vec3::new(-0.5, 1.0, -1.5));
        assert_eq!(a.min(b), Vec3::new(-1.0, -2.0, 0.5));
        assert_eq!(a.max(b), Vec3::new(1.0, 5.0, 3.0));
        assert!(!Vec3::new(1.0, f64::INFINITY, 0.0).is_finite());
        assert_eq!(a.xy().with_z(a.z), a);
    }
}
