//! Аффинные преобразования: [`Xform2`] на плоскости и [`Xform3`] в пространстве.
//!
//! Преобразование — это матрица `a` и сдвиг `t`: p' = a·p + t.
//! Порядок задаётся словами, а не умножением: `x.then(&y)` — сначала `x`,
//! потом `y`. Так нервюру сначала строят в своей плоскости, потом ставят
//! на место по размаху, и путаницы с порядком нет.

use crate::vec::{Vec2, Vec3};

/// Аффинное преобразование плоскости: p' = a·p + t.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Xform2 {
    /// Матрица 2×2 по строкам.
    pub a: [[f64; 2]; 2],
    /// Сдвиг.
    pub t: Vec2,
}

impl Xform2 {
    /// Тождественное преобразование.
    pub const IDENTITY: Self = Self {
        a: [[1.0, 0.0], [0.0, 1.0]],
        t: Vec2::ZERO,
    };

    /// Сдвиг на вектор `d`.
    pub fn translation(d: Vec2) -> Self {
        Self {
            t: d,
            ..Self::IDENTITY
        }
    }

    /// Поворот вокруг начала координат на угол в радианах, против часовой стрелки.
    pub fn rotation(angle: f64) -> Self {
        let (s, c) = angle.sin_cos();
        Self {
            a: [[c, -s], [s, c]],
            t: Vec2::ZERO,
        }
    }

    /// Масштаб по осям относительно начала координат.
    pub fn scaling(sx: f64, sy: f64) -> Self {
        Self {
            a: [[sx, 0.0], [0.0, sy]],
            t: Vec2::ZERO,
        }
    }

    /// Применить к направлению (без сдвига).
    pub fn apply_vector(&self, v: Vec2) -> Vec2 {
        let a = &self.a;
        Vec2::new(a[0][0] * v.x + a[0][1] * v.y, a[1][0] * v.x + a[1][1] * v.y)
    }

    /// Применить к точке.
    pub fn apply_point(&self, p: Vec2) -> Vec2 {
        self.apply_vector(p) + self.t
    }

    /// Сначала `self`, потом `next`.
    pub fn then(&self, next: &Self) -> Self {
        let (p, q) = (&self.a, &next.a);
        Self {
            a: [
                [
                    q[0][0] * p[0][0] + q[0][1] * p[1][0],
                    q[0][0] * p[0][1] + q[0][1] * p[1][1],
                ],
                [
                    q[1][0] * p[0][0] + q[1][1] * p[1][0],
                    q[1][0] * p[0][1] + q[1][1] * p[1][1],
                ],
            ],
            t: next.apply_point(self.t),
        }
    }

    /// Определитель линейной части: больше нуля — ориентация сохраняется.
    pub fn determinant(&self) -> f64 {
        self.a[0][0] * self.a[1][1] - self.a[0][1] * self.a[1][0]
    }

    /// Обратное преобразование; `None`, если оно вырождено или не число.
    pub fn inverse(&self) -> Option<Self> {
        let det = self.determinant();
        if det == 0.0 || !det.is_finite() {
            return None;
        }
        let [[a, b], [c, d]] = self.a;
        let linear = Self {
            a: [[d / det, -b / det], [-c / det, a / det]],
            t: Vec2::ZERO,
        };
        Some(Self {
            t: -linear.apply_vector(self.t),
            ..linear
        })
    }
}

/// Аффинное преобразование пространства: p' = a·p + t.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Xform3 {
    /// Матрица 3×3 по строкам.
    pub a: [[f64; 3]; 3],
    /// Сдвиг.
    pub t: Vec3,
}

fn mat3_mul(q: &[[f64; 3]; 3], p: &[[f64; 3]; 3]) -> [[f64; 3]; 3] {
    std::array::from_fn(|i| {
        std::array::from_fn(|j| q[i][0] * p[0][j] + q[i][1] * p[1][j] + q[i][2] * p[2][j])
    })
}

impl Xform3 {
    /// Тождественное преобразование.
    pub const IDENTITY: Self = Self {
        a: [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
        t: Vec3::ZERO,
    };

    /// Сдвиг на вектор `d`.
    pub fn translation(d: Vec3) -> Self {
        Self {
            t: d,
            ..Self::IDENTITY
        }
    }

    /// Поворот вокруг оси X (правило правой руки): Y переходит в Z.
    pub fn rotation_x(angle: f64) -> Self {
        let (s, c) = angle.sin_cos();
        Self::linear([[1.0, 0.0, 0.0], [0.0, c, -s], [0.0, s, c]])
    }

    /// Поворот вокруг оси Y (правило правой руки): Z переходит в X.
    pub fn rotation_y(angle: f64) -> Self {
        let (s, c) = angle.sin_cos();
        Self::linear([[c, 0.0, s], [0.0, 1.0, 0.0], [-s, 0.0, c]])
    }

    /// Поворот вокруг оси Z (правило правой руки): X переходит в Y.
    pub fn rotation_z(angle: f64) -> Self {
        let (s, c) = angle.sin_cos();
        Self::linear([[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]])
    }

    /// Поворот вокруг оси через начало координат (правило правой руки);
    /// `None` для нулевой оси.
    pub fn rotation_axis(axis: Vec3, angle: f64) -> Option<Self> {
        let k = axis.normalized()?;
        let (s, c) = angle.sin_cos();
        let v = 1.0 - c;
        Some(Self::linear([
            [
                c + k.x * k.x * v,
                k.x * k.y * v - k.z * s,
                k.x * k.z * v + k.y * s,
            ],
            [
                k.y * k.x * v + k.z * s,
                c + k.y * k.y * v,
                k.y * k.z * v - k.x * s,
            ],
            [
                k.z * k.x * v - k.y * s,
                k.z * k.y * v + k.x * s,
                c + k.z * k.z * v,
            ],
        ]))
    }

    /// Масштаб по осям относительно начала координат.
    pub fn scaling(sx: f64, sy: f64, sz: f64) -> Self {
        Self::linear([[sx, 0.0, 0.0], [0.0, sy, 0.0], [0.0, 0.0, sz]])
    }

    fn linear(a: [[f64; 3]; 3]) -> Self {
        Self { a, t: Vec3::ZERO }
    }

    /// Применить к направлению (без сдвига).
    pub fn apply_vector(&self, v: Vec3) -> Vec3 {
        let a = &self.a;
        Vec3::new(
            a[0][0] * v.x + a[0][1] * v.y + a[0][2] * v.z,
            a[1][0] * v.x + a[1][1] * v.y + a[1][2] * v.z,
            a[2][0] * v.x + a[2][1] * v.y + a[2][2] * v.z,
        )
    }

    /// Применить к точке.
    pub fn apply_point(&self, p: Vec3) -> Vec3 {
        self.apply_vector(p) + self.t
    }

    /// Сначала `self`, потом `next`.
    pub fn then(&self, next: &Self) -> Self {
        Self {
            a: mat3_mul(&next.a, &self.a),
            t: next.apply_point(self.t),
        }
    }

    /// Определитель линейной части: больше нуля — ориентация сохраняется.
    pub fn determinant(&self) -> f64 {
        let m = &self.a;
        m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
            + m[0][1] * (m[1][2] * m[2][0] - m[1][0] * m[2][2])
            + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0])
    }

    /// Обратное преобразование; `None`, если оно вырождено или не число.
    pub fn inverse(&self) -> Option<Self> {
        let det = self.determinant();
        if det == 0.0 || !det.is_finite() {
            return None;
        }
        let m = &self.a;
        let r = 1.0 / det;
        let linear = Self::linear([
            [
                (m[1][1] * m[2][2] - m[1][2] * m[2][1]) * r,
                (m[0][2] * m[2][1] - m[0][1] * m[2][2]) * r,
                (m[0][1] * m[1][2] - m[0][2] * m[1][1]) * r,
            ],
            [
                (m[1][2] * m[2][0] - m[1][0] * m[2][2]) * r,
                (m[0][0] * m[2][2] - m[0][2] * m[2][0]) * r,
                (m[0][2] * m[1][0] - m[0][0] * m[1][2]) * r,
            ],
            [
                (m[1][0] * m[2][1] - m[1][1] * m[2][0]) * r,
                (m[0][1] * m[2][0] - m[0][0] * m[2][1]) * r,
                (m[0][0] * m[1][1] - m[0][1] * m[1][0]) * r,
            ],
        ]);
        Some(Self {
            t: -linear.apply_vector(self.t),
            ..linear
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::testkit::{Rng, for_cases};
    use std::f64::consts::FRAC_PI_2;

    fn close2(a: Vec2, b: Vec2, tol: f64) -> bool {
        (a - b).length() <= tol
    }

    fn close3(a: Vec3, b: Vec3, tol: f64) -> bool {
        (a - b).length() <= tol
    }

    fn random_xform2(rng: &mut Rng) -> Xform2 {
        Xform2::scaling(rng.range_f64(0.1, 10.0), rng.range_f64(0.1, 10.0))
            .then(&Xform2::rotation(rng.range_f64(-4.0, 4.0)))
            .then(&Xform2::translation(Vec2::new(
                rng.range_f64(-1e3, 1e3),
                rng.range_f64(-1e3, 1e3),
            )))
    }

    fn random_unit(rng: &mut Rng) -> Vec3 {
        loop {
            let v = Vec3::new(
                rng.range_f64(-1.0, 1.0),
                rng.range_f64(-1.0, 1.0),
                rng.range_f64(-1.0, 1.0),
            );
            if let Some(u) = v.normalized().filter(|_| v.length() > 0.1) {
                return u;
            }
        }
    }

    fn random_xform3(rng: &mut Rng) -> Xform3 {
        let turn = Xform3::rotation_axis(random_unit(rng), rng.range_f64(-4.0, 4.0)).unwrap();
        Xform3::scaling(
            rng.range_f64(0.1, 10.0),
            rng.range_f64(0.1, 10.0),
            rng.range_f64(0.1, 10.0),
        )
        .then(&turn)
        .then(&Xform3::translation(Vec3::new(
            rng.range_f64(-1e3, 1e3),
            rng.range_f64(-1e3, 1e3),
            rng.range_f64(-1e3, 1e3),
        )))
    }

    #[test]
    fn order_of_operations_is_explicit() {
        let shift = Xform2::translation(Vec2::new(1.0, 0.0));
        let turn = Xform2::rotation(FRAC_PI_2);
        // Сначала сдвиг, потом поворот: (0, 0) → (1, 0) → (0, 1).
        let p = shift.then(&turn).apply_point(Vec2::ZERO);
        assert!(close2(p, Vec2::new(0.0, 1.0), 1e-15), "{p:?}");
        // Сначала поворот, потом сдвиг: (0, 0) → (0, 0) → (1, 0).
        let p = turn.then(&shift).apply_point(Vec2::ZERO);
        assert!(close2(p, Vec2::new(1.0, 0.0), 1e-15), "{p:?}");
    }

    #[test]
    fn named_rotations_follow_right_hand_rule() {
        let q = FRAC_PI_2;
        assert!(close3(
            Xform3::rotation_z(q).apply_vector(Vec3::X),
            Vec3::Y,
            1e-15
        ));
        assert!(close3(
            Xform3::rotation_x(q).apply_vector(Vec3::Y),
            Vec3::Z,
            1e-15
        ));
        assert!(close3(
            Xform3::rotation_y(q).apply_vector(Vec3::Z),
            Vec3::X,
            1e-15
        ));
        assert!(Xform3::rotation_axis(Vec3::ZERO, 1.0).is_none());
    }

    #[test]
    fn axis_rotation_matches_named_axes() {
        for angle in [0.3, -1.2, 2.9] {
            let pairs = [
                (Vec3::X, Xform3::rotation_x(angle)),
                (Vec3::Y, Xform3::rotation_y(angle)),
                (Vec3::Z, Xform3::rotation_z(angle)),
            ];
            for (axis, named) in pairs {
                let general = Xform3::rotation_axis(axis, angle).unwrap();
                let same = general
                    .a
                    .iter()
                    .flatten()
                    .zip(named.a.iter().flatten())
                    .all(|(x, y)| (x - y).abs() < 1e-15);
                assert!(same, "ось {axis:?}, угол {angle}");
            }
        }
    }

    #[test]
    fn inverse_undoes_2d() {
        for_cases("inverse_undoes_2d", 5000, |rng| {
            let x = random_xform2(rng);
            let inv = x.inverse().expect("преобразование невырождено");
            let p = Vec2::new(rng.range_f64(-1e3, 1e3), rng.range_f64(-1e3, 1e3));
            let back = inv.apply_point(x.apply_point(p));
            assert!(close2(back, p, 1e-8), "{x:?}: {p:?} → {back:?}");
            let id = x.then(&inv);
            assert!(close2(id.apply_point(p), p, 1e-8));
        });
    }

    #[test]
    fn inverse_undoes_3d() {
        for_cases("inverse_undoes_3d", 5000, |rng| {
            let x = random_xform3(rng);
            let inv = x.inverse().expect("преобразование невырождено");
            let p = Vec3::new(
                rng.range_f64(-1e3, 1e3),
                rng.range_f64(-1e3, 1e3),
                rng.range_f64(-1e3, 1e3),
            );
            let back = inv.apply_point(x.apply_point(p));
            assert!(close3(back, p, 1e-8), "{x:?}: {p:?} → {back:?}");
        });
    }

    #[test]
    fn rotations_keep_lengths() {
        for_cases("rotations_keep_lengths", 5000, |rng| {
            let turn = Xform3::rotation_axis(random_unit(rng), rng.range_f64(-10.0, 10.0)).unwrap();
            let v = Vec3::new(
                rng.range_f64(-1e3, 1e3),
                rng.range_f64(-1e3, 1e3),
                rng.range_f64(-1e3, 1e3),
            );
            let len = turn.apply_vector(v).length();
            assert!((len - v.length()).abs() <= 1e-12 * v.length().max(1.0));
            assert!((turn.determinant() - 1.0).abs() < 1e-14);
        });
    }

    #[test]
    fn degenerate_has_no_inverse() {
        assert!(Xform2::scaling(1.0, 0.0).inverse().is_none());
        assert!(Xform3::scaling(1.0, 2.0, 0.0).inverse().is_none());
        assert!(Xform3::scaling(f64::NAN, 1.0, 1.0).inverse().is_none());
    }
}
