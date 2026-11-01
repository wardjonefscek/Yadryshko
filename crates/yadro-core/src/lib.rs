//! # yadro — своё CAD-ядро
//!
//! Ступень 0, «Фундамент»: допуски, 2D-сетка в микрометрах, векторы и
//! преобразования, точные предикаты со своей длинной арифметикой,
//! отладочный вывод в SVG и STL, проверки сетки треугольников.
//!
//! Единицы — миллиметры. 2D-координаты хранятся целыми микрометрами
//! (см. [`grid`]); 3D — числами `f64`.
//!
//! Хрупкие места помечены в коде комментарием `РИСК:`.

pub mod bigint;
pub mod exact;
pub mod grid;
pub mod mesh;
pub mod predicates;
pub mod stl;
pub mod svg;
pub mod testkit;
pub mod tol;
pub mod trace;
pub mod vec;
pub mod xform;

pub use grid::GridPoint;
pub use mesh::Mesh;
pub use predicates::Sign;
pub use vec::{Vec2, Vec3};
pub use xform::{Xform2, Xform3};

/// Версия ядра.
pub const VERSION: &str = env!("CARGO_PKG_VERSION");
