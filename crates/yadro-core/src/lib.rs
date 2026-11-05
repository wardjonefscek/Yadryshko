//! # yadro — своё CAD-ядро
//!
//! Ступень 0, «Фундамент»: допуски, 2D-сетка в микрометрах, векторы и
//! преобразования, точные предикаты со своей длинной арифметикой,
//! отладочный вывод в SVG и STL, проверки сетки треугольников.
//!
//! Ступень 1, «2D-ядро» (идёт): точные операции с отрезками ([`seg`]),
//! контуры ([`ring`]), области с отверстиями ([`region`]) и булевы операции
//! над ними методом snap rounding ([`boolean`]); дуги ([`arc`]), контуры из
//! отрезков и дуг ([`path`]), эквидистанта и скругления ([`offset`]); профили
//! крыла из .dat и NACA ([`airfoil`]), контур из поля плотности TopOpt ([`density`]).
//!
//! Единицы — миллиметры. 2D-координаты хранятся целыми микрометрами
//! (см. [`grid`]); 3D — числами `f64`.
//!
//! Хрупкие места помечены в коде комментарием `РИСК:`.

pub mod airfoil;
pub mod arc;
pub mod bigint;
pub mod boolean;
pub mod density;
pub mod exact;
pub mod grid;
pub mod mesh;
pub mod offset;
pub mod path;
pub mod predicates;
pub mod region;
pub mod ring;
pub mod seg;
pub mod stl;
pub mod svg;
pub mod testkit;
pub mod tol;
pub mod trace;
pub mod vec;
pub mod xform;

pub use airfoil::Airfoil;
pub use arc::Arc;
pub use boolean::Op;
pub use density::{ContourOptions, DensityField};
pub use grid::GridPoint;
pub use mesh::Mesh;
pub use path::{Path, PathVertex, Segment};
pub use predicates::Sign;
pub use region::Region;
pub use ring::{Location, Ring};
pub use vec::{Vec2, Vec3};
pub use xform::{Xform2, Xform3};

/// Версия ядра.
pub const VERSION: &str = env!("CARGO_PKG_VERSION");
