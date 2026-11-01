//! Счётчики предикатов и журнал «решений на грани».
//!
//! Каждый вызов предиката отмечается: хватило ли быстрого расчёта или
//! пришлось считать точно. Точный путь значит, что точки почти вырождены, —
//! именно такие места отладчик будет подсвечивать. Счётчики свои у каждого
//! потока, поэтому параллельные тесты друг другу не мешают.

use crate::predicates::Sign;
use std::cell::{Cell, RefCell};

/// Какой предикат сработал.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum Predicate {
    /// Ориентация на плоскости (`f64`).
    Orient2d,
    /// Ориентация в пространстве (`f64`).
    Orient3d,
}

/// Счётчики вызовов в текущем потоке.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct Counters {
    /// orient2d: знак дал быстрый расчёт.
    pub orient2d_fast: u64,
    /// orient2d: понадобился точный путь.
    pub orient2d_exact: u64,
    /// orient3d: знак дал быстрый расчёт.
    pub orient3d_fast: u64,
    /// orient3d: понадобился точный путь.
    pub orient3d_exact: u64,
}

impl Counters {
    /// Все счётчики по нулям.
    pub const ZERO: Self = Self {
        orient2d_fast: 0,
        orient2d_exact: 0,
        orient3d_fast: 0,
        orient3d_exact: 0,
    };
}

/// Случай, когда понадобился точный путь.
#[derive(Clone, Debug, PartialEq)]
pub struct NearZero {
    /// Какой предикат.
    pub predicate: Predicate,
    /// Входные координаты подряд: x, y (и z) каждой точки.
    pub inputs: Vec<f64>,
    /// Точный ответ.
    pub sign: Sign,
}

thread_local! {
    static COUNTERS: Cell<Counters> = const { Cell::new(Counters::ZERO) };
    static CAPTURE: RefCell<Option<Vec<NearZero>>> = const { RefCell::new(None) };
}

/// Текущие счётчики потока.
pub fn counters() -> Counters {
    COUNTERS.with(Cell::get)
}

/// Обнулить счётчики потока.
pub fn reset_counters() {
    COUNTERS.with(|c| c.set(Counters::ZERO));
}

/// Начать запись случаев точного пути; прежняя запись стирается.
pub fn start_capture() {
    CAPTURE.with(|c| *c.borrow_mut() = Some(Vec::new()));
}

/// Забрать записанные случаи и остановить запись.
pub fn take_captured() -> Vec<NearZero> {
    CAPTURE.with(|c| c.borrow_mut().take().unwrap_or_default())
}

pub(crate) fn note_fast(p: Predicate) {
    COUNTERS.with(|c| {
        let mut v = c.get();
        match p {
            Predicate::Orient2d => v.orient2d_fast += 1,
            Predicate::Orient3d => v.orient3d_fast += 1,
        }
        c.set(v);
    });
}

pub(crate) fn note_exact(p: Predicate, inputs: &[f64], sign: Sign) {
    COUNTERS.with(|c| {
        let mut v = c.get();
        match p {
            Predicate::Orient2d => v.orient2d_exact += 1,
            Predicate::Orient3d => v.orient3d_exact += 1,
        }
        c.set(v);
    });
    CAPTURE.with(|c| {
        if let Some(log) = c.borrow_mut().as_mut() {
            log.push(NearZero {
                predicate: p,
                inputs: inputs.to_vec(),
                sign,
            });
        }
    });
}
