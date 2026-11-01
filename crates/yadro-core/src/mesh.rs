//! Сетка треугольников: объём, габарит, замкнутость, характеристика Эйлера.
//!
//! Внешняя сторона треугольника — та, откуда обход его вершин виден
//! против часовой стрелки. У замкнутого тела каждое ребро принадлежит
//! ровно двум треугольникам, и они проходят его в противоположных
//! направлениях.

use crate::vec::Vec3;
use std::collections::HashSet;
use std::fmt;

/// Сетка треугольников: вершины в мм и тройки индексов вершин.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Mesh {
    /// Вершины.
    pub vertices: Vec<Vec3>,
    /// Треугольники: индексы вершин, снаружи обход против часовой стрелки.
    pub triangles: Vec<[u32; 3]>,
}

/// Чем сетка не годится как замкнутое тело.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum MeshIssue {
    /// Индекс вершины вне списка вершин.
    BadIndex {
        /// Номер треугольника.
        triangle: usize,
    },
    /// В треугольнике повторяется вершина.
    DegenerateTriangle {
        /// Номер треугольника.
        triangle: usize,
    },
    /// У ребра нет пары: в теле дыра.
    OpenEdge {
        /// Начало ребра.
        from: u32,
        /// Конец ребра.
        to: u32,
    },
    /// Ребро пройдено в одну сторону дважды: соседние треугольники
    /// ориентированы несогласованно или ребро несут больше двух треугольников.
    RepeatedEdge {
        /// Начало ребра.
        from: u32,
        /// Конец ребра.
        to: u32,
    },
}

impl fmt::Display for MeshIssue {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match *self {
            Self::BadIndex { triangle } => {
                write!(f, "треугольник {triangle}: индекс вершины вне списка")
            }
            Self::DegenerateTriangle { triangle } => {
                write!(f, "треугольник {triangle}: вершина повторяется")
            }
            Self::OpenEdge { from, to } => write!(f, "у ребра {from}→{to} нет пары: в теле дыра"),
            Self::RepeatedEdge { from, to } => write!(
                f,
                "ребро {from}→{to} пройдено дважды в одну сторону: несогласованная ориентация"
            ),
        }
    }
}

impl std::error::Error for MeshIssue {}

impl Mesh {
    /// Прямоугольный брусок от `min` до `max` (12 треугольников).
    pub fn cuboid(min: Vec3, max: Vec3) -> Self {
        let pick = |bit: bool, lo: f64, hi: f64| if bit { hi } else { lo };
        let vertices = (0..8u32)
            .map(|i| {
                Vec3::new(
                    pick(i & 1 != 0, min.x, max.x),
                    pick(i & 2 != 0, min.y, max.y),
                    pick(i & 4 != 0, min.z, max.z),
                )
            })
            .collect();
        let triangles = vec![
            [0, 4, 6], // −X
            [0, 6, 2],
            [1, 3, 7], // +X
            [1, 7, 5],
            [0, 1, 5], // −Y
            [0, 5, 4],
            [2, 6, 7], // +Y
            [2, 7, 3],
            [0, 2, 3], // −Z
            [0, 3, 1],
            [4, 5, 7], // +Z
            [4, 7, 6],
        ];
        Self {
            vertices,
            triangles,
        }
    }

    /// Вершины треугольника `t`; `None`, если индекс вне списка.
    pub fn triangle(&self, t: usize) -> Option<[Vec3; 3]> {
        let [i, j, k] = *self.triangles.get(t)?;
        Some([
            *self.vertices.get(i as usize)?,
            *self.vertices.get(j as usize)?,
            *self.vertices.get(k as usize)?,
        ])
    }

    /// Единичная внешняя нормаль треугольника `t` (ноль для вырожденного).
    pub fn normal(&self, t: usize) -> Vec3 {
        self.triangle(t)
            .and_then(|[a, b, c]| (b - a).cross(c - a).normalized())
            .unwrap_or(Vec3::ZERO)
    }

    /// Ориентированный объём, мм³: у замкнутого тела с внешними нормалями — больше нуля.
    pub fn signed_volume(&self) -> f64 {
        (0..self.triangles.len())
            .filter_map(|t| self.triangle(t))
            .map(|[a, b, c]| a.dot(b.cross(c)))
            .sum::<f64>()
            / 6.0
    }

    /// Габарит: углы (минимум, максимум); `None` для пустой сетки.
    pub fn bbox(&self) -> Option<(Vec3, Vec3)> {
        let first = *self.vertices.first()?;
        Some(
            self.vertices
                .iter()
                .fold((first, first), |(lo, hi), &v| (lo.min(v), hi.max(v))),
        )
    }

    /// Проверить, что сетка — замкнутое тело с согласованной ориентацией.
    ///
    /// Ошибки ищутся в порядке треугольников, поэтому ответ не зависит
    /// от случая и его удобно сравнивать в тестах.
    pub fn check_closed(&self) -> Result<(), MeshIssue> {
        let n = self.vertices.len();
        let mut edges: HashSet<(u32, u32)> = HashSet::with_capacity(self.triangles.len() * 3);
        for (t, &[i, j, k]) in self.triangles.iter().enumerate() {
            if [i, j, k].iter().any(|&v| v as usize >= n) {
                return Err(MeshIssue::BadIndex { triangle: t });
            }
            if i == j || j == k || k == i {
                return Err(MeshIssue::DegenerateTriangle { triangle: t });
            }
            for (from, to) in [(i, j), (j, k), (k, i)] {
                if !edges.insert((from, to)) {
                    return Err(MeshIssue::RepeatedEdge { from, to });
                }
            }
        }
        for &[i, j, k] in &self.triangles {
            for (from, to) in [(i, j), (j, k), (k, i)] {
                if !edges.contains(&(to, from)) {
                    return Err(MeshIssue::OpenEdge { from, to });
                }
            }
        }
        Ok(())
    }

    /// Характеристика Эйлера V − E + F по вершинам и рёбрам, которые есть в треугольниках.
    ///
    /// У замкнутого тела без сквозных отверстий она равна 2, с g отверстиями —
    /// 2 − 2g (пластина с n отверстиями: 2 − 2n).
    pub fn euler_characteristic(&self) -> i64 {
        let mut vertices = HashSet::new();
        let mut edges = HashSet::new();
        for &[i, j, k] in &self.triangles {
            vertices.extend([i, j, k]);
            for (a, b) in [(i, j), (j, k), (k, i)] {
                edges.insert((a.min(b), a.max(b)));
            }
        }
        vertices.len() as i64 - edges.len() as i64 + self.triangles.len() as i64
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn brick() -> Mesh {
        Mesh::cuboid(Vec3::new(1.0, 2.0, 3.0), Vec3::new(21.0, 32.0, 6.0))
    }

    #[test]
    fn brick_is_closed_with_right_volume() {
        let m = brick();
        assert_eq!(m.check_closed(), Ok(()));
        assert_eq!(m.signed_volume(), 20.0 * 30.0 * 3.0);
        assert_eq!(m.euler_characteristic(), 2);
        assert_eq!(
            m.bbox(),
            Some((Vec3::new(1.0, 2.0, 3.0), Vec3::new(21.0, 32.0, 6.0)))
        );
        // Нормали смотрят наружу: у грани −X нормаль (−1, 0, 0).
        assert_eq!(m.normal(0), -Vec3::X);
        assert_eq!(m.normal(11), Vec3::Z);
    }

    #[test]
    fn hole_is_found() {
        let mut m = brick();
        m.triangles.pop();
        assert!(matches!(m.check_closed(), Err(MeshIssue::OpenEdge { .. })));
    }

    #[test]
    fn flipped_triangle_is_found() {
        let mut m = brick();
        m.triangles[3].swap(1, 2);
        assert!(matches!(
            m.check_closed(),
            Err(MeshIssue::RepeatedEdge { .. })
        ));
    }

    #[test]
    fn broken_triangles_are_found() {
        let mut m = brick();
        m.triangles[5] = [0, 0, 1];
        assert_eq!(
            m.check_closed(),
            Err(MeshIssue::DegenerateTriangle { triangle: 5 })
        );
        m.triangles[5] = [0, 99, 1];
        assert_eq!(m.check_closed(), Err(MeshIssue::BadIndex { triangle: 5 }));
        assert_eq!(m.triangle(5), None);
    }

    #[test]
    fn inside_out_brick_has_negative_volume() {
        let mut m = brick();
        for t in &mut m.triangles {
            t.swap(1, 2);
        }
        assert_eq!(m.check_closed(), Ok(()));
        assert_eq!(m.signed_volume(), -1800.0);
        assert!(Mesh::default().bbox().is_none());
    }
}
