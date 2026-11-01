//! STL: двоичный — компактный, для слайсера; текстовый — чтобы читать глазами.
//!
//! Координаты пишутся в миллиметрах, как ждут слайсеры. Точность `f32`
//! на детали размером в метр — около 0,06 мкм, для печати с запасом.

use crate::mesh::Mesh;
use crate::vec::Vec3;
use std::collections::HashMap;
use std::fs::File;
use std::io::{self, BufWriter, Write};
use std::path::Path;

fn bad_data(msg: &str) -> io::Error {
    io::Error::new(io::ErrorKind::InvalidData, msg)
}

/// Записать двоичный STL.
pub fn write_binary(mesh: &Mesh, mut out: impl Write) -> io::Result<()> {
    let count = u32::try_from(mesh.triangles.len())
        .map_err(|_| bad_data("слишком много треугольников для STL"))?;
    // Заголовок не должен начинаться со слова «solid»: так помечают текстовые файлы.
    let mut header = [0u8; 80];
    let tag = b"yadro binary STL";
    header[..tag.len()].copy_from_slice(tag);
    out.write_all(&header)?;
    out.write_all(&count.to_le_bytes())?;
    for t in 0..mesh.triangles.len() {
        let corners = mesh
            .triangle(t)
            .ok_or_else(|| bad_data("индекс вершины вне списка"))?;
        for v in std::iter::once(mesh.normal(t)).chain(corners) {
            for c in [v.x, v.y, v.z] {
                out.write_all(&(c as f32).to_le_bytes())?;
            }
        }
        out.write_all(&[0, 0])?;
    }
    Ok(())
}

/// Записать текстовый STL.
pub fn write_ascii(mesh: &Mesh, name: &str, mut out: impl Write) -> io::Result<()> {
    let name: String = name
        .chars()
        .map(|c| if c.is_whitespace() { '_' } else { c })
        .collect();
    writeln!(out, "solid {name}")?;
    for t in 0..mesh.triangles.len() {
        let corners = mesh
            .triangle(t)
            .ok_or_else(|| bad_data("индекс вершины вне списка"))?;
        let n = mesh.normal(t);
        writeln!(out, "  facet normal {} {} {}", n.x, n.y, n.z)?;
        writeln!(out, "    outer loop")?;
        for v in corners {
            writeln!(out, "      vertex {} {} {}", v.x, v.y, v.z)?;
        }
        writeln!(out, "    endloop")?;
        writeln!(out, "  endfacet")?;
    }
    writeln!(out, "endsolid {name}")
}

/// Сохранить двоичный STL в файл.
pub fn save_binary(mesh: &Mesh, path: impl AsRef<Path>) -> io::Result<()> {
    let mut w = BufWriter::new(File::create(path)?);
    write_binary(mesh, &mut w)?;
    w.flush()
}

/// Сохранить текстовый STL в файл.
pub fn save_ascii(mesh: &Mesh, name: &str, path: impl AsRef<Path>) -> io::Result<()> {
    let mut w = BufWriter::new(File::create(path)?);
    write_ascii(mesh, name, &mut w)?;
    w.flush()
}

/// Прочитать двоичный STL.
///
/// Одинаковые вершины (бит в бит, −0 считается равным +0) склеиваются,
/// поэтому замкнутость прочитанной сетки можно проверить сразу.
pub fn read_binary(bytes: &[u8]) -> io::Result<Mesh> {
    let count_bytes: [u8; 4] = bytes
        .get(80..84)
        .and_then(|s| s.try_into().ok())
        .ok_or_else(|| bad_data("файл короче заголовка STL"))?;
    let count = u32::from_le_bytes(count_bytes) as usize;
    let expected = count
        .checked_mul(50)
        .and_then(|n| n.checked_add(84))
        .ok_or_else(|| bad_data("неверное число треугольников"))?;
    if bytes.len() != expected {
        return Err(bad_data("длина файла не совпадает с числом треугольников"));
    }
    let word =
        |at: usize| u32::from_le_bytes([bytes[at], bytes[at + 1], bytes[at + 2], bytes[at + 3]]);
    let mut mesh = Mesh::default();
    let mut index: HashMap<[u32; 3], u32> = HashMap::new();
    for t in 0..count {
        let record = 84 + t * 50;
        let mut tri = [0u32; 3];
        for (k, slot) in tri.iter_mut().enumerate() {
            let at = record + 12 + k * 12;
            let bits: [u32; 3] = std::array::from_fn(|c| match word(at + c * 4) {
                0x8000_0000 => 0, // −0 → +0
                b => b,
            });
            let next = mesh.vertices.len() as u32;
            *slot = *index.entry(bits).or_insert_with(|| {
                let [x, y, z] = bits.map(|b| f64::from(f32::from_bits(b)));
                mesh.vertices.push(Vec3::new(x, y, z));
                next
            });
        }
        mesh.triangles.push(tri);
    }
    Ok(mesh)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn brick() -> Mesh {
        Mesh::cuboid(Vec3::new(0.0, 0.0, 0.0), Vec3::new(20.0, 30.0, 3.0))
    }

    #[test]
    fn binary_round_trip_keeps_a_closed_body() {
        let mut bytes = Vec::new();
        write_binary(&brick(), &mut bytes).unwrap();
        assert_eq!(bytes.len(), 84 + 12 * 50);
        assert!(!bytes.starts_with(b"solid"));
        let back = read_binary(&bytes).unwrap();
        assert_eq!(back.vertices.len(), 8);
        assert_eq!(back.check_closed(), Ok(()));
        assert_eq!(back.signed_volume(), 1800.0);
    }

    #[test]
    fn ascii_has_every_facet() {
        let mut out = Vec::new();
        write_ascii(&brick(), "брусок 20x30", &mut out).unwrap();
        let text = String::from_utf8(out).unwrap();
        assert!(text.starts_with("solid брусок_20x30\n"));
        assert_eq!(text.matches("facet normal").count(), 12);
        assert!(text.trim_end().ends_with("endsolid брусок_20x30"));
    }

    #[test]
    fn broken_input_is_rejected() {
        let mut m = brick();
        m.triangles[0] = [0, 1, 42];
        assert!(write_binary(&m, Vec::new()).is_err());
        assert!(read_binary(&[0u8; 10]).is_err());
        let mut bytes = Vec::new();
        write_binary(&brick(), &mut bytes).unwrap();
        bytes.pop();
        assert!(read_binary(&bytes).is_err());
    }
}
