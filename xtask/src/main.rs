//! Стенд проекта одной командой. Он написан на Rust, поэтому одинаково
//! работает на Windows, Linux и macOS: нужен только сам Rust.
//!
//! ```text
//! cargo xtask check           быстрый прогон — после каждой правки
//! cargo xtask check --full    полный прогон — перед выпуском версии
//! cargo xtask night [часы]    ночной поиск ошибок, по умолчанию 8 часов
//! ```
//!
//! Днём сид постоянный, чтобы прогоны повторялись. Ночью каждый раунд
//! получает новый сид от текущего времени, а сиды и найденные ошибки
//! пишутся в отчёт `reports/night-ГГГГ-ММ-ДД-ЧЧММ.md`.

use std::collections::BTreeMap;
use std::env;
use std::fs;
use std::path::{Path, PathBuf};
use std::process::{Command, ExitCode};
use std::time::{Duration, Instant, SystemTime, UNIX_EPOCH};

/// Случаев на тест в полном и ночном прогонах, если не задано YADRO_CASES.
const FULL_CASES: &str = "100000";
/// Сколько примеров одной ошибки сохранять в отчёте.
const EXAMPLES_PER_TEST: usize = 3;

fn main() -> ExitCode {
    let args: Vec<String> = env::args().skip(1).collect();
    let result = match args.first().map(String::as_str) {
        Some("check") => check(args.iter().any(|a| a == "--full")),
        Some("night") => match args.get(1).map(|h| h.replace(',', ".").parse::<f64>()) {
            None => night(8.0),
            Some(Ok(hours)) if hours > 0.0 => night(hours),
            Some(_) => Err("число часов — положительное число, например 7.5".into()),
        },
        _ => {
            eprintln!(
                "Как запускать:\n  cargo xtask check          быстрый прогон\n  \
                 cargo xtask check --full   полный прогон\n  \
                 cargo xtask night [часы]   ночной поиск ошибок"
            );
            return ExitCode::from(2);
        }
    };
    match result {
        Ok(()) => ExitCode::SUCCESS,
        Err(msg) => {
            eprintln!("\n✗ {msg}");
            ExitCode::FAILURE
        }
    }
}

/// Корень проекта — папка над xtask.
fn root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .expect("xtask лежит внутри проекта")
        .to_path_buf()
}

fn cargo(args: &[&str]) -> Command {
    let mut cmd = Command::new(env::var("CARGO").unwrap_or_else(|_| "cargo".into()));
    cmd.args(args).current_dir(root());
    cmd
}

fn run(title: &str, mut cmd: Command) -> Result<(), String> {
    println!("\n== {title}");
    let status = cmd
        .status()
        .map_err(|e| format!("не удалось запустить cargo: {e}"))?;
    if status.success() {
        Ok(())
    } else {
        Err(format!("не прошёл шаг «{title}»"))
    }
}

fn cases_or(default: &str) -> String {
    env::var("YADRO_CASES").unwrap_or_else(|_| default.into())
}

fn check(full: bool) -> Result<(), String> {
    let started = Instant::now();
    run(
        "Формат кода (rustfmt)",
        cargo(&["fmt", "--all", "--", "--check"]),
    )?;
    run(
        "Анализатор (clippy)",
        cargo(&[
            "clippy",
            "--workspace",
            "--all-targets",
            "--quiet",
            "--",
            "-D",
            "warnings",
        ]),
    )?;
    let mut test = cargo(&["test", "--workspace", "--quiet"]);
    let title = if full {
        let cases = cases_or(FULL_CASES);
        test.arg("--release").env("YADRO_CASES", &cases);
        format!("Тесты (полный прогон, случаев на тест: {cases})")
    } else {
        "Тесты (быстрый прогон)".to_string()
    };
    run(&title, test)?;
    let mode = if full {
        "полный"
    } else {
        "быстрый"
    };
    println!(
        "\nВсё зелёное: {mode} прогон за {}.",
        duration(started.elapsed())
    );
    Ok(())
}

/// SplitMix64: из одного числа — хорошо перемешанное другое.
fn mix(mut z: u64) -> u64 {
    z = z.wrapping_add(0x9E37_79B9_7F4A_7C15);
    z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    z ^ (z >> 31)
}

/// Дата и время UTC по секундам с 1970 года: (год, месяц, день, часы, минуты).
/// Алгоритм Хиннанта, без библиотек.
fn utc(secs: u64) -> (i64, u32, u32, u32, u32) {
    let z = (secs / 86_400) as i64 + 719_468;
    let rem = secs % 86_400;
    let era = z.div_euclid(146_097);
    let doe = z.rem_euclid(146_097);
    let yoe = (doe - doe / 1460 + doe / 36_524 - doe / 146_096) / 365;
    let doy = doe - (365 * yoe + yoe / 4 - yoe / 100);
    let mp = (5 * doy + 2) / 153;
    let day = (doy - (153 * mp + 2) / 5 + 1) as u32;
    let month = (if mp < 10 { mp + 3 } else { mp - 9 }) as u32;
    let year = yoe + era * 400 + i64::from(month <= 2);
    (
        year,
        month,
        day,
        (rem / 3600) as u32,
        (rem % 3600 / 60) as u32,
    )
}

fn duration(d: Duration) -> String {
    let s = d.as_secs();
    match s {
        0..60 => format!("{s} с"),
        60..3600 => format!("{} мин {} с", s / 60, s % 60),
        _ => format!("{} ч {} мин", s / 3600, s % 3600 / 60),
    }
}

/// Одна найденная ошибка: в каком раунде, с каким сидом и что напечатал тест.
struct Example {
    round: u64,
    round_seed: u64,
    lines: Vec<String>,
}

/// Ошибки одного теста за ночь: сколько раз упал и первые примеры.
#[derive(Default)]
struct TestFailures {
    rounds: u64,
    examples: Vec<Example>,
}

/// Разобрать вывод `cargo test`: какие тесты упали и что они напечатали.
///
/// Строки «упал случай» и «Повторить только его» сохраняются всегда;
/// остальные — первые 12, чтобы длинный вывод не раздувал отчёт.
fn failed_tests(output: &str) -> Vec<(String, Vec<String>)> {
    let mut found: Vec<(String, Vec<String>)> = Vec::new();
    for line in output.lines() {
        if let Some(name) = line
            .strip_prefix("---- ")
            .and_then(|l| l.strip_suffix(" stdout ----"))
        {
            found.push((name.to_string(), Vec::new()));
        } else if line.starts_with("failures:") || line.starts_with("test result:") {
            // Конец блока с выводом упавшего теста.
            found.push((String::new(), Vec::new()));
        } else if let Some((name, lines)) = found.last_mut()
            && !name.is_empty()
            && !line.trim().is_empty()
        {
            let key = line.contains("упал случай") || line.contains("Повторить только его");
            if key || lines.len() < 12 {
                lines.push(line.to_string());
            }
        }
    }
    found.retain(|(name, _)| !name.is_empty());
    found
}

struct Night {
    started_secs: u64,
    night_seed: u64,
    cases: String,
    rounds: u64,
    elapsed: Duration,
    failures: BTreeMap<String, TestFailures>,
}

impl Night {
    fn report(&self) -> String {
        let (y, mo, d, h, mi) = utc(self.started_secs);
        let verdict = if self.failures.is_empty() {
            "**всё зелёное**".to_string()
        } else {
            format!("**упали тесты: {}**", self.failures.len())
        };
        let mut s = format!(
            "# Ночной поиск ошибок — {y:04}-{mo:02}-{d:02} {h:02}:{mi:02} UTC\n\n\
             - Итог: {verdict}\n\
             - Раундов: {}, случаев на тест в раунде: {}, длительность: {}\n\
             - Сид ночи: {}; сид раунда r — это mix(сид ночи + r)\n",
            self.rounds,
            self.cases,
            duration(self.elapsed),
            self.night_seed
        );
        for (test, f) in &self.failures {
            s.push_str(&format!("\n## {test} — упал в раундах: {}\n", f.rounds));
            for ex in &f.examples {
                s.push_str(&format!(
                    "\nРаунд {}, сид раунда {}:\n\n```text\n{}\n```\n",
                    ex.round,
                    ex.round_seed,
                    ex.lines.join("\n")
                ));
            }
        }
        s
    }
}

fn night(hours: f64) -> Result<(), String> {
    let now = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .unwrap_or_default();
    let started = Instant::now();
    let deadline = started + Duration::from_secs_f64(hours * 3600.0);
    let mut state = Night {
        started_secs: now.as_secs(),
        night_seed: mix(now.as_nanos() as u64),
        cases: cases_or(FULL_CASES),
        rounds: 0,
        elapsed: Duration::ZERO,
        failures: BTreeMap::new(),
    };
    let (y, mo, d, h, mi) = utc(state.started_secs);
    let dir = root().join("reports");
    fs::create_dir_all(&dir).map_err(|e| format!("не удалось создать папку отчётов: {e}"))?;
    let path = dir.join(format!("night-{y:04}-{mo:02}-{d:02}-{h:02}{mi:02}.md"));

    run(
        "Сборка тестов (release)",
        cargo(&["test", "--workspace", "--release", "--no-run", "--quiet"]),
    )?;
    println!(
        "\n== Ночной поиск ошибок: {hours} ч, сид ночи {}\nОтчёт: {}",
        state.night_seed,
        path.display()
    );

    while Instant::now() < deadline {
        let round_seed = mix(state.night_seed.wrapping_add(state.rounds));
        let out = cargo(&["test", "--workspace", "--release", "--quiet"])
            .env("YADRO_SEED", round_seed.to_string())
            .env("YADRO_CASES", &state.cases)
            // Без трассировки стека: в отчёте нужны сообщение и команда повтора.
            .env("RUST_BACKTRACE", "0")
            .output()
            .map_err(|e| format!("не удалось запустить cargo: {e}"))?;
        if !out.status.success() {
            let text = format!(
                "{}\n{}",
                String::from_utf8_lossy(&out.stdout),
                String::from_utf8_lossy(&out.stderr)
            );
            let mut tests = failed_tests(&text);
            if tests.is_empty() {
                tests.push((
                    "(не разобрать вывод)".into(),
                    text.lines().take(12).map(String::from).collect(),
                ));
            }
            for (test, lines) in tests {
                println!("Раунд {}: упал {test}", state.rounds);
                let entry = state.failures.entry(test).or_default();
                entry.rounds += 1;
                if entry.examples.len() < EXAMPLES_PER_TEST {
                    entry.examples.push(Example {
                        round: state.rounds,
                        round_seed,
                        lines,
                    });
                }
            }
        } else if state.rounds.is_multiple_of(50) {
            println!("Раунд {}: всё зелёное", state.rounds);
        }
        state.rounds += 1;
        state.elapsed = started.elapsed();
        // Отчёт переписывается после каждого раунда: если прогон прервут,
        // всё найденное уже в файле.
        fs::write(&path, state.report()).map_err(|e| format!("не удалось записать отчёт: {e}"))?;
    }

    println!("\n{}", state.report());
    if state.failures.is_empty() {
        Ok(())
    } else {
        Err(format!(
            "найдены ошибки, подробности в отчёте {}",
            path.display()
        ))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn dates_are_right() {
        assert_eq!(utc(0), (1970, 1, 1, 0, 0));
        assert_eq!(utc(951_782_400), (2000, 2, 29, 0, 0));
        assert_eq!(utc(1_793_361_600 + 23 * 60), (2026, 10, 30, 12, 23));
    }

    #[test]
    fn failed_tests_are_parsed() {
        let out = "\
running 3 tests
..F
failures:

---- predicates::tests::orient2d_symmetry stdout ----

thread 'x' panicked at src/testkit.rs:163:13:
orient2d_symmetry: упал случай 21 из 100000.
Повторить только его: YADRO_SEED=42 YADRO_CASES=1 cargo test -p yadro-core orient2d_symmetry

failures:
    predicates::tests::orient2d_symmetry

test result: FAILED. 2 passed; 1 failed";
        let found = failed_tests(out);
        assert_eq!(found.len(), 1);
        assert_eq!(found[0].0, "predicates::tests::orient2d_symmetry");
        assert!(found[0].1.iter().any(|l| l.contains("YADRO_SEED=42")));
    }

    #[test]
    fn seeds_differ_between_rounds() {
        let base = mix(12345);
        assert_ne!(mix(base), mix(base.wrapping_add(1)));
    }
}
