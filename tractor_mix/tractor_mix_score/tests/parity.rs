use std::collections::HashMap;
use std::fs;
use std::path::{Path, PathBuf};
use std::process::Command;

use tempfile::TempDir;
use tractor_mix_score::{run_score, ScoreConfig};

fn crate_root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
}

fn rust_bin() -> PathBuf {
    std::env::var_os("CARGO_BIN_EXE_tractor-mix-score")
        .map(PathBuf::from)
        .unwrap_or_else(|| {
            crate_root()
                .join("target")
                .join(if cfg!(debug_assertions) {
                    "debug"
                } else {
                    "release"
                })
                .join("tractor-mix-score")
        })
}

fn rscript_available() -> bool {
    Command::new("Rscript")
        .arg("-e")
        .arg("quit(status=0)")
        .status()
        .map(|s| s.success())
        .unwrap_or(false)
}

fn gmmat_available() -> bool {
    Command::new("Rscript")
        .arg("-e")
        .arg("suppressPackageStartupMessages(library(GMMAT)); quit(status=0)")
        .status()
        .map(|s| s.success())
        .unwrap_or(false)
}

fn run_gen_fixture(out: &Path) {
    let status = Command::new("Rscript")
        .arg(crate_root().join("testdata/gen_oracle.R"))
        .arg(out)
        .status()
        .expect("Rscript gen_oracle");
    assert!(status.success(), "gen_oracle.R failed");
}

fn run_r_oracle(
    null_rds: &Path,
    dosage_files: &[PathBuf],
    out: &Path,
    ac_threshold: i32,
) {
    let mut cmd = Command::new("Rscript");
    cmd.arg(crate_root().join("testdata/run_oracle.R"))
        .arg("--null-rds")
        .arg(null_rds)
        .arg("--out")
        .arg(out)
        .arg("--ac-threshold")
        .arg(ac_threshold.to_string())
        .arg("--n-core")
        .arg("1")
        .arg("--chunk-size")
        .arg("16")
        .arg("--oracle-r")
        .arg(crate_root().join("oracle/TractorMix.score.R"))
        .arg("--dosage-files");
    for f in dosage_files {
        cmd.arg(f);
    }
    let status = cmd.status().expect("R oracle");
    assert!(status.success(), "R oracle failed");
}

fn run_rust_score(
    null_export: &Path,
    dosage_files: &[PathBuf],
    out: &Path,
    ac_threshold: i32,
    threads: usize,
) {
    let cfg = ScoreConfig {
        null_export: null_export.to_path_buf(),
        dosage_files: dosage_files.to_vec(),
        out_tsv: out.to_path_buf(),
        ac_threshold,
        chunk_size: 16,
        threads,
    };
    run_score(&cfg).expect("rust score");
}

fn parse_tsv(path: &Path) -> (Vec<String>, Vec<HashMap<String, String>>) {
    let text = fs::read_to_string(path).unwrap();
    let mut lines = text.lines();
    let header: Vec<String> = lines
        .next()
        .unwrap()
        .split('\t')
        .map(String::from)
        .collect();
    let mut rows = Vec::new();
    for line in lines {
        if line.trim().is_empty() {
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        let mut row = HashMap::new();
        for (i, col) in header.iter().enumerate() {
            row.insert(col.clone(), fields.get(i).unwrap_or(&"").to_string());
        }
        rows.push(row);
    }
    (header, rows)
}

fn r_round5(s: &str) -> String {
    if s == "NA" || s.is_empty() {
        return "NA".to_string();
    }
    let x: f64 = s.parse().unwrap();
    format!("{}", (x * 1e5).round() / 1e5)
}

fn r_signif5(s: &str) -> String {
    if s == "NA" || s.is_empty() {
        return "NA".to_string();
    }
    let x: f64 = s.parse().unwrap();
    if x == 0.0 {
        return "0".to_string();
    }
    let exp = x.abs().log10().floor() as i32;
    let factor = 10_f64.powi(5 - 1 - exp);
    format!("{}", (x * factor).round() / factor)
}

fn compare_tsv(r_path: &Path, rs_path: &Path) {
    let (h1, r_rows) = parse_tsv(r_path);
    let (h2, rs_rows) = parse_tsv(rs_path);
    assert_eq!(h1, h2, "header mismatch");
    assert_eq!(r_rows.len(), rs_rows.len(), "row count mismatch");

    let numeric_round5 = |c: &str| c == "Chi2" || c.starts_with("Eff_anc") || c.starts_with("SE_anc");
    let numeric_signif5 = |c: &str| c == "P" || c.starts_with("Pval_anc");

    for (i, (r, rs)) in r_rows.iter().zip(rs_rows.iter()).enumerate() {
        for col in &h1 {
            let a = r.get(col).unwrap();
            let b = rs.get(col).unwrap();
            if numeric_round5(col) {
                assert_eq!(r_round5(a), r_round5(b), "row {i} col {col}: {a} vs {b}");
            } else if numeric_signif5(col) {
                assert_eq!(r_signif5(a), r_signif5(b), "row {i} col {col}: {a} vs {b}");
            } else {
                assert_eq!(a, b, "row {i} col {col}: {a} vs {b}");
            }
        }
    }
}

#[test]
#[ignore = "requires Rscript + GMMAT"]
fn parity_oracle_fixture() {
    if !rscript_available() || !gmmat_available() {
        eprintln!("skip: Rscript/GMMAT not available");
        return;
    }
    let tmp = TempDir::new().unwrap();
    let fixture = tmp.path().join("fixture");
    run_gen_fixture(&fixture);

    let null_rds = fixture.join("null.rds");
    let null_export = fixture.join("null_export");
    let dosages: Vec<PathBuf> = (0..5)
        .map(|i| {
            fixture.join(format!(
                "dosages/anc_{:02}.dosage.txt.gz",
                i
            ))
        })
        .collect();

    let r_out = tmp.path().join("r.tsv");
    let rs_out = tmp.path().join("rust.tsv");
    run_r_oracle(&null_rds, &dosages, &r_out, 10);
    run_rust_score(&null_export, &dosages, &rs_out, 10, 1);
    compare_tsv(&r_out, &rs_out);
}

#[test]
#[ignore = "requires Rscript + GMMAT"]
fn parallel_matches_serial_on_fixture() {
    if !rscript_available() || !gmmat_available() {
        eprintln!("skip: Rscript/GMMAT not available");
        return;
    }
    let tmp = TempDir::new().unwrap();
    let fixture = tmp.path().join("fixture");
    run_gen_fixture(&fixture);

    let null_export = fixture.join("null_export");
    let dosages: Vec<PathBuf> = (0..5)
        .map(|i| fixture.join(format!("dosages/anc_{:02}.dosage.txt.gz", i)))
        .collect();

    let out1 = tmp.path().join("t1.tsv");
    let out8 = tmp.path().join("t8.tsv");
    run_rust_score(&null_export, &dosages, &out1, 10, 1);
    run_rust_score(&null_export, &dosages, &out8, 10, 8);
    compare_tsv(&out1, &out8);
    assert_eq!(
        fs::read(&out1).unwrap(),
        fs::read(&out8).unwrap(),
        "parallel output should be byte-identical"
    );
}

#[test]
fn cli_help() {
    let output = Command::new(rust_bin())
        .arg("--help")
        .output()
        .expect("run bin");
    assert!(output.status.success());
    let stdout = String::from_utf8_lossy(&output.stdout);
    assert!(stdout.contains("null-export"));
}
