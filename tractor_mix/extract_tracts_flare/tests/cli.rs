use std::fs;
use std::path::PathBuf;
use std::process::Command;

use tempfile::TempDir;

fn rust_bin() -> PathBuf {
    std::env::var_os("CARGO_BIN_EXE_extract-tracts-flare")
        .or_else(|| std::env::var_os("CARGO_BIN_EXE_extract_tracts_flare"))
        .map(PathBuf::from)
        .unwrap_or_else(|| {
            PathBuf::from(env!("CARGO_MANIFEST_DIR"))
                .join("target")
                .join(if cfg!(debug_assertions) {
                    "debug"
                } else {
                    "release"
                })
                .join("extract-tracts-flare")
        })
}

fn testdata(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("testdata")
        .join(name)
}

fn run(args: &[&str]) -> std::process::Output {
    Command::new(rust_bin()).args(args).output().unwrap()
}

#[test]
fn missing_output_dir() {
    let tmp = TempDir::new().unwrap();
    let missing = tmp.path().join("nope");
    let out = run(&[
        "--vcf",
        testdata("tiny.vcf").to_str().unwrap(),
        "--num-ancs",
        "5",
        "--output-dir",
        missing.to_str().unwrap(),
    ]);
    assert!(!out.status.success());
    let err = String::from_utf8_lossy(&out.stderr);
    assert!(err.contains("does not exist"), "{err}");
}

#[test]
fn bad_extension() {
    let tmp = TempDir::new().unwrap();
    let f = tmp.path().join("x.bcf");
    fs::write(&f, b"not a vcf").unwrap();
    let out_dir = tmp.path().join("out");
    fs::create_dir(&out_dir).unwrap();
    let out = run(&[
        "--vcf",
        f.to_str().unwrap(),
        "--num-ancs",
        "5",
        "--output-dir",
        out_dir.to_str().unwrap(),
    ]);
    assert!(!out.status.success());
    let err = String::from_utf8_lossy(&out.stderr);
    assert!(
        err.contains("Unexpected file extension") || err.contains("unexpected file extension"),
        "{err}"
    );
}

#[test]
fn unphased_errors_with_context() {
    let tmp = TempDir::new().unwrap();
    let out_dir = tmp.path().join("out");
    fs::create_dir(&out_dir).unwrap();
    let out = run(&[
        "--vcf",
        testdata("unphased.vcf").to_str().unwrap(),
        "--num-ancs",
        "2",
        "--output-dir",
        out_dir.to_str().unwrap(),
    ]);
    assert!(!out.status.success());
    let err = String::from_utf8_lossy(&out.stderr);
    assert!(err.contains("fewer than 4"), "{err}");
    assert!(err.contains("S1"), "{err}");
}

#[test]
fn truncated_genotype_errors() {
    let tmp = TempDir::new().unwrap();
    let out_dir = tmp.path().join("out");
    fs::create_dir(&out_dir).unwrap();
    let out = run(&[
        "--vcf",
        testdata("truncated.vcf").to_str().unwrap(),
        "--num-ancs",
        "2",
        "--output-dir",
        out_dir.to_str().unwrap(),
    ]);
    assert!(!out.status.success());
    let err = String::from_utf8_lossy(&out.stderr);
    assert!(err.contains("fewer than 4"), "{err}");
}

#[test]
fn missing_keep_ids() {
    let tmp = TempDir::new().unwrap();
    let out_dir = tmp.path().join("out");
    fs::create_dir(&out_dir).unwrap();
    let keep = tmp.path().join("keep.txt");
    fs::write(&keep, "NOT_IN_VCF\n").unwrap();
    let out = run(&[
        "--vcf",
        testdata("tiny.vcf").to_str().unwrap(),
        "--num-ancs",
        "5",
        "--output-dir",
        out_dir.to_str().unwrap(),
        "--samples",
        keep.to_str().unwrap(),
    ]);
    assert!(!out.status.success());
    let err = String::from_utf8_lossy(&out.stderr);
    assert!(err.contains("missing"), "{err}");
}

#[test]
fn missing_vcf() {
    let tmp = TempDir::new().unwrap();
    let out_dir = tmp.path().join("out");
    fs::create_dir(&out_dir).unwrap();
    let out = run(&[
        "--vcf",
        tmp.path().join("absent.vcf").to_str().unwrap(),
        "--num-ancs",
        "5",
        "--output-dir",
        out_dir.to_str().unwrap(),
    ]);
    assert!(!out.status.success());
}
