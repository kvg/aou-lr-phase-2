use std::fs;
use std::io::{BufRead, Read};
use std::path::{Path, PathBuf};
use std::process::Command;

use extract_tracts_flare::extract::{extract_tracts_flare, ExtractConfig};
use flate2::read::MultiGzDecoder;
use tempfile::TempDir;

fn crate_root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
}

fn testdata(name: &str) -> PathBuf {
    crate_root().join("testdata").join(name)
}

fn oracle() -> PathBuf {
    crate_root().join("oracle").join("extract_tracts_flare.py")
}

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

fn read_text(path: &Path) -> String {
    let name = path.file_name().unwrap().to_string_lossy();
    if name.ends_with(".gz") {
        let f = fs::File::open(path).unwrap();
        let mut s = String::new();
        MultiGzDecoder::new(f).read_to_string(&mut s).unwrap();
        s
    } else {
        fs::read_to_string(path).unwrap()
    }
}

fn list_outputs(dir: &Path, prefix: &str) -> Vec<PathBuf> {
    let mut v: Vec<_> = fs::read_dir(dir)
        .unwrap()
        .filter_map(|e| e.ok())
        .map(|e| e.path())
        .filter(|p| {
            let n = p.file_name().unwrap().to_string_lossy();
            n.starts_with(prefix)
                && (n.contains(".dosage.")
                    || n.contains(".hapcount.")
                    || n.contains(".anc") && n.contains(".vcf"))
                && n != "dosage_sample_order.txt"
        })
        .collect();
    v.sort();
    v
}

fn run_python(vcf: &Path, out: &Path, num_ancs: usize, compress: bool, output_vcf: bool) {
    fs::create_dir_all(out).unwrap();
    let mut cmd = Command::new("python3");
    cmd.arg(oracle())
        .arg("--vcf")
        .arg(vcf)
        .arg("--num-ancs")
        .arg(num_ancs.to_string())
        .arg("--output-dir")
        .arg(out);
    if compress {
        cmd.arg("--compress-output");
    }
    if output_vcf {
        cmd.arg("--output-vcf");
    }
    let status = cmd.status().expect("python3");
    assert!(status.success(), "python oracle failed: {status}");
}

fn run_rust_cli(
    vcf: &Path,
    out: &Path,
    num_ancs: usize,
    compress: bool,
    output_vcf: bool,
    samples: Option<&Path>,
) {
    fs::create_dir_all(out).unwrap();
    let mut cmd = Command::new(rust_bin());
    cmd.arg("--vcf")
        .arg(vcf)
        .arg("--num-ancs")
        .arg(num_ancs.to_string())
        .arg("--output-dir")
        .arg(out);
    if compress {
        cmd.arg("--compress-output");
    }
    if output_vcf {
        cmd.arg("--output-vcf");
    }
    if let Some(s) = samples {
        cmd.arg("--samples").arg(s);
    }
    let output = cmd.output().expect("rust bin");
    assert!(
        output.status.success(),
        "rust failed: {}\n{}",
        String::from_utf8_lossy(&output.stderr),
        String::from_utf8_lossy(&output.stdout)
    );
}

fn assert_dirs_match(py_dir: &Path, rs_dir: &Path, prefix: &str) {
    let py_files: Vec<_> = list_outputs(py_dir, prefix)
        .into_iter()
        .filter(|p| {
            let n = p.file_name().unwrap().to_string_lossy();
            n.contains(".dosage.")
                || n.contains(".hapcount.")
                || n.ends_with(".vcf")
                || n.contains(".vcf.")
        })
        .collect();
    let rs_files: Vec<_> = list_outputs(rs_dir, prefix)
        .into_iter()
        .filter(|p| {
            let n = p.file_name().unwrap().to_string_lossy();
            n.contains(".dosage.")
                || n.contains(".hapcount.")
                || n.ends_with(".vcf")
                || n.contains(".vcf.")
        })
        .collect();
    let py_names: Vec<_> = py_files
        .iter()
        .map(|p| p.file_name().unwrap().to_string_lossy().into_owned())
        .collect();
    let rs_names: Vec<_> = rs_files
        .iter()
        .map(|p| p.file_name().unwrap().to_string_lossy().into_owned())
        .collect();
    assert_eq!(py_names, rs_names, "output filename mismatch");
    for (p, r) in py_files.iter().zip(rs_files.iter()) {
        let a = read_text(p);
        let b = read_text(r);
        assert_eq!(
            a,
            b,
            "content mismatch for {}",
            p.file_name().unwrap().to_string_lossy()
        );
    }
}

#[test]
fn parity_tiny_plain() {
    let tmp = TempDir::new().unwrap();
    let py = tmp.path().join("py");
    let rs = tmp.path().join("rs");
    let vcf = testdata("tiny.vcf");
    run_python(&vcf, &py, 5, false, false);
    run_rust_cli(&vcf, &rs, 5, false, false, None);
    assert_dirs_match(&py, &rs, "tiny");
}

#[test]
fn parity_tiny_gzip_and_output_vcf() {
    let tmp = TempDir::new().unwrap();
    let py = tmp.path().join("py");
    let rs = tmp.path().join("rs");
    let vcf = testdata("tiny.vcf");
    run_python(&vcf, &py, 5, true, true);
    run_rust_cli(&vcf, &rs, 5, true, true, None);
    assert_dirs_match(&py, &rs, "tiny");
}

#[test]
fn parity_gzip_input() {
    let tmp = TempDir::new().unwrap();
    let gz = tmp.path().join("tiny.vcf.gz");
    {
        use flate2::write::GzEncoder;
        use flate2::Compression;
        use std::io::Write;
        let mut enc = GzEncoder::new(fs::File::create(&gz).unwrap(), Compression::default());
        enc.write_all(&fs::read(testdata("tiny.vcf")).unwrap())
            .unwrap();
        enc.finish().unwrap();
    }
    let py = tmp.path().join("py");
    let rs = tmp.path().join("rs");
    run_python(&gz, &py, 5, true, false);
    run_rust_cli(&gz, &rs, 5, true, false, None);
    assert_dirs_match(&py, &rs, "tiny");
}

#[test]
fn parity_generated_with_anp() {
    let tmp = TempDir::new().unwrap();
    let vcf = tmp.path().join("gen.vcf");
    let status = Command::new("python3")
        .arg(crate_root().join("testdata/gen_flare_vcf.py"))
        .arg("--out")
        .arg(&vcf)
        .arg("--sites")
        .arg("25")
        .arg("--samples")
        .arg("12")
        .arg("--num-ancs")
        .arg("5")
        .arg("--anp")
        .status()
        .unwrap();
    assert!(status.success());
    let py = tmp.path().join("py");
    let rs = tmp.path().join("rs");
    run_python(&vcf, &py, 5, false, true);
    run_rust_cli(&vcf, &rs, 5, false, true, None);
    assert_dirs_match(&py, &rs, "gen");
}

#[test]
fn samples_projection_reorders_columns() {
    let tmp = TempDir::new().unwrap();
    let keep = tmp.path().join("keep.txt");
    fs::write(&keep, "S3\nS1\n").unwrap();
    let out = tmp.path().join("rs");
    run_rust_cli(&testdata("tiny.vcf"), &out, 5, false, false, Some(&keep));
    let dos = fs::read_to_string(out.join("tiny.anc0.dosage.txt")).unwrap();
    let header = dos.lines().next().unwrap();
    let cols: Vec<_> = header.split('\t').collect();
    assert_eq!(&cols[..5], &["CHROM", "POS", "ID", "REF", "ALT"]);
    assert_eq!(&cols[5..], &["S3", "S1"]);
    let order = fs::read_to_string(out.join("dosage_sample_order.txt")).unwrap();
    assert_eq!(order, "S3\nS1\n");
}

#[test]
fn lib_extract_matches_cli() {
    let tmp = TempDir::new().unwrap();
    let out = tmp.path().join("lib");
    fs::create_dir_all(&out).unwrap();
    extract_tracts_flare(&ExtractConfig {
        vcf: testdata("tiny.vcf"),
        num_ancs: 5,
        output_dir: Some(out.clone()),
        output_vcf: false,
        compress_output: false,
        samples: None,
        threads: 1,
    })
    .unwrap();
    assert!(out.join("tiny.anc0.dosage.txt").exists());
}

#[test]
fn header_column_count() {
    let tmp = TempDir::new().unwrap();
    let out = tmp.path().join("rs");
    run_rust_cli(&testdata("tiny.vcf"), &out, 5, false, false, None);
    let dos = fs::read_to_string(out.join("tiny.anc2.dosage.txt")).unwrap();
    let header = dos.lines().next().unwrap();
    assert_eq!(header.split('\t').count(), 8); // 5 meta + 3 samples
}

#[test]
#[ignore]
fn scale_smoke_1k() {
    let tmp = TempDir::new().unwrap();
    let vcf = tmp.path().join("scale.vcf.gz");
    let status = Command::new("python3")
        .arg(crate_root().join("testdata/gen_flare_vcf.py"))
        .arg("--out")
        .arg(&vcf)
        .arg("--sites")
        .arg("1000")
        .arg("--samples")
        .arg("1000")
        .arg("--num-ancs")
        .arg("5")
        .status()
        .unwrap();
    assert!(status.success());
    let rs = tmp.path().join("rs");
    let start = std::time::Instant::now();
    run_rust_cli(&vcf, &rs, 5, true, false, None);
    let elapsed = start.elapsed();
    assert!(rs.join("scale.anc0.dosage.txt.gz").exists());
    eprintln!("scale 1000x1000 rust {:?}", elapsed);
    let mut first = String::new();
    let f = fs::File::open(rs.join("scale.anc0.dosage.txt.gz")).unwrap();
    std::io::BufReader::new(MultiGzDecoder::new(f))
        .read_line(&mut first)
        .unwrap();
    assert_eq!(first.split('\t').count(), 5 + 1000);
}
