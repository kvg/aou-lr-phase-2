use std::collections::HashMap;
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

fn parse_dosage(path: &PathBuf) -> (Vec<String>, HashMap<String, Vec<i64>>) {
    let text = fs::read_to_string(path).unwrap();
    let mut lines = text.lines();
    let header: Vec<String> = lines
        .next()
        .unwrap()
        .split('\t')
        .map(str::to_string)
        .collect();
    let mut rows = HashMap::new();
    for line in lines {
        if line.is_empty() {
            continue;
        }
        let f: Vec<&str> = line.split('\t').collect();
        let id = f[2].to_string();
        let vals: Vec<i64> = f[5..].iter().map(|s| s.parse().unwrap()).collect();
        rows.insert(id, vals);
    }
    (header, rows)
}

#[test]
fn ru_dosage_c_times_ancestry_and_comparison_encodings() {
    let tmp = TempDir::new().unwrap();
    let out = tmp.path().join("out");
    fs::create_dir(&out).unwrap();
    let status = Command::new(rust_bin())
        .args([
            "--vcf",
            testdata("ru_dosage.vcf").to_str().unwrap(),
            "--num-ancs",
            "3",
            "--output-dir",
            out.to_str().unwrap(),
        ])
        .status()
        .unwrap();
    assert!(status.success());

    let (_, a0) = parse_dosage(&out.join("ru_dosage.anc0.dosage.txt"));
    let (_, a1) = parse_dosage(&out.join("ru_dosage.anc1.dosage.txt"));
    let (_, a2) = parse_dosage(&out.join("ru_dosage.anc2.dosage.txt"));

    // Classic SNV (allele==1 only), one locus row.
    assert_eq!(a0["snv1"], vec![0, 0, 0]);
    assert_eq!(a1["snv1"], vec![1, 0, 2]);
    assert_eq!(a2["snv1"], vec![0, 0, 0]);

    // INS: C = 10 (REF) or 13 (ALT). Multi-allelic VNTR stays one row.
    assert_eq!(a0["vntr_ins"], vec![10, 26, 0]);
    assert_eq!(a1["vntr_ins"], vec![13, 0, 20]);
    assert_eq!(a2["vntr_ins"], vec![0, 0, 0]);

    // DEL: C = 20 (REF) or 15 (ALT).
    assert_eq!(a0["vntr_del"], vec![20, 0, 15]);
    assert_eq!(a1["vntr_del"], vec![15, 0, 15]);
    assert_eq!(a2["vntr_del"], vec![0, 40, 0]);

    // Multi-allelic: C = 6 / 7 / 9; still one locus in the dosage file.
    assert_eq!(a0["vntr_multi"], vec![7, 0, 12]);
    assert_eq!(a1["vntr_multi"], vec![9, 13, 0]);
    assert_eq!(a2["vntr_multi"], vec![0, 0, 0]);

    // extra-units-vs-REF (no CN_REF): REF=0, ALT=1.
    assert_eq!(a0["rel_ins"], vec![1, 0, 0]);
    assert_eq!(a1["rel_ins"], vec![0, 2, 0]);
    assert_eq!(a2["rel_ins"], vec![0, 0, 0]);

    let (coll_hdr, c0) = parse_dosage(&out.join("ru_dosage.collapse.anc0.dosage.txt"));
    let (_, c1) = parse_dosage(&out.join("ru_dosage.collapse.anc1.dosage.txt"));
    let (_, c2) = parse_dosage(&out.join("ru_dosage.collapse.anc2.dosage.txt"));
    assert_eq!(&coll_hdr[..5], &["CHROM", "POS", "ID", "REF", "ALT"]);
    assert!(!c0.contains_key("snv1"), "collapse is RU_TEST sites only");
    assert_eq!(c0["vntr_ins"], vec![0, 2, 0]);
    assert_eq!(c1["vntr_ins"], vec![1, 0, 0]);
    assert_eq!(c0["vntr_del"], vec![0, 0, 1]);
    assert_eq!(c1["vntr_del"], vec![1, 0, 1]);
    assert_eq!(c2["vntr_del"], vec![0, 0, 0]);
    assert_eq!(c0["vntr_multi"], vec![1, 0, 0]);
    assert_eq!(c1["vntr_multi"], vec![1, 1, 0]);
    assert_eq!(c0["rel_ins"], vec![1, 0, 0]);
    assert_eq!(c1["rel_ins"], vec![0, 2, 0]);
    assert_eq!(c2["rel_ins"], vec![0, 0, 0]);

    let (_, s0) = parse_dosage(&out.join("ru_dosage.split.anc0.dosage.txt"));
    let (_, s1) = parse_dosage(&out.join("ru_dosage.split.anc1.dosage.txt"));
    assert_eq!(s0["vntr_multi:1"], vec![1, 0, 0]);
    assert_eq!(s1["vntr_multi:1"], vec![0, 1, 0]);
    assert_eq!(s0["vntr_multi:2"], vec![0, 0, 0]);
    assert_eq!(s1["vntr_multi:2"], vec![1, 0, 0]);
    assert_eq!(s0["vntr_ins:1"], vec![0, 2, 0]);
    assert_eq!(s1["vntr_ins:1"], vec![1, 0, 0]);

    let hap0 = fs::read_to_string(out.join("ru_dosage.anc0.hapcount.txt")).unwrap();
    let hap_ins: Vec<&str> = hap0
        .lines()
        .find(|l| l.contains("vntr_ins"))
        .unwrap()
        .split('\t')
        .collect();
    // S1 0|1:0:1 → hap anc0=1; S2 1|1:0:0 → hap anc0=2; S3 0|0:1:1 → hap anc0=0
    assert_eq!(&hap_ins[5..], &["1", "2", "0"]);
}
