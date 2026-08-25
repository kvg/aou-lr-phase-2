use std::path::PathBuf;

use propagate_annotations::{propagate_annotations, PropagateConfig};

fn base_cfg(dir: &PathBuf, out_path: PathBuf) -> PropagateConfig {
    PropagateConfig {
        snv_indel_vcf: Some(dir.join("snv_indel.vcf")),
        sv_vcf: Some(dir.join("sv.vcf")),
        flare_vcf: Some(dir.join("flare.vcf")),
        output: out_path,
        compress_output: false,
        stats_json: None,
        stats_tsv: None,
        min_match_rate_snv_indel: None,
        min_match_rate_sv: None,
        min_match_rate_flare: None,
        threads: 2,
        dense_format: false,
        region: None,
        progress_every: 100_000,
    }
}

#[test]
fn unions_snv_sv_and_flare() {
    let dir = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("testdata");
    let out = tempfile::NamedTempFile::new().unwrap();
    let out_path = out.path().with_extension("vcf");
    let stats_json = out.path().with_extension("stats.json");
    let stats_tsv = out.path().with_extension("stats.tsv");

    let mut cfg = base_cfg(&dir, out_path.clone());
    cfg.stats_json = Some(stats_json.clone());
    cfg.stats_tsv = Some(stats_tsv.clone());

    let stats = propagate_annotations(&cfg).expect("propagate");

    assert_eq!(stats.mode, "union_stream");
    assert_eq!(stats.union_sites, stats.output_sites);
    assert_eq!(stats.union_sites, 4); // snv 100,200,300 + sv 500
    assert_eq!(stats.chromosome, "chr22");
    assert!((stats.snv_flare_overlap_rate.unwrap() - 1.0).abs() < 1e-9);

    let snv = stats
        .sources
        .iter()
        .find(|s| s.label == "snv_indel")
        .expect("snv stats");
    assert_eq!(snv.source_sites, 3);
    assert_eq!(snv.union_sites_with_source, 3);

    let sv = stats.sources.iter().find(|s| s.label == "sv").expect("sv stats");
    assert_eq!(sv.source_sites, 1);
    assert_eq!(sv.union_sites_with_source, 1);

    assert_eq!(stats.sites_by_n_sources.get(&1), Some(&1)); // 500 sv-only
    assert_eq!(stats.sites_by_n_sources.get(&2), Some(&3)); // snv+flare

    let text = std::fs::read_to_string(&out_path).unwrap();
    assert!(text.contains("##FORMAT=<ID=GQ,"));
    assert!(text.contains("##FORMAT=<ID=AN1,"));
    assert!(text.contains("chr22\t500\t"));
    assert!(!text.lines().any(|l| l.starts_with("chr22\t400\t")));

    // Sparse FORMAT: SV-only site should not carry FLARE columns.
    let sv_line = text
        .lines()
        .find(|l| l.starts_with("chr22\t500\t"))
        .expect("sv line");
    let fmt = sv_line.split('\t').nth(8).unwrap();
    assert!(fmt.contains("GQ") || fmt.contains("DP") || fmt.contains("AD"));
    assert!(!fmt.contains("AN1"));

    // SNV+FLARE site should carry GQ/DP/AD + AN* .
    let snv_line = text
        .lines()
        .find(|l| l.starts_with("chr22\t100\t"))
        .expect("snv line");
    let fmt = snv_line.split('\t').nth(8).unwrap();
    assert!(fmt.contains("GQ"));
    assert!(fmt.contains("AN1"));
    assert!(!fmt.split(':').any(|t| t == "GT"));
    assert!(!fmt.split(':').any(|t| t == "PL" || t == "RNC" || t == "FT"));

    assert!(stats_json.exists());
    assert!(stats_tsv.exists());
}

#[test]
fn min_flare_overlap_fails_when_too_low() {
    let dir = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("testdata");
    let out = tempfile::NamedTempFile::new().unwrap();
    let err = propagate_annotations(&PropagateConfig {
        snv_indel_vcf: Some(dir.join("snv_indel.vcf")),
        sv_vcf: None,
        flare_vcf: Some(dir.join("flare_partial.vcf")),
        output: out.path().with_extension("vcf"),
        compress_output: false,
        stats_json: None,
        stats_tsv: None,
        min_match_rate_snv_indel: None,
        min_match_rate_sv: None,
        min_match_rate_flare: Some(0.95),
        threads: 1,
        dense_format: false,
        region: None,
        progress_every: 100_000,
    })
    .unwrap_err();
    assert!(err.to_string().contains("flare_on_snv"));
}

#[test]
fn writes_bcf_via_bcftools_when_available() {
    if std::process::Command::new("bcftools")
        .arg("--version")
        .output()
        .map(|o| o.status.success())
        .unwrap_or(false)
        == false
    {
        eprintln!("skipping BCF test: bcftools not on PATH");
        return;
    }
    let dir = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("testdata");
    let tmp = tempfile::tempdir().unwrap();
    let out_path = tmp.path().join("out.bcf");
    let mut cfg = base_cfg(&dir, out_path.clone());
    cfg.sv_vcf = None;
    propagate_annotations(&cfg).expect("propagate bcf");
    assert!(out_path.exists());
    let probe = std::process::Command::new("bcftools")
        .args(["view", "-H", out_path.to_str().unwrap()])
        .output()
        .expect("bcftools view");
    assert!(probe.status.success());
    let n = String::from_utf8_lossy(&probe.stdout).lines().count();
    assert_eq!(n, 3);
}
