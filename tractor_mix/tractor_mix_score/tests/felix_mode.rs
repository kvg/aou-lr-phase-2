use std::fs::File;
use std::io::Write;
use std::path::Path;

use flate2::write::GzEncoder;
use flate2::Compression;
use tempfile::TempDir;
use tractor_mix_score::null::{write_sigma_i_csc, CscMatrix};
use tractor_mix_score::{run_score, ScoreConfig, ScoreModeArg};

fn write_i32(w: &mut impl Write, x: i32) {
    w.write_all(&x.to_le_bytes()).unwrap();
}
fn write_f64(w: &mut impl Write, x: f64) {
    w.write_all(&x.to_le_bytes()).unwrap();
}

fn write_felix_export(dir: &Path, n: usize, vr: f64) {
    let mut colptr = vec![0usize];
    let mut rowidx = Vec::new();
    let mut values = Vec::new();
    for i in 0..n {
        rowidx.push(i);
        values.push(1.0);
        colptr.push(rowidx.len());
    }
    let sigma = CscMatrix {
        n,
        colptr,
        rowidx,
        values,
    };
    std::fs::create_dir_all(dir).unwrap();
    std::fs::write(
        dir.join("meta.json"),
        format!(
            r#"{{"n":{n},"p":1,"nnz":{n},"family":"gaussian","source":"felix","variance_ratio":{vr},"tau0":1.0,"trait_type":"quantitative"}}"#
        ),
    )
    .unwrap();
    let ids: String = (0..n).map(|i| format!("S{i:03}\n")).collect();
    std::fs::write(dir.join("id_include.txt"), ids).unwrap();
    write_sigma_i_csc(&dir.join("sigma_i.csc.bin"), &sigma).unwrap();

    let mut f = File::create(dir.join("sigma_i_x.bin")).unwrap();
    write_i32(&mut f, n as i32);
    write_i32(&mut f, 1);
    for _ in 0..n {
        write_f64(&mut f, 0.0);
    }
    let mut f = File::create(dir.join("cov.bin")).unwrap();
    write_i32(&mut f, 1);
    write_f64(&mut f, 1.0);
    let mut f = File::create(dir.join("residuals.bin")).unwrap();
    write_i32(&mut f, n as i32);
    for i in 0..n {
        write_f64(&mut f, ((i as f64) * 0.11).sin() * 0.2);
    }
    std::fs::write(dir.join("variance_ratio.txt"), format!("{vr}\n")).unwrap();
}

fn write_f64_dosages(path: &Path, n: usize, n_sites: usize, anc: usize) {
    let mut enc = GzEncoder::new(File::create(path).unwrap(), Compression::default());
    let samples: Vec<String> = (0..n).map(|i| format!("S{i:03}")).collect();
    write!(enc, "CHROM\tPOS\tID\tREF\tALT\t{}\n", samples.join("\t")).unwrap();
    for s in 0..n_sites {
        let mut row = vec![
            "22".into(),
            ((s + 1) * 100).to_string(),
            format!("v{s}"),
            "A".into(),
            "G".into(),
        ];
        for i in 0..n {
            let g = ((s + anc + i) % 5) as f64 * 0.25;
            row.push(format!("{g}"));
        }
        write!(enc, "{}\n", row.join("\t")).unwrap();
    }
    enc.finish().unwrap();
}

#[test]
fn felix_mode_writes_cct_columns_and_parses_f64() {
    let tmp = TempDir::new().unwrap();
    let export = tmp.path().join("null_export");
    write_felix_export(&export, 40, 1.25);
    let d0 = tmp.path().join("anc_00.dosage.txt.gz");
    let d1 = tmp.path().join("anc_01.dosage.txt.gz");
    write_f64_dosages(&d0, 40, 8, 0);
    write_f64_dosages(&d1, 40, 8, 1);
    let out = tmp.path().join("out.tsv");
    run_score(&ScoreConfig {
        null_export: export,
        dosage_files: vec![d0, d1],
        out_tsv: out.clone(),
        ac_threshold: 50,
        chunk_size: 4,
        threads: 2,
        mode: ScoreModeArg::Auto,
        variance_ratio: None,
        min_copy_var: 1e-6,
    })
    .unwrap();
    let text = std::fs::read_to_string(&out).unwrap();
    let mut lines = text.lines();
    let header = lines.next().unwrap();
    assert!(header.contains("P_cct_admixed"));
    assert!(header.contains("P_hom_admixed"));
    assert!(header.contains("BETA_ancALL"));
    assert!(!header.contains("Chi2"));
    let rows: Vec<_> = lines.filter(|l| !l.is_empty()).collect();
    assert_eq!(rows.len(), 8);
    let cols: Vec<_> = rows[0].split('\t').collect();
    assert_eq!(cols.len(), header.split('\t').count());
    assert!(cols.iter().any(|c| *c != "NA"));
}

#[test]
fn auto_mode_gmmat_export_keeps_legacy_columns() {
    let tmp = TempDir::new().unwrap();
    let export = tmp.path().join("null_export");
    let n = 20usize;
    let mut colptr = vec![0usize];
    let mut rowidx = Vec::new();
    let mut values = Vec::new();
    for i in 0..n {
        rowidx.push(i);
        values.push(1.0);
        colptr.push(rowidx.len());
    }
    std::fs::create_dir_all(&export).unwrap();
    std::fs::write(
        export.join("meta.json"),
        format!(r#"{{"n":{n},"p":1,"nnz":{n},"family":"binomial"}}"#),
    )
    .unwrap();
    std::fs::write(
        export.join("id_include.txt"),
        (0..n).map(|i| format!("S{i:03}\n")).collect::<String>(),
    )
    .unwrap();
    write_sigma_i_csc(
        &export.join("sigma_i.csc.bin"),
        &CscMatrix {
            n,
            colptr,
            rowidx,
            values,
        },
    )
    .unwrap();
    let mut f = File::create(export.join("sigma_i_x.bin")).unwrap();
    write_i32(&mut f, n as i32);
    write_i32(&mut f, 1);
    for _ in 0..n {
        write_f64(&mut f, 0.0);
    }
    let mut f = File::create(export.join("cov.bin")).unwrap();
    write_i32(&mut f, 1);
    write_f64(&mut f, 1.0);
    let mut f = File::create(export.join("residuals.bin")).unwrap();
    write_i32(&mut f, n as i32);
    for i in 0..n {
        write_f64(&mut f, (i as f64 - 10.0) * 0.01);
    }

    let d0 = tmp.path().join("anc_00.dosage.txt.gz");
    let d1 = tmp.path().join("anc_01.dosage.txt.gz");
    for (anc, path) in [d0.clone(), d1.clone()].into_iter().enumerate() {
        let mut enc = GzEncoder::new(File::create(path).unwrap(), Compression::default());
        let samples: Vec<String> = (0..n).map(|i| format!("S{i:03}")).collect();
        write!(enc, "CHROM\tPOS\tID\tREF\tALT\t{}\n", samples.join("\t")).unwrap();
        let mut row = vec!["22".into(), "100".into(), "v0".into(), "A".into(), "G".into()];
        for i in 0..n {
            row.push((((i + anc) % 3) as u8).to_string());
        }
        write!(enc, "{}\n", row.join("\t")).unwrap();
        enc.finish().unwrap();
    }

    let out = tmp.path().join("legacy.tsv");
    run_score(&ScoreConfig {
        null_export: export,
        dosage_files: vec![d0, d1],
        out_tsv: out.clone(),
        ac_threshold: 5,
        chunk_size: 8,
        threads: 1,
        mode: ScoreModeArg::Auto,
        variance_ratio: None,
        min_copy_var: 1e-6,
    })
    .unwrap();
    let header = std::fs::read_to_string(&out)
        .unwrap()
        .lines()
        .next()
        .unwrap()
        .to_string();
    assert!(header.contains("Chi2"));
    assert!(header.contains("Eff_anc0"));
    assert!(!header.contains("P_cct_admixed"));
}
