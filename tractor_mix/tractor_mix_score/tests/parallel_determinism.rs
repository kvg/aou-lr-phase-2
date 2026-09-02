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

fn write_tiny_fixture(base: &Path, n: usize, n_sites: usize) {
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
    let export = base.join("null_export");
    std::fs::create_dir_all(&export).unwrap();
    std::fs::write(
        export.join("meta.json"),
        format!(r#"{{"n":{n},"p":2,"nnz":{n},"family":"binomial"}}"#),
    )
    .unwrap();
    let ids: String = (0..n).map(|i| format!("S{i:03}\n")).collect();
    std::fs::write(export.join("id_include.txt"), ids).unwrap();
    write_sigma_i_csc(&export.join("sigma_i.csc.bin"), &sigma).unwrap();

    let mut f = File::create(export.join("sigma_i_x.bin")).unwrap();
    write_i32(&mut f, n as i32);
    write_i32(&mut f, 2);
    for _ in 0..(n * 2) {
        write_f64(&mut f, 0.01);
    }
    let mut f = File::create(export.join("cov.bin")).unwrap();
    write_i32(&mut f, 2);
    for i in 0..4 {
        write_f64(&mut f, if i % 3 == 0 { 1.0 } else { 0.0 });
    }
    let mut f = File::create(export.join("residuals.bin")).unwrap();
    write_i32(&mut f, n as i32);
    for i in 0..n {
        write_f64(&mut f, ((i as f64) * 0.07).sin() * 0.1);
    }

    let samples: Vec<String> = (0..n).map(|i| format!("S{i:03}")).collect();
    for anc in 0..2 {
        let path = base.join(format!("anc_{anc:02}.dosage.txt.gz"));
        let mut enc = GzEncoder::new(File::create(path).unwrap(), Compression::default());
        write!(
            enc,
            "CHROM\tPOS\tID\tREF\tALT\t{}\n",
            samples.join("\t")
        )
        .unwrap();
        for s in 0..n_sites {
            let mut row = vec![
                "22".into(),
                ((s + 1) * 100).to_string(),
                format!("v{s}"),
                "A".into(),
                "G".into(),
            ];
            for i in 0..n {
                row.push((((s + anc + i) % 3) as u8).to_string());
            }
            write!(enc, "{}\n", row.join("\t")).unwrap();
        }
        enc.finish().unwrap();
    }
}

#[test]
fn parallel_matches_serial_on_synthetic_fixture() {
    let tmp = TempDir::new().unwrap();
    let base = tmp.path();
    write_tiny_fixture(base, 40, 64);

    let dosages = vec![
        base.join("anc_00.dosage.txt.gz"),
        base.join("anc_01.dosage.txt.gz"),
    ];
    let out1 = base.join("t1.tsv");
    let out4 = base.join("t4.tsv");

    run_score(&ScoreConfig {
        null_export: base.join("null_export"),
        dosage_files: dosages.clone(),
        out_tsv: out1.clone(),
        ac_threshold: 5,
        chunk_size: 16,
        threads: 1,
        mode: ScoreModeArg::Legacy,
        variance_ratio: None,
        min_copy_var: 1e-6,
    })
    .unwrap();
    run_score(&ScoreConfig {
        null_export: base.join("null_export"),
        dosage_files: dosages,
        out_tsv: out4.clone(),
        ac_threshold: 5,
        chunk_size: 16,
        threads: 4,
        mode: ScoreModeArg::Legacy,
        variance_ratio: None,
        min_copy_var: 1e-6,
    })
    .unwrap();

    assert_eq!(std::fs::read(out1).unwrap(), std::fs::read(out4).unwrap());
}
