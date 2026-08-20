use std::path::PathBuf;
use std::process::Command;

use tractor_mix_score::null::{csc_matmul_cols, CscMatrix};

fn crate_root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
}

fn gmmat_available() -> bool {
    Command::new("Rscript")
        .arg("-e")
        .arg("suppressPackageStartupMessages(library(GMMAT)); quit(status=0)")
        .status()
        .map(|s| s.success())
        .unwrap_or(false)
}

#[test]
fn csc_matmul_matches_naive_dense() {
    let n = 20;
    let mut dense = vec![0.0; n * n];
    for i in 0..n {
        dense[i + i * n] = 1.0 + 0.01 * i as f64;
        if i + 3 < n {
            dense[i + (i + 3) * n] = 0.05;
            dense[(i + 3) + i * n] = 0.05;
        }
    }
    let mut colptr = vec![0usize];
    let mut rowidx = Vec::new();
    let mut values = Vec::new();
    for col in 0..n {
        for row in 0..n {
            let v = dense[row + col * n];
            if v != 0.0 {
                rowidx.push(row);
                values.push(v);
            }
        }
        colptr.push(rowidx.len());
    }
    let csc = CscMatrix {
        n,
        colptr,
        rowidx,
        values,
    };

    let k = 3;
    let g: Vec<f64> = (0..n * k).map(|i| (i as f64 * 0.13).sin()).collect();
    let mut y_csc = vec![0.0; n * k];
    csc_matmul_cols(&csc, &g, k, &mut y_csc);

    let mut y_dense = vec![0.0; n * k];
    for j in 0..k {
        for i in 0..n {
            let mut s = 0.0;
            for t in 0..n {
                s += dense[i + t * n] * g[t + j * n];
            }
            y_dense[i + j * n] = s;
        }
    }
    for (a, b) in y_csc.iter().zip(y_dense.iter()) {
        assert!((a - b).abs() < 1e-10, "{a} vs {b}");
    }
}

#[test]
#[ignore = "requires Rscript + GMMAT"]
fn generate_fixture_smoke() {
    if !gmmat_available() {
        eprintln!("skip: Rscript/GMMAT not available");
        return;
    }
    let fixture = crate_root().join("testdata/fixture");
    let status = Command::new("Rscript")
        .arg(crate_root().join("testdata/gen_oracle.R"))
        .arg(&fixture)
        .status()
        .unwrap();
    assert!(status.success());
    assert!(fixture.join("null_export/meta.json").exists());
}
