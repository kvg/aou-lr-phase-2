//! Quick local scale benchmark (no R/GMMAT). Run:
//!   cargo run --release --example micro_bench -- 9000 5000 5

use std::fs::File;
use std::io::Write;
use std::path::PathBuf;
use std::time::Instant;

use flate2::write::GzEncoder;
use flate2::Compression;
use tractor_mix_score::null::{write_sigma_i_csc, CscMatrix};
use tractor_mix_score::{run_score, ScoreConfig};

fn write_i32(w: &mut impl Write, x: i32) {
    w.write_all(&x.to_le_bytes()).unwrap();
}
fn write_f64(w: &mut impl Write, x: f64) {
    w.write_all(&x.to_le_bytes()).unwrap();
}

fn make_sparse_pd(n: usize, target_nnz: usize) -> CscMatrix {
    use std::collections::HashSet;
    let mut colptr = vec![0usize];
    let mut rowidx = Vec::new();
    let mut values = Vec::new();
    let mut seen = HashSet::new();
    // Diagonal (PD ridge)
    for i in 0..n {
        rowidx.push(i);
        values.push(1.0);
        colptr.push(rowidx.len());
    }
    seen.extend((0..n).map(|i| (i.min(i), i.max(i)))); // placeholder wrong

    let mut rng_state = 0xDEADBEEF_u64;
    while rowidx.len() < target_nnz {
        rng_state = rng_state.wrapping_mul(6364136223846793005).wrapping_add(1);
        let row = (rng_state as usize) % n;
        rng_state = rng_state.wrapping_mul(6364136223846793005).wrapping_add(1);
        let col = (rng_state as usize) % n;
        let (a, b) = if row <= col { (row, col) } else { (col, row) };
        if seen.insert((a, b)) && a != b {
            rowidx.push(row);
            values.push(0.03);
            // append to col's list — rebuild CSC properly below
        }
    }

    // Build symmetric CSC with ~target_nnz off-diagonal + diagonal
    let mut cols: Vec<Vec<(usize, f64)>> = vec![Vec::new(); n];
    for i in 0..n {
        cols[i].push((i, 1.0));
    }
    let mut rng_state = 0xCAFEBABE_u64;
    let mut off = 0usize;
    while off < target_nnz.saturating_sub(n) {
        rng_state = rng_state.wrapping_mul(6364136223846793005).wrapping_add(1);
        let row = (rng_state as usize) % n;
        rng_state = rng_state.wrapping_mul(6364136223846793005).wrapping_add(1);
        let col = (rng_state as usize) % n;
        if row == col {
            continue;
        }
        cols[col].push((row, 0.03));
        off += 1;
    }
    colptr = vec![0];
    rowidx.clear();
    values.clear();
    for col in 0..n {
        cols[col].sort_by_key(|(r, _)| *r);
        for (r, v) in &cols[col] {
            rowidx.push(*r);
            values.push(*v);
        }
        colptr.push(rowidx.len());
    }
    CscMatrix {
        n,
        colptr,
        rowidx,
        values,
    }
}

fn write_null_export(dir: &PathBuf, n: usize, p: usize, sigma: &CscMatrix) {
    std::fs::create_dir_all(dir).unwrap();
    std::fs::write(
        dir.join("meta.json"),
        format!(
            r#"{{"n":{n},"p":{p},"nnz":{},"family":"binomial","tractor_mix_score_sha":"bench"}}"#,
            sigma.values.len()
        ),
    )
    .unwrap();
    let ids: String = (0..n).map(|i| format!("S{i:05}\n")).collect();
    std::fs::write(dir.join("id_include.txt"), ids).unwrap();
    write_sigma_i_csc(&dir.join("sigma_i.csc.bin"), sigma).unwrap();
    let mut f = File::create(dir.join("sigma_i_x.bin")).unwrap();
    write_i32(&mut f, n as i32);
    write_i32(&mut f, p as i32);
    for _ in 0..(n * p) {
        write_f64(&mut f, 0.01);
    }
    let mut f = File::create(dir.join("cov.bin")).unwrap();
    write_i32(&mut f, p as i32);
    for i in 0..(p * p) {
        write_f64(&mut f, if i % (p + 1) == 0 { 1.0 } else { 0.0 });
    }
    let mut f = File::create(dir.join("residuals.bin")).unwrap();
    write_i32(&mut f, n as i32);
    for i in 0..n {
        write_f64(&mut f, ((i as f64 * 0.13).sin()) * 0.1);
    }
}

fn write_dosage(path: &PathBuf, n: usize, n_sites: usize, _n_anc: usize, anc: usize) {
    let mut enc = GzEncoder::new(File::create(path).unwrap(), Compression::default());
    let samples: Vec<String> = (0..n).map(|i| format!("S{i:05}")).collect();
    write!(
        enc,
        "CHROM\tPOS\tID\tREF\tALT\t{}\n",
        samples.join("\t")
    )
    .unwrap();
    for s in 0..n_sites {
        let mut row = vec![
            "22".to_string(),
            ((s + 1) * 100).to_string(),
            format!("v{s}"),
            "A".to_string(),
            "G".to_string(),
        ];
        for i in 0..n {
            let g = ((s + anc + i) % 3) as u8;
            row.push(g.to_string());
        }
        write!(enc, "{}\n", row.join("\t")).unwrap();
    }
    enc.finish().unwrap();
}

fn main() {
    let n: usize = std::env::args()
        .nth(1)
        .and_then(|s| s.parse().ok())
        .unwrap_or(9000);
    let n_sites: usize = std::env::args()
        .nth(2)
        .and_then(|s| s.parse().ok())
        .unwrap_or(5000);
    let n_anc: usize = std::env::args()
        .nth(3)
        .and_then(|s| s.parse().ok())
        .unwrap_or(5);

    let tmp = std::env::temp_dir().join("tractor_mix_score_bench");
    let _ = std::fs::remove_dir_all(&tmp);
    std::fs::create_dir_all(&tmp).unwrap();

    let target_nnz = std::env::args()
        .nth(4)
        .and_then(|s| s.parse().ok())
        .unwrap_or(6_500_000);
    let p = 3;
    println!("Building synthetic null (n={n}, target_nnz={target_nnz})...");
    let t0 = Instant::now();
    let sigma = make_sparse_pd(n, target_nnz);
    let nnz = sigma.values.len();
    write_null_export(&tmp.join("null_export"), n, p, &sigma);
    println!("  null export: {:.2?} (nnz={nnz})", t0.elapsed());

    let mut dosage_files = Vec::new();
    let t1 = Instant::now();
    for a in 0..n_anc {
        let pth = tmp.join(format!("anc_{a:02}.dosage.txt.gz"));
        write_dosage(&pth, n, n_sites, n_anc, a);
        dosage_files.push(pth);
    }
    println!(
        "  {n_anc} dosage files × {n_sites} sites: {:.2?}",
        t1.elapsed()
    );

    let out = tmp.join("results.tsv");
    let cfg = ScoreConfig {
        null_export: tmp.join("null_export"),
        dosage_files,
        out_tsv: out.clone(),
        ac_threshold: 50,
        chunk_size: 2048,
        threads: std::env::args()
            .nth(5)
            .and_then(|s| s.parse().ok())
            .unwrap_or(8),
    };

    let t2 = Instant::now();
    run_score(&cfg).expect("score");
    let elapsed = t2.elapsed();
    let lines = std::fs::read_to_string(&out).unwrap().lines().count() - 1;
    println!(
        "  tractor-mix-score (release): {:.2?} for {lines} variants ({:.0} var/s)",
        elapsed,
        lines as f64 / elapsed.as_secs_f64()
    );
    println!("Output: {}", out.display());
}
