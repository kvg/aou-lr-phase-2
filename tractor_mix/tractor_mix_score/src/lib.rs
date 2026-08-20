pub mod dosage;
pub mod error;
pub mod format;
pub mod null;
pub mod output;
pub mod score;

use std::io::Write;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::time::Instant;

use rayon::prelude::*;

use dosage::{count_variants, read_dosage_header, DosageReader, VariantMeta};
use error::{Result, ScoreError};
use null::load_null_export;
use output::{open_output, write_header, write_variant_row};
use score::{score_variant, VariantStats};

pub struct ScoreConfig {
    pub null_export: std::path::PathBuf,
    pub dosage_files: Vec<std::path::PathBuf>,
    pub out_tsv: std::path::PathBuf,
    pub ac_threshold: i32,
    pub chunk_size: usize,
    pub threads: usize,
}

struct Progress {
    total: usize,
    done: AtomicUsize,
    start: Instant,
    last_log: std::sync::Mutex<Instant>,
}

impl Progress {
    fn new(total: usize) -> Self {
        let now = Instant::now();
        Self {
            total,
            done: AtomicUsize::new(0),
            start: now,
            last_log: std::sync::Mutex::new(now),
        }
    }

    fn add(&self, n: usize) {
        let processed = self.done.fetch_add(n, Ordering::Relaxed) + n;
        let mut last = self.last_log.lock().expect("progress lock");
        if last.elapsed().as_secs() >= 30 || processed >= self.total {
            let elapsed = self.start.elapsed().as_secs_f64().max(1e-9);
            let rate = processed as f64 / elapsed;
            let remaining = self.total.saturating_sub(processed);
            let eta = remaining as f64 / rate.max(1e-9);
            eprintln!(
                "tractor-mix-score: {processed}/{} variants ({rate:.0} var/s, ETA {eta:.0}s)",
                self.total
            );
            *last = Instant::now();
        }
    }
}

fn score_chunk_variants(
    null: &null::NullModel,
    chunk_geno: &[Vec<u8>],
    n: usize,
    n_anc: usize,
    n_chunk: usize,
    ac_threshold: i32,
    threads: usize,
) -> Vec<VariantStats> {
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(threads.max(1))
        .build()
        .expect("thread pool");

    pool.install(|| {
        (0..n_chunk)
            .into_par_iter()
            .map(|vi| {
                let mut geno_all = vec![0u8; n * n_anc];
                for anc in 0..n_anc {
                    let src = &chunk_geno[anc][vi * n..(vi + 1) * n];
                    geno_all[anc * n..(anc + 1) * n].copy_from_slice(src);
                }
                score_variant(null, &geno_all, n_anc, ac_threshold)
            })
            .collect()
    })
}

pub fn run_score(cfg: &ScoreConfig) -> Result<()> {
    if cfg.threads == 0 {
        return Err(ScoreError::msg("--threads must be >= 1"));
    }

    let null = load_null_export(&cfg.null_export)?;
    let n_anc = cfg.dosage_files.len();
    if n_anc == 0 {
        return Err(ScoreError::msg("at least one --dosage-file required"));
    }

    let header0 = read_dosage_header(&cfg.dosage_files[0])?;
    for path in cfg.dosage_files.iter().skip(1) {
        let h = read_dosage_header(path)?;
        if h != header0 {
            return Err(ScoreError::msg(format!(
                "sample order mismatch in {}",
                path.display()
            )));
        }
    }

    let n_snps = count_variants(&cfg.dosage_files[0])?;
    eprintln!(
        "tractor-mix-score: n={} variants={} ancestries={} threads={} chunk_size={}",
        null.meta.n, n_snps, n_anc, cfg.threads, cfg.chunk_size
    );

    let mut readers = Vec::with_capacity(n_anc);
    for path in &cfg.dosage_files {
        let mut r = DosageReader::open(path, &null.id_include, &header0)?;
        r.skip_header()?;
        readers.push(r);
    }

    let mut out = open_output(&cfg.out_tsv)?;
    write_header(&mut out, n_anc)?;

    let n = null.meta.n;
    let mut chunk_meta: Vec<VariantMeta> = Vec::new();
    let mut chunk_geno: Vec<Vec<u8>> = vec![Vec::new(); n_anc];
    let mut processed = 0usize;
    let progress = Progress::new(n_snps);
    let run_start = Instant::now();

    while processed < n_snps {
        chunk_meta.clear();
        let mut n_chunk = 0usize;
        for (anc, reader) in readers.iter_mut().enumerate() {
            let nread = reader.read_chunk(cfg.chunk_size, &mut chunk_geno[anc], &mut chunk_meta)?;
            if anc == 0 {
                n_chunk = nread;
            } else if nread != n_chunk {
                return Err(ScoreError::msg("dosage files differ in length"));
            }
        }
        if n_chunk == 0 {
            break;
        }

        let stats = score_chunk_variants(
            &null,
            &chunk_geno,
            n,
            n_anc,
            n_chunk,
            cfg.ac_threshold,
            cfg.threads,
        );

        for (vi, st) in stats.into_iter().enumerate() {
            write_variant_row(&mut out, &chunk_meta[vi], &st, n_anc)?;
        }

        processed += n_chunk;
        progress.add(n_chunk);
    }

    out.flush()?;
    let elapsed = run_start.elapsed().as_secs_f64();
    eprintln!(
        "tractor-mix-score: done {processed} variants in {elapsed:.1}s ({:.0} var/s)",
        processed as f64 / elapsed.max(1e-9)
    );
    Ok(())
}
