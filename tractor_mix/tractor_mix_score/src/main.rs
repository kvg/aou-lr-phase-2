use std::path::PathBuf;
use std::process::ExitCode;

use clap::Parser;
use tractor_mix_score::{run_score, ScoreConfig};

#[derive(Parser, Debug)]
#[command(
    name = "tractor-mix-score",
    version,
    about = "Streaming TractorMix.score replacement (sparse GRM / binomial null)."
)]
struct Cli {
    /// Directory written by fit_null.R (--out-null-export)
    #[arg(long)]
    null_export: PathBuf,

    /// Ancestry dosage files in order (anc0, anc1, ...)
    #[arg(long = "dosage-files", num_args = 1..)]
    dosage_files: Vec<PathBuf>,

    /// Output TSV path
    #[arg(long)]
    out: PathBuf,

    /// Minimum ancestry allele count (strict > threshold)
    #[arg(long, default_value_t = 50)]
    ac_threshold: i32,

    /// Variants per read chunk
    #[arg(long, default_value_t = 2048)]
    chunk_size: usize,

    /// Worker threads for parallel variant scoring within each chunk
    #[arg(long, default_value_t = 1)]
    threads: usize,
}

fn main() -> ExitCode {
    let cli = Cli::parse();
    if cli.ac_threshold < 0 {
        eprintln!("error: --ac-threshold must be non-negative");
        return ExitCode::FAILURE;
    }
    if cli.threads == 0 {
        eprintln!("error: --threads must be >= 1");
        return ExitCode::FAILURE;
    }
    let cfg = ScoreConfig {
        null_export: cli.null_export,
        dosage_files: cli.dosage_files,
        out_tsv: cli.out,
        ac_threshold: cli.ac_threshold,
        chunk_size: cli.chunk_size,
        threads: cli.threads,
    };
    match run_score(&cfg) {
        Ok(()) => ExitCode::SUCCESS,
        Err(e) => {
            eprintln!("error: {e}");
            ExitCode::FAILURE
        }
    }
}
