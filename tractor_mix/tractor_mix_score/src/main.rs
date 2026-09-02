use std::path::PathBuf;
use std::process::ExitCode;

use clap::{Parser, ValueEnum};
use tractor_mix_score::{run_score, ScoreConfig, ScoreModeArg};

#[derive(Copy, Clone, Debug, Eq, PartialEq, ValueEnum)]
enum ModeArg {
    Auto,
    Legacy,
    Felix,
}

impl From<ModeArg> for ScoreModeArg {
    fn from(m: ModeArg) -> Self {
        match m {
            ModeArg::Auto => ScoreModeArg::Auto,
            ModeArg::Legacy => ScoreModeArg::Legacy,
            ModeArg::Felix => ScoreModeArg::Felix,
        }
    }
}

#[derive(Parser, Debug)]
#[command(
    name = "tractor-mix-score",
    version,
    about = "Streaming Tractor-Mix / FELIX score tests (sparse GRM null)."
)]
struct Cli {
    /// Directory written by fit_null.R or export_felix_null.R (--out-null-export)
    #[arg(long)]
    null_export: PathBuf,

    /// Ancestry dosage files in order (anc0, anc1, ...). Integers or f64.
    #[arg(long = "dosage-files", num_args = 1..)]
    dosage_files: Vec<PathBuf>,

    /// Output TSV path
    #[arg(long)]
    out: PathBuf,

    /// Minimum ancestry allele count for legacy GMMAT mode (strict > threshold).
    /// Ignored in --mode felix (zero-variance ancestries are dropped instead).
    #[arg(long, default_value_t = 50)]
    ac_threshold: i32,

    /// Scoring mode: auto detects FELIX vs GMMAT from meta.json
    #[arg(long, value_enum, default_value_t = ModeArg::Auto)]
    mode: ModeArg,

    /// Override the scalar variance ratio from the null export (FELIX mode)
    #[arg(long)]
    variance_ratio: Option<f64>,

    /// Skip a variant in FELIX mode when total copy-number variance is below this
    #[arg(long, default_value_t = 1e-6)]
    min_copy_var: f64,

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
    if let Some(vr) = cli.variance_ratio {
        if !(vr.is_finite() && vr > 0.0) {
            eprintln!("error: --variance-ratio must be positive");
            return ExitCode::FAILURE;
        }
    }
    let cfg = ScoreConfig {
        null_export: cli.null_export,
        dosage_files: cli.dosage_files,
        out_tsv: cli.out,
        ac_threshold: cli.ac_threshold,
        chunk_size: cli.chunk_size,
        threads: cli.threads,
        mode: cli.mode.into(),
        variance_ratio: cli.variance_ratio,
        min_copy_var: cli.min_copy_var,
    };
    match run_score(&cfg) {
        Ok(()) => ExitCode::SUCCESS,
        Err(e) => {
            eprintln!("error: {e}");
            ExitCode::FAILURE
        }
    }
}
