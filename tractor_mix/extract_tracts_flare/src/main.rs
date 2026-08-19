use std::path::PathBuf;
use std::process::ExitCode;

use clap::Parser;
use extract_tracts_flare::{extract_tracts_flare, ExtractConfig};

#[derive(Parser, Debug)]
#[command(
    name = "extract-tracts-flare",
    version,
    about = "Extract ancestry-specific dosage and hapcount files from a FLARE VCF (Tractor extract_tracts_flare.py compatible)."
)]
struct Cli {
    /// Path to VCF file (*.vcf or *.vcf.gz)
    #[arg(long)]
    vcf: PathBuf,

    /// Number of ancestral populations within the VCF file
    #[arg(long)]
    num_ancs: usize,

    /// Directory for output files. Directory must already exist.
    #[arg(long)]
    output_dir: Option<PathBuf>,

    /// Write ancestry-specific VCF files
    #[arg(long, default_value_t = false)]
    output_vcf: bool,

    /// gzip (not bgzip) all output files
    #[arg(long, default_value_t = false)]
    compress_output: bool,

    /// Subset and reorder sample columns to this ID list (one ID per line)
    #[arg(long, alias = "keep")]
    samples: Option<PathBuf>,

    /// Reserved for future parallel parsing; extract is currently single-threaded
    #[arg(long, default_value_t = 1)]
    threads: usize,
}

fn main() -> ExitCode {
    let cli = Cli::parse();
    let cfg = ExtractConfig {
        vcf: cli.vcf,
        num_ancs: cli.num_ancs,
        output_dir: cli.output_dir,
        output_vcf: cli.output_vcf,
        compress_output: cli.compress_output,
        samples: cli.samples,
        threads: cli.threads,
    };
    match extract_tracts_flare(&cfg) {
        Ok(()) => ExitCode::SUCCESS,
        Err(e) => {
            eprintln!("error: {e}");
            ExitCode::FAILURE
        }
    }
}
