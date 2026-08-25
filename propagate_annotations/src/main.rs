use std::path::PathBuf;
use std::process::ExitCode;

use clap::Parser;
use propagate_annotations::region::Region;
use propagate_annotations::{propagate_annotations, PropagateConfig};

#[derive(Parser, Debug)]
#[command(
    name = "propagate-annotations",
    version,
    about = "Union SNV/indel, SV, and FLARE VCFs into one fully annotated callset."
)]
struct Cli {
    /// SNV/indel joint callset (.vcf or .vcf.gz)
    #[arg(long)]
    snv_indel_vcf: Option<PathBuf>,

    /// SV integrated partition (.vcf or .vcf.gz; convert BCF with bcftools view first)
    #[arg(long)]
    sv_vcf: Option<PathBuf>,

    /// FLARE local-ancestry VCF (.vcf or .vcf.gz)
    #[arg(long)]
    flare_vcf: Option<PathBuf>,

    /// Output path (.bcf preferred; also .vcf / .vcf.gz)
    #[arg(long)]
    output: PathBuf,

    /// JSON sanity-check / match-rate statistics
    #[arg(long)]
    stats_json: Option<PathBuf>,

    /// TSV summary of match statistics
    #[arg(long)]
    stats_tsv: Option<PathBuf>,

    /// When SNV/indel and FLARE are both provided, fail if FLARE covers fewer than
    /// this fraction of SNV/indel sites (0–1)
    #[arg(long)]
    min_match_rate_flare: Option<f64>,

    /// Deprecated no-op (kept for WDL compatibility)
    #[arg(long, hide = true)]
    min_match_rate_snv_indel: Option<f64>,

    /// Deprecated no-op (kept for WDL compatibility)
    #[arg(long, hide = true)]
    min_match_rate_sv: Option<f64>,

    /// Write gzip-compressed VCF when the path ends with .gz (ignored for .bcf)
    #[arg(long)]
    compress_output: Option<bool>,

    /// Rayon threads for per-site sample assembly (0 = rayon default / all CPUs)
    #[arg(long, default_value_t = 0)]
    threads: usize,

    /// Emit every FORMAT tag on every site with '.' fillers (larger / slower)
    #[arg(long, default_value_t = false)]
    dense_format: bool,

    /// Optional region filter CHR or CHR:START-END
    #[arg(long)]
    region: Option<String>,

    /// Progress log every N union sites
    #[arg(long, default_value_t = 100_000)]
    progress_every: u64,
}

fn main() -> ExitCode {
    let cli = Cli::parse();
    let compress = cli.compress_output.unwrap_or_else(|| {
        let s = cli.output.to_string_lossy();
        (s.ends_with(".gz") || s.ends_with(".bgz")) && !s.ends_with(".bcf.gz")
    });
    let region = match cli.region.as_deref().map(Region::parse).transpose() {
        Ok(r) => r,
        Err(e) => {
            eprintln!("error: {e}");
            return ExitCode::FAILURE;
        }
    };
    let cfg = PropagateConfig {
        snv_indel_vcf: cli.snv_indel_vcf,
        sv_vcf: cli.sv_vcf,
        flare_vcf: cli.flare_vcf,
        output: cli.output,
        compress_output: compress,
        stats_json: cli.stats_json,
        stats_tsv: cli.stats_tsv,
        min_match_rate_snv_indel: cli.min_match_rate_snv_indel,
        min_match_rate_sv: cli.min_match_rate_sv,
        min_match_rate_flare: cli.min_match_rate_flare,
        threads: cli.threads,
        dense_format: cli.dense_format,
        region,
        progress_every: cli.progress_every,
    };
    match propagate_annotations(&cfg) {
        Ok(stats) => {
            eprintln!(
                "OK: {} union sites -> {} output sites on {} ({} samples)",
                stats.union_sites,
                stats.output_sites,
                stats.chromosome,
                stats.n_samples
            );
            if let Some(rate) = stats.snv_flare_overlap_rate {
                eprintln!("  flare covers {:.2}% of SNV/indel sites", rate * 100.0);
            }
            for src in &stats.sources {
                eprintln!(
                    "  {}: annotated {}/{} union sites ({:.2}%)",
                    src.label,
                    src.union_sites_with_source,
                    stats.union_sites,
                    src.union_annotation_rate * 100.0,
                );
            }
            ExitCode::SUCCESS
        }
        Err(e) => {
            eprintln!("error: {e}");
            ExitCode::FAILURE
        }
    }
}
