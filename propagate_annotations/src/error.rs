use std::path::PathBuf;

use thiserror::Error;

#[derive(Debug, Error)]
pub enum PropagateError {
    #[error("path not found: {0}")]
    PathMissing(PathBuf),
    #[error("invalid VCF extension (expected .vcf or .vcf.gz): {0}")]
    BadExtension(String),
    #[error("no annotation source VCF provided (need at least one of snv-indel, sv, flare)")]
    NoAnnotationSource,
    #[error("base VCF has no sample columns")]
    NoSamples,
    #[error("{path}: line {line}: {message}")]
    ParseLine {
        path: PathBuf,
        line: usize,
        message: String,
    },
    #[error("{path}: too few fields on data line {line}")]
    TooFewFields { path: PathBuf, line: usize },
    #[error("{path}: duplicate variant key at line {line}: {key}")]
    DuplicateKey {
        path: PathBuf,
        line: usize,
        key: String,
    },
    #[error("sample {sample} missing from {path}")]
    SampleMissing { sample: String, path: PathBuf },
    #[error("match rate too low for {annotation_source}: {rate:.4} < threshold {threshold:.4}")]
    MatchRateTooLow {
        annotation_source: String,
        rate: f64,
        threshold: f64,
    },
    #[error("stats error: {0}")]
    Stats(String),
    #[error("{path}: line {line}: input is not sorted by (CHROM, POS, REF, ALT): {message}")]
    Unsorted {
        path: PathBuf,
        line: usize,
        message: String,
    },
    #[error("I/O error: {0}")]
    Io(#[from] std::io::Error),
    #[error("bcftools failed: {0}")]
    Bcftools(String),
    #[error("invalid region {region}: {message}")]
    BadRegion { region: String, message: String },
}

pub type Result<T> = std::result::Result<T, PropagateError>;
