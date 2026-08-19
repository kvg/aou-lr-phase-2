use std::path::PathBuf;

use thiserror::Error;

#[derive(Debug, Error)]
pub enum ExtractError {
    #[error("file not found: {0}")]
    NotFound(PathBuf),
    #[error("the path '{0}' does not exist")]
    PathMissing(PathBuf),
    #[error("unexpected file extension: {0}. VCF file must end in .vcf or .vcf.gz")]
    BadExtension(String),
    #[error("the output directory '{0}' does not exist")]
    OutputDirMissing(PathBuf),
    #[error("num-ancs must be a positive integer, got {0}")]
    BadNumAncs(usize),
    #[error("{0}: expected VCF #CHROM header line")]
    MissingHeader(PathBuf),
    #[error("{path}:{line}: expected at least 10 tab-separated VCF fields")]
    TooFewFields { path: PathBuf, line: usize },
    #[error("{path}:{line}: sample {sample} (column {column}): genotype field has fewer than 4 '|'/':' tokens: {field:?}")]
    TruncatedGenotype {
        path: PathBuf,
        line: usize,
        sample: String,
        column: usize,
        field: String,
    },
    #[error("{path}: {count} analysis samples missing from VCF (e.g. {examples:?})")]
    MissingSamples {
        path: PathBuf,
        count: usize,
        examples: Vec<String>,
    },
    #[error("{path}: duplicate sample ID in VCF header: {sample}")]
    DuplicateSample { path: PathBuf, sample: String },
    #[error("{0}: no sample IDs in --samples file")]
    EmptySamples(PathBuf),
    #[error(transparent)]
    Io(#[from] std::io::Error),
}

pub type Result<T> = std::result::Result<T, ExtractError>;
