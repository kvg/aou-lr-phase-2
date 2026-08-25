pub mod error;
pub mod fields;
pub mod index;
pub mod merge;
pub mod region;
pub mod stats;
pub mod vcf_io;

pub use error::{PropagateError, Result};
pub use merge::{propagate_annotations, PropagateConfig};
pub use stats::PropagateStats;
