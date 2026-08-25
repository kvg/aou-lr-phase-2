use crate::error::{PropagateError, Result};

/// 1-based inclusive genomic interval (VCF / bcftools style).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Region {
    pub chrom: String,
    pub start: u64,
    pub end: u64,
}

impl Region {
    pub fn parse(spec: &str) -> Result<Self> {
        let (chrom, rest) = match spec.split_once(':') {
            Some(pair) => pair,
            None => {
                return Ok(Self {
                    chrom: spec.to_string(),
                    start: 1,
                    end: u64::MAX,
                });
            }
        };
        if chrom.is_empty() {
            return Err(PropagateError::BadRegion {
                region: spec.to_string(),
                message: "empty chromosome".into(),
            });
        }
        let (start_s, end_s) = rest.split_once('-').ok_or_else(|| PropagateError::BadRegion {
            region: spec.to_string(),
            message: "expected CHR:START-END".into(),
        })?;
        let start: u64 = start_s.parse().map_err(|_| PropagateError::BadRegion {
            region: spec.to_string(),
            message: format!("bad start: {start_s}"),
        })?;
        let end: u64 = if end_s.is_empty() {
            u64::MAX
        } else {
            end_s.parse().map_err(|_| PropagateError::BadRegion {
                region: spec.to_string(),
                message: format!("bad end: {end_s}"),
            })?
        };
        if start == 0 || end < start {
            return Err(PropagateError::BadRegion {
                region: spec.to_string(),
                message: "start must be >= 1 and end >= start".into(),
            });
        }
        Ok(Self {
            chrom: chrom.to_string(),
            start,
            end,
        })
    }

    pub fn contains(&self, chrom: &str, pos: u64) -> bool {
        chrom == self.chrom && pos >= self.start && pos <= self.end
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_full_chrom_and_interval() {
        let r = Region::parse("chr22").unwrap();
        assert!(r.contains("chr22", 1));
        assert!(r.contains("chr22", u64::MAX));
        assert!(!r.contains("chr21", 1));

        let r = Region::parse("chr22:100-200").unwrap();
        assert!(r.contains("chr22", 100));
        assert!(r.contains("chr22", 200));
        assert!(!r.contains("chr22", 99));
        assert!(!r.contains("chr22", 201));
    }
}
