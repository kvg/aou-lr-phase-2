use std::collections::BTreeSet;

use serde::Serialize;

use crate::fields::SourceSpec;

#[derive(Debug, Clone, Serialize)]
pub struct PropagateStats {
    pub mode: String,
    pub output_vcf: String,
    pub chromosome: String,
    pub n_samples: usize,
    pub union_sites: u64,
    pub output_sites: u64,
    pub snv_sites: u64,
    pub flare_on_snv: u64,
    pub sites_by_n_sources: std::collections::BTreeMap<u8, u64>,
    pub snv_flare_overlap_rate: Option<f64>,
    pub sources: Vec<SourceStats>,
}

#[derive(Debug, Clone, Serialize)]
pub struct SourceStats {
    pub label: String,
    pub path: String,
    pub source_sites: u64,
    pub union_sites_with_source: u64,
    pub union_sites_without_source: u64,
    pub union_annotation_rate: f64,
    pub sample_gt_fill_rate: f64,
    pub sample_gt_filled: u64,
    pub sample_gt_total: u64,
}

#[derive(Debug, Default)]
pub struct StatsCollector {
    pub output_vcf: String,
    pub n_samples: usize,
    pub chromosomes: BTreeSet<String>,
    pub union_sites: u64,
    pub output_sites: u64,
    pub sites_by_n_sources: std::collections::BTreeMap<u8, u64>,
    pub sources: Vec<SourceStatsCollector>,
    pub snv_sites: u64,
    pub flare_on_snv: u64,
}

#[derive(Debug, Default)]
pub struct SourceStatsCollector {
    pub label: String,
    pub path: String,
    pub union_sites_with_source: u64,
    pub sample_gt_filled: u64,
    pub sample_gt_total: u64,
}

impl StatsCollector {
    pub fn new(
        output_vcf: &str,
        n_samples: usize,
        specs: &[(String, SourceSpec, String)],
    ) -> Self {
        Self {
            output_vcf: output_vcf.to_string(),
            n_samples,
            sources: specs
                .iter()
                .map(|(label, _spec, path)| SourceStatsCollector {
                    label: label.to_string(),
                    path: path.to_string(),
                    ..Default::default()
                })
                .collect(),
            ..Default::default()
        }
    }

    pub fn record_site(
        &mut self,
        chrom: &str,
        n_sources_matched: u8,
        source_matches: &[(usize, bool, u64, u64)],
        has_snv: bool,
        has_flare: bool,
    ) {
        self.union_sites += 1;
        self.output_sites += 1;
        self.chromosomes.insert(chrom.to_string());
        *self
            .sites_by_n_sources
            .entry(n_sources_matched)
            .or_default() += 1;
        if has_snv {
            self.snv_sites += 1;
            if has_flare {
                self.flare_on_snv += 1;
            }
        }
        for (idx, matched, gt_filled, gt_total) in source_matches {
            if *matched {
                self.sources[*idx].union_sites_with_source += 1;
                self.sources[*idx].sample_gt_filled += gt_filled;
                self.sources[*idx].sample_gt_total += gt_total;
            }
        }
    }

    pub fn finalize(self, source_site_counts: &[u64]) -> PropagateStats {
        let chromosome = if self.chromosomes.len() == 1 {
            self.chromosomes.iter().next().cloned().unwrap_or_default()
        } else if self.chromosomes.is_empty() {
            String::new()
        } else {
            format!("multi[{}]", self.chromosomes.len())
        };

        let union_sites = self.union_sites;
        let sources = self
            .sources
            .into_iter()
            .zip(source_site_counts.iter().copied())
            .map(|(s, source_sites)| {
                let matched = s.union_sites_with_source;
                SourceStats {
                    label: s.label,
                    path: s.path,
                    source_sites,
                    union_sites_with_source: matched,
                    union_sites_without_source: union_sites.saturating_sub(matched),
                    union_annotation_rate: rate(matched, union_sites),
                    sample_gt_fill_rate: rate(s.sample_gt_filled, s.sample_gt_total),
                    sample_gt_filled: s.sample_gt_filled,
                    sample_gt_total: s.sample_gt_total,
                }
            })
            .collect();

        let snv_flare_overlap_rate = if self.snv_sites == 0 {
            None
        } else {
            Some(rate(self.flare_on_snv, self.snv_sites))
        };

        PropagateStats {
            mode: "union_stream".to_string(),
            output_vcf: self.output_vcf,
            chromosome,
            n_samples: self.n_samples,
            union_sites,
            output_sites: self.output_sites,
            snv_sites: self.snv_sites,
            flare_on_snv: self.flare_on_snv,
            sites_by_n_sources: self.sites_by_n_sources,
            snv_flare_overlap_rate,
            sources,
        }
    }
}

impl PropagateStats {
    pub fn write_json(&self, path: &std::path::Path) -> crate::error::Result<()> {
        let text = serde_json::to_string_pretty(self).map_err(|e| {
            crate::error::PropagateError::Stats(format!("serialize json: {e}"))
        })?;
        std::fs::write(path, text)?;
        Ok(())
    }

    pub fn write_tsv(&self, path: &std::path::Path) -> crate::error::Result<()> {
        let mut lines = vec![
            "metric\tvalue".to_string(),
            format!("mode\t{}", self.mode),
            format!("chromosome\t{}", self.chromosome),
            format!("n_samples\t{}", self.n_samples),
            format!("union_sites\t{}", self.union_sites),
            format!("output_sites\t{}", self.output_sites),
        ];
        if let Some(rate) = self.snv_flare_overlap_rate {
            lines.push(format!("snv_flare_overlap_rate\t{rate:.6}"));
        }
        for (n, count) in &self.sites_by_n_sources {
            lines.push(format!("union_sites_with_{n}_sources\t{count}"));
        }
        for src in &self.sources {
            let p = &src.label;
            lines.push(format!("{p}_source_sites\t{}", src.source_sites));
            lines.push(format!(
                "{p}_union_sites_with_source\t{}",
                src.union_sites_with_source
            ));
            lines.push(format!(
                "{p}_union_sites_without_source\t{}",
                src.union_sites_without_source
            ));
            lines.push(format!(
                "{p}_union_annotation_rate\t{:.6}",
                src.union_annotation_rate
            ));
            lines.push(format!(
                "{p}_sample_gt_fill_rate\t{:.6}",
                src.sample_gt_fill_rate
            ));
        }
        std::fs::write(path, lines.join("\n") + "\n")?;
        Ok(())
    }

    pub fn check_min_match_rates(
        &self,
        _snv_min: Option<f64>,
        _sv_min: Option<f64>,
        flare_min: Option<f64>,
    ) -> crate::error::Result<()> {
        if self.union_sites != self.output_sites {
            return Err(crate::error::PropagateError::Stats(format!(
                "internal error: union_sites ({}) != output_sites ({})",
                self.union_sites, self.output_sites
            )));
        }

        if let Some(threshold) = flare_min {
            if let Some(rate) = self.snv_flare_overlap_rate {
                if rate < threshold {
                    return Err(crate::error::PropagateError::MatchRateTooLow {
                        annotation_source: "flare_on_snv".to_string(),
                        rate,
                        threshold,
                    });
                }
            }
        }

        Ok(())
    }
}

fn rate(num: u64, den: u64) -> f64 {
    if den == 0 {
        0.0
    } else {
        num as f64 / den as f64
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rate_handles_zero() {
        assert_eq!(rate(0, 0), 0.0);
        assert_eq!(rate(5, 10), 0.5);
    }
}
