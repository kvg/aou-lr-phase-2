use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};
use std::process::{Child, ChildStdin, Command, Stdio};
use std::time::Instant;

use rayon::prelude::*;

use crate::error::{PropagateError, Result};
use crate::fields::{
    all_output_format_ids, format_description, format_ids_for_labels, format_number, format_type,
    FLARE, SNV_INDEL, SV,
};
use crate::index::{union_sample_lists, AnnotRecord, AnnotationStream};
use crate::region::Region;
use crate::stats::{PropagateStats, StatsCollector};
use crate::vcf_io::{bytes_to_string, merge_info, parse_info, VcfWriter};

#[derive(Debug, Clone)]
pub struct PropagateConfig {
    pub snv_indel_vcf: Option<PathBuf>,
    pub sv_vcf: Option<PathBuf>,
    pub flare_vcf: Option<PathBuf>,
    pub output: PathBuf,
    pub compress_output: bool,
    pub stats_json: Option<PathBuf>,
    pub stats_tsv: Option<PathBuf>,
    pub min_match_rate_snv_indel: Option<f64>,
    pub min_match_rate_sv: Option<f64>,
    pub min_match_rate_flare: Option<f64>,
    /// Rayon worker threads for per-site sample assembly (0 = rayon default).
    pub threads: usize,
    /// Emit every declared FORMAT tag with `.` fillers (legacy / denser output).
    pub dense_format: bool,
    /// Optional CHR or CHR:START-END filter (inputs should already be subset when sharding).
    pub region: Option<Region>,
    /// Progress log every N sites.
    pub progress_every: u64,
}

impl Default for PropagateConfig {
    fn default() -> Self {
        Self {
            snv_indel_vcf: None,
            sv_vcf: None,
            flare_vcf: None,
            output: PathBuf::from("annotated.bcf"),
            compress_output: false,
            stats_json: None,
            stats_tsv: None,
            min_match_rate_snv_indel: None,
            min_match_rate_sv: None,
            min_match_rate_flare: None,
            threads: 0,
            dense_format: false,
            region: None,
            progress_every: 100_000,
        }
    }
}

struct SourceInput {
    path: PathBuf,
    spec: crate::fields::SourceSpec,
}

enum OutputSink {
    Vcf(VcfWriter),
    Bcf {
        child: Child,
        stdin: BufWriter<ChildStdin>,
    },
}

impl OutputSink {
    fn open(path: &Path, compress_vcf_gz: bool) -> Result<Self> {
        let name = path
            .file_name()
            .and_then(|s| s.to_str())
            .unwrap_or_default();
        if name.ends_with(".bcf") || name.ends_with(".bcf.gz") {
            let mut cmd = Command::new("bcftools");
            cmd.arg("view").arg("-Ob").arg("-o").arg(path).arg("-");
            let mut child = cmd
                .stdin(Stdio::piped())
                .stdout(Stdio::null())
                .stderr(Stdio::piped())
                .spawn()
                .map_err(|e| {
                    PropagateError::Bcftools(format!(
                        "failed to spawn bcftools view (is bcftools on PATH?): {e}"
                    ))
                })?;
            let stdin = child.stdin.take().ok_or_else(|| {
                PropagateError::Bcftools("bcftools stdin unavailable".into())
            })?;
            Ok(Self::Bcf {
                child,
                stdin: BufWriter::with_capacity(1 << 20, stdin),
            })
        } else {
            Ok(Self::Vcf(VcfWriter::create(path, compress_vcf_gz)?))
        }
    }

    fn write_line(&mut self, line: &[u8]) -> Result<()> {
        match self {
            Self::Vcf(w) => w.write_line(line),
            Self::Bcf { stdin, .. } => {
                stdin.write_all(line)?;
                if !line.ends_with(b"\n") {
                    stdin.write_all(b"\n")?;
                }
                Ok(())
            }
        }
    }

    fn finish(self) -> Result<()> {
        match self {
            Self::Vcf(w) => w.finish(),
            Self::Bcf { child, stdin } => {
                stdin.into_inner().map_err(|e| PropagateError::Io(e.into_error()))?;
                let output = child.wait_with_output()?;
                if !output.status.success() {
                    let err = String::from_utf8_lossy(&output.stderr);
                    return Err(PropagateError::Bcftools(format!(
                        "bcftools view -Ob failed ({}): {err}",
                        output.status
                    )));
                }
                Ok(())
            }
        }
    }
}

/// Sorted streaming union of SNV/indel + SV + FLARE.
///
/// Inputs must be sorted by (CHROM, POS, REF, ALT). Memory stays O(samples ×
/// FORMAT fields × open streams), not O(sites).
pub fn propagate_annotations(cfg: &PropagateConfig) -> Result<PropagateStats> {
    if cfg.threads > 0 {
        rayon::ThreadPoolBuilder::new()
            .num_threads(cfg.threads)
            .build_global()
            .ok(); // ignore if pool already built (tests / re-entry)
    }

    let inputs = collect_inputs(cfg)?;
    if inputs.is_empty() {
        return Err(PropagateError::NoAnnotationSource);
    }

    let sample_lists: Vec<Vec<String>> = inputs
        .iter()
        .map(|i| AnnotationStream::read_sample_ids(&i.path))
        .collect::<Result<_>>()?;
    let output_samples = union_sample_lists(&sample_lists)?;

    let mut streams: Vec<AnnotationStream> = inputs
        .iter()
        .map(|i| AnnotationStream::open(i.path.clone(), i.spec, &output_samples))
        .collect::<Result<_>>()?;

    let source_specs: Vec<(String, crate::fields::SourceSpec, String)> = inputs
        .iter()
        .map(|i| {
            (
                i.spec.label.to_string(),
                i.spec,
                i.path.display().to_string(),
            )
        })
        .collect();

    let mut stats = StatsCollector::new(
        &cfg.output.display().to_string(),
        output_samples.len(),
        &source_specs,
    );

    let header_lines = choose_header_lines(&streams);
    let mut writer = OutputSink::open(&cfg.output, cfg.compress_output)?;
    write_output_header(&mut writer, &header_lines, &output_samples)?;

    let n_samples = output_samples.len();
    let progress_every = cfg.progress_every.max(1);
    let t0 = Instant::now();
    let mut sites_done: u64 = 0;
    let mut line_buf = String::with_capacity(1 << 20);

    loop {
        let min_key = streams
            .iter()
            .filter_map(|s| s.peek_key().map(|(c, k)| (c.to_string(), k.clone())))
            .min();
        let Some((chrom, key)) = min_key else {
            break;
        };

        let mut matching: Vec<Option<AnnotRecord>> = Vec::with_capacity(streams.len());
        for stream in streams.iter_mut() {
            let take = matches!(
                stream.peek_key(),
                Some((c, k)) if c == chrom && *k == key
            );
            if take {
                matching.push(stream.current.take());
                stream.advance()?;
            } else {
                matching.push(None);
            }
        }

        if let Some(region) = &cfg.region {
            // Inputs are sorted; once we're past the window, stop.
            if chrom != region.chrom || key.pos > region.end {
                break;
            }
            if key.pos < region.start {
                continue;
            }
        }

        let shell = matching
            .iter()
            .flatten()
            .next()
            .expect("at least one matching record");
        let id = if shell.id.is_empty() { "." } else { shell.id.as_str() };
        let qual = if shell.qual.is_empty() {
            "."
        } else {
            shell.qual.as_str()
        };
        let filter = if shell.filter.is_empty() {
            "."
        } else {
            shell.filter.as_str()
        };
        let (merged_info, src_labels) = merge_info_from_records(&matching, &streams);
        let (out_format, merged_samples, match_meta) =
            merge_sample_fields(n_samples, &matching, &streams, cfg.dense_format);

        let n_matched = src_labels.len() as u8;
        let has_snv = src_labels.iter().any(|l| *l == SNV_INDEL.label);
        let has_flare = src_labels.iter().any(|l| *l == FLARE.label);
        let source_matches: Vec<(usize, bool, u64, u64)> = match_meta
            .iter()
            .enumerate()
            .map(|(i, m)| (i, m.matched, m.gt_filled, m.gt_total))
            .collect();
        stats.record_site(&chrom, n_matched, &source_matches, has_snv, has_flare);

        line_buf.clear();
        use std::fmt::Write as _;
        let _ = write!(
            line_buf,
            "{chrom}\t{pos}\t{id}\t{ref_}\t{alt}\t{qual}\t{filter}\t{info}\t{fmt}",
            pos = key.pos,
            ref_ = key.ref_allele,
            alt = key.alt_allele,
            info = merged_info,
            fmt = out_format,
        );
        for sample in &merged_samples {
            line_buf.push('\t');
            line_buf.push_str(sample);
        }
        line_buf.push('\n');
        writer.write_line(line_buf.as_bytes())?;

        sites_done += 1;
        if sites_done % progress_every == 0 {
            let secs = t0.elapsed().as_secs_f64().max(1e-6);
            eprintln!(
                "propagate-annotations: {} union sites ({:.1} sites/s, {:.1} min, {})",
                sites_done,
                sites_done as f64 / secs,
                secs / 60.0,
                chrom
            );
        }
    }

    writer.finish()?;

    let source_site_counts: Vec<u64> = streams.iter().map(|s| s.sites_read).collect();
    let stats = stats.finalize(&source_site_counts);

    stats.check_min_match_rates(
        cfg.min_match_rate_snv_indel,
        cfg.min_match_rate_sv,
        cfg.min_match_rate_flare,
    )?;

    if let Some(path) = &cfg.stats_json {
        stats.write_json(path)?;
    }
    if let Some(path) = &cfg.stats_tsv {
        stats.write_tsv(path)?;
    }

    let secs = t0.elapsed().as_secs_f64().max(1e-6);
    eprintln!(
        "propagate-annotations: finished {} sites in {:.1} min ({:.1} sites/s)",
        sites_done,
        secs / 60.0,
        sites_done as f64 / secs
    );

    Ok(stats)
}

fn collect_inputs(cfg: &PropagateConfig) -> Result<Vec<SourceInput>> {
    let mut out = Vec::new();
    if let Some(p) = &cfg.snv_indel_vcf {
        out.push(SourceInput {
            path: p.clone(),
            spec: SNV_INDEL,
        });
    }
    if let Some(p) = &cfg.sv_vcf {
        out.push(SourceInput {
            path: p.clone(),
            spec: SV,
        });
    }
    if let Some(p) = &cfg.flare_vcf {
        out.push(SourceInput {
            path: p.clone(),
            spec: FLARE,
        });
    }
    Ok(out)
}

fn choose_header_lines(streams: &[AnnotationStream]) -> Vec<Vec<u8>> {
    for label in [SNV_INDEL.label, SV.label, FLARE.label] {
        if let Some(src) = streams.iter().find(|s| s.spec.label == label) {
            if !src.header_lines.is_empty() {
                return src.header_lines.clone();
            }
        }
    }
    Vec::new()
}

fn write_output_header(
    writer: &mut OutputSink,
    header_lines: &[Vec<u8>],
    samples: &[String],
) -> Result<()> {
    use std::collections::HashSet;
    let existing_fmt: HashSet<String> = header_lines
        .iter()
        .filter_map(|line| {
            let s = bytes_to_string(line);
            if let Some(rest) = s.strip_prefix("##FORMAT=<ID=") {
                rest.split(',').next().map(|id| id.to_string())
            } else {
                None
            }
        })
        .collect();
    let existing_info: HashSet<String> = header_lines
        .iter()
        .filter_map(|line| {
            let s = bytes_to_string(line);
            if let Some(rest) = s.strip_prefix("##INFO=<ID=") {
                rest.split(',').next().map(|id| id.to_string())
            } else {
                None
            }
        })
        .collect();

    for line in header_lines {
        if line.starts_with(b"#CHROM") {
            continue;
        }
        if line.starts_with(b"##INFO=<ID=SRC,") {
            continue;
        }
        writer.write_line(line)?;
    }

    writer.write_line(
        b"##INFO=<ID=SRC,Number=.,Type=String,Description=\"Annotation sources matched for this record\">\n",
    )?;

    for spec in [SNV_INDEL, SV, FLARE] {
        for &(src_key, dst_key) in spec.info_fields {
            if existing_info.contains(dst_key) {
                continue;
            }
            let line = format!(
                "##INFO=<ID={dst_key},Number=.,Type=String,Description=\"Propagated from {label} INFO/{src_key}\">\n",
                label = spec.label,
            );
            writer.write_line(line.as_bytes())?;
        }
    }

    for id in all_output_format_ids() {
        if existing_fmt.contains(id) {
            continue;
        }
        let line = format!(
            "##FORMAT=<ID={id},Number={num},Type={ty},Description=\"{desc}\">\n",
            id = id,
            num = format_number(id),
            ty = format_type(id),
            desc = format_description(id),
        );
        writer.write_line(line.as_bytes())?;
    }

    let chrom_line = format!(
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{}\n",
        samples.join("\t"),
    );
    writer.write_line(chrom_line.as_bytes())?;
    Ok(())
}

fn merge_info_from_records(
    matching: &[Option<AnnotRecord>],
    streams: &[AnnotationStream],
) -> (String, Vec<&'static str>) {
    use std::collections::HashMap;
    let mut extras = Vec::new();
    let mut labels = Vec::new();
    for (rec, stream) in matching.iter().zip(streams.iter()) {
        let Some(rec) = rec else {
            continue;
        };
        labels.push(stream.spec.label);
        let info_map: HashMap<String, String> = parse_info(&rec.info).into_iter().collect();
        for &(src_key, dst_key) in stream.spec.info_fields {
            if let Some(val) = info_map.get(src_key) {
                extras.push((dst_key.to_string(), val.clone()));
            }
        }
    }
    let mut merged = merge_info(".", &extras);
    if !labels.is_empty() {
        merged = merge_info(&merged, &[("SRC".to_string(), labels.join(","))]);
    }
    (merged, labels)
}

struct SourceMatchMeta {
    matched: bool,
    gt_filled: u64,
    gt_total: u64,
}

/// Pre-parsed colon fields for one matched source at the current site.
struct SourceSampleView<'a> {
    tag_idx: Vec<Option<usize>>,
    samples: &'a [String],
    n_src_tags: usize,
}

fn merge_sample_fields(
    n_samples: usize,
    matching: &[Option<AnnotRecord>],
    streams: &[AnnotationStream],
    dense: bool,
) -> (String, Vec<String>, Vec<SourceMatchMeta>) {
    let mut meta: Vec<SourceMatchMeta> = (0..streams.len())
        .map(|_| SourceMatchMeta {
            matched: false,
            gt_filled: 0,
            gt_total: n_samples as u64,
        })
        .collect();

    let matched_labels: Vec<&str> = matching
        .iter()
        .zip(streams.iter())
        .filter_map(|(rec, stream)| rec.as_ref().map(|_| stream.spec.label))
        .collect();

    let out_ids: Vec<&'static str> = if dense {
        all_output_format_ids()
    } else {
        format_ids_for_labels(&matched_labels)
    };
    let out_format = out_ids.join(":");

    let mut slot_src: Vec<Option<(usize, usize)>> = vec![None; out_ids.len()];
    let mut views: Vec<Option<SourceSampleView<'_>>> = Vec::with_capacity(streams.len());

    for (si, (rec, stream)) in matching.iter().zip(streams.iter()).enumerate() {
        let Some(rec) = rec else {
            views.push(None);
            continue;
        };
        meta[si].matched = true;
        let tags: Vec<&str> = if rec.format.is_empty() {
            Vec::new()
        } else {
            rec.format.split(':').collect()
        };
        let n_src_tags = tags.len();
        let tag_idx: Vec<Option<usize>> = out_ids
            .iter()
            .map(|dst| {
                stream
                    .spec
                    .format_fields
                    .iter()
                    .find(|&&(_, d)| d == *dst)
                    .and_then(|&(src, _)| tags.iter().position(|t| *t == src))
            })
            .collect();
        for (oi, dst) in out_ids.iter().enumerate() {
            if slot_src[oi].is_some() {
                continue; // first matched source wins on shared IDs (GQ/DP/AD)
            }
            if stream
                .spec
                .format_fields
                .iter()
                .any(|&(_, d)| d == *dst)
            {
                slot_src[oi] = Some((si, oi));
            }
        }
        // Count GT fill rate from primary FORMAT tag.
        let primary = stream.spec.format_fields[0].0;
        if let Some(gt_i) = tags.iter().position(|t| *t == primary) {
            for sample in &rec.samples {
                let val = nth_field(sample, gt_i, n_src_tags);
                if val != "." {
                    meta[si].gt_filled += 1;
                }
            }
        }
        views.push(Some(SourceSampleView {
            tag_idx,
            samples: &rec.samples,
            n_src_tags,
        }));
    }

    let samples = if n_samples >= 2_000 {
        (0..n_samples)
            .into_par_iter()
            .map(|sample_i| {
                build_sample_genotype(&out_ids, &slot_src, &views, sample_i)
            })
            .collect()
    } else {
        (0..n_samples)
            .map(|sample_i| build_sample_genotype(&out_ids, &slot_src, &views, sample_i))
            .collect()
    };

    (out_format, samples, meta)
}

fn build_sample_genotype(
    out_ids: &[&str],
    slot_src: &[Option<(usize, usize)>],
    views: &[Option<SourceSampleView<'_>>],
    sample_i: usize,
) -> String {
    if out_ids.is_empty() {
        return ".".to_string();
    }
    let mut out = String::with_capacity(out_ids.len() * 8);
    for (oi, _) in out_ids.iter().enumerate() {
        if oi > 0 {
            out.push(':');
        }
        let val = match slot_src[oi] {
            Some((si, out_slot)) => {
                if let Some(view) = &views[si] {
                    match view.tag_idx[out_slot] {
                        Some(ti) => {
                            let raw = view.samples.get(sample_i).map(|s| s.as_str()).unwrap_or(".");
                            nth_field(raw, ti, view.n_src_tags)
                        }
                        None => ".",
                    }
                } else {
                    "."
                }
            }
            None => ".",
        };
        out.push_str(val);
    }
    out
}

fn nth_field(sample: &str, idx: usize, n_tags: usize) -> &str {
    if sample == "." || n_tags == 0 {
        return ".";
    }
    let mut start = 0usize;
    let bytes = sample.as_bytes();
    let mut seen = 0usize;
    for (i, &b) in bytes.iter().enumerate() {
        if b == b':' {
            if seen == idx {
                let field = &sample[start..i];
                return if field.is_empty() { "." } else { field };
            }
            seen += 1;
            start = i + 1;
        }
    }
    if seen == idx {
        let field = &sample[start..];
        return if field.is_empty() { "." } else { field };
    }
    "."
}
