use std::fs::File;
use std::io::{BufRead, BufReader, Read};
use std::path::{Path, PathBuf};

use flate2::read::MultiGzDecoder;

use crate::error::{ExtractError, Result};
use crate::io_out::OutWriter;
use crate::parse::{
    apply_haplotype, first_four_tokens, parse_nonneg_int, parse_vcf_filename, split_tab_max9,
    split_tabs, strip_ascii_ws, strip_hash_space, trim_crlf, write_ancestry_gt,
};
use crate::ru::{
    count_alts, haplotype_c, nth_alt, parse_ru_site, push_number, split_locus_id, RuSite,
};

const READ_BUF: usize = 1 << 20;

#[derive(Debug, Clone)]
pub struct ExtractConfig {
    pub vcf: PathBuf,
    pub num_ancs: usize,
    pub output_dir: Option<PathBuf>,
    pub output_vcf: bool,
    pub compress_output: bool,
    pub samples: Option<PathBuf>,
    /// Accepted for CLI compatibility; the extract loop is currently single-threaded.
    pub threads: usize,
}

pub fn extract_tracts_flare(cfg: &ExtractConfig) -> Result<()> {
    if !cfg.vcf.exists() {
        return Err(ExtractError::PathMissing(cfg.vcf.clone()));
    }
    if cfg.num_ancs == 0 {
        return Err(ExtractError::BadNumAncs(cfg.num_ancs));
    }
    let _ = cfg.threads;
    if let Some(dir) = &cfg.output_dir {
        if !dir.exists() {
            return Err(ExtractError::OutputDirMissing(dir.clone()));
        }
    }

    let filename = cfg
        .vcf
        .file_name()
        .and_then(|s| s.to_str())
        .ok_or_else(|| ExtractError::BadExtension(cfg.vcf.display().to_string()))?;
    let (prefix, zipped) = parse_vcf_filename(filename).map_err(ExtractError::BadExtension)?;

    let out_dir: &Path = cfg
        .output_dir
        .as_deref()
        .or_else(|| cfg.vcf.parent())
        .unwrap_or_else(|| Path::new("."));
    let output_path = out_dir.join(&prefix);
    let ext = if cfg.compress_output { ".gz" } else { "" };

    let keep = match &cfg.samples {
        Some(p) => Some(read_sample_ids(p)?),
        None => None,
    };

    let file = File::open(&cfg.vcf)?;
    let mut reader: Box<dyn BufRead> = if zipped {
        Box::new(BufReader::with_capacity(
            READ_BUF,
            MultiGzDecoder::new(file),
        ))
    } else {
        Box::new(BufReader::with_capacity(READ_BUF, file))
    };

    let mut line_buf = Vec::with_capacity(1 << 16);
    let mut line_no = 0usize;
    let mut vcf_header = Vec::new();
    let mut seen_chrom_header = false;

    let mut dosage_writers: Vec<OutWriter> = Vec::new();
    let mut hap_writers: Vec<OutWriter> = Vec::new();
    let mut collapse_writers: Vec<OutWriter> = Vec::new();
    let mut split_writers: Vec<OutWriter> = Vec::new();
    let mut vcf_writers: Vec<OutWriter> = Vec::new();

    let mut vcf_sample_ids: Vec<String> = Vec::new();
    let mut out_indices: Vec<usize> = Vec::new();

    let mut dosage_row = vec![Vec::new(); cfg.num_ancs];
    let mut hap_row = vec![Vec::new(); cfg.num_ancs];
    let mut collapse_row = vec![Vec::new(); cfg.num_ancs];
    let mut vcf_row = vec![Vec::new(); cfg.num_ancs];
    let mut counts_dos = vec![0u8; cfg.num_ancs];
    let mut counts_hap = vec![0u8; cfg.num_ancs];
    let mut ru_dos = vec![0.0f64; cfg.num_ancs];
    let mut coll_dos = vec![0u8; cfg.num_ancs];

    loop {
        line_buf.clear();
        let n = reader.read_until(b'\n', &mut line_buf)?;
        if n == 0 {
            break;
        }
        line_no += 1;

        if line_buf.starts_with(b"##") {
            if cfg.output_vcf && !seen_chrom_header {
                vcf_header.extend_from_slice(&line_buf);
                if !line_buf.ends_with(b"\n") {
                    vcf_header.push(b'\n');
                }
            }
            continue;
        }

        if line_buf.starts_with(b"#") {
            seen_chrom_header = true;
            let stripped = strip_hash_space(&line_buf);
            let parts = split_tab_max9(stripped).ok_or_else(|| ExtractError::TooFewFields {
                path: cfg.vcf.clone(),
                line: line_no,
            })?;
            let sample_blob = trim_crlf(parts[9]);
            vcf_sample_ids = if sample_blob.is_empty() {
                Vec::new()
            } else {
                split_tabs(sample_blob)
                    .into_iter()
                    .map(|s| String::from_utf8_lossy(s).into_owned())
                    .collect()
            };

            let mut seen = std::collections::HashMap::with_capacity(vcf_sample_ids.len());
            for (i, id) in vcf_sample_ids.iter().enumerate() {
                if seen.insert(id.clone(), i).is_some() {
                    return Err(ExtractError::DuplicateSample {
                        path: cfg.vcf.clone(),
                        sample: id.clone(),
                    });
                }
            }

            let out_sample_ids = if let Some(keep_ids) = &keep {
                if keep_ids.is_empty() {
                    return Err(ExtractError::EmptySamples(cfg.samples.clone().unwrap()));
                }
                let mut missing = Vec::new();
                out_indices.clear();
                for id in keep_ids {
                    match seen.get(id) {
                        Some(&i) => out_indices.push(i),
                        None => missing.push(id.clone()),
                    }
                }
                if !missing.is_empty() {
                    return Err(ExtractError::MissingSamples {
                        path: cfg.vcf.clone(),
                        count: missing.len(),
                        examples: missing.into_iter().take(5).collect(),
                    });
                }
                keep_ids.clone()
            } else {
                out_indices = (0..vcf_sample_ids.len()).collect();
                vcf_sample_ids.clone()
            };

            dosage_writers.reserve(cfg.num_ancs);
            hap_writers.reserve(cfg.num_ancs);
            collapse_writers.reserve(cfg.num_ancs);
            split_writers.reserve(cfg.num_ancs);
            for i in 0..cfg.num_ancs {
                let dos = format!("{}.anc{i}.dosage.txt{ext}", output_path.display());
                let hap = format!("{}.anc{i}.hapcount.txt{ext}", output_path.display());
                let coll = format!("{}.collapse.anc{i}.dosage.txt{ext}", output_path.display());
                let split = format!("{}.split.anc{i}.dosage.txt{ext}", output_path.display());
                dosage_writers.push(OutWriter::create(Path::new(&dos), cfg.compress_output)?);
                hap_writers.push(OutWriter::create(Path::new(&hap), cfg.compress_output)?);
                collapse_writers.push(OutWriter::create(Path::new(&coll), cfg.compress_output)?);
                split_writers.push(OutWriter::create(Path::new(&split), cfg.compress_output)?);
                if cfg.output_vcf {
                    let v = format!("{}.anc{i}.vcf{ext}", output_path.display());
                    vcf_writers.push(OutWriter::create(Path::new(&v), cfg.compress_output)?);
                }
            }

            write_dosage_header(
                &mut dosage_writers,
                &mut hap_writers,
                parts,
                &out_sample_ids,
            )?;
            write_sample_header(&mut collapse_writers, parts, &out_sample_ids)?;
            write_sample_header(&mut split_writers, parts, &out_sample_ids)?;
            if cfg.output_vcf {
                let projected = keep.is_some();
                for w in &mut vcf_writers {
                    w.write_all(&vcf_header)?;
                    if projected {
                        write_chrom_header_vcf(w, &line_buf, &vcf_sample_ids, &out_indices)?;
                    } else {
                        w.write_all(trim_crlf(&line_buf))?;
                        w.write_all(b"\n")?;
                    }
                }
            }

            write_sample_order(out_dir, &out_sample_ids)?;
            continue;
        }

        if !seen_chrom_header {
            continue;
        }

        process_data_line(
            cfg,
            line_no,
            &line_buf,
            &vcf_sample_ids,
            &out_indices,
            &mut dosage_writers,
            &mut hap_writers,
            &mut collapse_writers,
            &mut split_writers,
            &mut vcf_writers,
            &mut dosage_row,
            &mut hap_row,
            &mut collapse_row,
            &mut vcf_row,
            &mut counts_dos,
            &mut counts_hap,
            &mut ru_dos,
            &mut coll_dos,
        )?;
    }

    if !seen_chrom_header {
        return Err(ExtractError::MissingHeader(cfg.vcf.clone()));
    }

    for w in dosage_writers {
        w.finish()?;
    }
    for w in hap_writers {
        w.finish()?;
    }
    for w in collapse_writers {
        w.finish()?;
    }
    for w in split_writers {
        w.finish()?;
    }
    for w in vcf_writers {
        w.finish()?;
    }
    Ok(())
}

fn write_sample_order(out_dir: &Path, ids: &[String]) -> Result<()> {
    let path = out_dir.join("dosage_sample_order.txt");
    let mut f = File::create(path)?;
    use std::io::Write;
    for id in ids {
        writeln!(f, "{id}")?;
    }
    Ok(())
}

fn write_dosage_header(
    dosage: &mut [OutWriter],
    hap: &mut [OutWriter],
    parts: [&[u8]; 10],
    out_samples: &[String],
) -> Result<()> {
    // CHROM POS ID REF ALT + selected sample IDs + newline
    let mut header = Vec::with_capacity(64 + out_samples.len() * 16);
    header.extend_from_slice(trim_crlf(parts[0]));
    for idx in 1..5 {
        header.push(b'\t');
        header.extend_from_slice(trim_crlf(parts[idx]));
    }
    for id in out_samples {
        header.push(b'\t');
        header.extend_from_slice(id.as_bytes());
    }
    header.push(b'\n');
    for w in dosage.iter_mut().chain(hap.iter_mut()) {
        w.write_all(&header)?;
    }
    Ok(())
}

fn write_sample_header(
    writers: &mut [OutWriter],
    parts: [&[u8]; 10],
    out_samples: &[String],
) -> Result<()> {
    let mut header = Vec::with_capacity(64 + out_samples.len() * 16);
    header.extend_from_slice(trim_crlf(parts[0]));
    for idx in 1..5 {
        header.push(b'\t');
        header.extend_from_slice(trim_crlf(parts[idx]));
    }
    for id in out_samples {
        header.push(b'\t');
        header.extend_from_slice(id.as_bytes());
    }
    header.push(b'\n');
    for w in writers.iter_mut() {
        w.write_all(&header)?;
    }
    Ok(())
}

fn write_chrom_header_vcf(
    w: &mut OutWriter,
    original_line: &[u8],
    vcf_samples: &[String],
    out_indices: &[usize],
) -> Result<()> {
    // Copy ## already written; rewrite #CHROM with FORMAT as in the original
    // line, projecting sample columns.
    let trimmed = trim_crlf(original_line);
    let parts = split_tab_max9(trimmed).ok_or_else(|| ExtractError::TooFewFields {
        path: PathBuf::from("<header>"),
        line: 0,
    })?;
    let mut buf = Vec::new();
    for i in 0..9 {
        if i > 0 {
            buf.push(b'\t');
        }
        buf.extend_from_slice(parts[i]);
    }
    for &idx in out_indices {
        buf.push(b'\t');
        buf.extend_from_slice(vcf_samples[idx].as_bytes());
    }
    buf.push(b'\n');
    w.write_all(&buf)?;
    Ok(())
}

#[allow(clippy::too_many_arguments)]
fn process_data_line(
    cfg: &ExtractConfig,
    line_no: usize,
    line_buf: &[u8],
    vcf_sample_ids: &[String],
    out_indices: &[usize],
    dosage_writers: &mut [OutWriter],
    hap_writers: &mut [OutWriter],
    collapse_writers: &mut [OutWriter],
    split_writers: &mut [OutWriter],
    vcf_writers: &mut [OutWriter],
    dosage_row: &mut [Vec<u8>],
    hap_row: &mut [Vec<u8>],
    collapse_row: &mut [Vec<u8>],
    vcf_row: &mut [Vec<u8>],
    counts_dos: &mut [u8],
    counts_hap: &mut [u8],
    ru_dos: &mut [f64],
    coll_dos: &mut [u8],
) -> Result<()> {
    let stripped = strip_ascii_ws(line_buf);
    if stripped.is_empty() {
        return Ok(());
    }
    let parts = split_tab_max9(stripped).ok_or_else(|| ExtractError::TooFewFields {
        path: cfg.vcf.clone(),
        line: line_no,
    })?;
    let genos = split_tabs(parts[9]);
    if genos.len() != vcf_sample_ids.len() {
        return Err(ExtractError::TooFewFields {
            path: cfg.vcf.clone(),
            line: line_no,
        });
    }

    if let Some(ru) = parse_ru_site(parts[7], parts[4]) {
        return process_ru_line(
            cfg,
            line_no,
            parts,
            &genos,
            vcf_sample_ids,
            out_indices,
            &ru,
            dosage_writers,
            hap_writers,
            collapse_writers,
            split_writers,
            vcf_writers,
            dosage_row,
            hap_row,
            collapse_row,
            vcf_row,
            counts_hap,
            ru_dos,
            coll_dos,
        );
    }

    for j in 0..cfg.num_ancs {
        dosage_row[j].clear();
        hap_row[j].clear();
        for i in 0..5 {
            if i > 0 {
                dosage_row[j].push(b'\t');
                hap_row[j].push(b'\t');
            }
            dosage_row[j].extend_from_slice(parts[i]);
            hap_row[j].extend_from_slice(parts[i]);
        }
        if cfg.output_vcf {
            vcf_row[j].clear();
            for i in 0..8 {
                if i > 0 {
                    vcf_row[j].push(b'\t');
                }
                vcf_row[j].extend_from_slice(parts[i]);
            }
            vcf_row[j].extend_from_slice(b"\tGT");
        }
    }

    for &si in out_indices {
        let field = genos[si];
        let toks = first_four_tokens(field).ok_or_else(|| ExtractError::TruncatedGenotype {
            path: cfg.vcf.clone(),
            line: line_no,
            sample: vcf_sample_ids
                .get(si)
                .cloned()
                .unwrap_or_else(|| format!("index {si}")),
            column: si + 10,
            field: String::from_utf8_lossy(field).into_owned(),
        })?;
        let (geno_a, geno_b, call_a, call_b) = (toks[0], toks[1], toks[2], toks[3]);
        counts_dos.fill(0);
        counts_hap.fill(0);
        apply_haplotype(geno_a, call_a, cfg.num_ancs, counts_dos, counts_hap);
        apply_haplotype(geno_b, call_b, cfg.num_ancs, counts_dos, counts_hap);
        for j in 0..cfg.num_ancs {
            dosage_row[j].push(b'\t');
            dosage_row[j].push(b'0' + counts_dos[j]);
            hap_row[j].push(b'\t');
            hap_row[j].push(b'0' + counts_hap[j]);
            if cfg.output_vcf {
                vcf_row[j].push(b'\t');
                write_ancestry_gt(&mut vcf_row[j], geno_a, geno_b, call_a, call_b, j);
            }
        }
    }

    for j in 0..cfg.num_ancs {
        dosage_row[j].push(b'\n');
        hap_row[j].push(b'\n');
        dosage_writers[j].write_all(&dosage_row[j])?;
        hap_writers[j].write_all(&hap_row[j])?;
        if cfg.output_vcf {
            vcf_row[j].push(b'\n');
            vcf_writers[j].write_all(&vcf_row[j])?;
        }
    }
    Ok(())
}

fn write_locus_prefix(buf: &mut Vec<u8>, parts: [&[u8]; 10]) {
    buf.clear();
    for i in 0..5 {
        if i > 0 {
            buf.push(b'\t');
        }
        buf.extend_from_slice(parts[i]);
    }
}

#[allow(clippy::too_many_arguments)]
fn apply_ru_haplotype(
    allele_tok: &[u8],
    ancestry_tok: &[u8],
    num_ancs: usize,
    site: &RuSite,
    n_alts: usize,
    ru_dos: &mut [f64],
    hap: &mut [u8],
    collapse: &mut [u8],
    split: &mut [Vec<u8>],
) {
    let Some(j) = parse_nonneg_int(ancestry_tok) else {
        return;
    };
    if j >= num_ancs {
        return;
    }
    hap[j] = hap[j].saturating_add(1);
    let Some(allele) = parse_nonneg_int(allele_tok) else {
        return;
    };
    ru_dos[j] += haplotype_c(allele, site);
    if allele >= 1 && allele <= n_alts {
        collapse[j] = collapse[j].saturating_add(1);
        split[allele - 1][j] = split[allele - 1][j].saturating_add(1);
    }
}

#[allow(clippy::too_many_arguments)]
fn process_ru_line(
    cfg: &ExtractConfig,
    line_no: usize,
    parts: [&[u8]; 10],
    genos: &[&[u8]],
    vcf_sample_ids: &[String],
    out_indices: &[usize],
    site: &RuSite,
    dosage_writers: &mut [OutWriter],
    hap_writers: &mut [OutWriter],
    collapse_writers: &mut [OutWriter],
    split_writers: &mut [OutWriter],
    vcf_writers: &mut [OutWriter],
    dosage_row: &mut [Vec<u8>],
    hap_row: &mut [Vec<u8>],
    collapse_row: &mut [Vec<u8>],
    vcf_row: &mut [Vec<u8>],
    counts_hap: &mut [u8],
    ru_dos: &mut [f64],
    coll_dos: &mut [u8],
) -> Result<()> {
    let n_alts = count_alts(parts[4]).max(site.ru.len());
    let mut split_counts: Vec<Vec<u8>> = (0..n_alts).map(|_| vec![0u8; cfg.num_ancs]).collect();
    let mut split_rows: Vec<Vec<Vec<u8>>> = (0..n_alts)
        .map(|ai| {
            let id = split_locus_id(parts[2], parts[0], parts[1], ai + 1);
            let alt = nth_alt(parts[4], ai).unwrap_or(b".");
            (0..cfg.num_ancs)
                .map(|_| {
                    let mut buf = Vec::new();
                    buf.extend_from_slice(parts[0]);
                    buf.push(b'\t');
                    buf.extend_from_slice(parts[1]);
                    buf.push(b'\t');
                    buf.extend_from_slice(&id);
                    buf.push(b'\t');
                    buf.extend_from_slice(parts[3]);
                    buf.push(b'\t');
                    buf.extend_from_slice(alt);
                    buf
                })
                .collect()
        })
        .collect();

    for j in 0..cfg.num_ancs {
        write_locus_prefix(&mut dosage_row[j], parts);
        write_locus_prefix(&mut hap_row[j], parts);
        write_locus_prefix(&mut collapse_row[j], parts);
        if cfg.output_vcf {
            vcf_row[j].clear();
            for i in 0..8 {
                if i > 0 {
                    vcf_row[j].push(b'\t');
                }
                vcf_row[j].extend_from_slice(parts[i]);
            }
            vcf_row[j].extend_from_slice(b"\tGT");
        }
    }

    for &si in out_indices {
        let field = genos[si];
        let toks = first_four_tokens(field).ok_or_else(|| ExtractError::TruncatedGenotype {
            path: cfg.vcf.clone(),
            line: line_no,
            sample: vcf_sample_ids
                .get(si)
                .cloned()
                .unwrap_or_else(|| format!("index {si}")),
            column: si + 10,
            field: String::from_utf8_lossy(field).into_owned(),
        })?;
        let (geno_a, geno_b, call_a, call_b) = (toks[0], toks[1], toks[2], toks[3]);
        ru_dos.fill(0.0);
        counts_hap.fill(0);
        coll_dos.fill(0);
        for row in split_counts.iter_mut() {
            row.fill(0);
        }
        apply_ru_haplotype(
            geno_a,
            call_a,
            cfg.num_ancs,
            site,
            n_alts,
            ru_dos,
            counts_hap,
            coll_dos,
            &mut split_counts,
        );
        apply_ru_haplotype(
            geno_b,
            call_b,
            cfg.num_ancs,
            site,
            n_alts,
            ru_dos,
            counts_hap,
            coll_dos,
            &mut split_counts,
        );
        for j in 0..cfg.num_ancs {
            dosage_row[j].push(b'\t');
            push_number(&mut dosage_row[j], ru_dos[j]);
            hap_row[j].push(b'\t');
            hap_row[j].push(b'0' + counts_hap[j]);
            collapse_row[j].push(b'\t');
            collapse_row[j].push(b'0' + coll_dos[j]);
            if cfg.output_vcf {
                vcf_row[j].push(b'\t');
                write_ancestry_gt(&mut vcf_row[j], geno_a, geno_b, call_a, call_b, j);
            }
            for ai in 0..n_alts {
                split_rows[ai][j].push(b'\t');
                split_rows[ai][j].push(b'0' + split_counts[ai][j]);
            }
        }
    }

    for j in 0..cfg.num_ancs {
        dosage_row[j].push(b'\n');
        hap_row[j].push(b'\n');
        collapse_row[j].push(b'\n');
        dosage_writers[j].write_all(&dosage_row[j])?;
        hap_writers[j].write_all(&hap_row[j])?;
        collapse_writers[j].write_all(&collapse_row[j])?;
        for ai in 0..n_alts {
            split_rows[ai][j].push(b'\n');
            split_writers[j].write_all(&split_rows[ai][j])?;
        }
        if cfg.output_vcf {
            vcf_row[j].push(b'\n');
            vcf_writers[j].write_all(&vcf_row[j])?;
        }
    }
    Ok(())
}

pub fn read_sample_ids(path: &Path) -> Result<Vec<String>> {
    if !path.exists() {
        return Err(ExtractError::NotFound(path.to_path_buf()));
    }
    let mut raw = String::new();
    File::open(path)?.read_to_string(&mut raw)?;
    let mut ids = Vec::new();
    for line in raw.lines() {
        let s = line.trim();
        if s.is_empty() {
            continue;
        }
        let id = s.split_whitespace().next().unwrap();
        ids.push(id.to_string());
    }
    Ok(ids)
}
