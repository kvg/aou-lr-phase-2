use std::collections::{HashMap, VecDeque};
use std::io::BufRead;
use std::path::{Path, PathBuf};

use crate::error::{PropagateError, Result};
use crate::fields::SourceSpec;
use crate::vcf_io::{
    bytes_to_string, open_vcf_reader, parse_data_line, split_tab_max9, trim_crlf, VariantKey,
};

/// One parsed VCF record with samples remapped to the output sample order.
#[derive(Debug, Clone)]
pub struct AnnotRecord {
    pub chrom: String,
    pub key: VariantKey,
    pub id: String,
    pub qual: String,
    pub filter: String,
    pub format: String,
    pub samples: Vec<String>,
    pub info: String,
}

/// Streaming annotation source.
///
/// Requires inputs sorted by `(CHROM, POS)` only (bcftools / TBI order). Records
/// that share a position are buffered and emitted in `(REF, ALT)` order so the
/// multi-way merge key `(CHROM, POS, REF, ALT)` stays consistent — SV callsets
/// often have several ALTs at one POS without REF/ALT sort order.
pub struct AnnotationStream {
    pub path: PathBuf,
    pub spec: SourceSpec,
    pub sample_ids: Vec<String>,
    pub sample_map: Vec<Option<usize>>,
    pub header_lines: Vec<Vec<u8>>,
    reader: Box<dyn BufRead>,
    line_buf: Vec<u8>,
    line_no: usize,
    pub current: Option<AnnotRecord>,
    pub sites_read: u64,
    /// Next records already read (same POS batch), sorted by VariantKey.
    pending: VecDeque<AnnotRecord>,
    /// First record of the next POS batch (read-ahead), if any.
    lookahead: Option<AnnotRecord>,
    last_chrom_pos: Option<(String, u64)>,
}

impl AnnotationStream {
    pub fn open(
        path: PathBuf,
        spec: SourceSpec,
        output_samples: &[String],
    ) -> Result<Self> {
        if !path.exists() {
            return Err(PropagateError::PathMissing(path.clone()));
        }
        let (sample_ids, header_lines) = read_header(&path)?;
        let sample_map = map_samples(output_samples, &sample_ids);
        let (mut reader, _) = open_vcf_reader(&path)?;
        skip_to_data(&mut reader)?;
        let mut stream = Self {
            path,
            spec,
            sample_ids,
            sample_map,
            header_lines,
            reader,
            line_buf: Vec::with_capacity(1 << 16),
            line_no: 0,
            current: None,
            sites_read: 0,
            pending: VecDeque::new(),
            lookahead: None,
            last_chrom_pos: None,
        };
        stream.advance()?;
        Ok(stream)
    }

    pub fn read_sample_ids(path: &Path) -> Result<Vec<String>> {
        Ok(read_header(path)?.0)
    }

    /// Advance to the next data record. Returns Ok(None) at EOF via `current`.
    pub fn advance(&mut self) -> Result<()> {
        if let Some(rec) = self.pending.pop_front() {
            self.current = Some(rec);
            return Ok(());
        }
        self.fill_pos_batch()?;
        self.current = self.pending.pop_front();
        Ok(())
    }

    /// Read one `(chrom, pos)` group into `pending`, sorted by `(REF, ALT)`.
    fn fill_pos_batch(&mut self) -> Result<()> {
        debug_assert!(self.pending.is_empty());
        let first = match self.lookahead.take() {
            Some(r) => r,
            None => match self.read_one()? {
                Some(r) => r,
                None => return Ok(()),
            },
        };
        let chrom = first.chrom.clone();
        let pos = first.key.pos;
        self.check_chrom_pos_order(&chrom, pos)?;

        let mut batch = vec![first];
        loop {
            match self.read_one()? {
                None => break,
                Some(rec) => {
                    if rec.chrom == chrom && rec.key.pos == pos {
                        batch.push(rec);
                    } else {
                        self.lookahead = Some(rec);
                        break;
                    }
                }
            }
        }
        batch.sort_by(|a, b| a.key.cmp(&b.key));
        self.sites_read += batch.len() as u64;
        self.pending.extend(batch);
        Ok(())
    }

    fn check_chrom_pos_order(&mut self, chrom: &str, pos: u64) -> Result<()> {
        if let Some((prev_chrom, prev_pos)) = &self.last_chrom_pos {
            let went_back = chrom < prev_chrom.as_str()
                || (chrom == prev_chrom.as_str() && pos < *prev_pos);
            if went_back {
                return Err(PropagateError::Unsorted {
                    path: self.path.clone(),
                    line: self.line_no,
                    message: format!(
                        "record {chrom}:{pos} goes before previous {prev_chrom}:{prev_pos}"
                    ),
                });
            }
        }
        self.last_chrom_pos = Some((chrom.to_string(), pos));
        Ok(())
    }

    fn read_one(&mut self) -> Result<Option<AnnotRecord>> {
        loop {
            self.line_buf.clear();
            let n = self.reader.read_until(b'\n', &mut self.line_buf)?;
            if n == 0 {
                return Ok(None);
            }
            self.line_no += 1;
            if self.line_buf.starts_with(b"#") {
                continue;
            }
            let (chrom, key, id, qual, filter, info, format, samples) =
                parse_data_line(&self.line_buf, &self.path, self.line_no)?;
            let remapped = remap_samples(&samples, &self.sample_map);
            return Ok(Some(AnnotRecord {
                chrom,
                key,
                id,
                qual,
                filter,
                format,
                samples: remapped,
                info,
            }));
        }
    }

    pub fn peek_key(&self) -> Option<(&str, &VariantKey)> {
        self.current
            .as_ref()
            .map(|r| (r.chrom.as_str(), &r.key))
    }
}

fn skip_to_data(reader: &mut Box<dyn BufRead>) -> Result<()> {
    let mut line_buf = Vec::with_capacity(4096);
    loop {
        line_buf.clear();
        let n = reader.read_until(b'\n', &mut line_buf)?;
        if n == 0 {
            return Ok(());
        }
        if line_buf.starts_with(b"#CHROM") {
            return Ok(());
        }
    }
}

fn read_header(path: &Path) -> Result<(Vec<String>, Vec<Vec<u8>>)> {
    let (mut reader, _) = open_vcf_reader(path)?;
    let mut header_lines = Vec::new();
    let mut line_buf = Vec::with_capacity(1 << 16);
    let mut sample_ids = Vec::new();

    loop {
        line_buf.clear();
        let n = reader.read_until(b'\n', &mut line_buf)?;
        if n == 0 {
            break;
        }
        if line_buf.starts_with(b"##") {
            header_lines.push(line_buf.clone());
            continue;
        }
        if line_buf.starts_with(b"#CHROM") {
            header_lines.push(line_buf.clone());
            let stripped = trim_crlf(&line_buf);
            let parts = split_tab_max9(stripped).ok_or_else(|| PropagateError::TooFewFields {
                path: path.to_path_buf(),
                line: header_lines.len(),
            })?;
            let blob = trim_crlf(parts[9]);
            sample_ids = if blob.is_empty() {
                Vec::new()
            } else {
                crate::vcf_io::split_tabs(blob)
                    .into_iter()
                    .map(bytes_to_string)
                    .collect()
            };
            break;
        }
    }
    Ok((sample_ids, header_lines))
}

pub fn union_sample_lists(lists: &[Vec<String>]) -> Result<Vec<String>> {
    let mut seen = HashMap::new();
    let mut out = Vec::new();
    for samples in lists {
        for s in samples {
            if seen.insert(s.clone(), ()).is_none() {
                out.push(s.clone());
            }
        }
    }
    if out.is_empty() {
        return Err(PropagateError::NoSamples);
    }
    Ok(out)
}

fn map_samples(output_samples: &[String], source_samples: &[String]) -> Vec<Option<usize>> {
    let index: HashMap<&str, usize> = source_samples
        .iter()
        .enumerate()
        .map(|(i, s)| (s.as_str(), i))
        .collect();
    output_samples
        .iter()
        .map(|s| index.get(s.as_str()).copied())
        .collect()
}

fn remap_samples(samples: &[String], sample_map: &[Option<usize>]) -> Vec<String> {
    sample_map
        .iter()
        .map(|idx| match idx {
            Some(i) => samples.get(*i).cloned().unwrap_or_else(|| ".".to_string()),
            None => ".".to_string(),
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fields::SV;
    use std::io::Write;

    #[test]
    fn remap_samples_maps_indices() {
        let src = vec!["0/0".into(), "0/1".into(), "1/1".into()];
        let map = vec![Some(2), Some(0)];
        let out = remap_samples(&src, &map);
        assert_eq!(out, vec!["1/1", "0/0"]);
    }

    #[test]
    fn union_sample_lists_preserves_first_seen_order() {
        let a = vec!["S1".into(), "S2".into()];
        let b = vec!["S2".into(), "S3".into()];
        let out = union_sample_lists(&[a, b]).unwrap();
        assert_eq!(out, vec!["S1", "S2", "S3"]);
    }

    #[test]
    fn same_pos_alts_emitted_in_ref_alt_order() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("sv.vcf");
        {
            let mut f = std::fs::File::create(&path).unwrap();
            write!(
                f,
                "##fileformat=VCFv4.2\n\
                 #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n\
                 chr22\t100\t.\tN\t<DUP>\t.\t.\t.\tGT\t0/1\n\
                 chr22\t100\t.\tN\t<DEL>\t.\t.\t.\tGT\t0/1\n\
                 chr22\t200\t.\tN\t<INS>\t.\t.\t.\tGT\t0/1\n"
            )
            .unwrap();
        }
        let mut stream = AnnotationStream::open(path, SV, &["S1".into()]).unwrap();
        let k1 = stream.peek_key().unwrap().1.clone();
        assert_eq!(k1.alt_allele, "<DEL>"); // DEL < DUP lexicographically
        stream.advance().unwrap();
        let k2 = stream.peek_key().unwrap().1.clone();
        assert_eq!(k2.alt_allele, "<DUP>");
        stream.advance().unwrap();
        let k3 = stream.peek_key().unwrap().1.clone();
        assert_eq!(k3.pos, 200);
    }
}
