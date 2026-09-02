use std::fs::File;
use std::io::{BufRead, BufReader, Read};
use std::path::Path;

use flate2::read::MultiGzDecoder;

use crate::error::{Result, ScoreError};

#[derive(Debug, Clone)]
pub struct VariantMeta {
    pub chrom: String,
    pub pos: String,
    pub id: String,
    pub r#ref: String,
    pub alt: String,
}

pub struct DosageReader {
    reader: Box<dyn BufRead>,
    sample_cols: Vec<usize>,
    n_samples_file: usize,
}

impl DosageReader {
    pub fn open(path: &Path, id_include: &[String], file_sample_ids: &[String]) -> Result<Self> {
        let sample_cols = map_sample_columns(id_include, file_sample_ids)?;
        let reader: Box<dyn BufRead> = if path
            .extension()
            .is_some_and(|e| e == "gz")
        {
            Box::new(BufReader::with_capacity(
                1 << 20,
                MultiGzDecoder::new(File::open(path)?),
            ))
        } else {
            Box::new(BufReader::with_capacity(1 << 20, File::open(path)?))
        };
        Ok(Self {
            reader,
            sample_cols,
            n_samples_file: file_sample_ids.len(),
        })
    }

    pub fn skip_header(&mut self) -> Result<()> {
        let mut line = String::new();
        self.reader.read_line(&mut line)?;
        Ok(())
    }

    pub fn read_chunk(
        &mut self,
        nrows: usize,
        out_genotypes: &mut Vec<f64>,
        out_meta: &mut Vec<VariantMeta>,
    ) -> Result<usize> {
        out_meta.clear();
        out_genotypes.clear();
        let n_out = self.sample_cols.len();
        let mut read = 0usize;
        while read < nrows {
            let mut line = String::new();
            let nbytes = self.reader.read_line(&mut line)?;
            if nbytes == 0 {
                break;
            }
            let line = line.trim_end_matches(['\r', '\n']);
            if line.is_empty() {
                continue;
            }
            let fields: Vec<&str> = line.split('\t').collect();
            if fields.len() < 5 + self.n_samples_file {
                return Err(ScoreError::msg(format!(
                    "malformed dosage line with {} fields",
                    fields.len()
                )));
            }
            out_meta.push(VariantMeta {
                chrom: fields[0].to_string(),
                pos: fields[1].to_string(),
                id: fields[2].to_string(),
                r#ref: fields[3].to_string(),
                alt: fields[4].to_string(),
            });
            let base = out_genotypes.len();
            out_genotypes.resize(base + n_out, 0.0);
            for (j, &col) in self.sample_cols.iter().enumerate() {
                let raw = fields[5 + col];
                let v: f64 = if raw == "NA" || raw == "." || raw.is_empty() {
                    0.0
                } else {
                    raw.parse().map_err(|_| {
                        ScoreError::msg(format!("invalid genotype {raw}"))
                    })?
                };
                out_genotypes[base + j] = v;
            }
            read += 1;
        }
        Ok(read)
    }
}

pub fn read_dosage_header(path: &Path) -> Result<Vec<String>> {
    let mut reader: Box<dyn BufRead> = if path.extension().is_some_and(|e| e == "gz") {
        Box::new(BufReader::new(MultiGzDecoder::new(File::open(path)?)))
    } else {
        Box::new(BufReader::new(File::open(path)?))
    };
    let mut line = String::new();
    reader.read_line(&mut line)?;
    let fields: Vec<&str> = line.trim_end().split('\t').collect();
    if fields.len() < 6 {
        return Err(ScoreError::msg("dosage header missing sample columns"));
    }
    Ok(fields[5..].iter().map(|s| s.to_string()).collect())
}

fn map_sample_columns(id_include: &[String], file_sample_ids: &[String]) -> Result<Vec<usize>> {
    let mut cols = Vec::with_capacity(id_include.len());
    for id in id_include {
        let pos = file_sample_ids
            .iter()
            .position(|s| s == id)
            .ok_or_else(|| ScoreError::msg(format!("sample {id} missing from dosage header")))?;
        cols.push(pos);
    }
    Ok(cols)
}

pub fn count_variants(path: &Path) -> Result<usize> {
    let mut reader: Box<dyn Read> = if path.extension().is_some_and(|e| e == "gz") {
        Box::new(MultiGzDecoder::new(File::open(path)?))
    } else {
        Box::new(File::open(path)?)
    };
    let mut buf = [0u8; 8192];
    let mut lines = 0usize;
    let mut partial = false;
    loop {
        let n = reader.read(&mut buf)?;
        if n == 0 {
            break;
        }
        for &b in &buf[..n] {
            if b == b'\n' {
                lines += 1;
                partial = false;
            } else {
                partial = true;
            }
        }
    }
    if partial {
        lines += 1;
    }
    Ok(lines.saturating_sub(1))
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    use tempfile::NamedTempFile;

    #[test]
    fn parses_float_and_integer_dosages() {
        let mut f = NamedTempFile::new().unwrap();
        writeln!(f, "CHROM\tPOS\tID\tREF\tALT\ts1\ts2").unwrap();
        writeln!(f, "1\t100\tv1\tA\tG\t0.5\t2").unwrap();
        f.flush().unwrap();
        let ids = vec!["s1".into(), "s2".into()];
        let include = vec!["s1".into(), "s2".into()];
        let mut r = DosageReader::open(f.path(), &include, &ids).unwrap();
        r.skip_header().unwrap();
        let mut g = Vec::new();
        let mut m = Vec::new();
        let n = r.read_chunk(10, &mut g, &mut m).unwrap();
        assert_eq!(n, 1);
        assert!((g[0] - 0.5).abs() < 1e-9);
        assert!((g[1] - 2.0).abs() < 1e-9);
    }
}
