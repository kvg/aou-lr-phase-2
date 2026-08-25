use std::fs::File;
use std::io::{BufRead, BufReader, Write};
use std::path::Path;

use flate2::read::MultiGzDecoder;
use flate2::write::GzEncoder;
use flate2::Compression;
use memchr::memchr;

use crate::error::{PropagateError, Result};

const READ_BUF: usize = 1 << 20;

#[derive(Debug, Clone, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct VariantKey {
    pub pos: u64,
    pub ref_allele: String,
    pub alt_allele: String,
}

impl VariantKey {
    pub fn display(&self) -> String {
        format!("{}:{}:{}>{}", self.pos, self.ref_allele, self.alt_allele, "")
    }
}

pub fn parse_vcf_filename(filename: &str) -> Result<(String, bool)> {
    if filename.ends_with(".vcf.gz") {
        Ok((filename.trim_end_matches(".vcf.gz").to_string(), true))
    } else if filename.ends_with(".vcf") {
        Ok((filename.trim_end_matches(".vcf").to_string(), false))
    } else if filename.ends_with(".vcf.bgz") {
        Ok((filename.trim_end_matches(".vcf.bgz").to_string(), true))
    } else {
        Err(PropagateError::BadExtension(filename.to_string()))
    }
}

pub fn open_vcf_reader(path: &Path) -> Result<(Box<dyn BufRead>, bool)> {
    if !path.exists() {
        return Err(PropagateError::PathMissing(path.to_path_buf()));
    }
    let filename = path
        .file_name()
        .and_then(|s| s.to_str())
        .ok_or_else(|| PropagateError::BadExtension(path.display().to_string()))?;
    let (_, gz) = parse_vcf_filename(filename)?;
    let file = File::open(path)?;
    let reader: Box<dyn BufRead> = if gz {
        Box::new(BufReader::with_capacity(
            READ_BUF,
            MultiGzDecoder::new(file),
        ))
    } else {
        Box::new(BufReader::with_capacity(READ_BUF, file))
    };
    Ok((reader, gz))
}

pub struct VcfWriter {
    inner: WriterKind,
}

enum WriterKind {
    Plain(File),
    Gz(GzEncoder<File>),
}

impl VcfWriter {
    pub fn create(path: &Path, gzip: bool) -> Result<Self> {
        let file = File::create(path)?;
        let inner = if gzip {
            WriterKind::Gz(GzEncoder::new(file, Compression::default()))
        } else {
            WriterKind::Plain(file)
        };
        Ok(Self { inner })
    }

    pub fn write_line(&mut self, line: &[u8]) -> Result<()> {
        let w: &mut dyn Write = match &mut self.inner {
            WriterKind::Plain(f) => f,
            WriterKind::Gz(e) => e,
        };
        w.write_all(line)?;
        if !line.ends_with(b"\n") {
            w.write_all(b"\n")?;
        }
        Ok(())
    }

    pub fn finish(self) -> Result<()> {
        match self.inner {
            WriterKind::Plain(_) => Ok(()),
            WriterKind::Gz(encoder) => encoder.finish().map(|_| ()).map_err(PropagateError::Io),
        }
    }
}

pub fn split_tab_max9(line: &[u8]) -> Option<[&[u8]; 10]> {
    let mut parts = [&b""[..]; 10];
    let mut start = 0usize;
    for i in 0..9 {
        match memchr(b'\t', &line[start..]) {
            Some(rel) => {
                parts[i] = &line[start..start + rel];
                start += rel + 1;
            }
            None => return None,
        }
    }
    parts[9] = &line[start..];
    Some(parts)
}

pub fn split_tabs(s: &[u8]) -> Vec<&[u8]> {
    let mut out = Vec::new();
    let mut start = 0usize;
    while let Some(rel) = memchr(b'\t', &s[start..]) {
        out.push(&s[start..start + rel]);
        start += rel + 1;
    }
    out.push(&s[start..]);
    out
}

pub fn trim_crlf(s: &[u8]) -> &[u8] {
    let mut end = s.len();
    if end > 0 && s[end - 1] == b'\n' {
        end -= 1;
    }
    if end > 0 && s[end - 1] == b'\r' {
        end -= 1;
    }
    &s[..end]
}

pub fn bytes_to_string(b: &[u8]) -> String {
    String::from_utf8_lossy(b).into_owned()
}

pub fn parse_data_line(
    line: &[u8],
    path: &Path,
    line_no: usize,
) -> Result<(
    String,
    VariantKey,
    String,
    String,
    String,
    String,
    String,
    Vec<String>,
)> {
    let stripped = trim_crlf(line);
    let parts = split_tab_max9(stripped).ok_or_else(|| PropagateError::TooFewFields {
        path: path.to_path_buf(),
        line: line_no,
    })?;
    let chrom = bytes_to_string(parts[0]);
    let pos = bytes_to_string(parts[1])
        .parse::<u64>()
        .map_err(|_| PropagateError::ParseLine {
            path: path.to_path_buf(),
            line: line_no,
            message: format!("bad POS: {}", bytes_to_string(parts[1])),
        })?;
    let key = VariantKey {
        pos,
        ref_allele: bytes_to_string(parts[3]),
        alt_allele: bytes_to_string(parts[4]),
    };
    let id = bytes_to_string(parts[2]);
    let qual = bytes_to_string(parts[5]);
    let filter = bytes_to_string(parts[6]);
    let info = bytes_to_string(parts[7]);
    let format = bytes_to_string(parts[8]);
    let sample_blob = trim_crlf(parts[9]);
    let samples = if sample_blob.is_empty() {
        Vec::new()
    } else {
        split_tabs(sample_blob)
            .into_iter()
            .map(bytes_to_string)
            .collect()
    };
    Ok((chrom, key, id, qual, filter, info, format, samples))
}

#[derive(Debug, Clone)]
pub struct ParsedFormat {
    pub tags: Vec<String>,
    pub values: Vec<Vec<String>>,
}

pub fn parse_format(format: &str, sample_fields: &[String]) -> ParsedFormat {
    let tags: Vec<String> = if format.is_empty() {
        Vec::new()
    } else {
        format.split(':').map(|s| s.to_string()).collect()
    };
    let n_tags = tags.len();
    let values: Vec<Vec<String>> = sample_fields
        .iter()
        .map(|field| {
            if n_tags == 0 {
                return Vec::new();
            }
            let mut vals: Vec<String> = field.split(':').map(|s| s.to_string()).collect();
            if vals.len() < n_tags {
                vals.resize(n_tags, ".".to_string());
            } else if vals.len() > n_tags {
                vals.truncate(n_tags);
            }
            vals
        })
        .collect();
    ParsedFormat { tags, values }
}

pub fn get_tag_values(parsed: &ParsedFormat, tag: &str) -> Vec<String> {
    let Some(idx) = parsed.tags.iter().position(|t| t == tag) else {
        return vec![".".to_string(); parsed.values.len()];
    };
    parsed
        .values
        .iter()
        .map(|v| {
            v.get(idx)
                .cloned()
                .unwrap_or_else(|| ".".to_string())
        })
        .collect()
}

pub fn parse_info(info: &str) -> Vec<(String, String)> {
    if info.is_empty() || info == "." {
        return Vec::new();
    }
    info.split(';')
        .filter_map(|part| {
            if part.is_empty() {
                return None;
            }
            if let Some((k, v)) = part.split_once('=') {
                Some((k.to_string(), v.to_string()))
            } else {
                Some((part.to_string(), String::new()))
            }
        })
        .collect()
}

pub fn merge_info(base: &str, extras: &[(String, String)]) -> String {
    let mut items = parse_info(base);
    for (k, v) in extras {
        if v.is_empty() {
            if !items.iter().any(|(ek, _)| ek == k) {
                items.push((k.clone(), String::new()));
            }
        } else if !items.iter().any(|(ek, _)| ek == k) {
            items.push((k.clone(), v.clone()));
        }
    }
    if items.is_empty() {
        return ".".to_string();
    }
    items
        .into_iter()
        .map(|(k, v)| {
            if v.is_empty() {
                k
            } else {
                format!("{k}={v}")
            }
        })
        .collect::<Vec<_>>()
        .join(";")
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parse_format_roundtrip() {
        let pf = parse_format("GT:DP:GQ", &["0/1:10:30".to_string(), "./.:.:.".to_string()]);
        assert_eq!(get_tag_values(&pf, "DP"), vec!["10", "."]);
    }

    #[test]
    fn merge_info_adds_keys() {
        let out = merge_info("AF=0.1", &[("SRC".into(), "snv".into())]);
        assert!(out.contains("AF=0.1"));
        assert!(out.contains("SRC=snv"));
    }
}
