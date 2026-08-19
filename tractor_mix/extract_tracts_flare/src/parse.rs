use memchr::memchr;

/// Python `str.split("\\t", 9)`: at most 10 parts (9 splits).
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

/// Python `str.strip("# ")` — strip `#` and space from both ends, keep newlines.
pub fn strip_hash_space(s: &[u8]) -> &[u8] {
    let mut start = 0usize;
    let mut end = s.len();
    while start < end && (s[start] == b'#' || s[start] == b' ') {
        start += 1;
    }
    while end > start && (s[end - 1] == b'#' || s[end - 1] == b' ') {
        end -= 1;
    }
    &s[start..end]
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

/// Python `str.strip()` whitespace (ASCII) used on data lines.
pub fn strip_ascii_ws(s: &[u8]) -> &[u8] {
    let mut start = 0usize;
    let mut end = s.len();
    while start < end && s[start].is_ascii_whitespace() {
        start += 1;
    }
    while end > start && s[end - 1].is_ascii_whitespace() {
        end -= 1;
    }
    &s[start..end]
}

/// Split on tabs without allocating strings. Empty input → one empty field.
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

/// First four tokens of a FLARE sample field split on `|` or `:`.
/// Extra tokens (ANP1/ANP2, …) are ignored, matching `re.split('[|:]', geno)[:4]`.
pub fn first_four_tokens(field: &[u8]) -> Option<[&[u8]; 4]> {
    let mut toks = [&b""[..]; 4];
    let mut start = 0usize;
    let mut n = 0usize;
    for i in 0..field.len() {
        if field[i] == b'|' || field[i] == b':' {
            if n == 4 {
                return Some(toks);
            }
            toks[n] = &field[start..i];
            n += 1;
            start = i + 1;
            if n == 4 {
                return Some(toks);
            }
        }
    }
    if n < 4 {
        toks[n] = &field[start..];
        n += 1;
    }
    if n == 4 {
        Some(toks)
    } else {
        None
    }
}

pub fn parse_nonneg_int(s: &[u8]) -> Option<usize> {
    if s.is_empty() {
        return None;
    }
    let mut n = 0usize;
    for &b in s {
        if !b.is_ascii_digit() {
            return None;
        }
        n = n.checked_mul(10)?.checked_add((b - b'0') as usize)?;
    }
    Some(n)
}

/// Dosage (alt `"1"` copies on this ancestry) and hapcount for one sample.
/// Ancestry tokens are matched as integers against `0..num_ancs`, equivalent to
/// Python `call == str(j)` for typical FLARE IDs (no leading zeros).
pub fn sample_counts(
    geno_a: &[u8],
    geno_b: &[u8],
    call_a: &[u8],
    call_b: &[u8],
    num_ancs: usize,
) -> (Vec<u8>, Vec<u8>) {
    let mut dosage = vec![0u8; num_ancs];
    let mut hapcount = vec![0u8; num_ancs];
    apply_haplotype(geno_a, call_a, num_ancs, &mut dosage, &mut hapcount);
    apply_haplotype(geno_b, call_b, num_ancs, &mut dosage, &mut hapcount);
    (dosage, hapcount)
}

pub fn apply_haplotype(
    allele: &[u8],
    ancestry: &[u8],
    num_ancs: usize,
    dosage: &mut [u8],
    hapcount: &mut [u8],
) {
    if let Some(j) = parse_nonneg_int(ancestry) {
        if j < num_ancs {
            hapcount[j] = hapcount[j].saturating_add(1);
            if allele == b"1" {
                dosage[j] = dosage[j].saturating_add(1);
            }
        }
    }
}

pub fn write_ancestry_gt(
    buf: &mut Vec<u8>,
    geno_a: &[u8],
    geno_b: &[u8],
    call_a: &[u8],
    call_b: &[u8],
    anc: usize,
) {
    if parse_nonneg_int(call_a) == Some(anc) {
        buf.extend_from_slice(geno_a);
    } else {
        buf.push(b'.');
    }
    buf.push(b'|');
    if parse_nonneg_int(call_b) == Some(anc) {
        buf.extend_from_slice(geno_b);
    } else {
        buf.push(b'.');
    }
}

pub fn parse_vcf_filename(filename: &str) -> Result<(String, bool), String> {
    if filename.ends_with(".vcf.gz") {
        let stem = filename.trim_end_matches(".vcf.gz");
        Ok((stem.to_string(), true))
    } else if filename.ends_with(".vcf") {
        let stem = filename.trim_end_matches(".vcf");
        Ok((stem.to_string(), false))
    } else {
        Err(filename.to_string())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn four_tokens_basic_flare() {
        let t = first_four_tokens(b"0|1:2:3").unwrap();
        assert_eq!(t, [b"0".as_slice(), b"1", b"2", b"3"]);
    }

    #[test]
    fn four_tokens_ignores_anp() {
        let t = first_four_tokens(b"0|1:2:3:0.1,0.8,0.1:0.2,0.7,0.1").unwrap();
        assert_eq!(t, [b"0".as_slice(), b"1", b"2", b"3"]);
    }

    #[test]
    fn four_tokens_unphased_too_few() {
        assert!(first_four_tokens(b"0/1:2:3").is_none());
    }

    #[test]
    fn four_tokens_truncated() {
        assert!(first_four_tokens(b"0|1:2").is_none());
    }

    #[test]
    fn counts_mixed_ancestry_het() {
        // 0|1:0:1 → dos0=1, dos1=1; hap 1/1
        let (d, h) = sample_counts(b"0", b"1", b"0", b"1", 2);
        assert_eq!(d, vec![0, 1]);
        assert_eq!(h, vec![1, 1]);
    }

    #[test]
    fn counts_both_alt_same_anc() {
        let (d, h) = sample_counts(b"1", b"1", b"1", b"1", 3);
        assert_eq!(d, vec![0, 2, 0]);
        assert_eq!(h, vec![0, 2, 0]);
    }

    #[test]
    fn counts_allele_two_not_dosage() {
        let (d, h) = sample_counts(b"2", b"1", b"0", b"0", 2);
        assert_eq!(d, vec![1, 0]);
        assert_eq!(h, vec![2, 0]);
    }

    #[test]
    fn counts_missing_and_oor_ancestry() {
        let (d, h) = sample_counts(b"1", b"1", b".", b"5", 5);
        assert_eq!(d, vec![0, 0, 0, 0, 0]);
        assert_eq!(h, vec![0, 0, 0, 0, 0]);
    }

    #[test]
    fn strip_hash_keeps_newline() {
        let s = strip_hash_space(b"#CHROM\tPOS\n");
        assert_eq!(s, b"CHROM\tPOS\n");
    }

    #[test]
    fn split_tab_max9_header() {
        let line = b"CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\n";
        let p = split_tab_max9(line).unwrap();
        assert_eq!(p[0], b"CHROM");
        assert_eq!(p[4], b"ALT");
        assert_eq!(p[9], b"S1\tS2\n");
    }

    #[test]
    fn ancestry_gt_layout() {
        let mut buf = Vec::new();
        write_ancestry_gt(&mut buf, b"0", b"1", b"0", b"1", 0);
        assert_eq!(buf, b"0|.");
        buf.clear();
        write_ancestry_gt(&mut buf, b"0", b"1", b"0", b"1", 1);
        assert_eq!(buf, b".|1");
    }
}
