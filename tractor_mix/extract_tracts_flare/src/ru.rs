use memchr::memchr;

/// Repeat-unit dosage site parsed from VCF INFO.
#[derive(Debug, Clone, PartialEq)]
pub struct RuSite {
    pub cn_ref: f64,
    /// +1 for INS/DUP, -1 for DEL.
    pub sign: i8,
    /// Per-ALT repeat-unit counts (Number=A).
    pub ru: Vec<f64>,
}

pub fn has_ru_test(info: &[u8]) -> bool {
    for (key, val) in info_items(info) {
        if key == b"RU_TEST" {
            return match val {
                None | Some(b"") | Some(b"1") | Some(b"True") | Some(b"true") => true,
                Some(b"0") | Some(b"False") | Some(b"false") => false,
                Some(_) => true,
            };
        }
    }
    false
}

pub fn parse_ru_site(info: &[u8], alt: &[u8]) -> Option<RuSite> {
    if !has_ru_test(info) {
        return None;
    }
    let ru = parse_ru_list(info_get(info, b"RU")?)?;
    if ru.is_empty() {
        return None;
    }
    let cn_ref = info_get(info, b"CN_REF").and_then(parse_f64).unwrap_or(0.0);
    Some(RuSite {
        cn_ref,
        sign: svtype_sign(info, alt),
        ru,
    })
}

/// Copy-number contributed by one haplotype allele (`0` = REF, `1..` = ALT index).
pub fn haplotype_c(allele: usize, site: &RuSite) -> f64 {
    if allele == 0 {
        return site.cn_ref;
    }
    let extra = site.ru.get(allele - 1).copied().unwrap_or(0.0);
    site.cn_ref + f64::from(site.sign) * extra
}

pub fn count_alts(alt: &[u8]) -> usize {
    if alt.is_empty() {
        return 0;
    }
    1 + alt.iter().filter(|&&b| b == b',').count()
}

pub fn nth_alt(alt: &[u8], idx: usize) -> Option<&[u8]> {
    let mut start = 0usize;
    let mut n = 0usize;
    for i in 0..=alt.len() {
        if i == alt.len() || alt[i] == b',' {
            if n == idx {
                return Some(&alt[start..i]);
            }
            n += 1;
            start = i + 1;
        }
    }
    None
}

pub fn split_locus_id<'a>(
    vid: &'a [u8],
    chrom: &'a [u8],
    pos: &'a [u8],
    alt_idx1: usize,
) -> Vec<u8> {
    let mut out = Vec::with_capacity(vid.len() + 8);
    if vid == b"." || vid.is_empty() {
        out.extend_from_slice(chrom);
        out.push(b':');
        out.extend_from_slice(pos);
    } else {
        out.extend_from_slice(vid);
    }
    out.push(b':');
    out.extend_from_slice(alt_idx1.to_string().as_bytes());
    out
}

pub fn push_number(buf: &mut Vec<u8>, v: f64) {
    if !v.is_finite() {
        buf.push(b'0');
        return;
    }
    if (v - v.round()).abs() < 1e-9 && v.abs() < 1e15 {
        buf.extend_from_slice((v.round() as i64).to_string().as_bytes());
        return;
    }
    buf.extend_from_slice(v.to_string().as_bytes());
}

pub fn info_get<'a>(info: &'a [u8], key: &[u8]) -> Option<&'a [u8]> {
    for (k, v) in info_items(info) {
        if k == key {
            return v;
        }
    }
    None
}

fn parse_ru_list(raw: &[u8]) -> Option<Vec<f64>> {
    if raw.is_empty() {
        return None;
    }
    let mut out = Vec::new();
    let mut start = 0usize;
    for i in 0..=raw.len() {
        if i == raw.len() || raw[i] == b',' {
            let tok = &raw[start..i];
            out.push(parse_f64(tok)?);
            start = i + 1;
        }
    }
    Some(out)
}

fn parse_f64(s: &[u8]) -> Option<f64> {
    if s.is_empty() || s == b"." {
        return None;
    }
    std::str::from_utf8(s).ok()?.parse().ok()
}

fn svtype_sign(info: &[u8], alt: &[u8]) -> i8 {
    if let Some(v) = info_get(info, b"SVTYPE") {
        let end = v
            .iter()
            .position(|&b| b == b',' || b == b':')
            .unwrap_or(v.len());
        if v[..end].eq_ignore_ascii_case(b"DEL") {
            return -1;
        }
        return 1;
    }
    if alt.len() >= 4 && alt[..4].eq_ignore_ascii_case(b"<DEL") {
        return -1;
    }
    1
}

fn info_items(info: &[u8]) -> InfoIter<'_> {
    InfoIter { rest: info }
}

struct InfoIter<'a> {
    rest: &'a [u8],
}

impl<'a> Iterator for InfoIter<'a> {
    type Item = (&'a [u8], Option<&'a [u8]>);

    fn next(&mut self) -> Option<Self::Item> {
        if self.rest.is_empty() {
            return None;
        }
        let (item, rest) = match memchr(b';', self.rest) {
            Some(i) => (&self.rest[..i], &self.rest[i + 1..]),
            None => {
                let item = self.rest;
                self.rest = b"";
                (item, &b""[..])
            }
        };
        self.rest = rest;
        if item.is_empty() {
            return self.next();
        }
        Some(match memchr(b'=', item) {
            Some(eq) => (&item[..eq], Some(&item[eq + 1..])),
            None => (item, None),
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn flag_and_equals_forms() {
        assert!(has_ru_test(b"SVTYPE=INS;RU_TEST;RU=3"));
        assert!(has_ru_test(b"RU_TEST=1;RU=3"));
        assert!(!has_ru_test(b"RU_TEST=0;RU=3"));
        assert!(!has_ru_test(b"SVTYPE=INS;RU=3"));
    }

    #[test]
    fn c_ref_and_alt_ins_del() {
        let ins = parse_ru_site(b"RU_TEST;SVTYPE=INS;CN_REF=10;RU=3", b"ACAGCAGCAG").unwrap();
        assert_eq!(haplotype_c(0, &ins), 10.0);
        assert_eq!(haplotype_c(1, &ins), 13.0);
        let del = parse_ru_site(b"RU_TEST;SVTYPE=DEL;CN_REF=20;RU=5", b"<DEL>").unwrap();
        assert_eq!(haplotype_c(0, &del), 20.0);
        assert_eq!(haplotype_c(1, &del), 15.0);
        let rel = parse_ru_site(b"RU_TEST;SVTYPE=INS;RU=1", b"ACAG").unwrap();
        assert_eq!(haplotype_c(0, &rel), 0.0);
        assert_eq!(haplotype_c(1, &rel), 1.0);
    }

    #[test]
    fn multi_allelic_ru() {
        let site =
            parse_ru_site(b"RU_TEST;SVTYPE=INS;CN_REF=6;RU=1,3", b"ACAG,ACAGCAGCAG").unwrap();
        assert_eq!(haplotype_c(1, &site), 7.0);
        assert_eq!(haplotype_c(2, &site), 9.0);
        assert_eq!(count_alts(b"ACAG,ACAGCAGCAG"), 2);
        assert_eq!(nth_alt(b"ACAG,ACAGCAGCAG", 1).unwrap(), b"ACAGCAGCAG");
    }

    #[test]
    fn push_int_and_float() {
        let mut b = Vec::new();
        push_number(&mut b, 13.0);
        assert_eq!(b, b"13");
        b.clear();
        push_number(&mut b, -5.0);
        assert_eq!(b, b"-5");
        b.clear();
        push_number(&mut b, 10.5);
        assert_eq!(b, b"10.5");
    }
}
