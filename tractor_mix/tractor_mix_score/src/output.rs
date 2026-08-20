use std::fs::File;
use std::io::{BufWriter, Write};

use crate::dosage::VariantMeta;
use crate::error::Result;
use crate::score::VariantStats;

pub fn write_header(w: &mut impl Write, n_anc: usize) -> Result<()> {
    let mut cols = vec![
        "CHR".to_string(),
        "POS".to_string(),
        "ID".to_string(),
        "REF".to_string(),
        "ALT".to_string(),
        "Chi2".to_string(),
        "P".to_string(),
    ];
    for i in 0..n_anc {
        cols.push(format!("Eff_anc{i}"));
    }
    for i in 0..n_anc {
        cols.push(format!("SE_anc{i}"));
    }
    for i in 0..n_anc {
        cols.push(format!("Pval_anc{i}"));
    }
    for i in 0..n_anc {
        cols.push(format!("AC_count{i}"));
    }
    for i in 0..n_anc {
        cols.push(format!("include_anc{i}"));
    }
    writeln!(w, "{}", cols.join("\t"))?;
    Ok(())
}

fn fmt_num(x: f64) -> String {
    if x.is_nan() {
        "NA".to_string()
    } else {
        format!("{x}")
    }
}

pub fn write_variant_row(
    w: &mut impl Write,
    meta: &VariantMeta,
    stats: &VariantStats,
    n_anc: usize,
) -> Result<()> {
    let mut fields = vec![
        meta.chrom.clone(),
        meta.pos.clone(),
        meta.id.clone(),
        meta.r#ref.clone(),
        meta.alt.clone(),
    ];
    if stats.ac.is_empty() {
        let n_stat = 2 + 3 * n_anc + 2 * n_anc;
        for _ in 0..n_stat {
            fields.push("NA".to_string());
        }
    } else {
        fields.push(fmt_num(stats.chi2));
        fields.push(fmt_num(stats.p));
        for i in 0..n_anc {
            fields.push(fmt_num(stats.eff[i]));
        }
        for i in 0..n_anc {
            fields.push(fmt_num(stats.se[i]));
        }
        for i in 0..n_anc {
            fields.push(fmt_num(stats.pval[i]));
        }
        for i in 0..n_anc {
            fields.push(stats.ac[i].to_string());
        }
        for i in 0..n_anc {
            fields.push(stats.include[i].to_string());
        }
    }
    writeln!(w, "{}", fields.join("\t"))?;
    Ok(())
}

pub fn open_output(path: &std::path::Path) -> Result<BufWriter<File>> {
    Ok(BufWriter::with_capacity(1 << 20, File::create(path)?))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn header_has_expected_column_count() {
        let mut buf = Vec::<u8>::new();
        write_header(&mut buf, 5).unwrap();
        let line = String::from_utf8(buf).unwrap();
        let cols: Vec<_> = line.trim_end().split('\t').collect();
        assert_eq!(cols.len(), 7 + 5 * 5);
        assert_eq!(cols[0], "CHR");
        assert_eq!(cols[6], "P");
    }

    #[test]
    fn all_na_row_emits_expected_field_count() {
        let meta = VariantMeta {
            chrom: "22".into(),
            pos: "100".into(),
            id: "v1".into(),
            r#ref: "A".into(),
            alt: "G".into(),
        };
        let stats = VariantStats::default();
        let mut buf = Vec::<u8>::new();
        write_variant_row(&mut buf, &meta, &stats, 3).unwrap();
        let line = String::from_utf8(buf).unwrap();
        let cols: Vec<_> = line.trim_end().split('\t').collect();
        assert_eq!(cols.len(), 5 + (2 + 3 * 3 + 2 * 3));
        assert_eq!(cols[0], "22");
        assert_eq!(cols[5], "NA");
    }
}
