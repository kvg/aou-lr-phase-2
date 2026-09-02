use std::fs::File;
use std::io::{BufWriter, Write};

use crate::dosage::VariantMeta;
use crate::error::Result;
use crate::score::{ScoreMode, VariantStats};

pub fn write_header(w: &mut impl Write, n_anc: usize, mode: ScoreMode) -> Result<()> {
    match mode {
        ScoreMode::Felix => write_felix_header(w, n_anc),
        ScoreMode::Legacy | ScoreMode::Auto => write_legacy_header(w, n_anc),
    }
}

fn write_legacy_header(w: &mut impl Write, n_anc: usize) -> Result<()> {
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

fn write_felix_header(w: &mut impl Write, n_anc: usize) -> Result<()> {
    let mut cols = vec![
        "CHR".to_string(),
        "POS".to_string(),
        "MarkerID".to_string(),
        "Allele1".to_string(),
        "Allele2".to_string(),
    ];
    for i in 1..=n_anc {
        cols.push(format!("AC_Allele2_anc{i}"));
        cols.push(format!("AF_Allele2_anc{i}"));
        cols.push(format!("BETA_anc{i}"));
        cols.push(format!("SE_anc{i}"));
        cols.push(format!("Tstat_anc{i}"));
        cols.push(format!("var_anc{i}"));
        cols.push(format!("p.value_anc{i}"));
    }
    cols.extend([
        "AC_Allele2_ancALL".into(),
        "AF_Allele2_ancALL".into(),
        "BETA_ancALL".into(),
        "SE_ancALL".into(),
        "Tstat_ancALL".into(),
        "var_ancALL".into(),
        "p.value_ancALL".into(),
        "P_het_admixed".into(),
        "P_hom_admixed".into(),
        "P_cct_admixed".into(),
        "P_het_admixed_c".into(),
        "P_hom_admixed_c".into(),
        "P_cct_admixed_c".into(),
    ]);
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

fn fmt_ac_legacy(x: f64) -> String {
    if !x.is_finite() {
        "NA".to_string()
    } else if (x - x.round()).abs() < 1e-9 {
        format!("{}", x.round() as i64)
    } else {
        format!("{x}")
    }
}

pub fn write_variant_row(
    w: &mut impl Write,
    meta: &VariantMeta,
    stats: &VariantStats,
    n_anc: usize,
    mode: ScoreMode,
) -> Result<()> {
    match mode {
        ScoreMode::Felix => write_felix_row(w, meta, stats, n_anc),
        ScoreMode::Legacy | ScoreMode::Auto => write_legacy_row(w, meta, stats, n_anc),
    }
}

fn write_legacy_row(
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
            fields.push(fmt_ac_legacy(stats.ac[i]));
        }
        for i in 0..n_anc {
            fields.push(stats.include[i].to_string());
        }
    }
    writeln!(w, "{}", fields.join("\t"))?;
    Ok(())
}

fn write_felix_row(
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
        let n_stat = 7 * (n_anc + 1) + 6;
        for _ in 0..n_stat {
            fields.push("NA".to_string());
        }
    } else {
        for i in 0..n_anc {
            fields.push(fmt_num(stats.ac[i]));
            fields.push(fmt_num(stats.af[i]));
            fields.push(fmt_num(stats.eff[i]));
            fields.push(fmt_num(stats.se[i]));
            fields.push(fmt_num(stats.tstat[i]));
            fields.push(fmt_num(stats.var[i]));
            fields.push(fmt_num(stats.pval[i]));
        }
        fields.push(fmt_num(stats.ac_all));
        fields.push(fmt_num(stats.af_all));
        fields.push(fmt_num(stats.beta_all));
        fields.push(fmt_num(stats.se_all));
        fields.push(fmt_num(stats.tstat_all));
        fields.push(fmt_num(stats.var_all));
        fields.push(fmt_num(stats.p_all));
        fields.push(fmt_num(stats.p_het));
        fields.push(fmt_num(stats.p_hom));
        fields.push(fmt_num(stats.p_cct));
        fields.push(fmt_num(stats.p_het));
        fields.push(fmt_num(stats.p_hom));
        fields.push(fmt_num(stats.p_cct));
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
        write_header(&mut buf, 5, ScoreMode::Legacy).unwrap();
        let line = String::from_utf8(buf).unwrap();
        let cols: Vec<_> = line.trim_end().split('\t').collect();
        assert_eq!(cols.len(), 7 + 5 * 5);
        assert_eq!(cols[0], "CHR");
        assert_eq!(cols[6], "P");
    }

    #[test]
    fn felix_header_has_cct_columns() {
        let mut buf = Vec::<u8>::new();
        write_header(&mut buf, 2, ScoreMode::Felix).unwrap();
        let line = String::from_utf8(buf).unwrap();
        let cols: Vec<_> = line.trim_end().split('\t').collect();
        assert_eq!(cols.len(), 5 + 7 * 3 + 6);
        assert_eq!(cols[2], "MarkerID");
        assert!(cols.contains(&"P_cct_admixed"));
        assert!(cols.contains(&"P_cct_admixed_c"));
        assert!(cols.contains(&"BETA_ancALL"));
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
        write_variant_row(&mut buf, &meta, &stats, 3, ScoreMode::Legacy).unwrap();
        let line = String::from_utf8(buf).unwrap();
        let cols: Vec<_> = line.trim_end().split('\t').collect();
        assert_eq!(cols.len(), 5 + (2 + 3 * 3 + 2 * 3));
        assert_eq!(cols[0], "22");
        assert_eq!(cols[5], "NA");
    }
}
