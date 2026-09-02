use crate::format::{cct, pchisq_upper, round_digits, signif_digits};
use crate::null::{csc_matmul_cols, NullModel, ScoringKind};

const VAR_EPS: f64 = 1e-12;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ScoreMode {
    Auto,
    Legacy,
    Felix,
}

pub fn resolve_mode(requested: ScoreMode, null: &NullModel) -> ScoreMode {
    match requested {
        ScoreMode::Auto => match null.scoring {
            ScoringKind::Felix => ScoreMode::Felix,
            ScoringKind::Gmmat => ScoreMode::Legacy,
        },
        other => other,
    }
}

#[derive(Debug, Clone, Default)]
pub struct VariantStats {
    pub chi2: f64,
    pub p: f64,
    pub eff: Vec<f64>,
    pub se: Vec<f64>,
    pub pval: Vec<f64>,
    pub ac: Vec<f64>,
    pub include: Vec<i32>,
    pub af: Vec<f64>,
    pub tstat: Vec<f64>,
    pub var: Vec<f64>,
    pub p_het: f64,
    pub p_hom: f64,
    pub p_cct: f64,
    pub beta_all: f64,
    pub se_all: f64,
    pub tstat_all: f64,
    pub var_all: f64,
    pub p_all: f64,
    pub ac_all: f64,
    pub af_all: f64,
}

pub fn score_variant(
    null: &NullModel,
    genotypes: &[f64],
    n_anc: usize,
    ac_threshold: i32,
    mode: ScoreMode,
) -> VariantStats {
    match resolve_mode(mode, null) {
        ScoreMode::Felix => score_variant_felix(null, genotypes, n_anc),
        ScoreMode::Legacy => score_variant_legacy(null, genotypes, n_anc, ac_threshold),
        ScoreMode::Auto => unreachable!("resolve_mode never returns Auto"),
    }
}

fn col_ac(g: &[f64], n: usize, anc: usize) -> f64 {
    let start = anc * n;
    g[start..start + n].iter().sum()
}

fn col_var(g: &[f64], n: usize, anc: usize) -> f64 {
    let start = anc * n;
    let slice = &g[start..start + n];
    let mean = slice.iter().sum::<f64>() / n as f64;
    slice.iter().map(|x| (x - mean) * (x - mean)).sum::<f64>() / n as f64
}

fn score_variant_legacy(
    null: &NullModel,
    genotypes: &[f64],
    n_anc: usize,
    ac_threshold: i32,
) -> VariantStats {
    let n = null.meta.n;
    let mut ac = vec![0.0; n_anc];
    for anc in 0..n_anc {
        ac[anc] = col_ac(genotypes, n, anc);
    }

    let filter_mask: Vec<bool> = ac.iter().map(|&a| a > ac_threshold as f64).collect();
    if !filter_mask.iter().any(|&x| x) {
        return VariantStats::default();
    }

    let k = filter_mask.iter().filter(|&&x| x).count();
    let iffilter = k < n_anc;

    let mut g = vec![0.0; n * k];
    let mut col = 0usize;
    for anc in 0..n_anc {
        if !filter_mask[anc] {
            continue;
        }
        g[col * n..(col + 1) * n].copy_from_slice(&genotypes[anc * n..(anc + 1) * n]);
        col += 1;
    }

    let Some((score, inv)) = sandwich_inverse(null, &g, k, false) else {
        return VariantStats::default();
    };

    let mut joint_chi2 = f64::NAN;
    if !iffilter {
        joint_chi2 = quadratic_form(&score, &inv, k);
    }

    let mut anc_eff = vec![f64::NAN; n_anc];
    let mut anc_se = vec![f64::NAN; n_anc];
    let mut anc_pval = vec![f64::NAN; n_anc];

    for j in 0..k {
        let mut eff = 0.0;
        for i in 0..k {
            eff += inv[j + i * k] * score[i];
        }
        let se = inv[j + j * k].sqrt();
        let pval = pchisq_upper((eff / se).powi(2), 1);
        let anc_idx = filter_mask
            .iter()
            .enumerate()
            .filter(|(_, &keep)| keep)
            .nth(j)
            .map(|(i, _)| i)
            .unwrap();
        anc_eff[anc_idx] = eff;
        anc_se[anc_idx] = se;
        anc_pval[anc_idx] = pval;
    }

    let joint_p = if joint_chi2.is_finite() {
        pchisq_upper(joint_chi2, n_anc)
    } else {
        f64::NAN
    };

    let include: Vec<i32> = filter_mask.iter().map(|&b| i32::from(b)).collect();

    VariantStats {
        chi2: round_digits(joint_chi2, 5),
        p: signif_digits(joint_p, 5),
        eff: anc_eff.iter().map(|&x| round_digits(x, 5)).collect(),
        se: anc_se.iter().map(|&x| round_digits(x, 5)).collect(),
        pval: anc_pval.iter().map(|&x| signif_digits(x, 5)).collect(),
        ac,
        include,
        ..VariantStats::default()
    }
}

fn score_variant_felix(null: &NullModel, genotypes: &[f64], n_anc: usize) -> VariantStats {
    let n = null.meta.n;
    let mut ac = vec![0.0; n_anc];
    let mut af = vec![0.0; n_anc];
    let mut include = vec![0i32; n_anc];
    for anc in 0..n_anc {
        ac[anc] = col_ac(genotypes, n, anc);
        af[anc] = ac[anc] / (2.0 * n as f64);
        if col_var(genotypes, n, anc) > VAR_EPS {
            include[anc] = 1;
        }
    }

    let kept: Vec<usize> = (0..n_anc).filter(|&a| include[a] == 1).collect();
    let ac_all: f64 = ac.iter().sum();
    let af_all = ac_all / (2.0 * n as f64);

    let mut g_sum = vec![0.0; n];
    for anc in 0..n_anc {
        for i in 0..n {
            g_sum[i] += genotypes[i + anc * n];
        }
    }

    let (beta_all, se_all, tstat_all, var_all, p_all) = single_df_test(null, &g_sum);
    let p_hom = p_all;

    let mut anc_eff = vec![f64::NAN; n_anc];
    let mut anc_se = vec![f64::NAN; n_anc];
    let mut anc_pval = vec![f64::NAN; n_anc];
    let mut anc_tstat = vec![f64::NAN; n_anc];
    let mut anc_var = vec![f64::NAN; n_anc];

    for &anc in &kept {
        let col = &genotypes[anc * n..(anc + 1) * n];
        let (b, se, t, v, p) = single_df_test(null, col);
        anc_eff[anc] = b;
        anc_se[anc] = se;
        anc_tstat[anc] = t;
        anc_var[anc] = v;
        anc_pval[anc] = p;
    }

    let p_het = if kept.is_empty() {
        f64::NAN
    } else if kept.len() == 1 {
        anc_pval[kept[0]]
    } else {
        let k = kept.len();
        let mut g = vec![0.0; n * k];
        for (j, &anc) in kept.iter().enumerate() {
            g[j * n..(j + 1) * n].copy_from_slice(&genotypes[anc * n..(anc + 1) * n]);
        }
        match sandwich_inverse(null, &g, k, true) {
            Some((score, inv)) => {
                let chi2 = quadratic_form(&score, &inv, k);
                pchisq_upper(chi2, k)
            }
            None => f64::NAN,
        }
    };

    let p_cct = combine_cct(p_het, p_hom);

    VariantStats {
        chi2: f64::NAN,
        p: f64::NAN,
        eff: anc_eff,
        se: anc_se,
        pval: anc_pval,
        ac,
        include,
        af,
        tstat: anc_tstat,
        var: anc_var,
        p_het,
        p_hom,
        p_cct,
        beta_all,
        se_all,
        tstat_all,
        var_all,
        p_all,
        ac_all,
        af_all,
    }
}

fn combine_cct(p_het: f64, p_hom: f64) -> f64 {
    let p = cct(&[p_het, p_hom]);
    if p == 0.0 {
        if p_hom == 0.0 {
            p_hom
        } else {
            p_het
        }
    } else {
        p
    }
}

fn single_df_test(null: &NullModel, g: &[f64]) -> (f64, f64, f64, f64, f64) {
    match sandwich_inverse(null, g, 1, true) {
        Some((score, inv)) => {
            let s = score[0];
            let var1 = if inv[0].abs() > 0.0 {
                1.0 / inv[0]
            } else {
                return (f64::NAN, f64::NAN, f64::NAN, f64::NAN, f64::NAN);
            };
            if !(var1.is_finite() && var1 > 0.0) {
                return (f64::NAN, f64::NAN, s, var1, 1.0);
            }
            let beta = s / var1;
            let se = 1.0 / var1.sqrt();
            let p = pchisq_upper(s * s / var1, 1);
            (beta, se, s, var1, p)
        }
        None => (f64::NAN, f64::NAN, f64::NAN, f64::NAN, f64::NAN),
    }
}

fn project_columns(null: &NullModel, g: &mut [f64], k: usize) {
    let n = null.meta.n;
    let p = null.meta.p;
    let (Some(xv), Some(xxvx)) = (&null.xv, &null.xxvx_inv) else {
        return;
    };
    for j in 0..k {
        let col = j * n;
        let mut z = vec![0.0; p];
        for i in 0..n {
            let gi = g[col + i];
            if gi == 0.0 {
                continue;
            }
            for r in 0..p {
                z[r] += xv[r + i * p] * gi;
            }
        }
        for i in 0..n {
            let mut s = g[col + i];
            for r in 0..p {
                s -= xxvx[i + r * n] * z[r];
            }
            g[col + i] = s;
        }
    }
}

/// Score = G' residuals; Var = varRatio * (G' Sigma_i G - sandwich).
/// When `project` is true and XV/XXVX_inv are present, G is replaced by getadjG.
fn sandwich_inverse(
    null: &NullModel,
    g_in: &[f64],
    k: usize,
    project: bool,
) -> Option<(Vec<f64>, Vec<f64>)> {
    if k == 0 {
        return None;
    }
    let n = null.meta.n;
    let p = null.meta.p;
    let mut g = g_in.to_vec();
    if project {
        project_columns(null, &mut g, k);
    }

    let mut sigma_g = vec![0.0; n * k];
    csc_matmul_cols(&null.sigma_i, &g, k, &mut sigma_g);

    let mut score = vec![0.0; k];
    for j in 0..k {
        let mut s = 0.0;
        for i in 0..n {
            s += g[i + j * n] * null.residuals[i];
        }
        score[j] = s;
    }

    let mut xsigma_g = vec![0.0; p * k];
    for j in 0..k {
        for row in 0..p {
            let mut s = 0.0;
            for i in 0..n {
                s += null.sigma_i_x[i + row * n] * g[i + j * n];
            }
            xsigma_g[row + j * p] = s;
        }
    }

    let mut var_score = vec![0.0; k * k];
    for a in 0..k {
        for b in 0..k {
            let mut v = 0.0;
            for i in 0..n {
                v += g[i + a * n] * sigma_g[i + b * n];
            }
            for r in 0..p {
                for c in 0..p {
                    v -= xsigma_g[r + a * p] * null.cov[r + c * p] * xsigma_g[c + b * p];
                }
            }
            var_score[a + b * k] = v * null.variance_ratio;
        }
    }

    let inv = invert_kxk(&var_score, k)?;
    Some((score, inv))
}

fn quadratic_form(v: &[f64], inv: &[f64], k: usize) -> f64 {
    let mut tmp = vec![0.0; k];
    for i in 0..k {
        let mut s = 0.0;
        for j in 0..k {
            s += inv[i + j * k] * v[j];
        }
        tmp[i] = s;
    }
    let mut q = 0.0;
    for i in 0..k {
        q += v[i] * tmp[i];
    }
    q
}

fn invert_kxk(a: &[f64], k: usize) -> Option<Vec<f64>> {
    if k == 0 {
        return None;
    }
    let mut aug = vec![0.0; k * 2 * k];
    for i in 0..k {
        for j in 0..k {
            aug[i + j * k] = a[i + j * k];
        }
        aug[i + (k + i) * k] = 1.0;
    }
    for col in 0..k {
        let mut pivot = col;
        let mut best = aug[pivot + col * k].abs();
        for row in (col + 1)..k {
            let v = aug[row + col * k].abs();
            if v > best {
                best = v;
                pivot = row;
            }
        }
        if best < 1e-12 {
            return None;
        }
        if pivot != col {
            for c in 0..(2 * k) {
                aug.swap(col + c * k, pivot + c * k);
            }
        }
        let div = aug[col + col * k];
        for c in 0..(2 * k) {
            aug[col + c * k] /= div;
        }
        for row in 0..k {
            if row == col {
                continue;
            }
            let factor = aug[row + col * k];
            if factor == 0.0 {
                continue;
            }
            for c in 0..(2 * k) {
                aug[row + c * k] -= factor * aug[col + c * k];
            }
        }
    }
    let mut inv = vec![0.0; k * k];
    for i in 0..k {
        for j in 0..k {
            inv[i + j * k] = aug[i + (k + j) * k];
        }
    }
    Some(inv)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::null::{CscMatrix, NullMeta, ScoringKind};

    fn toy_null() -> NullModel {
        let n = 4;
        let p = 2;
        NullModel {
            meta: NullMeta {
                n,
                p,
                nnz: n,
                family: "binomial".into(),
                source: None,
                variance_ratio: None,
                tau0: None,
                trait_type: None,
                tractor_mix_score_sha: None,
            },
            id_include: (1..=n).map(|i| format!("S{i}")).collect(),
            sigma_i: CscMatrix {
                n,
                colptr: (0..=n).map(|i| i).collect(),
                rowidx: (0..n).collect(),
                values: vec![1.0; n],
            },
            sigma_i_x: vec![0.1; n * p],
            cov: vec![1.0, 0.0, 0.0, 1.0],
            residuals: vec![0.5, -0.5, 0.25, -0.25],
            variance_ratio: 1.0,
            xv: None,
            xxvx_inv: None,
            scoring: ScoringKind::Gmmat,
        }
    }

    fn identity_null(n: usize) -> NullModel {
        NullModel {
            meta: NullMeta {
                n,
                p: 1,
                nnz: n,
                family: "binomial".into(),
                source: None,
                variance_ratio: None,
                tau0: None,
                trait_type: None,
                tractor_mix_score_sha: None,
            },
            id_include: (1..=n).map(|i| format!("S{i}")).collect(),
            sigma_i: CscMatrix {
                n,
                colptr: (0..=n).collect(),
                rowidx: (0..n).collect(),
                values: vec![1.0; n],
            },
            sigma_i_x: vec![0.0; n],
            cov: vec![1.0],
            residuals: vec![1.0, -1.0, 0.5, -0.5][..n].to_vec(),
            variance_ratio: 1.0,
            xv: None,
            xxvx_inv: None,
            scoring: ScoringKind::Gmmat,
        }
    }

    fn felix_identity_null(n: usize, vr: f64) -> NullModel {
        let mut m = identity_null(n);
        m.meta.family = "gaussian".into();
        m.meta.source = Some("felix".into());
        m.variance_ratio = vr;
        m.scoring = ScoringKind::Felix;
        m
    }

    #[test]
    fn invert_kxk_roundtrip() {
        let a = vec![2.0, 1.0, 1.0, 3.0];
        let inv = invert_kxk(&a, 2).unwrap();
        for i in 0..2 {
            for j in 0..2 {
                let mut s = 0.0;
                for t in 0..2 {
                    s += a[i + t * 2] * inv[t + j * 2];
                }
                let expect = if i == j { 1.0 } else { 0.0 };
                assert!((s - expect).abs() < 1e-10, "({i},{j})={s}");
            }
        }
    }

    #[test]
    fn ac_filter_all_fail_returns_default() {
        let null = toy_null();
        let geno = vec![0.0; 4 * 2];
        let stats = score_variant(&null, &geno, 2, 50, ScoreMode::Legacy);
        assert!(stats.ac.is_empty());
    }

    #[test]
    fn partial_ac_filter_sets_joint_to_na_and_include_mask() {
        let null = identity_null(4);
        let geno = vec![
            0.0, 1.0, 0.0, 1.0, //
            0.0, 0.0, 0.0, 0.0, //
            1.0, 0.0, 1.0, 0.0,
        ];
        let stats = score_variant(&null, &geno, 3, 1, ScoreMode::Legacy);
        assert!(stats.p.is_nan(), "joint p should be NA when ancestries dropped");
        assert!(stats.chi2.is_nan(), "joint chi2 should be NA when ancestries dropped");
        assert_eq!(stats.include, vec![1, 0, 1], "expected anc1 dropped only");
        assert_eq!(stats.ac, vec![2.0, 0.0, 2.0]);
        assert!(!stats.eff[0].is_nan());
        assert!(stats.eff[1].is_nan());
        assert!(!stats.eff[2].is_nan());
    }

    #[test]
    fn all_ancestries_kept_computes_joint_stats() {
        let null = identity_null(4);
        let geno = vec![
            0.0, 1.0, 0.0, 1.0, //
            0.0, 0.0, 1.0, 1.0, //
            1.0, 0.0, 1.0, 0.0,
        ];
        let stats = score_variant(&null, &geno, 3, 1, ScoreMode::Legacy);
        assert_eq!(stats.include, vec![1, 1, 1]);
        assert!(!stats.chi2.is_nan(), "joint chi2 should be present");
        assert!(!stats.p.is_nan(), "joint p should be present");
    }

    #[test]
    fn singular_varscore_returns_default_stats() {
        let null = identity_null(4);
        let geno = vec![
            1.0, 0.0, 1.0, 0.0, //
            1.0, 0.0, 1.0, 0.0,
        ];
        let stats = score_variant(&null, &geno, 2, 0, ScoreMode::Legacy);
        assert!(stats.ac.is_empty(), "singular solve should map to all-NA row behavior");
    }

    #[test]
    fn felix_drops_zero_variance_ancestry() {
        let null = felix_identity_null(4, 1.0);
        let geno = vec![
            0.0, 1.0, 0.0, 1.0, //
            0.0, 0.0, 0.0, 0.0, //
            1.0, 0.0, 1.0, 0.0,
        ];
        let stats = score_variant(&null, &geno, 3, 50, ScoreMode::Felix);
        assert_eq!(stats.include, vec![1, 0, 1]);
        assert!(stats.pval[1].is_nan());
        assert!(!stats.pval[0].is_nan());
        assert!(stats.p_het.is_finite());
        assert!(stats.p_hom.is_finite());
        assert!(stats.p_cct.is_finite());
        assert_eq!(stats.ac_all, 4.0);
    }

    #[test]
    fn felix_homogeneous_is_sum_of_ancestries() {
        let null = felix_identity_null(4, 1.0);
        let g0 = vec![1.0, 0.0, 1.0, 0.0];
        let g1 = vec![0.0, 1.0, 0.0, 1.0];
        let mut geno = Vec::new();
        geno.extend_from_slice(&g0);
        geno.extend_from_slice(&g1);
        let stats = score_variant(&null, &geno, 2, 0, ScoreMode::Felix);
        let gsum: Vec<f64> = (0..4).map(|i| g0[i] + g1[i]).collect();
        let (_, _, _, _, p_sum) = single_df_test(&null, &gsum);
        assert!((stats.p_hom - p_sum).abs() < 1e-12);
    }

    #[test]
    fn felix_variance_ratio_scales_var() {
        let g = vec![0.0, 1.0, 0.0, 1.0];
        let n1 = felix_identity_null(4, 1.0);
        let n4 = felix_identity_null(4, 4.0);
        let (_, _, _, v1, _) = single_df_test(&n1, &g);
        let (_, _, _, v4, _) = single_df_test(&n4, &g);
        assert!((v4 - 4.0 * v1).abs() < 1e-10, "v1={v1} v4={v4}");
    }

    #[test]
    fn felix_f64_dosages_are_used() {
        let null = felix_identity_null(4, 1.0);
        let geno = vec![0.5, 1.25, 0.0, 0.75, 0.1, 0.2, 0.3, 0.4];
        let stats = score_variant(&null, &geno, 2, 50, ScoreMode::Felix);
        assert!((stats.ac[0] - 2.5).abs() < 1e-12);
        assert!((stats.ac[1] - 1.0).abs() < 1e-12);
        assert!(stats.p_hom.is_finite());
    }

    #[test]
    fn getadjg_projects_out_x() {
        let n = 4;
        let p = 1;
        let mut null = identity_null(n);
        // XV = 1 x n of ones; XXVX_inv = n x 1 of 0.25 so gtilde = G - mean(G)
        null.xv = Some(vec![1.0; n]);
        null.xxvx_inv = Some(vec![0.25; n]);
        let mut g = vec![1.0, 3.0, 5.0, 7.0];
        project_columns(&null, &mut g, 1);
        let mean = 4.0;
        for (i, &v) in g.iter().enumerate() {
            assert!((v - ([1.0, 3.0, 5.0, 7.0][i] - mean)).abs() < 1e-12);
        }
        let _ = p;
    }
}
