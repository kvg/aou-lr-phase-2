use crate::format::{pchisq_upper, round_digits, signif_digits};
use crate::null::{csc_matmul_cols, NullModel};

#[derive(Debug, Clone, Default)]
pub struct VariantStats {
    pub chi2: f64,
    pub p: f64,
    pub eff: Vec<f64>,
    pub se: Vec<f64>,
    pub pval: Vec<f64>,
    pub ac: Vec<i32>,
    pub include: Vec<i32>,
}

pub fn score_variant(
    null: &NullModel,
    genotypes: &[u8],
    n_anc: usize,
    ac_threshold: i32,
) -> VariantStats {
    let n = null.meta.n;
    let p = null.meta.p;
    let mut ac = vec![0i32; n_anc];
    for anc in 0..n_anc {
        let mut sum = 0i32;
        for i in 0..n {
            sum += genotypes[i + anc * n] as i32;
        }
        ac[anc] = sum;
    }

    let filter_mask: Vec<bool> = ac.iter().map(|&a| a > ac_threshold).collect();
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
        for i in 0..n {
            g[i + col * n] = genotypes[i + anc * n] as f64;
        }
        col += 1;
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
            var_score[a + b * k] = v;
        }
    }

    let inv = match invert_kxk(&var_score, k) {
        Some(m) => m,
        None => return VariantStats::default(),
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
    }
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
    // k x (2k) augmented matrix, column-major (stride k).
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
    use crate::null::CscMatrix;

    fn toy_null() -> NullModel {
        let n = 4;
        let p = 2;
        NullModel {
            meta: crate::null::NullMeta {
                n,
                p,
                nnz: n,
                family: "binomial".into(),
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
        }
    }

    fn identity_null(n: usize) -> NullModel {
        NullModel {
            meta: crate::null::NullMeta {
                n,
                p: 1,
                nnz: n,
                family: "binomial".into(),
                tractor_mix_score_sha: None,
            },
            id_include: (1..=n).map(|i| format!("S{i}")).collect(),
            sigma_i: CscMatrix {
                n,
                colptr: (0..=n).collect(),
                rowidx: (0..n).collect(),
                values: vec![1.0; n],
            },
            // Make the X*Sigma term zero so VarScore is just G'G in tests.
            sigma_i_x: vec![0.0; n],
            cov: vec![1.0],
            residuals: vec![1.0, -1.0, 0.5, -0.5][..n].to_vec(),
        }
    }

    #[test]
    fn invert_kxk_roundtrip() {
        let a = vec![2.0, 1.0, 1.0, 3.0];
        let inv = invert_kxk(&a, 2).unwrap();
        // A * inv ≈ I
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
        let geno = vec![0u8; 4 * 2];
        let stats = score_variant(&null, &geno, 2, 50);
        assert!(stats.ac.is_empty());
    }

    #[test]
    fn partial_ac_filter_sets_joint_to_na_and_include_mask() {
        let null = identity_null(4);
        // anc0: AC=2, anc1: AC=0, anc2: AC=2 (non-collinear kept ancestries)
        let geno = vec![
            // anc0
            0, 1, 0, 1, //
            // anc1
            0, 0, 0, 0, //
            // anc2
            1, 0, 1, 0,
        ];
        let stats = score_variant(&null, &geno, 3, 1);
        assert!(stats.p.is_nan(), "joint p should be NA when ancestries dropped");
        assert!(stats.chi2.is_nan(), "joint chi2 should be NA when ancestries dropped");
        assert_eq!(stats.include, vec![1, 0, 1], "expected anc1 dropped only");
        assert_eq!(stats.ac, vec![2, 0, 2]);
        assert!(!stats.eff[0].is_nan());
        assert!(stats.eff[1].is_nan());
        assert!(!stats.eff[2].is_nan());
    }

    #[test]
    fn all_ancestries_kept_computes_joint_stats() {
        let null = identity_null(4);
        // all AC > 1
        let geno = vec![
            // anc0
            0, 1, 0, 1, //
            // anc1
            0, 0, 1, 1, //
            // anc2
            1, 0, 1, 0,
        ];
        let stats = score_variant(&null, &geno, 3, 1);
        assert_eq!(stats.include, vec![1, 1, 1]);
        assert!(!stats.chi2.is_nan(), "joint chi2 should be present");
        assert!(!stats.p.is_nan(), "joint p should be present");
    }

    #[test]
    fn singular_varscore_returns_default_stats() {
        let null = identity_null(4);
        // anc0 and anc1 are identical -> singular VarScore for k=2.
        let geno = vec![
            // anc0
            1, 0, 1, 0, //
            // anc1 (identical)
            1, 0, 1, 0,
        ];
        let stats = score_variant(&null, &geno, 2, 0);
        assert!(stats.ac.is_empty(), "singular solve should map to all-NA row behavior");
    }
}
