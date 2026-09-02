use std::fs::File;
use std::io::{Read, Write};
use std::path::Path;

use serde::Deserialize;

use crate::error::{Result, ScoreError};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ScoringKind {
    Gmmat,
    Felix,
}

#[derive(Debug, Clone, Deserialize)]
pub struct NullMeta {
    pub n: usize,
    pub p: usize,
    pub nnz: usize,
    pub family: String,
    #[serde(default)]
    pub source: Option<String>,
    #[serde(default)]
    pub variance_ratio: Option<f64>,
    #[serde(default)]
    pub tau0: Option<f64>,
    #[serde(default)]
    pub trait_type: Option<String>,
    pub tractor_mix_score_sha: Option<String>,
}

impl NullMeta {
    pub fn is_quantitative(&self) -> bool {
        let fam = self.family.to_ascii_lowercase();
        fam == "gaussian" || fam == "quantitative" || fam == "gaussian()"
    }

    pub fn is_binomial(&self) -> bool {
        let fam = self.family.to_ascii_lowercase();
        fam == "binomial" || fam == "binary"
    }

    pub fn scoring_kind(&self) -> ScoringKind {
        let src = self
            .source
            .as_deref()
            .unwrap_or("")
            .to_ascii_lowercase();
        if src == "felix" || src == "saige" || self.is_quantitative() {
            ScoringKind::Felix
        } else {
            ScoringKind::Gmmat
        }
    }
}

#[derive(Debug, Clone)]
pub struct NullModel {
    pub meta: NullMeta,
    pub id_include: Vec<String>,
    pub sigma_i: CscMatrix,
    pub sigma_i_x: Vec<f64>,
    pub cov: Vec<f64>,
    pub residuals: Vec<f64>,
    pub variance_ratio: f64,
    /// Optional p x n column-major XV for FELIX getadjG.
    pub xv: Option<Vec<f64>>,
    /// Optional n x p column-major XXVX_inv for FELIX getadjG.
    pub xxvx_inv: Option<Vec<f64>>,
    pub scoring: ScoringKind,
}

impl NullModel {
    pub fn has_getadj(&self) -> bool {
        self.xv.is_some() && self.xxvx_inv.is_some()
    }
}

#[derive(Debug, Clone)]
pub struct CscMatrix {
    pub n: usize,
    pub colptr: Vec<usize>,
    pub rowidx: Vec<usize>,
    pub values: Vec<f64>,
}

fn read_i32_le(r: &mut impl Read) -> Result<i32> {
    let mut buf = [0u8; 4];
    r.read_exact(&mut buf)?;
    Ok(i32::from_le_bytes(buf))
}

fn read_f64_le(r: &mut impl Read) -> Result<f64> {
    let mut buf = [0u8; 8];
    r.read_exact(&mut buf)?;
    Ok(f64::from_le_bytes(buf))
}

pub fn parse_variance_ratio_file_public(path: &Path) -> Option<f64> {
    parse_variance_ratio_file(path)
}

fn parse_variance_ratio_file(path: &Path) -> Option<f64> {
    let text = std::fs::read_to_string(path).ok()?;
    for line in text.lines() {
        let t = line.trim();
        if t.is_empty() {
            continue;
        }
        let first = t.split_whitespace().next()?;
        if let Ok(v) = first.parse::<f64>() {
            if v.is_finite() && v > 0.0 {
                return Some(v);
            }
        }
    }
    None
}

pub fn load_null_export(dir: &Path) -> Result<NullModel> {
    let meta: NullMeta = serde_json::from_reader(File::open(dir.join("meta.json"))?)?;
    if !meta.is_binomial() && !meta.is_quantitative() {
        return Err(ScoreError::msg(format!(
            "unsupported family {} (binomial or quantitative/gaussian)",
            meta.family
        )));
    }

    let id_include = std::fs::read_to_string(dir.join("id_include.txt"))?
        .lines()
        .map(str::trim)
        .filter(|s| !s.is_empty())
        .map(String::from)
        .collect::<Vec<_>>();
    if id_include.len() != meta.n {
        return Err(ScoreError::msg(format!(
            "id_include count {} != meta.n {}",
            id_include.len(),
            meta.n
        )));
    }

    let sigma_i = read_sigma_i_csc(&dir.join("sigma_i.csc.bin"))?;
    let sigma_i_x = read_dense_np(&dir.join("sigma_i_x.bin"), meta.n, meta.p)?;
    let cov = read_dense_p(&dir.join("cov.bin"), meta.p)?;
    let residuals = read_residuals(&dir.join("residuals.bin"), meta.n)?;

    let mut variance_ratio = meta.variance_ratio.unwrap_or(1.0);
    if let Some(vr) = parse_variance_ratio_file(&dir.join("variance_ratio.txt")) {
        if meta.variance_ratio.is_none() {
            variance_ratio = vr;
        }
    }
    if !(variance_ratio.is_finite() && variance_ratio > 0.0) {
        return Err(ScoreError::msg(format!(
            "invalid variance_ratio {variance_ratio}"
        )));
    }

    let xv_path = dir.join("xv.bin");
    let xx_path = dir.join("xxvx_inv.bin");
    let (xv, xxvx_inv) = if xv_path.exists() && xx_path.exists() {
        (
            Some(read_dense_np(&xv_path, meta.p, meta.n)?),
            Some(read_dense_np(&xx_path, meta.n, meta.p)?),
        )
    } else {
        (None, None)
    };

    let scoring = meta.scoring_kind();

    Ok(NullModel {
        meta,
        id_include,
        sigma_i,
        sigma_i_x,
        cov,
        residuals,
        variance_ratio,
        xv,
        xxvx_inv,
        scoring,
    })
}

fn read_sigma_i_csc(path: &Path) -> Result<CscMatrix> {
    let mut f = File::open(path)?;
    let n = read_i32_le(&mut f)? as usize;
    let nnz = read_i32_le(&mut f)? as usize;
    let mut colptr = Vec::with_capacity(n + 1);
    for _ in 0..=n {
        colptr.push(read_i32_le(&mut f)? as usize);
    }
    let mut rowidx = Vec::with_capacity(nnz);
    for _ in 0..nnz {
        rowidx.push(read_i32_le(&mut f)? as usize);
    }
    let mut values = Vec::with_capacity(nnz);
    for _ in 0..nnz {
        values.push(read_f64_le(&mut f)?);
    }
    Ok(CscMatrix {
        n,
        colptr,
        rowidx,
        values,
    })
}

fn read_dense_np(path: &Path, n: usize, p: usize) -> Result<Vec<f64>> {
    let mut f = File::open(path)?;
    let rn = read_i32_le(&mut f)? as usize;
    let rp = read_i32_le(&mut f)? as usize;
    if rn != n || rp != p {
        return Err(ScoreError::msg(format!(
            "matrix dims {rn}x{rp} != expected {n}x{p} ({})",
            path.display()
        )));
    }
    let mut out = vec![0.0; n * p];
    for v in &mut out {
        *v = read_f64_le(&mut f)?;
    }
    Ok(out)
}

fn read_dense_p(path: &Path, p: usize) -> Result<Vec<f64>> {
    let mut f = File::open(path)?;
    let rp = read_i32_le(&mut f)? as usize;
    if rp != p {
        return Err(ScoreError::msg(format!("cov dim {rp} != expected {p}")));
    }
    let mut out = vec![0.0; p * p];
    for v in &mut out {
        *v = read_f64_le(&mut f)?;
    }
    Ok(out)
}

fn read_residuals(path: &Path, n: usize) -> Result<Vec<f64>> {
    let mut f = File::open(path)?;
    let rn = read_i32_le(&mut f)? as usize;
    if rn != n {
        return Err(ScoreError::msg(format!("residuals len {rn} != expected {n}")));
    }
    let mut out = vec![0.0; n];
    for v in &mut out {
        *v = read_f64_le(&mut f)?;
    }
    Ok(out)
}

pub fn write_sigma_i_csc(path: &Path, mat: &CscMatrix) -> Result<()> {
    let mut f = File::create(path)?;
    write_i32(&mut f, mat.n as i32)?;
    write_i32(&mut f, mat.values.len() as i32)?;
    for &cp in &mat.colptr {
        write_i32(&mut f, cp as i32)?;
    }
    for &ri in &mat.rowidx {
        write_i32(&mut f, ri as i32)?;
    }
    for &v in &mat.values {
        write_f64(&mut f, v)?;
    }
    Ok(())
}

fn write_i32(w: &mut impl Write, x: i32) -> Result<()> {
    w.write_all(&x.to_le_bytes())?;
    Ok(())
}

fn write_f64(w: &mut impl Write, x: f64) -> Result<()> {
    w.write_all(&x.to_le_bytes())?;
    Ok(())
}

/// y = A * x where A is n x n CSC and x is length n (column vector).
pub fn csc_matvec(mat: &CscMatrix, x: &[f64], y: &mut [f64]) {
    y.fill(0.0);
    for col in 0..mat.n {
        let start = mat.colptr[col];
        let end = mat.colptr[col + 1];
        let xv = x[col];
        if xv == 0.0 {
            continue;
        }
        for idx in start..end {
            y[mat.rowidx[idx]] += mat.values[idx] * xv;
        }
    }
}

/// Y = A * G where A is n x n CSC, G is n x k column-major, Y is n x k column-major.
pub fn csc_matmul_cols(mat: &CscMatrix, g: &[f64], k: usize, y: &mut [f64]) {
    y.fill(0.0);
    let n = mat.n;
    for j in 0..k {
        for col in 0..n {
            let start = mat.colptr[col];
            let end = mat.colptr[col + 1];
            let gv = g[col + j * n];
            if gv == 0.0 {
                continue;
            }
            for idx in start..end {
                y[mat.rowidx[idx] + j * n] += mat.values[idx] * gv;
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn csc_matvec_matches_dense() {
        let mat = CscMatrix {
            n: 3,
            colptr: vec![0, 1, 2, 3],
            rowidx: vec![0, 1, 2],
            values: vec![2.0, 3.0, 4.0],
        };
        let x = vec![1.0, 2.0, 3.0];
        let mut y = vec![0.0; 3];
        csc_matvec(&mat, &x, &mut y);
        assert_eq!(y, vec![2.0, 6.0, 12.0]);
    }

    #[test]
    fn family_detection() {
        let mut m = NullMeta {
            n: 1,
            p: 1,
            nnz: 1,
            family: "binomial".into(),
            source: None,
            variance_ratio: None,
            tau0: None,
            trait_type: None,
            tractor_mix_score_sha: None,
        };
        assert_eq!(m.scoring_kind(), ScoringKind::Gmmat);
        m.family = "gaussian".into();
        assert_eq!(m.scoring_kind(), ScoringKind::Felix);
        m.family = "binomial".into();
        m.source = Some("felix".into());
        assert_eq!(m.scoring_kind(), ScoringKind::Felix);
    }
}
