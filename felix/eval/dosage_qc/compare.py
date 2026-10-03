#!/usr/bin/env python3
"""Compare FELIX step-2 outputs for the M6 Option A spike."""

import sys
from pathlib import Path

import numpy as np
import pandas as pd

OUT = Path(__file__).resolve().parent / "out"
PCOLS = ["p.value_anc1", "p.value_anc2", "p.value_ancALL", "P_het_admixed", "P_hom_admixed", "P_cct_admixed"]


def load(name):
    path = OUT / f"{name}.txt"
    if not path.is_file():
        return None
    df = pd.read_csv(path, sep="\t", na_values=["NA"])
    df["key"] = df["POS"].astype(int)
    return df.set_index("key")


def max_rel_diff(a, b, cols):
    rows = []
    for c in cols:
        x = pd.to_numeric(a[c], errors="coerce")
        y = pd.to_numeric(b[c], errors="coerce")
        both = x.notna() & y.notna()
        nan_mismatch = int((x.isna() != y.isna()).sum())
        with np.errstate(divide="ignore", invalid="ignore"):
            # compare on -log10 scale; tiny p-values matter most
            lx, ly = -np.log10(x[both]), -np.log10(y[both])
            d = float(np.max(np.abs(lx - ly))) if both.any() else float("nan")
        rows.append({"column": c, "n": int(both.sum()), "nan_mismatch": nan_mismatch, "max_abs_diff_log10p": d})
    return pd.DataFrame(rows)


def lam(p):
    from scipy.stats import chi2

    p = pd.to_numeric(p, errors="coerce").dropna()
    p = p[(p > 0) & (p <= 1)]
    return float(np.median(chi2.isf(p, 1)) / chi2.ppf(0.5, 1)) if len(p) else float("nan")


for ph in ("y_bin10", "y_bin2", "y_q"):
    print(f"\n===== {ph} =====")
    a, b = load(f"felixla_{ph}"), load(f"snv_admix_{ph}")
    if a is not None and b is not None:
        common = a.index.intersection(b.index)
        print(f"A-1 SNV equivalence: felixla {len(a)} vs dosage-VCF {len(b)} sites, {len(common)} shared")
        print(max_rel_diff(a.loc[common], b.loc[common], PCOLS).to_string(index=False))
        spa = [c for c in a.columns if c.startswith("Is.SPA")]
        if spa:
            print("SPA used (felixla / dosage-VCF):",
                  {c: (int((a[c].astype(str) == "true").sum()), int((b[c].astype(str) == "true").sum())) for c in spa})
    r1, r2, rs = load(f"rep_s1_{ph}"), load(f"rep_s2_{ph}"), load(f"rep_shift_{ph}")
    if r1 is not None:
        print(f"repeat loci tested: {len(r1)} (of 40)")
        print("lambda P_cct_admixed (40 null loci; noisy):", round(lam(r1["P_cct_admixed"]), 3))
        for name, other in (("A-2 scale s vs 2s", r2), ("shift (x-min)/(max-min) vs scale s", rs)):
            if other is None:
                print(name, ": missing output")
                continue
            common = r1.index.intersection(other.index)
            print(f"{name}: {len(common)} shared loci")
            print(max_rel_diff(r1.loc[common], other.loc[common], PCOLS).to_string(index=False))
        b1 = pd.to_numeric(r1["BETA_ancALL"], errors="coerce")
        b2 = pd.to_numeric(r2["BETA_ancALL"], errors="coerce") if r2 is not None else None
        if b2 is not None:
            ratio = (b2 / b1).replace([np.inf, -np.inf], np.nan).dropna()
            print("BETA ratio (2s / s), expect 2:", round(float(ratio.median()), 4))
        print("AC_Allele2_ancALL range (s):", float(pd.to_numeric(r1["AC_Allele2_ancALL"]).min()),
              float(pd.to_numeric(r1["AC_Allele2_ancALL"]).max()))
sys.exit(0)
