#!/usr/bin/env python3
"""Null calibration of FELIX dosage-VCF tests on 2,000 repeat loci."""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import chi2

OUT = Path(__file__).resolve().parent / "out"
COLS = ["p.value_anc1", "p.value_anc2", "p.value_ancALL", "P_het_admixed", "P_hom_admixed", "P_cct_admixed"]

for ph in sys.argv[1:] or ["y_bin2", "y_bin10", "y_q"]:
    path = OUT / f"rep_null_{ph}.txt"
    if not path.is_file():
        print(ph, "pending")
        continue
    d = pd.read_csv(path, sep="\t", na_values=["NA"])
    print(f"\n== {ph}: {len(d)} loci ==")
    rows = []
    for c in COLS:
        p = pd.to_numeric(d[c], errors="coerce")
        n_na = int(p.isna().sum())
        p = p.dropna()
        p = p[(p > 0) & (p <= 1)]
        lam = np.median(chi2.isf(p, 1)) / chi2.ppf(0.5, 1)
        rows.append({
            "test": c,
            "n": len(p),
            "dropped": n_na,
            "lambda": round(float(lam), 3),
            "rate@0.05": round(float((p < 0.05).mean()), 4),
            "rate@0.01": round(float((p < 0.01).mean()), 4),
            "rate@0.001": round(float((p < 0.001).mean()), 4),
        })
    print(pd.DataFrame(rows).to_string(index=False))
    spa = d["Is.SPA_ancALL"].astype(str).str.lower().eq("true").mean()
    print(f"SPA applied (ancALL): {spa:.3f} of loci")
    ac = pd.to_numeric(d["AC_Allele2_anc1"], errors="coerce")
    print("ancestry-1 tests dropped with negative summed dosage:",
          int(((ac < 0) & d["p.value_anc1"].isna()).sum()), "of", int(d["p.value_anc1"].isna().sum()))
