#!/usr/bin/env python3
"""Compare patched-FELIX validation outputs (val/) with stock outputs (out/)."""

import filecmp
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import chi2

HERE = Path(__file__).resolve().parent
OUT, VAL = HERE / "out", HERE / "val"
PCOLS = ["p.value_anc1", "p.value_anc2", "p.value_ancALL", "P_het_admixed", "P_hom_admixed", "P_cct_admixed"]


def load(path):
    d = pd.read_csv(path, sep="\t", na_values=["NA"])
    return d.set_index(d["POS"].astype(int))


print("== V1: switch off, patched image vs stock FELIX (byte-identical TSV) ==")
for name in ("felixla_y_bin10", "snv_admix_y_bin10", "rep_s1_y_bin10"):
    a, b = OUT / f"{name}.txt", VAL / f"off_{name}.txt"
    same = b.is_file() and filecmp.cmp(a, b, shallow=False)
    print(f"  {name}: {'IDENTICAL' if same else 'DIFFERENT or missing'}")

print("\n== V2: switch on, SNVs vs stock dosage-VCF run ==")
a, b = load(OUT / "snv_admix_y_bin10.txt"), load(VAL / "on_snv_admix_y_bin10.txt")
common = a.index.intersection(b.index)
print(f"  sites: stock {len(a)}, carrier-QC {len(b)}, shared {len(common)}")
for c in PCOLS:
    x, y = pd.to_numeric(a.loc[common, c]), pd.to_numeric(b.loc[common, c])
    both = x.notna() & y.notna()
    diff = float(np.max(np.abs(np.log10(x[both]) - np.log10(y[both])))) if both.any() else float("nan")
    print(f"  {c:15s} tested stock={int(x.notna().sum()):3d} carrier={int(y.notna().sum()):3d} "
          f"max|dlog10p| on both={diff:.2e}")

print("\n== V3: switch on, 2,000 null repeat loci (minMAC=20 carriers) ==")
for ph in ("y_bin2", "y_bin10", "y_q"):
    path = VAL / f"on_rep_null_{ph}.txt"
    if not path.is_file():
        print(f"  {ph}: missing")
        continue
    d = pd.read_csv(path, sep="\t", na_values=["NA"])
    stock = pd.read_csv(OUT / f"rep_null_{ph}.txt", sep="\t", na_values=["NA"])
    rows = []
    for c in PCOLS:
        p = pd.to_numeric(d[c], errors="coerce")
        dropped = int(p.isna().sum())
        p = p.dropna()
        p = p[(p > 0) & (p <= 1)]
        rows.append({
            "test": c, "n": len(p), "dropped": dropped,
            "stock_dropped": int(pd.to_numeric(stock[c], errors="coerce").isna().sum()),
            "lambda": round(float(np.median(chi2.isf(p, 1)) / chi2.ppf(0.5, 1)), 3),
            "rate@0.05": round(float((p < 0.05).mean()), 4),
            "rate@0.01": round(float((p < 0.01).mean()), 4),
            "rate@0.001": round(float((p < 0.001).mean()), 4),
        })
    print(f"  {ph}: loci out {len(d)} (stock {len(stock)})")
    print(pd.DataFrame(rows).to_string(index=False))

print("\n== V4: wrapper end to end ==")
w = VAL / "wrapper_rep_null_y_bin2.tsv"
if w.is_file():
    d = pd.read_csv(w, sep="\t", na_values=["NA"])
    print(f"  wrapper output: {len(d)} loci; P_cct_admixed_c non-NA {int(d['P_cct_admixed_c'].notna().sum())}")
else:
    print("  wrapper output missing (see val/wrapper.log)")
