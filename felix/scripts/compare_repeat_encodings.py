#!/usr/bin/env python3
"""Compare repeat-unit encodings (dosage vs collapse vs split) on the same null.

Expects result TSVs from tractor-mix-score --mode felix with filenames containing
the encoding tag, e.g. phenotype.chr22.ru_dosage.felix.tsv.

Reports λGC on P_cct_admixed_c (or best available joint column) and genome-wide
hit counts per encoding.
"""

from __future__ import annotations

import argparse
import math
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ENCODINGS = ("dosage", "collapse", "split")
ENCODING_COLORS = {
    "dosage": "#1f4e79",
    "collapse": "#c0392b",
    "split": "#27ae60",
}
JOINT_P_COLS = [
    "P_cct_admixed_c",
    "P_cct_admixed",
    "P_hom_admixed_c",
    "P_het_admixed_c",
    "p.value_ancALL",
]


def erfcinv_scalar(p: float) -> float:
    lo, hi = 0.0, 10.0
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        if math.erfc(mid) > p:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def lambda_gc_from_p(pvals: np.ndarray) -> float:
    p = pvals[np.isfinite(pvals) & (pvals > 0) & (pvals <= 1)]
    if p.size < 10:
        return float("nan")
    x = 2.0 * (np.array([erfcinv_scalar(float(v)) for v in p]) ** 2)
    return float(np.median(x) / 0.4549364)


def pick_joint_p(df: pd.DataFrame) -> str | None:
    for c in JOINT_P_COLS:
        if c in df.columns:
            return c
    return None


def infer_encoding(path: Path) -> str | None:
    name = path.name.lower()
    for enc in ENCODINGS:
        if f".ru_{enc}." in name or f"_ru_{enc}." in name or f"_{enc}." in name:
            return enc
    return None


def discover_by_encoding(directory: Path, phenotype: str | None) -> dict[str, Path]:
    out: dict[str, Path] = {}
    for path in sorted(directory.glob("*.tsv")):
        if path.name.endswith(".raw.txt"):
            continue
        enc = infer_encoding(path)
        if enc is None:
            continue
        if phenotype is not None:
            safe = re.sub(r"[^A-Za-z0-9._-]+", "_", phenotype.strip())
            if safe not in path.name and phenotype not in path.name:
                continue
        out[enc] = path
    return out


def qq_plot(pvals: np.ndarray, title: str, out_png: Path, color: str) -> float:
    p = np.sort(pvals[np.isfinite(pvals) & (pvals > 0) & (pvals <= 1)])
    lam = lambda_gc_from_p(p)
    fig, ax = plt.subplots(figsize=(5, 5))
    if p.size:
        n = p.size
        exp = -np.log10(np.arange(1, n + 1) / (n + 1))
        obs = -np.log10(p)
        ax.scatter(exp, obs, s=4, alpha=0.45, c=color, rasterized=True)
        m = max(float(exp.max()), float(obs.max()))
        ax.plot([0, m], [0, m], ls="--", c="gray", lw=1)
        ax.set_title(f"{title}\nλGC≈{lam:.3f} n={n:,}")
    else:
        ax.set_title(title + " (empty)")
    ax.set_xlabel(r"Expected $-\log_{10}(p)$")
    ax.set_ylabel(r"Observed $-\log_{10}(p)$")
    fig.savefig(out_png, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return lam


def summarize_encoding(path: Path, pcol: str, p_threshold: float) -> dict:
    df = pd.read_csv(path, sep="\t", low_memory=False)
    p = pd.to_numeric(df[pcol], errors="coerce").to_numpy(dtype=float)
    n_tested = int(np.isfinite(p).sum())
    n_hits = int(np.sum(np.isfinite(p) & (p < p_threshold)))
    return {
        "n_variants": len(df),
        "n_tested": n_tested,
        "lambda_gc": lambda_gc_from_p(p),
        "n_gw_hits": n_hits,
        "p_column": pcol,
    }


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--results-dir",
        type=Path,
        required=True,
        help="Directory with *.ru_{dosage,collapse,split}.felix.tsv",
    )
    ap.add_argument("--phenotype", default=None, help="Filter to one phenotype name")
    ap.add_argument("--p-threshold", type=float, default=5e-8)
    ap.add_argument("--out-dir", type=Path, required=True)
    args = ap.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    found = discover_by_encoding(args.results_dir, args.phenotype)
    if not found:
        raise SystemExit(f"No encoding-tagged TSVs in {args.results_dir}")

    rows = []
    for enc in ENCODINGS:
        if enc not in found:
            continue
        path = found[enc]
        df_head = pd.read_csv(path, sep="\t", nrows=5)
        pcol = pick_joint_p(df_head)
        if pcol is None:
            full = pd.read_csv(path, sep="\t", nrows=0)
            pcol = pick_joint_p(full)
        if pcol is None:
            raise SystemExit(f"No joint p column in {path}")
        stats = summarize_encoding(path, pcol, args.p_threshold)
        rows.append({"encoding": enc, "result_file": str(path), **stats})

        p = pd.read_csv(path, sep="\t", usecols=[pcol], low_memory=False)[pcol].to_numpy(
            dtype=float
        )
        qq_plot(
            p,
            f"{enc} encoding",
            args.out_dir / f"qq_{enc}.png",
            ENCODING_COLORS[enc],
        )

    summary = pd.DataFrame(rows)
    summary.to_csv(args.out_dir / "encoding_comparison.tsv", sep="\t", index=False)

    cal = summary[["encoding", "lambda_gc", "n_tested", "n_gw_hits", "p_column"]]
    cal.to_csv(args.out_dir / "calibration_summary.tsv", sep="\t", index=False)

    lines = [
        "# Repeat encoding comparison (same FELIX null)",
        "",
        "Compare λGC and genome-wide hit counts across dosage / collapse / split.",
        "Do not interpret hit-count differences without checking calibration.",
        "",
    ]
    try:
        lines.append(cal.to_markdown(index=False))
    except Exception:
        lines.append("```\n" + cal.to_string(index=False) + "\n```")
    (args.out_dir / "encoding_comparison.md").write_text("\n".join(lines) + "\n")
    print(f"Wrote {args.out_dir / 'calibration_summary.tsv'} ({len(rows)} encodings)")


if __name__ == "__main__":
    main()
