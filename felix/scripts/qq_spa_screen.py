#!/usr/bin/env python3
"""Side-by-side QQ plot: FELIX results as reported, and after removing tests
whose SPA correction silently failed.

This is a POST-PROCESSING QC FIGURE. It is not part of FelixGenome.wdl and it
changes nothing about how a run executes.

What the two panels mean
------------------------
LEFT ("as reported") is the QQ of the p-value column exactly as FELIX wrote it.
Tests flagged by the `spa_fallback` screen are drawn in the alarm colour.

RIGHT ("SPA-fallback tests removed") drops those tests and recomputes the
expected quantiles on the retained set.

The screen is the one defined in flag_joint_test_artifacts.py and is imported
from there, so there is a single definition of it:
(Is.SPA == false) AND (|Tstat|/sqrt(var) > --spa-cutoff), in any ancestry
stratum or the pooled test. SPA was attempted and bailed, so the reported
p-value is the uncorrected normal-approximation score p.

READ THE TWO NUMBERS SEPARATELY. The screen removes a tail spike; it does not
change lambda, because lambda is a median over ~1.6M tests and the spike is a
few dozen of them. On BI_174 (63 cases) the screen took the repeat-dosage arm
from 52 genome-wide hits to 5 while lambda moved 0.6419 -> 0.6418. A panel that
still shows lambda = 0.64 after screening has NOT been fixed -- the deflation
is a separate problem (CCT conservatism at low case counts) and needs the case
floor, not this screen.

Usage
-----
    ./qq_spa_screen.py --results merged.felix.tsv.gz \
        --out-png qq_spa_screen.png --out-summary-json qq_spa_screen.json
"""
from __future__ import annotations

import argparse
import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

_spec = importlib.util.spec_from_file_location(
    "fjta", Path(__file__).with_name("flag_joint_test_artifacts.py"))
_fjta = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_fjta)

CHI2_MEDIAN = 0.4549364231195733  # qchisq(0.5, df=1)
ALARM = "#d1495b"


def lambda_gc(p: np.ndarray) -> float:
    """Median-chi2 genomic inflation factor, 1 df."""
    from scipy.stats import chi2
    p = np.asarray(p, dtype=float)
    p = p[np.isfinite(p) & (p > 0) & (p <= 1)]
    if p.size == 0:
        return float("nan")
    return float(np.median(chi2.isf(p, 1)) / CHI2_MEDIAN)


def qq_xy(p: np.ndarray):
    """Observed and expected -log10(p), sorted most-significant first."""
    p = np.asarray(p, dtype=float)
    p = np.sort(p[np.isfinite(p) & (p > 0) & (p <= 1)])
    n = p.size
    exp = -np.log10((np.arange(1, n + 1) - 0.5) / n)
    return exp, -np.log10(p)


def thin(exp, obs, keep_above=2.0, max_points=20000):
    """Keep every point in the tail; subsample the bulk in rank space.

    The tail is what the figure is about, so it is never thinned. Uses a
    ceiling step so the last bulk point is always retained.
    """
    tail = obs > keep_above
    idx_tail = np.flatnonzero(tail)
    idx_bulk = np.flatnonzero(~tail)
    budget = max(max_points - idx_tail.size, 1)
    if idx_bulk.size > budget:
        step = -(-idx_bulk.size // budget)      # ceiling division
        sel = idx_bulk[::step]
        if sel[-1] != idx_bulk[-1]:
            # a stride that does not land on the last element silently drops
            # the least-significant point and shortens the QQ tail-to-origin
            sel = np.append(sel, idx_bulk[-1])
        idx_bulk = sel
    idx = np.union1d(idx_tail, idx_bulk)
    return exp[idx], obs[idx]


def panel(ax, exp, obs, flagged_mask, title, lam, n_tested, n_removed=None,
          lim=None, annotate=True):
    if lim is None:
        lim = max(exp.max(), obs.max()) * 1.05
    ax.plot([0, lim], [0, lim], ls="--", lw=0.8, color="#888888", zorder=1)
    if flagged_mask is None:
        e, o = thin(exp, obs)
        ax.scatter(e, o, s=5, color="#3b6ea5", linewidths=0, zorder=2)
    else:
        e_k, o_k = thin(exp[~flagged_mask], obs[~flagged_mask])
        ax.scatter(e_k, o_k, s=5, color="#3b6ea5", linewidths=0, zorder=2)
        ax.scatter(exp[flagged_mask], obs[flagged_mask], s=18, facecolors="none",
                   edgecolors=ALARM, linewidths=0.9, zorder=3,
                   label="SPA correction failed")
        ax.legend(frameon=False, loc="upper left", handletextpad=0.3)
    ax.set_xlabel(r"Expected $-\log_{10}(p)$")
    if title:
        ax.set_title(title)
    if annotate:
        note = f"$\\lambda$ = {lam:.3f}\nn = {n_tested:,}"
        if n_removed is not None:
            note += f"\n{n_removed:,} removed"
        # upper-left, dropped clear of the legend when one is present
        y = 0.74 if flagged_mask is not None else 0.97
        ax.text(0.03, y, note, transform=ax.transAxes, ha="left", va="top")
    ax.set_xlim(0, max(exp.max(), 1.0) * 1.08)
    ax.set_ylim(0, lim)


def build(df: pd.DataFrame, pcol: str, num_ancs: int, spa_cutoff: float,
          p_threshold: float, out_png: Path, label: str = ""):
    p = pd.to_numeric(df[pcol], errors="coerce")
    ok = np.isfinite(p) & (p > 0) & (p <= 1)
    mask, _, _ = _fjta.spa_fallback_mask(df, num_ancs, spa_cutoff)

    p_ok = p[ok].to_numpy(dtype=float)
    m_ok = mask[ok].to_numpy(dtype=bool)
    order = np.argsort(p_ok)
    obs = -np.log10(p_ok[order])
    flag_sorted = m_ok[order]
    n = obs.size
    exp = -np.log10((np.arange(1, n + 1) - 0.5) / n)

    kept = p_ok[~m_ok]
    exp_k, obs_k = qq_xy(kept)

    lam_all, lam_kept = lambda_gc(p_ok), lambda_gc(kept)
    summary = dict(
        p_col=pcol, spa_cutoff=spa_cutoff, p_threshold=p_threshold,
        n_tested=int(n), n_spa_fallback=int(m_ok.sum()),
        lambda_as_reported=lam_all, lambda_screened=lam_kept,
        n_sig_as_reported=int((p_ok < p_threshold).sum()),
        n_sig_screened=int((kept < p_threshold).sum()),
        max_neglog10p_as_reported=float(obs.max()) if n else float("nan"),
        max_neglog10p_screened=float(obs_k.max()) if obs_k.size else float("nan"),
    )

    # Top row spans the full observed range so the left panel's spike is shown
    # in full (clipping it would hide the very thing the figure is about);
    # the bottom row repeats both panels zoomed to the bulk, where the
    # as-reported and screened distributions are indistinguishable.
    full = max(obs.max(), exp.max()) * 1.05
    zoom = max(np.percentile(obs, 99.999), obs_k.max(), exp.max()) * 1.10

    fig, axes = plt.subplots(2, 2, figsize=(7.2, 6.4))
    panel(axes[0][0], exp, obs, flag_sorted, "As reported by FELIX",
          lam_all, n, lim=full)
    panel(axes[0][1], exp_k, obs_k, None, "SPA-fallback tests removed",
          lam_kept, int(kept.size), n_removed=int(m_ok.sum()), lim=full)
    panel(axes[1][0], exp, obs, flag_sorted, "", lam_all, n,
          lim=zoom, annotate=False)
    panel(axes[1][1], exp_k, obs_k, None, "", lam_kept, int(kept.size),
          lim=zoom, annotate=False)
    axes[1][0].get_legend().remove()
    for ax in (axes[0][0], axes[0][1]):
        ax.set_xlabel("")
    for ax in (axes[0][0], axes[1][0]):
        ax.set_ylabel(r"Observed $-\log_{10}(p)$")
    axes[0][0].text(-0.28, 0.5, "full range", transform=axes[0][0].transAxes,
                    rotation=90, va="center", ha="center")
    axes[1][0].text(-0.28, 0.5, "bulk, zoomed", transform=axes[1][0].transAxes,
                    rotation=90, va="center", ha="center")
    if label:
        fig.suptitle(label, x=0.01, ha="left")
    fig.tight_layout()
    fig.savefig(out_png, dpi=300)
    plt.close(fig)
    return summary


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--results", type=Path, required=True)
    ap.add_argument("--out-png", type=Path, required=True)
    ap.add_argument("--out-summary-json", type=Path, default=None)
    ap.add_argument("--p-col", default=None, help="default: the joint CCT column")
    ap.add_argument("--p-threshold", type=float, default=5e-8)
    ap.add_argument("--num-ancs", type=int, default=5)
    ap.add_argument("--spa-cutoff", type=float, default=_fjta.SPA_CUTOFF_DEFAULT)
    ap.add_argument("--label", default="", help="figure suptitle, e.g. the phenotype")
    a = ap.parse_args()

    df = pd.read_csv(a.results, sep="\t", low_memory=False)
    pcol = a.p_col or _fjta.pick_joint_col(df)
    if pcol not in df.columns:
        raise SystemExit(f"p-value column {pcol!r} not in {a.results}")
    summary = build(df, pcol, a.num_ancs, a.spa_cutoff, a.p_threshold,
                    a.out_png, a.label)
    if a.out_summary_json:
        a.out_summary_json.write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))
    print(f"Wrote {a.out_png}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
