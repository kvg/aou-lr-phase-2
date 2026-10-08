#!/usr/bin/env python3
"""Figure: the Mendelian violation rate is not comparable across nanc."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib as mpl
import numpy as np
import pandas as pd

K5 = "#1b5e9c"
K6 = "#c9792b"
PROJ = "#4d4d4d"


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--in-dir", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    agg = pd.read_csv(Path(a.in_dir) / "k_bias_summary.tsv", sep="\t")

    import matplotlib.pyplot as plt
    fig, (axa, axb) = plt.subplots(1, 2, figsize=(9.6, 3.8),
                                   gridspec_kw={"width_ratios": [1.0, 1.0], "wspace": 0.34})

    # ---- (a) rate vs LAI error, one line per K ---------------------------
    u = agg[(agg.regime == "uniform") & (agg.scoring == "raw")]
    ks = sorted(u.K.unique())
    cmap = plt.cm.viridis(np.linspace(0.08, 0.88, len(ks)))
    for col, k in zip(cmap, ks):
        s = u[u.K == k].sort_values("eps")
        lw, z = (2.1, 5) if k in (5, 6) else (1.0, 3)
        axa.plot(s.eps * 100, s.rate_mean * 100, "-o", color=col, lw=lw, ms=3.4, zorder=z)
        # K=4 and K=7 end within a label height of their neighbours; the series
        # are monotone in K, so labelling the extremes and the two candidates
        # identifies every line without stacking text.
        if k in (3, 5, 6, 8):
            axa.annotate(f"K={k}", xy=(s.eps.iloc[-1] * 100, s.rate_mean.iloc[-1] * 100),
                         xytext=(3, 0), textcoords="offset points", fontsize=6.5,
                         color=col, va="center",
                         fontweight="bold" if k in (5, 6) else "normal")
    axa.set_xlabel("LAI misclassification rate per haplotype-locus (%)")
    axa.set_ylabel("Mendelian violations per\ninformative locus (%)")
    axa.set_title("More ancestry labels look worse at equal accuracy", loc="left")
    axa.set_xlim(-0.4, 11.6)
    axa.margins(y=0.08)

    # ---- (b) same truth, one extra noisy label ---------------------------
    eps = 0.02
    labels, k5v, k6v, k6p = [], [], [], []
    for purity in (0.85, 0.95):
        reg = f"split_purity{purity}"
        sel = agg[(agg.regime == reg) & (np.isclose(agg.eps, eps))]
        labels.append(f"split purity\n{int(purity * 100)}%")
        k5v.append(float(sel[(sel.K == 5) & (sel.scoring == "raw")].rate_mean.iloc[0]) * 100)
        k6v.append(float(sel[(sel.K == 6) & (sel.scoring == "raw")].rate_mean.iloc[0]) * 100)
        k6p.append(float(sel[(sel.K == 6) & (sel.scoring == "projected")].rate_mean.iloc[0]) * 100)

    x = np.arange(len(labels))
    w = 0.26
    axb.bar(x - w, k5v, w, color=K5, label="K=5 (true partition)")
    axb.bar(x, k6v, w, color=K6, label="K=6, scored as-is")
    axb.bar(x + w, k6p, w, color=PROJ, label="K=6, projected to panels")
    for xi, vals in zip(x, zip(k5v, k6v, k6p)):
        for dx, v in zip((-w, 0, w), vals):
            axb.text(xi + dx, v + 0.25, f"{v:.1f}", ha="center", fontsize=6.3)
    axb.set_xticks(x)
    axb.set_xticklabels(labels)
    axb.set_ylabel("Mendelian violations per\ninformative locus (%)")
    axb.set_title("Identical truth and accuracy; only the labels differ", loc="left")
    axb.set_ylim(0, max(k6v) * 1.30)
    axb.legend(frameon=False, fontsize=6.8, loc="upper left")

    for ax, letter in ((axa, "a"), (axb, "b")):
        ax.text(-0.02, 1.11, letter, transform=ax.transAxes, fontweight="bold",
                fontsize=10, va="top", ha="right")

    fig.savefig(a.out, dpi=300, bbox_inches="tight")

    r = fig.canvas.get_renderer()
    texts = [(t, t.get_window_extent(r)) for t in fig.findobj(mpl.text.Text)
             if t.get_text().strip() and t.get_visible()]
    tickl = {ax: set(ax.get_xticklabels(which="both") + ax.get_yticklabels(which="both"))
             for ax in fig.axes}
    ov = [(p[0].get_text(), q[0].get_text())
          for i, p in enumerate(texts) for q in texts[i + 1:] if p[1].overlaps(q[1])]
    ov += [(t.get_text(), "spine") for t, bt in texts for ax in fig.axes
           for s in ax.spines.values()
           if s.get_visible() and bt.overlaps(s.get_window_extent(r)) and t not in tickl[ax]]
    print("overlaps:", ov)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
