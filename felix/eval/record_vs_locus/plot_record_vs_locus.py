#!/usr/bin/env python3
"""Figure: what each association input delivers to FELIX on the defect fixture."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib as mpl
import numpy as np
import pandas as pd

LOCUS_COL = "#1b5e9c"   # locus path (focal)
RECORD_COL = "#c9792b"  # record path (comparator)

DEFECT_LABEL = {
    "1_split_across_records": "Locus split across\ntwo callset records",
    "2_impure_allele": "Allele length not a whole\nnumber of repeat units",
    "3_seqres_deletion_no_svtype": "Sequence-resolved deletion\nwith no SVTYPE",
    "4_mixed_symbolic_record": "One record with both a\ndeletion and an insertion ALT",
    "4b_mixed_seqresolved_record": "Same, sequence-resolved\n(no SVTYPE)",
}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--in-dir", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    d = Path(a.in_dir)

    cmp = pd.read_csv(d / "per_sample_ancestry_dosage.tsv", sep="\t")
    summary = pd.read_csv(d / "defect_summary.tsv", sep="\t")
    nz = cmp[cmp["expected"] != 0].copy()
    nz["record_path"] = nz["record_path"].fillna(0.0)

    import matplotlib.pyplot as plt
    fig, (axa, axb) = plt.subplots(
        1, 2, figsize=(10.0, 3.9), gridspec_kw={"width_ratios": [1.15, 1.0], "wspace": 0.42}
    )

    # ---- (a) share of per-ancestry dosages delivered correctly ------------
    s = summary.sort_values("defect", ascending=False).reset_index(drop=True)
    y = np.arange(len(s))
    frac_loc = s["locus_path_correct"] / s["nonzero_dik"]
    frac_rec = s["record_path_correct"] / s["nonzero_dik"]
    for yi, fr, fl in zip(y, frac_rec, frac_loc):
        axa.plot([fr, fl], [yi, yi], color="0.75", lw=1.2, zorder=1)
    axa.scatter(frac_rec, y, s=46, facecolors="white", edgecolors=RECORD_COL,
                linewidths=1.6, zorder=3, label="record-level input")
    axa.scatter(frac_loc, y, s=46, color=LOCUS_COL, zorder=3, label="locus-level input")
    axa.set_yticks(y)
    axa.set_yticklabels([DEFECT_LABEL[v] for v in s["defect"]])
    axa.set_xlim(-0.08, 1.14)
    axa.set_xticks([0, 0.5, 1.0])
    axa.set_xticklabels(["0", "50%", "100%"])
    axa.set_xlabel("Per-ancestry dosages delivered correctly")
    axa.set_title("Locus-level coding recovers every dosage", loc="left")
    for yi, n in zip(y, s["nonzero_dik"]):
        axa.text(1.10, yi, f"n={n}", va="center", ha="right", fontsize=6, color="0.4")
    axa.legend(frameon=False, loc="upper left", bbox_to_anchor=(0.04, 0.90), fontsize=7)
    axa.margins(y=0.14)

    # ---- (b) produced against expected dosage -----------------------------
    lim = 24
    axb.plot([-lim, lim], [-lim, lim], ls=(0, (4, 3)), lw=0.9, color="0.6", zorder=1)
    axb.axhline(0, lw=0.7, color="0.85", zorder=0)
    untested = ~nz["record_tested"]
    axb.scatter(nz.loc[untested, "expected"], nz.loc[untested, "record_path"],
                s=52, marker="s", facecolors="white", edgecolors=RECORD_COL,
                linewidths=1.6, zorder=3, label="record-level: locus never tested")
    axb.scatter(nz.loc[~untested, "expected"], nz.loc[~untested, "record_path"],
                s=38, marker="s", color=RECORD_COL, alpha=0.9, zorder=3,
                label="record-level: locus tested")
    axb.scatter(nz["expected"], nz["locus_path"], s=26, color=LOCUS_COL, zorder=4,
                label="locus-level")
    axb.set_xlabel("Expected per-ancestry dosage $D_{ik}$ (repeat units)")
    axb.set_ylabel("Dosage handed to FELIX")
    axb.set_title("Dropped loci enter as zero; mixed records flip sign", loc="left")
    axb.set_xlim(-6.0, 24.5)
    axb.set_ylim(-6.0, 24.5)
    axb.legend(frameon=False, loc="upper left", bbox_to_anchor=(0.16, 1.0), fontsize=7)

    flip = nz[(~untested) & (nz["record_path"] * nz["expected"] < 0)]
    if len(flip):
        r = flip.iloc[0]
        axb.annotate("insertion coded\nas a deletion",
                     xy=(r["expected"], r["record_path"]), xytext=(7.5, -5.2),
                     fontsize=6.5, color=RECORD_COL, ha="left", va="bottom",
                     arrowprops=dict(arrowstyle="-", lw=0.7, color=RECORD_COL,
                                     shrinkA=2, shrinkB=3))
    drop = nz[untested & (nz["expected"] < 0)]
    if len(drop):
        r = drop.sort_values("expected").iloc[0]
        axb.annotate("record dropped by the\nwhole-repeat-unit check",
                     xy=(r["expected"], 0.0), xytext=(-5.4, 7.0),
                     fontsize=6.5, color=RECORD_COL, ha="left", va="bottom",
                     arrowprops=dict(arrowstyle="-", lw=0.7, color=RECORD_COL,
                                     shrinkA=2, shrinkB=3))

    for ax, letter in ((axa, "a"), (axb, "b")):
        ax.text(-0.015, 1.1, letter, transform=ax.transAxes, fontweight="bold",
                fontsize=10, va="top", ha="right")

    fig.savefig(a.out, dpi=300, bbox_inches="tight")

    r = fig.canvas.get_renderer()
    texts = [(t, t.get_window_extent(r)) for t in fig.findobj(mpl.text.Text)
             if t.get_text().strip() and t.get_visible()]
    tickl = {ax: set(ax.get_xticklabels(which="both") + ax.get_yticklabels(which="both"))
             for ax in fig.axes}
    ov = [(x[0].get_text(), y[0].get_text())
          for i, x in enumerate(texts) for y in texts[i + 1:] if x[1].overlaps(y[1])]
    ov += [(t.get_text(), "spine")
           for t, bt in texts for ax in fig.axes for s in ax.spines.values()
           if s.get_visible() and bt.overlaps(s.get_window_extent(r)) and t not in tickl[ax]]
    print("overlaps:", ov)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
