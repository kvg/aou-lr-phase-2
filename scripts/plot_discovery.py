#!/usr/bin/env python3
"""Plot Ebert-style cumulative discovery curves from ebert_discovery.py output."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


FREQ_ORDER = ["shared", "major", "polymorphic", "singleton"]
FREQ_COLORS = {
    "shared": "#d62728",
    "major": "#9467bd",
    "polymorphic": "#1f77b4",
    "singleton": "#17becf",
}


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--discovery-tsv", required=True)
    p.add_argument("--out-png", required=True)
    p.add_argument("--title", default="Cumulative SV discovery")
    args = p.parse_args()

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    df = pd.read_csv(args.discovery_tsv, sep="\t")
    strata = sorted(df["stratum"].unique())
    n = len(strata)
    fig, axes = plt.subplots(1, n, figsize=(5 * n, 4), sharey=True)
    if n == 1:
        axes = [axes]

    for ax, stratum in zip(axes, strata):
        sub = df[df["stratum"] == stratum]
        # pivot to stacked areas
        wide = (
            sub.pivot_table(
                index="sample_index",
                columns="freq_class",
                values="cumulative_count",
                aggfunc="first",
            )
            .reindex(columns=[c for c in FREQ_ORDER if c in sub["freq_class"].unique()])
            .fillna(0)
        )
        x = wide.index.values
        bottom = None
        for fc in wide.columns:
            y = wide[fc].values
            if bottom is None:
                ax.fill_between(x, 0, y, color=FREQ_COLORS.get(fc, "#999999"), label=fc, alpha=0.9)
                bottom = y.copy()
            else:
                ax.fill_between(
                    x,
                    bottom,
                    bottom + y,
                    color=FREQ_COLORS.get(fc, "#999999"),
                    label=fc,
                    alpha=0.9,
                )
                bottom = bottom + y
        # ancestry boundary: first African sample
        anc = sub.drop_duplicates("sample_index").sort_values("sample_index")
        afr = anc[anc["ancestry"].str.lower().isin(["afr", "african"])]
        if not afr.empty:
            ax.axvline(afr["sample_index"].iloc[0] - 0.5, color="grey", ls="--", lw=1)
        ax.set_title(str(stratum))
        ax.set_xlabel("Discovery samples (ancestry-ordered)")
        ax.set_ylabel("Variant count (cumulative)")

    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper left", frameon=False)
    fig.suptitle(args.title)
    fig.tight_layout()
    Path(args.out_png).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.out_png, dpi=150)
    print(f"Wrote {args.out_png}")


if __name__ == "__main__":
    main()
