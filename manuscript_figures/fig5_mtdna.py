"""Figure 5: mtDNA haplogroups, heteroplasmy, and phenome-wide associations."""

from __future__ import annotations

from pathlib import Path

import numpy as np
from matplotlib.gridspec import GridSpec
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

from .style import (
    ANC_COLORS,
    ANC_ORDER,
    COL_IN,
    OKABE,
    P2_ANC_N,
    ax_panel,
    despine,
    new_figure,
    save,
    stamp_mockup,
)

# Macro-haplogroup palette.
HG_ORDER = ["L", "M", "N", "R", "H", "U", "J/T", "A/B/C/D", "other"]
HG_COLORS = {
    "L": OKABE["blue"],
    "M": OKABE["orange"],
    "N": OKABE["green"],
    "R": OKABE["pink"],
    "H": OKABE["sky"],
    "U": OKABE["vermillion"],
    "J/T": "#8E44AD",
    "A/B/C/D": "#2C7A7B",
    "other": OKABE["grey"],
}


def _haplogroup_fractions() -> dict[str, dict[str, float]]:
    # Plausible ancestry-stratified macro-haplogroup mix (placeholder).
    return {
        "AFR": {"L": 0.78, "M": 0.04, "N": 0.03, "R": 0.03, "H": 0.04, "U": 0.03, "J/T": 0.02, "A/B/C/D": 0.01, "other": 0.02},
        "AMR": {"L": 0.12, "M": 0.06, "N": 0.05, "R": 0.08, "H": 0.14, "U": 0.07, "J/T": 0.05, "A/B/C/D": 0.38, "other": 0.05},
        "EAS": {"L": 0.01, "M": 0.18, "N": 0.08, "R": 0.06, "H": 0.02, "U": 0.02, "J/T": 0.01, "A/B/C/D": 0.58, "other": 0.04},
        "EUR": {"L": 0.02, "M": 0.03, "N": 0.04, "R": 0.08, "H": 0.42, "U": 0.18, "J/T": 0.16, "A/B/C/D": 0.02, "other": 0.05},
        "MID": {"L": 0.06, "M": 0.10, "N": 0.08, "R": 0.14, "H": 0.18, "U": 0.16, "J/T": 0.18, "A/B/C/D": 0.04, "other": 0.06},
        "SAS": {"L": 0.03, "M": 0.42, "N": 0.08, "R": 0.18, "H": 0.08, "U": 0.08, "J/T": 0.05, "A/B/C/D": 0.03, "other": 0.05},
        "OTH": {"L": 0.18, "M": 0.12, "N": 0.08, "R": 0.10, "H": 0.16, "U": 0.10, "J/T": 0.08, "A/B/C/D": 0.12, "other": 0.06},
    }


def _panel_haplogroups(ax) -> None:
    fr = _haplogroup_fractions()
    ancs = [a for a in ANC_ORDER]
    x = np.arange(len(ancs))
    bottom = np.zeros(len(ancs))
    for hg in HG_ORDER:
        vals = np.array([fr[a][hg] * P2_ANC_N[a] for a in ancs])
        ax.bar(x, vals, bottom=bottom, color=HG_COLORS[hg], width=0.78, label=hg, zorder=2)
        bottom += vals
    ax.set_xticks(x)
    ax.set_xticklabels(ancs)
    ax.set_ylabel("Participants")
    ax.set_title("Mitochondrial macro-haplogroup composition", loc="left", pad=2)
    despine(ax)
    ax.legend(loc="upper left", ncol=3, fontsize=5.6, frameon=False, bbox_to_anchor=(0.62, 1.02))


def _panel_vaf(ax, rng: np.random.Generator) -> None:
    # Heteroplasmy VAF: HiFi floor ~3–5%; short-read NUMT-confounded floor ~10%.
    vaf = np.concatenate(
        [
            rng.beta(1.4, 40, 4000) * 0.25,  # low-level heteroplasmy
            rng.uniform(0.03, 0.15, 600),
        ]
    )
    vaf = vaf[(vaf > 0.005) & (vaf < 0.35)]
    bins = np.linspace(0, 0.30, 40)
    ax.hist(vaf, bins=bins, color=OKABE["blue"], alpha=0.85, zorder=2)
    ax.axvline(0.03, color=OKABE["green"], ls="--", lw=1.0, label="HiFi floor ~3%")
    ax.axvline(0.10, color=OKABE["vermillion"], ls="--", lw=1.0, label="srWGS NUMT floor ~10%")
    ax.set_xlabel("Heteroplasmic variant allele fraction")
    ax.set_ylabel("Sites × samples")
    ax.set_title("Heteroplasmy VAF (placeholder)", loc="left", pad=2)
    despine(ax)
    ax.legend(loc="upper right")
    ax.set_xlim(0, 0.30)


def _panel_age(ax, rng: np.random.Generator) -> None:
    n = 2200
    age = rng.uniform(20, 88, n)
    # Burden increases with age; slight ancestry offset.
    burden = np.clip(np.exp(-0.15 + 0.028 * (age - 40) + rng.normal(0, 0.32, n)), 0.2, 18)
    burden = burden * rng.uniform(0.7, 1.4, n)
    ax.scatter(age, burden, s=6, c="#9AA0A6", alpha=0.25, linewidths=0, rasterized=True)
    # Binned mean
    edges = np.arange(20, 90, 8)
    xs, ys, lo, hi = [], [], [], []
    for a, b in zip(edges[:-1], edges[1:]):
        m = (age >= a) & (age < b)
        if m.sum() < 20:
            continue
        xs.append((a + b) / 2)
        ys.append(np.mean(burden[m]))
        lo.append(np.percentile(burden[m], 25))
        hi.append(np.percentile(burden[m], 75))
    ax.fill_between(xs, lo, hi, color=OKABE["blue"], alpha=0.18, linewidth=0)
    ax.plot(xs, ys, color=OKABE["blue"], lw=1.4, marker="o", ms=3.2)
    ax.set_xlabel("Age (years)")
    ax.set_ylabel("Heteroplasmic sites per person")
    ax.set_title("Heteroplasmy burden vs age", loc="left", pad=2)
    despine(ax)
    ax.text(0.05, 0.92, r"$r \approx 0.34$ (placeholder)", transform=ax.transAxes, fontsize=6, color="#444444")


def _panel_circular(ax, rng: np.random.Generator) -> None:
    cats = [
        "Metabolic",
        "Cardio",
        "Neuro",
        "Renal",
        "Immune",
        "Neoplasm",
        "Heme",
        "Sensory",
        "GI",
        "Repro",
        "Musculo",
        "Other",
    ]
    n_cat = len(cats)
    theta_edges = np.linspace(0, 2 * np.pi, n_cat + 1)
    # background points
    planted = {
        "Metabolic": 11.4,
        "Cardio": 10.6,
        "Neuro": 9.2,
        "Heme": 8.8,
        "Immune": 8.1,
        "Neoplasm": 7.9,
    }
    for i, cat in enumerate(cats):
        th0, th1 = theta_edges[i], theta_edges[i + 1]
        n = 70
        th = rng.uniform(th0 + 0.05, th1 - 0.05, n)
        r = np.clip(rng.exponential(0.65, n), 0, 4.2)
        ax.scatter(th, r, s=5, c=OKABE["blue"] if i % 2 == 0 else OKABE["sky"], alpha=0.7, linewidths=0, rasterized=True)
        if cat in planted:
            hit_th = 0.5 * (th0 + th1)
            ax.scatter([hit_th], [planted[cat]], s=18, c=OKABE["vermillion"], zorder=4, linewidths=0)
            short = {"Metabolic": "T2D", "Cardio": "CAD", "Neuro": "PD", "Heme": "anemia", "Immune": "RA", "Neoplasm": "CLL"}
            ax.text(
                hit_th,
                planted[cat] + 0.55,
                short[cat],
                ha="center",
                va="bottom",
                fontsize=5.2,
                color=OKABE["vermillion"],
                fontweight="medium",
            )
        mid = 0.5 * (th0 + th1)
        ax.text(mid, 15.4, cat, ha="center", va="center", fontsize=5.6)
    ax.set_theta_offset(np.pi / 2)
    ax.set_theta_direction(-1)
    ax.set_ylim(0, 16.8)
    ax.set_yticks([5, 7.3, 10])
    ax.set_yticklabels(["5", "7.3", "10"], fontsize=5.5)
    ax.set_xticks(theta_edges[:-1])
    ax.set_xticklabels([])
    ax.plot(np.linspace(0, 2 * np.pi, 200), np.full(200, 7.3), color=OKABE["vermillion"], ls="--", lw=0.7)
    ax.set_title("Haplogroup PheWAS (circular Manhattan, placeholder)", loc="left", pad=10, y=1.08)
    ax.spines["polar"].set_linewidth(0.5)


def render() -> tuple[Path, Path]:
    rng = np.random.default_rng(5)
    ax_panel("fig5_haplogroups", 2.35, _panel_haplogroups, left=0.08, right=0.98, top=0.88, bottom=0.16)
    ax_panel("fig5_heteroplasmy_vaf", 2.85, lambda ax: _panel_vaf(ax, rng), width=COL_IN * 1.12, left=0.16, right=0.96, top=0.88, bottom=0.16)
    ax_panel("fig5_heteroplasmy_age", 2.85, lambda ax: _panel_age(ax, rng), width=COL_IN * 1.12, left=0.16, right=0.96, top=0.88, bottom=0.16)
    ax_panel(
        "fig5_circular_manhattan",
        3.55,
        lambda ax: _panel_circular(ax, rng),
        projection="polar",
        left=0.10,
        right=0.90,
        top=0.88,
        bottom=0.10,
    )
    fig = new_figure(7.55)
    gs = GridSpec(
        3,
        2,
        figure=fig,
        left=0.08,
        right=0.98,
        top=0.96,
        bottom=0.08,
        wspace=0.30,
        hspace=0.46,
        height_ratios=[0.95, 0.95, 1.45],
    )
    ax_a = fig.add_subplot(gs[0, :])
    ax_b = fig.add_subplot(gs[1, 0])
    ax_c = fig.add_subplot(gs[1, 1])
    ax_d = fig.add_subplot(gs[2, :], projection="polar")

    _panel_haplogroups(ax_a)
    _panel_vaf(ax_b, rng)
    _panel_age(ax_c, rng)
    _panel_circular(ax_d, rng)

    stamp_mockup(fig)
    return save(fig, "fig5_mtdna")
