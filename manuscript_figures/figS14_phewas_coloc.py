"""Figure S14: One SV, one highlighted disease, multi-omic colocalization.

A classic PheWAS (one SV × the phenome) with a single disease called out,
then stacked regional tracks for methylation, RNA, proteomics, and that
disease — the molecular reason the phecode lights up.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
from matplotlib.gridspec import GridSpec
from matplotlib.patches import Rectangle

from .style import (
    COL_IN,
    OKABE,
    P2_N,
    ax_panel,
    despine,
    new_figure,
    save,
    save_panel,
    stamp_mockup,
)

SV_X = 0.58  # position in the 1.2 Mb window
WINDOW_MB = 1.20

# Phecode groups in display order (classic PheWAS x-axis).
PHE_GROUPS = [
    ("Infectious", "Infectious", 8, "#B8B8B8"),
    ("Neoplasms", "Neoplasms", 9, "#A0A0A0"),
    ("Endocrine", "Endocrine", 12, OKABE["orange"]),
    ("Hematopoietic", "Heme", 8, OKABE["pink"]),
    ("Mental", "Mental", 8, "#C4C4C4"),
    ("Neurological", "Neuro", 8, "#B0B0B0"),
    ("Circulatory", "Circulatory", 14, OKABE["vermillion"]),
    ("Respiratory", "Resp.", 8, "#A8A8A8"),
    ("Digestive", "Digestive", 10, OKABE["green"]),
    ("Genitourinary", "GU", 14, OKABE["blue"]),
    ("Dermatologic", "Skin", 6, "#C8C8C8"),
    ("Musculoskeletal", "MSK", 10, OKABE["sky"]),
    ("Congenital", "Cong.", 5, "#D0D0D0"),
    ("Symptoms", "Symptoms", 8, "#BEBEBE"),
]

# Highlighted disease (index within Genitourinary group).
HIGHLIGHT_GROUP = "Genitourinary"
HIGHLIGHT_NAME = "Chronic kidney disease"
HIGHLIGHT_OFFSET = 4  # 0-based within the group


def _phewas_values(rng: np.random.Generator) -> tuple[np.ndarray, np.ndarray, np.ndarray, list[str], np.ndarray]:
    """Return x, logp, color, group label per point, and highlight index."""
    xs, logps, cols, groups = [], [], [], []
    x = 0.0
    highlight_i = None
    for name, _short, n, color in PHE_GROUPS:
        u = rng.uniform(1e-3, 1.0, n)
        lp = np.clip(-np.log10(u), 0.05, 3.6)
        # Plant APOL1-adjacent clinical signal in genitourinary + a little circulatory.
        if name == "Genitourinary":
            planted = np.array([4.2, 6.8, 8.1, 7.6, 12.4, 9.3, 8.6, 5.1, 3.4, 2.2, 1.8, 1.4, 1.1, 0.7])
            lp = planted[:n]
        elif name == "Circulatory":
            lp[2] = max(lp[2], 4.8)  # hypertension
            lp[5] = max(lp[5], 3.9)
        elif name == "Endocrine":
            lp[1] = max(lp[1], 3.1)
        for j in range(n):
            if name == HIGHLIGHT_GROUP and j == HIGHLIGHT_OFFSET:
                highlight_i = len(xs)
            xs.append(x + j)
            logps.append(lp[j])
            cols.append(color)
            groups.append(name)
        x += n + 1.6  # gap between groups
    return (
        np.asarray(xs),
        np.asarray(logps),
        np.asarray(cols),
        groups,
        np.array([highlight_i]),
    )


def _panel_phewas(ax, rng: np.random.Generator) -> None:
    x, logp, cols, groups, hi = _phewas_values(rng)
    hi = int(hi[0])

    ax.scatter(x, logp, s=11, c=cols, linewidths=0, zorder=2, rasterized=True)
    ax.axhline(7.5, color="#444444", ls="--", lw=0.6, zorder=1)

    # Highlight one disease.
    ax.scatter(
        [x[hi]],
        [logp[hi]],
        s=86,
        marker="D",
        c=OKABE["blue"],
        edgecolors="white",
        linewidths=0.7,
        zorder=5,
    )
    ax.annotate(
        HIGHLIGHT_NAME,
        xy=(x[hi], logp[hi]),
        xytext=(x[hi] - 18, logp[hi] + 2.4),
        fontsize=7.0,
        fontweight="bold",
        color=OKABE["blue"],
        ha="left",
        arrowprops=dict(arrowstyle="-|>", color=OKABE["blue"], lw=0.9, shrinkA=0, shrinkB=3),
    )
    ax.text(
        x[hi] - 18,
        logp[hi] + 1.55,
        "APOL1 tag SV  ·  AFR-stratified",
        fontsize=5.6,
        color="#444444",
        ha="left",
    )

    # Category strips under the axis.
    ymin = -1.85
    ax.set_ylim(ymin, 16.2)
    x_cursor = 0.0
    for name, short, n, color in PHE_GROUPS:
        ax.add_patch(
            Rectangle((x_cursor - 0.4, ymin), n - 0.2, 1.55, facecolor=color, alpha=0.35, edgecolor="none", zorder=0, clip_on=False)
        )
        ax.text(x_cursor + (n - 1) / 2, ymin + 0.72, short, ha="center", va="center", fontsize=5.3, color="#222222")
        x_cursor += n + 1.6

    ax.set_xlim(-1.2, x.max() + 2.0)
    ax.set_xticks([])
    ax.set_ylabel(r"$-\log_{10}(P)$")
    ax.set_title(
        f"PheWAS of one SV (APOL1 tag) in n={P2_N:,}  ·  one disease highlighted",
        loc="left",
        pad=2,
    )
    despine(ax, bottom=False)
    ax.spines["bottom"].set_visible(False)
    ax.tick_params(bottom=False)


def _shared_snps(rng: np.random.Generator, n: int = 220) -> tuple[np.ndarray, np.ndarray]:
    """Same variants in every track; LD-shaped association with the SV."""
    pos = np.sort(rng.uniform(0.02, 0.98, n))
    dist = np.abs(pos - SV_X)
    ld = np.exp(-0.5 * (dist / 0.07) ** 2)
    # A secondary LD block so it does not look like a single Gaussian.
    ld = np.maximum(ld, 0.35 * np.exp(-0.5 * ((pos - 0.41) / 0.045) ** 2))
    return pos, ld


def _track_logp(ld: np.ndarray, rng: np.random.Generator, height: float) -> np.ndarray:
    noise = rng.exponential(0.18, ld.size) * 0.55
    bg = np.clip(-np.log10(rng.uniform(1e-3, 1.0, ld.size)), 0, 2.8)
    return np.maximum(bg, height * ld + noise)


def _draw_track(ax, pos, logp, color, sv_h, ylabel, pp4, show_x=False) -> None:
    ax.scatter(pos, logp, s=5.5, c=color, linewidths=0, rasterized=True, zorder=2, alpha=0.9)
    ax.scatter([SV_X], [sv_h], marker="D", s=28, c="#111111", edgecolors="white", linewidths=0.4, zorder=5)
    ax.axvline(SV_X, color="#222222", ls=":", lw=0.6, zorder=1)
    ax.axhline(7.5, color=OKABE["vermillion"], ls="--", lw=0.45, alpha=0.7, zorder=1)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, max(sv_h, float(logp.max())) * 1.18)
    ax.set_ylabel(ylabel, fontsize=6.4)
    ax.set_yticks([0, 5, 10] if sv_h > 8 else [0, 4, 8])
    ax.text(0.99, 0.88, rf"PP$_4$={pp4:.2f}", transform=ax.transAxes, ha="right", va="top", fontsize=5.6, color="#333333")
    if show_x:
        ax.set_xticks([0, 0.25, 0.5, 0.75, 1.0])
        ax.set_xticklabels(["36.0", "36.3", "36.6", "36.9", "37.2"])
        ax.set_xlabel("chr22 position (Mb, placeholder window)")
    else:
        ax.set_xticks([])
        ax.spines["bottom"].set_visible(False)
    despine(ax)


def _gene_track(ax) -> None:
    """Simple gene models around APOL1."""
    ax.set_xlim(0, 1)
    ax.set_ylim(-0.35, 1.15)
    ax.axis("off")
    genes = [
        (0.08, 0.22, "MYH9", False),
        (0.36, 0.18, "APOL3", True),
        (0.52, 0.22, "APOL1", True),
        (0.78, 0.16, "APOL2", False),
    ]
    ax.plot([0.02, 0.98], [0.45, 0.45], color="#CCCCCC", lw=0.6, zorder=1)
    for start, width, name, on_sv in genes:
        y = 0.28 if on_sv else 0.62
        ax.add_patch(Rectangle((start, y), width, 0.28, facecolor=OKABE["blue"] if on_sv else "#888888", edgecolor="none", zorder=2))
        ax.text(start + width / 2, y + 0.42, name, ha="center", va="bottom", fontsize=5.6, color="#111111", fontstyle="italic")
    ax.scatter([SV_X], [0.08], marker="v", s=22, c="#111111", zorder=3)
    ax.text(SV_X + 0.015, 0.08, "SV", fontsize=5.4, va="center", color="#111111")


def _panel_effects(ax) -> None:
    """Same SV, four layers: direction of effect (placeholder)."""
    rows = [
        ("Haplotype\nmethylation", -0.31, 0.07, OKABE["orange"]),
        ("APOL1 RNA", 0.42, 0.08, OKABE["sky"]),
        ("APOL1 protein", 0.51, 0.07, OKABE["green"]),
        (HIGHLIGHT_NAME.replace(" ", "\n"), 0.38, 0.09, OKABE["blue"]),
    ]
    y = np.arange(len(rows))
    for i, (lab, beta, se, color) in enumerate(rows):
        ax.plot([beta - 1.96 * se, beta + 1.96 * se], [i, i], color=color, lw=1.35, zorder=2)
        ax.scatter([beta], [i], s=22, color=color, zorder=3, edgecolors="white", linewidths=0.3)
    ax.axvline(0, color="#888888", lw=0.6, ls="--")
    ax.set_yticks(y)
    ax.set_yticklabels([r[0] for r in rows], fontsize=5.8)
    ax.set_xlabel("Effect of SV (s.d.)")
    ax.set_xlim(-0.62, 0.78)
    ax.set_ylim(-0.55, len(rows) - 0.45)
    ax.set_title("Colocalizing effects", loc="left", pad=2)
    despine(ax)


def render() -> tuple[Path, Path]:
    rng = np.random.default_rng(31)
    ax_panel("figS14_phewas", 3.15, lambda ax: _panel_phewas(ax, rng), left=0.10, right=0.98, top=0.88, bottom=0.16)

    pos, ld = _shared_snps(rng)
    tracks = [
        (_track_logp(ld, rng, 9.4), OKABE["orange"], 10.6, "Methylation\n(meQTL)", 0.88),
        (_track_logp(ld, rng, 12.8), OKABE["sky"], 14.2, "RNA\n(eQTL)", 0.96),
        (_track_logp(ld, rng, 11.5), OKABE["green"], 13.1, "Proteomics\n(pQTL)", 0.93),
        (_track_logp(ld, rng, 10.2), OKABE["blue"], 12.4, "EHR\nCKD phecode", 0.91),
    ]
    fig_b = new_figure(4.85)
    gs_tr = GridSpec(
        5,
        1,
        figure=fig_b,
        left=0.16,
        right=0.98,
        top=0.90,
        bottom=0.08,
        hspace=0.08,
        height_ratios=[1.0, 1.0, 1.0, 1.0, 0.32],
    )
    for i, (logp, color, sv_h, ylab, pp4) in enumerate(tracks):
        ax = fig_b.add_subplot(gs_tr[i])
        _draw_track(ax, pos, logp, color, sv_h, ylab, pp4, show_x=(i == 3))
        if i == 0:
            ax.set_title(
                "Colocalization at APOL1  ·  same SNVs in every track, diamond is the SV",
                loc="left",
                pad=2,
            )
        if i < 3:
            ax.tick_params(labelbottom=False)
    ax_g = fig_b.add_subplot(gs_tr[4])
    _gene_track(ax_g)
    save_panel(fig_b, "figS14_coloc_tracks")
    ax_panel("figS14_effects", 3.15, _panel_effects, width=COL_IN * 1.15, left=0.22, right=0.96, top=0.88, bottom=0.12)

    fig = new_figure(8.55)
    gs = GridSpec(
        2,
        1,
        figure=fig,
        left=0.10,
        right=0.98,
        top=0.965,
        bottom=0.045,
        hspace=0.22,
        height_ratios=[1.05, 1.55],
    )
    ax_a = fig.add_subplot(gs[0, 0])
    gs_b = gs[1, 0].subgridspec(6, 2, height_ratios=[1.0, 1.0, 1.0, 1.0, 0.28, 0.02], width_ratios=[1.0, 0.38], wspace=0.08, hspace=0.08)

    _panel_phewas(ax_a, rng)

    pos, ld = _shared_snps(rng)
    tracks = [
        (_track_logp(ld, rng, 9.4), OKABE["orange"], 10.6, "Methylation\n(meQTL)", 0.88),
        (_track_logp(ld, rng, 12.8), OKABE["sky"], 14.2, "RNA\n(eQTL)", 0.96),
        (_track_logp(ld, rng, 11.5), OKABE["green"], 13.1, "Proteomics\n(pQTL)", 0.93),
        (_track_logp(ld, rng, 10.2), OKABE["blue"], 12.4, "EHR\nCKD phecode", 0.91),
    ]
    axes_tracks = []
    for i, (logp, color, sv_h, ylab, pp4) in enumerate(tracks):
        ax = fig.add_subplot(gs_b[i, 0], sharex=None)
        _draw_track(ax, pos, logp, color, sv_h, ylab, pp4, show_x=(i == 3))
        axes_tracks.append(ax)
        if i < 3:
            ax.tick_params(labelbottom=False)
    axes_tracks[0].set_title(
        "Colocalization at APOL1  ·  same SNVs in every track, diamond is the SV",
        loc="left",
        pad=2,
    )

    ax_g = fig.add_subplot(gs_b[4, 0])
    _gene_track(ax_g)

    ax_c = fig.add_subplot(gs_b[0:4, 1])
    _panel_effects(ax_c)

    # Hide the unused last row cell.
    ax_pad = fig.add_subplot(gs_b[5, :])
    ax_pad.axis("off")

    stamp_mockup(fig)
    return save(fig, "figS14_phewas_coloc")
