"""Figure 3, interim version: only data that exists. A panel without data is left blank (asterisk, explained in the footer).

A  Testable phecodes by participant group and PhecodeX chapter. Cell colour is the number printed in the cell.
   Observed: data/fig3_phecode_case_counts_current.tsv and fig3_cohort_summary_current.json, written by
   scripts/make_fig3_case_counts.py from the select_phenotypes.py cohort the FELIX scan is fitted on (9,775 participants).
B  Phecodes testable as a function of the case cutoff (same table). The floor used, 50 cases, is marked.
C  -log10 P of the heterogeneity test against the pooled test (P_het vs P_hom) at repeat loci, as a density heatmap.
   Needs the FelixGenome interim export interim_het_vs_pooled.tsv.gz (felix_04_interim_scan.ipynb) in data/, as
   fig3_het_vs_pooled.tsv.gz or under its own name; FIG3_HET_FILE overrides. Blank until then.
D  Squared correlation between repeat-dosage and carrier coding at length-polymorphic repeat loci
   (data/fig3_dosage_corr2_bands.json, band counts from felix_04_interim_scan.ipynb cell 8).
E  The card key only; the locus cards wait for association results.
"""

from __future__ import annotations

import csv
import gzip
import json
import os
from pathlib import Path

import matplotlib as mpl
import numpy as np
from matplotlib.colors import LinearSegmentedColormap, LogNorm, PowerNorm
from matplotlib.lines import Line2D

from . import fig3_scan_merged as M
from . import fig3_sim as S
from .fig3_scan import CARD_H, CARD_W, GREY, SLATE_D, SLATE_L, W, add, txt
from .fig3_scan_chapter_cards import panel_matrix_h
from .style import ANC_COLORS, DATADIR, new_figure, save

COUNTS = DATADIR / "fig3_phecode_case_counts_current.tsv"
SUMMARY = DATADIR / "fig3_cohort_summary_current.json"
BANDS = DATADIR / "fig3_dosage_corr2_bands.json"
HET_FILES = [DATADIR / "fig3_het_vs_pooled.tsv.gz", DATADIR / "interim_het_vs_pooled.tsv.gz"]
HET_COLS = (("p.value_ancALL", "p_pooled"), ("P_het_admixed_c", "P_het_admixed", "p_het"))  # pooled, het: accepted names

ROW_H = 1.50  # height of the axes in the B-D row (inches)
XLAB_Y = -0.135  # x-label position, axes fraction, shared by B, C and D
HOT = LinearSegmentedColormap.from_list("slate_density", ["#EEF1F4", "#A9BCC8", "#6F8A99", SLATE_D, "#16242D"])


def count_style(cov: S.Coverage):
    """Cell style for panel_matrix_h: colour is the testable-phecode count itself (square-root scale, 0 = white)."""
    vmax = max(cov.testable[ch][a] for ch in cov.chapters for a in cov.keys + ["ALL"])
    norm = PowerNorm(0.5, vmin=0, vmax=vmax, clip=True)

    def style(ch, a, k):
        if k == 0:
            return M.NOT_TESTABLE, "#B0B0B0", "0"
        fc = M.TAN_GREEN(norm(k))
        return fc, ("#FFFFFF" if M._lum(fc) < 0.30 else "#111111"), f"{k}"
    return style


# ------------------------------------------------------------------------------------------------ B: cutoff staircase
def panel_cutoff(ax, counts: dict):
    rows = [counts[r["phecode"]] for r in S.load_phecodes() if r["phecode"] in counts]  # the set panel A counts over
    thr = np.arange(1, 1001)

    def curve(vals):
        return (np.asarray(vals, float)[None, :] >= thr[:, None]).sum(axis=1)
    pooled = curve([r["n_cases"] for r in rows])
    ymax = float(np.ceil(pooled[9] * 1.03 / 250) * 250)
    ax.set_xscale("log")
    ax.set_xlim(10, 1000)
    ax.set_ylim(0, ymax)
    ax.set_xticks([10, 20, 50, 100, 200, 500, 1000])
    ax.set_xticklabels(["10", "20", "50", "100", "200", "500", "1000"])
    ax.minorticks_off()
    ax.set_yticks(np.arange(0, ymax + 1, 500))
    ax.axvline(S.MIN_CASES, color="#999999", lw=0.6, ls=(0, (2.2, 2.2)), zorder=1)
    for g in S.GROUP_ANC:
        ax.step(thr, curve([r["cases"][g] for r in rows]), where="post", color=ANC_COLORS[g], lw=0.9, zorder=3)
    ax.step(thr, pooled, where="post", color="#111111", lw=1.5, zorder=5)
    i = S.MIN_CASES - 1
    ax.plot([S.MIN_CASES], [pooled[i]], marker="o", ms=3.0, mfc="white", mec="#111111", mew=0.9, zorder=6)
    ax.annotate(f"{pooled[i]:,}", (S.MIN_CASES, pooled[i]), xytext=(5, 6), textcoords="offset points", fontsize=7,
                color="#111111", ha="left", va="center", zorder=7)
    ax.set_xlabel("Minimum cases", labelpad=2)
    ax.set_ylabel("Testable phecodes", labelpad=2)
    entries = [("Joint", "#111111", 1.5)] + [(g, ANC_COLORS[g], 0.9) for g in S.GROUP_ANC]
    ax.legend(handles=[Line2D([], [], color=c, lw=w, label=n) for n, c, w in entries], loc="upper right",
              bbox_to_anchor=(1.01, 1.02), ncol=2, fontsize=7, frameon=False, handlelength=1.3, handletextpad=0.4,
              columnspacing=0.8, labelspacing=0.25, borderaxespad=0.0)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    return pooled


# --------------------------------------------------------------------------------------- C: P_het vs P_hom heatmap
def het_path() -> Path | None:
    env = os.environ.get("FIG3_HET_FILE", "")
    for p in ([Path(env)] if env else []) + HET_FILES:
        if p.is_file():
            return p
    return None


def load_het(path: Path):
    """-log10 P arrays (pooled, heterogeneity) from the FELIX interim export; rows with a missing P are dropped."""
    opener = gzip.open if str(path).endswith(".gz") else open
    xs, ys = [], []
    with opener(path, "rt", newline="") as fh:
        rd = csv.reader(fh, delimiter="\t")
        head = next(rd)
        ix = next((head.index(c) for c in HET_COLS[0] if c in head), None)
        iy = next((head.index(c) for c in HET_COLS[1] if c in head), None)
        if ix is None or iy is None:
            raise KeyError(f"{path.name}: need a pooled column {HET_COLS[0]} and a heterogeneity column {HET_COLS[1]}; got {head}")
        for r in rd:
            try:
                a, b = float(r[ix]), float(r[iy])
            except (ValueError, IndexError):
                continue
            if a == a and b == b:
                xs.append(a)
                ys.append(b)
    nl = lambda v: -np.log10(np.clip(np.asarray(v, float), 1e-300, 1.0))  # noqa: E731
    return nl(xs), nl(ys)


def panel_het(ax, het):
    """Density heatmap of -log10 P_het (y) against -log10 P_hom (x); blank axes when there is no export yet."""
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    if het is None:
        ax.set_xticks([])
        ax.set_yticks([])
        for s in ("left", "bottom"):
            ax.spines[s].set_color("#C8C8C8")
            ax.spines[s].set_linewidth(0.5)
        return None
    x, y = het
    lim = float(min(7.5, np.ceil(max(x.max(), y.max()) * 2) / 2))
    step = 0.05
    edges = np.arange(0, lim + step / 2, step)
    h, _, _ = np.histogram2d(np.clip(x, 0, lim - 1e-9), np.clip(y, 0, lim - 1e-9), bins=[edges, edges])
    top = 10 ** int(np.ceil(np.log10(h.max())))
    mesh = np.ma.masked_less(h.T, 1)
    ax.pcolormesh(edges, edges, mesh, cmap=HOT, norm=LogNorm(1, top), rasterized=True, linewidth=0)
    ax.plot([0, lim], [0, lim], color="#8A8A8A", lw=0.6, zorder=3)
    ax.set_xlim(0, lim)
    ax.set_ylim(0, lim)
    ax.set_xticks(np.arange(0, lim + 0.1, 2))
    ax.set_yticks(np.arange(0, lim + 0.1, 2))
    ax.set_xlabel(r"$-\log_{10}P_{hom}$", labelpad=2)
    ax.set_ylabel(r"$-\log_{10}P_{het}$", labelpad=2)
    d = y - x
    share = dict(het=float((d > 1).mean()), pooled=float((d < -1).mean()), n=int(x.size))
    box = dict(boxstyle="square,pad=0.12", fc="white", ec="none", alpha=0.85)
    ax.text(0.03, 0.97, f"{share['het']:.1%}", transform=ax.transAxes, fontsize=7, ha="left", va="top", bbox=box, zorder=6)
    ax.text(0.97, 0.20, f"{share['pooled']:.1%}", transform=ax.transAxes, fontsize=7, ha="right", va="bottom", bbox=box, zorder=6)
    # colour key: a gradient strip inside the plot, lower right, with its end values
    cax = ax.inset_axes([0.50, 0.075, 0.30, 0.04])
    cax.imshow(np.linspace(0, 1, 128)[None, :], aspect="auto", cmap=HOT, extent=[0, 1, 0, 1])
    cax.set_xticks([])
    cax.set_yticks([])
    for sp in cax.spines.values():
        sp.set_linewidth(0.3)
    ax.text(0.47, 0.095, "1", transform=ax.transAxes, fontsize=7, ha="right", va="center", zorder=6)
    ax.text(0.83, 0.095, f"$10^{{{int(np.log10(top))}}}$", transform=ax.transAxes, fontsize=7, ha="left", va="center", zorder=6)
    return share


# ------------------------------------------------------------------------------- D: repeat dosage vs carrier coding
def panel_corr2(ax, rec: dict):
    bands = rec["bands"]
    tot = sum(b["n_loci"] for b in bands)
    pct = [100.0 * b["n_loci"] / tot for b in bands]
    xs = np.arange(len(bands))
    ax.bar(xs, pct, width=0.72, color=SLATE_D, lw=0, zorder=3)
    for x_, p in zip(xs, pct):
        ax.text(x_, p + 1.2, f"{p:.1f}%", fontsize=7, ha="center", va="bottom")
    ax.set_xticks(xs)
    ax.set_xticklabels([f"{b['lo']:g}\u2013{b['hi']:g}" for b in bands])
    ax.set_xlim(-0.6, len(bands) - 0.4)
    ax.set_ylim(0, 56)
    ax.set_yticks([0, 10, 20, 30, 40, 50])
    # no tick marks on categories, but the labels sit where the numeric axes of B and C put theirs
    ax.tick_params(axis="x", length=0, pad=mpl.rcParams["xtick.major.size"] + mpl.rcParams["xtick.major.pad"])
    ax.set_xlabel(r"corr$^2$(dosage, carrier)", labelpad=2)
    ax.set_ylabel("Repeat loci (%)", labelpad=2)
    ax.text(0.02, 0.98, f"{tot / 1e6:.2f} M loci", transform=ax.transAxes, fontsize=7, ha="left", va="top", color="#555555")
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    return pct


# ------------------------------------------------------------------------------------------------------------ render
def render():
    cov = S.build_coverage("testable", counts_path=COUNTS, summary_path=SUMMARY)
    assert not cov.simulated and cov.mode == "group"
    assert cov.n_total == int(sum(cov.totals.values())), (cov.n_total, cov.totals)
    pooled = [cov.testable[ch]["ALL"] for ch in cov.chapters]
    assert all(a >= b for a, b in zip(pooled, pooled[1:])), "chapter columns must be sorted by the Joint cell"
    counts = S.load_case_counts(COUNTS)
    corr2 = json.loads(BANDS.read_text())
    het_file = het_path()
    het = load_het(het_file) if het_file else None
    pending = [] if het is not None else ["panel C awaits the P_het export (interim_het_vs_pooled.tsv.gz)"]
    pending.append("panel E awaits the locus cards")

    th_expect = 0.150 * (len(cov.keys) + 2) + 0.04
    b_top = M.A_TOP + th_expect + 0.62
    f_top = b_top + ROW_H + 0.75
    fig = new_figure(f_top + CARD_H + 0.23)

    # ---- A: the table, coloured by its own counts (no legend: the number in each cell is the scale)
    M._head(fig, "A", "Testable phecodes by ancestry group and PhecodeX chapter", 0.02, M.A_HEAD_Y)
    th = panel_matrix_h(fig, cov, left=1.30, top=M.A_TOP, style=count_style(cov), pooled_label="Joint", cw=0.294,
                        x_pop=-1.22, pooled_gap=0.0, count_labels=("Participants with\nPacBio data", "In analysis"))
    assert abs(th - th_expect) < 1e-9, (th, th_expect)

    # ---- B-D: how the testing works
    sw = (W - 0.10) / 3
    ax_l, ax_w = 0.55, sw - 0.55 - 0.12
    head_y = b_top - 0.27
    M._head(fig, "B", "Phecodes testable at each case cutoff", 0.02, head_y)
    M._head(fig, "C", "Heterogeneity vs pooled test" + ("" if het is not None else "*"), 0.02 + sw, head_y)
    M._head(fig, "D", "Repeat dosage vs carrier coding", 0.02 + 2 * sw, head_y)
    ax_b = add(fig, 0.05 + ax_l, b_top, ax_w, ROW_H)
    ax_c = add(fig, 0.05 + sw + ax_l, b_top, ROW_H, ROW_H)  # square, so the diagonal is at 45 degrees
    ax_d = add(fig, 0.05 + 2 * sw + ax_l, b_top, ax_w, ROW_H)
    pooled_curve = panel_cutoff(ax_b, counts)
    share = panel_het(ax_c, het)
    pct = panel_corr2(ax_d, corr2)
    for a_ in (ax_b, ax_c, ax_d):
        a_.xaxis.set_label_coords(0.5, XLAB_Y)

    # ---- E: locus cards suppressed; the key stays
    M._head(fig, "E", "Representative loci: regional association and effect in each ancestry*", 0.02, f_top - 0.24)
    gap = (W - 0.10 - 4 * CARD_W) / 3
    M.key_card_annotated(fig, 0.05 + 3 * (CARD_W + gap), f_top)
    fig.text(0.995, 0.004, "*Not populated yet: " + "; ".join(pending) + ".", ha="right", va="bottom", fontsize=5.5,
             color="#8A8A8A", style="italic", transform=fig.transFigure)
    render.info = dict(pooled_at_floor=int(pooled_curve[S.MIN_CASES - 1]), corr2_pct=pct, het_share=share,
                       het_file=str(het_file) if het_file else None, n_analysis=cov.n_total)
    return save(fig, "fig3_interim")
