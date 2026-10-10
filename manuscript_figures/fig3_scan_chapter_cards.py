"""Figure 3 mock-up, layout 2: the chapter-card design.

A  Testable phecodes per ancestry (rows) and PheCodeX chapter (columns).
B  One card per chapter, in the column order of A: a mini genome-wide scan of
   the best phecode in that chapter at every SV.
C  Selected loci, as in layout 1 (shares card()).

Real: PheCodeX structure, phenotype-table denominators. Simulated: everything else.
"""

from __future__ import annotations

import numpy as np
from matplotlib.colors import PowerNorm
from matplotlib.patches import Rectangle

from . import fig3_sim as S
from .fig3_scan import (
    CARD_H, CARD_W, CMAP, GREY, SLATE_D, SLATE_L, STAMP_REAL, STAMP_SIM, SUBTITLE, add, card, key_card, letter, title, txt,
)
from .style import ANC_TEXT, FIG_W, HG38_MB, P2_ANC_N, P2_N, chrom_offsets, new_figure, save, stamp_mockup

H2 = 8.55


def share_style(cov: S.Coverage):
    """Default cell style: fill = share of the chapter's phecodes that are testable in the row."""
    norm = PowerNorm(0.5, vmin=0, vmax=0.55)

    def style(ch, a, k):
        frac = norm(k / cov.n_in_chapter[ch])
        return CMAP(frac), ("white" if frac > 0.5 else ("#B0B0B0" if k == 0 else "#111111")), f"{k}"
    return style


def panel_matrix_h(fig, cov: S.Coverage, left=1.22, top=1.34, style=None, pooled_label="All", cw=0.298, x_pop=-1.14, pooled_gap=None, count_labels=None):
    """Rows: phecodes per chapter, ancestry groups, pooled. Columns: chapters, then the total across chapters.

    ``style(chapter, row key, n testable) -> (fill, text colour, text | None)`` sets each cell; the default shades by
    the share of the chapter that is testable, and also grey-fills the chapter-count row and outlines the pooled row.
    """
    default = style is None
    style = style or share_style(cov)
    rh, gap, tw, tgap = 0.150, 0.07, 0.42, 0.05
    ncol = len(cov.chapters)
    ancs = cov.keys
    rows = ["CODES"] + ancs + ["ALL"]
    ytop, y = {}, 0.0
    for a in rows:
        if a == ancs[0]:
            y += 0.04
        if a == "ALL":
            y += gap if pooled_gap is None else pooled_gap
        ytop[a] = y
        y += rh
    total_h = y
    x_tot = cw * ncol + tgap
    total_w = x_tot + tw
    ax = add(fig, left, top, total_w, total_h)
    ax.set_xlim(0, total_w)
    ax.set_ylim(total_h, 0)
    ax.axis("off")
    group = cov.mode == "group"  # participants per group, or haplotypes per local ancestry
    for j, ch in enumerate(cov.chapters):
        ax.text(j * cw + cw / 2, -0.04, ch, rotation=90, ha="center", va="bottom", fontsize=7, clip_on=False)
        ax.add_patch(Rectangle((j * cw, ytop["CODES"]), cw, rh, fc="#F1F1F1" if default else "none", ec="white", lw=0.8))
        ax.text(j * cw + cw / 2, ytop["CODES"] + rh / 2, f"{cov.n_in_chapter[ch]}", ha="center", va="center", fontsize=7,
                color="#444444")
        for a in ancs + ["ALL"]:
            fc, tc, s = style(ch, a, cov.testable[ch][a])
            ax.add_patch(Rectangle((j * cw, ytop[a]), cw, rh, fc=fc, ec="white", lw=0.8))
            if s:
                ax.text(j * cw + cw / 2, ytop[a] + rh / 2, s, ha="center", va="center", fontsize=7, color=tc)
    # total column: all phecodes, and the testable phecodes of each row summed over chapters
    ax.text(x_tot + tw / 2, -0.04, "Total", rotation=90, ha="center", va="bottom", fontsize=7, fontweight="bold", clip_on=False)
    ax.plot([x_tot - tgap / 2] * 2, [0, total_h], color="#C8C8C8", lw=0.5)
    for a in rows:
        tot = sum(cov.n_in_chapter.values()) if a == "CODES" else sum(cov.testable[ch][a] for ch in cov.chapters)
        ax.text(x_tot + tw / 2, ytop[a] + rh / 2, f"{tot:,}", ha="center", va="center", fontsize=7,
                fontweight="bold" if a != "CODES" else "normal", color="#111111" if a != "CODES" else "#444444")
    if default:
        ax.add_patch(Rectangle((0, ytop["ALL"]), total_w, rh, fc="none", ec="#111111", lw=0.7))
    # left block, in reading order: population, then the sample counts (all samples, samples with phecodes)
    x_all, x_ph = -0.46, -0.04
    for a in rows[1:]:  # the grey chapter-count row carries no label (explained in the subtitle)
        yc = ytop[a] + rh / 2
        ax.text(x_pop, yc, pooled_label if a == "ALL" else a, ha="left", va="center", fontsize=7, fontweight="bold",
                color="#111111" if a == "ALL" else ANC_TEXT[a], clip_on=False)
        if group:
            n_all = P2_N if a == "ALL" else P2_ANC_N[a]
            n_ph = cov.n_total if a == "ALL" else cov.totals[a]
            ax.text(x_all, yc, f"{n_all:,}", ha="right", va="center", fontsize=7, color="#555555", clip_on=False)
            ax.text(x_ph, yc, f"{n_ph:,}", ha="right", va="center", fontsize=7, color="#555555", clip_on=False)
        else:
            n = 2 * cov.n_total if a == "ALL" else cov.totals[a]
            ax.text(x_ph, yc, f"{n / 1000:.1f}k", ha="right", va="center", fontsize=7, color="#555555", clip_on=False)
    labs = count_labels or ("All samples", "With phecodes")  # headers of the two sample-count columns (may span two lines)
    for x, lab in ((x_all, labs[0]), (x_ph, labs[1])) if group else ((x_ph, "Haplotypes"),):
        ax.text(x - 0.14, -0.04, lab, rotation=90, ha="center", va="bottom", multialignment="left", fontsize=7,
                color="#555555", clip_on=False)
    return total_h


def panel_chapter_cards(fig, cov, hits, left, top):
    rng = np.random.default_rng(S.SEED + 11)
    off = chrom_offsets()
    xmax = off[22] + HG38_MB[22] + 4
    ncol, nrow = 6, 3
    cwid, gap, cht = 1.115, 0.085, 0.80
    for i, ch in enumerate(cov.chapters):
        r, c = divmod(i, ncol)
        x0 = left + c * (cwid + gap)
        y0 = top + r * (cht + 0.20)
        k = cov.testable[ch]["ALL"]
        mine = [h for h in hits if h.chapter == ch]
        txt(fig, x0 + 0.02, y0, ch, fontsize=7, va="top", ha="left")
        txt(fig, x0 + cwid, y0, f"{len(mine)}", fontsize=7, va="top", ha="right", color="#555555", fontweight="bold")
        ax = add(fig, x0 + 0.17, y0 + 0.15, cwid - 0.17, cht - 0.15)
        if k > 0:
            n = 1600
            chroms = rng.choice(list(HG38_MB), size=n, p=np.array(list(HG38_MB.values())) / sum(HG38_MB.values()))
            pos = np.array([off[c_] + rng.uniform(0, HG38_MB[c_]) for c_ in chroms])
            lp = -np.log10(rng.beta(1, max(k, 1), n))
            lp = np.minimum(lp, S.GW_LOGP - 0.4)
            par = (chroms % 2).astype(bool)
            ax.scatter(pos[par], lp[par], s=1.6, color="#C3C7CA", lw=0, rasterized=True)
            ax.scatter(pos[~par], lp[~par], s=1.6, color="#9DA3A8", lw=0, rasterized=True)
            for h in mine:
                x = off[h.chrom] + h.pos_mb
                ax.scatter([x], [min(h.logp, 15.5)], s=7 if h.repeat else 6, fc=SLATE_D if h.repeat else SLATE_L,
                           ec="#16242D" if h.repeat else "#6F8A99", lw=0.3, zorder=4)
                if h.num:
                    edge = x > xmax - 220
                    ax.text(x + (-90 if edge else 90), min(h.logp, 15.5), str(h.num), fontsize=7, fontweight="bold",
                            va="center", ha="right" if edge else "left", zorder=5)
        else:
            ax.text(0.5, 0.2, "nothing testable", transform=ax.transAxes, ha="center", va="center", fontsize=7, color="#9A9A9A")
        ax.axhline(S.GW_LOGP, color=GREY, lw=0.5, ls=(0, (3, 2)), zorder=1)
        ax.set_xlim(-8, xmax + 6)
        ax.set_ylim(0, 16.5)
        ax.set_xticks([])
        ax.set_yticks([0, 8, 16])
        ax.tick_params(axis="y", length=2, pad=1)
        ax.spines["bottom"].set_linewidth(0.5)


def render():
    cov = S.build_coverage("testable")
    hits = S.build_hits(cov)
    fig = new_figure(H2)
    rng = np.random.default_rng(S.SEED + 9)
    letter(fig, "A", 0.02, 0.0)
    a_real = not cov.simulated
    title(fig, "Testable phecodes by %s and PheCodeX chapter" % ("ancestry group" if cov.mode == "group" else "local ancestry"), 0.30, 0.02,
          tag="observed counts" if a_real else "simulated", real=a_real)
    txt(fig, 0.30, 0.20, SUBTITLE[cov.mode].format(c=S.MIN_CASES, h=S.MIN_HAP)
        + "; grey row, phecodes per chapter; Total, summed across chapters",
        fontsize=7, color="#555555", va="top", ha="left")
    th = panel_matrix_h(fig, cov)
    b_top = 1.34 + th + 0.52
    letter(fig, "B", 0.02, b_top - 0.38)
    title(fig, "Best phecode per chapter at every SV; numbers mark the loci on the cards", 0.30, b_top - 0.36, tag="simulated")
    panel_chapter_cards(fig, cov, hits, 0.05, b_top)
    c_top = b_top + 3 * (0.80 + 0.20) + 0.24
    txt(fig, 3.6, c_top - 0.24, "Genome position, chr 1\u201322", fontsize=7, ha="center", va="top", color="#555555")
    letter(fig, "C", 0.02, c_top - 0.02)
    title(fig, "Selected loci", 0.30, c_top, tag="simulated")
    gap = (FIG_W - 0.10 - 4 * CARD_W) / 3
    for i, num in enumerate((1, 2, 4, None)):
        left, top = 0.05 + i * (CARD_W + gap), c_top + 0.30
        if num:
            card(fig, left, top, next(l for l in S.LOCI if l.num == num), rng)
        else:
            key_card(fig, left, top)
    stamp_mockup(fig, STAMP_REAL[cov.mode] if a_real else STAMP_SIM)
    return save(fig, "fig3_scan_chapter_cards")
