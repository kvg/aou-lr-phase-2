"""Figure 3 mock-up, layout 4 (results first): the table, a genome-wide best-phecode scan, and the loci cards.

A  Testable phecodes by ancestry group and PhecodeX chapter (observed counts; cell colour = simulated number of
   genome-wide-significant associations), as in fig3_scan_merged.
B  Best phecode association at every SV (simulated). The "how" is encoded as what each feature adds: light = significant
   with standard testing (carrier coding, ancestry-agnostic); dark = significant only with repeat-dosage coding (a repeat
   locus whose biallelic test is not); ring = significant only with the ancestry-aware test (the pooled test is not).
C  Three representative loci and the annotated card key.
The calibration and method-comparison panels moved to figSx_testing_validation.

Real: PhecodeX structure, case counts, cohort sizes. Simulated: cell colours, the scan, the cards.
"""

from __future__ import annotations

import dataclasses

import numpy as np

from . import fig3_scan_merged as M
from . import fig3_sim as S
from .fig3_scan import CARD_H, CARD_W, GREY, SLATE_D, SLATE_L, W, add, card, txt
from .fig3_scan_chapter_cards import panel_matrix_h
from .style import HG38_MB, chrom_offsets, new_figure, save, stamp_mockup

B_H = 1.10  # height of the scan's axes (inches)
EDGE_L, EDGE_D = "#6F8A99", "#16242D"
GAIN_DOSAGE = {1, 4}  # card numbers (full set) whose repeat-dosage test is significant and biallelic test is not (validation panel B)
GAIN_ANC = {2, 6}  # ... whose ancestry-aware test is significant and pooled test is not (validation panel C)
NUM_SIDE = {1: "r", 2: "a", 3: "a"}  # which side of its mark each card number sits on (3 has neighbours both sides)
N_PHECODES = 943  # phecodes tested per SV; sets the null distribution of the best P


def scan_loci(hits):
    """One mark per SV locus: its best phecode (largest -log10 P), with its repeat flag and card number."""
    d = {}
    for h in hits:
        e = d.setdefault((h.chrom, round(h.pos_mb, 2)), dict(chrom=h.chrom, pos=h.pos_mb, logp=-1.0, repeat=False, num=0))
        e["logp"] = max(e["logp"], h.logp)
        e["repeat"] |= h.repeat
        e["num"] = e["num"] or h.num
    loci = sorted(d.values(), key=lambda e: (e["chrom"], e["pos"]))
    rng = np.random.default_rng(S.SEED + 22)
    for e in loci:
        if e["num"]:
            e["dosage"], e["anc"] = e["num"] in GAIN_DOSAGE, e["num"] in GAIN_ANC
        else:
            e["dosage"] = bool(e["repeat"] and rng.random() < 0.6)
            e["anc"] = bool((not e["repeat"]) and rng.random() < 0.22)
    return loci


def panel_scan(fig, left, top, w, h, hits):
    ax = add(fig, left, top, w, h)
    off = chrom_offsets()
    xmax = off[22] + HG38_MB[22] + 4
    rng = np.random.default_rng(S.SEED + 21)
    for c in HG38_MB:
        if c % 2 == 1:
            ax.axvspan(off[c] - 4, off[c] + HG38_MB[c] + 4, color="#F1F3F5", lw=0, zorder=0)
    # background: best -log10 P of SVs with no signal = minimum of N_PHECODES null p-values (thinned for display)
    chroms = np.array(list(HG38_MB))
    wts = np.array([HG38_MB[c] for c in chroms], float)
    ch = rng.choice(chroms, size=3500, p=wts / wts.sum())
    xb = np.array([off[int(c)] + rng.uniform(1, HG38_MB[int(c)] - 1) for c in ch])
    yb = -np.log10(1 - (1 - rng.random(ch.size)) ** (1 / N_PHECODES))
    keep = yb < 6.8
    ax.scatter(xb[keep], yb[keep], s=2.2, color="#C3C7CA", lw=0, zorder=2, rasterized=True)
    loci = scan_loci(hits)
    for e in loci:  # flanking SVs in LD with the lead give each peak its tower
        x0 = off[e["chrom"]] + e["pos"]
        k = int(2 + max(e["logp"] - 9, 0))
        ys = np.maximum(e["logp"] - rng.exponential(1.8, k), 2.5)
        ax.scatter(x0 + rng.uniform(-0.8, 0.8, k) * 3, np.minimum(ys, 16), s=2.2, color="#C3C7CA", lw=0, zorder=2)
    u_in = xmax / w
    for e in loci:
        x, y = off[e["chrom"]] + e["pos"], min(e["logp"], 16)
        fc, ec = (SLATE_D, EDGE_D) if e["dosage"] else (SLATE_L, EDGE_L)
        ax.scatter([x], [y], s=16, fc=fc, ec=ec, lw=0.4, zorder=4, clip_on=False)
        if e["anc"]:
            ax.scatter([x], [y], s=52, fc="none", ec=SLATE_D, lw=0.9, zorder=5, clip_on=False)
        if e["num"] in M.SHOWN:
            n = M.SHOWN[e["num"]]
            dx = (0.10 if e["anc"] else 0.06) * u_in
            px, py, ha = {"r": (x + dx, y, "left"), "l": (x - dx, y, "right"), "a": (x, y + 1.7, "center")}[NUM_SIDE.get(n, "r")]
            ax.text(px, py, str(n), fontsize=7, fontweight="bold", va="center", ha=ha, zorder=6,
                    bbox=dict(boxstyle="square,pad=0.05", fc="white", ec="none", alpha=0.75))
    # genome-wide threshold (5e-8, stated in the caption): dark and heavier than a hairline, above the null cloud and
    # below the hit marks
    ax.axhline(S.GW_LOGP, color="#4A4A4A", lw=0.8, ls=(0, (4, 2.5)), zorder=3)
    ax.set_xlim(-8, xmax + 6)
    ax.set_ylim(0, 16.5)
    ax.set_yticks([0, 5, 10, 15])
    # every chromosome is labelled on one row. Labels sit under their ticks; where two would come within 0.03 in (only
    # 21 and 22) they are nudged apart, never more than 0.02 in from the tick
    ax.set_xticks([off[c] + HG38_MB[c] / 2 for c in HG38_MB])
    ax.set_xticklabels([])
    ax.tick_params(axis="x", length=2.5, width=0.5, pad=1.5)
    ax.tick_params(axis="y", length=2.5, pad=1.5)
    x_lo, x_hi = -8, xmax + 6
    cs = list(HG38_MB)
    ctr = np.array([(off[k] + HG38_MB[k] / 2 - x_lo) / (x_hi - x_lo) * w for k in cs])  # tick positions, inches from axes left
    half = np.array([0.027 * len(str(k)) for k in cs])  # a 7 pt digit is about 0.054 in wide
    pos = ctr.copy()
    for _ in range(200):
        moved = False
        for i in range(1, len(cs)):
            deficit = 0.03 - ((pos[i] - half[i]) - (pos[i - 1] + half[i - 1]))
            if deficit > 1e-6:
                pos[i - 1] -= deficit / 2
                pos[i] += deficit / 2
                moved = True
        pos = np.clip(pos, ctr - 0.02, ctr + 0.02)
        if not moved:
            break
    clear = (pos[1:] - half[1:]) - (pos[:-1] + half[:-1])
    assert clear.min() > 0.025, ("chromosome labels still touch", clear.min())
    for k, xc, p, c0 in zip(cs, [off[k] + HG38_MB[k] / 2 for k in cs], pos, ctr):
        ax.annotate(str(k), xy=(xc, 0), xytext=((p - c0) * 72, -4), textcoords="offset points", ha="center", va="top",
                    fontsize=7, annotation_clip=False)
    ax.annotate("Chromosome", xy=(0.5, 0), xycoords="axes fraction", xytext=(0, -4 - 8.6), textcoords="offset points",
                ha="center", va="top", fontsize=7)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.set_ylabel("Best $-\\log_{10}P$", labelpad=1)
    return loci


def scan_legend(fig, x_right, y_base):
    """Marker legend on the panel title's baseline, right-aligned: what each state of a mark means."""
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    if not hasattr(fig.canvas, "get_renderer"):
        FigureCanvasAgg(fig)
    rd = fig.canvas.get_renderer()
    fw, fh = fig.get_size_inches()

    def width(s):
        t = txt(fig, 0, 0, s, fontsize=7)
        wd = t.get_window_extent(rd).width / fig.dpi
        t.remove()
        return wd
    items = [("Standard testing", "light"), ("Dosage coding only", "dark"), ("Ancestry-aware only", "ring")]
    mk, pad, gap = 0.13, 0.05, 0.17
    total = sum(mk + pad + width(s) for s, _ in items) + gap * (len(items) - 1)
    x, yc = x_right - total, y_base - 0.035
    x0, y0, ht = x - 0.08, y_base - 0.16, 0.24  # overlay only as big as the legend row, so the saved crop is not widened
    ov = add(fig, x0, y0, x_right - x0 + 0.02, ht, zorder=20)
    ov.set_xlim(x0, x_right + 0.02)
    ov.set_ylim(y0 + ht, y0)  # data units are figure inches
    ov.axis("off")
    for s, kind in items:
        xc = x + mk / 2
        fc, ec = (SLATE_D, EDGE_D) if kind == "dark" else (SLATE_L, EDGE_L)
        ov.scatter([xc], [yc], s=16, fc=fc, ec=ec, lw=0.4)
        if kind == "ring":
            ov.scatter([xc], [yc], s=52, fc="none", ec=SLATE_D, lw=0.9)
        txt(fig, x + mk + pad, y_base, s, fontsize=7, color="#555555", va="baseline", ha="left")
        x += mk + pad + width(s) + gap


def render():
    cov = S.build_coverage("testable")
    pooled = [cov.testable[ch]["ALL"] for ch in cov.chapters]
    assert all(a >= b for a, b in zip(pooled, pooled[1:])), "chapter columns must be sorted by the Pooled cell"
    assert cov.mode == "group", "significance colours are defined for the 7-group view"
    hits = S.build_hits(cov)
    sig = M.sig_counts(cov, hits)
    rng = np.random.default_rng(S.SEED + 9)

    th_expect = 0.150 * (len(cov.keys) + 2) + 0.04
    b_top = M.A_TOP + th_expect + 0.62
    f_top = b_top + B_H + 0.68
    fig = new_figure(f_top + CARD_H + 0.23)

    # ---- A: table coloured by significance (as in fig3_scan_merged)
    M._head(fig, "A", "Testable phecodes by ancestry group and PhecodeX chapter*", 0.02, M.A_HEAD_Y)
    th = panel_matrix_h(fig, cov, left=1.30, top=M.A_TOP, style=M.sig_style(cov, sig), pooled_label="Pooled", cw=0.294,
                        x_pop=-1.22, pooled_gap=0.0, count_labels=("Participants with\nPacBio data", "In analysis"))
    assert abs(th - th_expect) < 1e-9, (th, th_expect)
    x_right = fig.axes[0].get_position().x1 * fig.get_size_inches()[0]  # the table's right edge
    M._legend_strip(fig, x_right, M.A_HEAD_Y, "Significant associations")

    # ---- B: best phecode at every SV, with the testing features shown as what they add
    M._head(fig, "B", "Best phecode association at every SV*", 0.02, b_top - 0.27)
    scan_legend(fig, x_right, b_top - 0.27)
    panel_scan(fig, 0.30, b_top, x_right - 0.30, B_H, hits)

    # ---- C: three loci and the card key
    M._head(fig, "C", "Representative loci: regional association and effect in each ancestry*", 0.02, f_top - 0.24)
    gap = (W - 0.10 - 4 * CARD_W) / 3
    by_num = {l.num: l for l in S.LOCI}
    for i, (orig, new) in enumerate(M.SHOWN.items()):
        card(fig, 0.05 + i * (CARD_W + gap), f_top, dataclasses.replace(by_num[orig], num=new), rng)
    M.key_card_annotated(fig, 0.05 + len(M.SHOWN) * (CARD_W + gap), f_top)
    stamp_mockup(fig, "*Not real data yet: panel A shows observed case counts with simulated colours; panels B and C are simulated.")
    return save(fig, "fig3_scan_results")
