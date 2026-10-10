"""Figure 3 mock-up, layout 3: layout 2's table, coloured by significance, replacing the hit map.

A  Testable phecodes (observed counts) by ancestry group and PhecodeX chapter, as in layout 2. Cell colour = number of
   genome-wide-significant SV-phecode associations in the chapter (pooled test for the Pooled row, within-group test
   for the groups; simulated); white = nothing testable, light grey = testable without a significant association.
   This folds layout 1's panel B (hit map) into the table.
B  Calibration of the FELIX headline test on null (GRM-simulated) phenotypes, one QQ curve per case-count stratum.
C-D  Dosage vs biallelic, ancestry-aware vs pooled (layout 1's D-E).
E  Three representative loci and the card key, in one row (layout 1's F, cut from seven loci to three).

Real: PhecodeX structure, case counts, cohort sizes. Simulated: the cell colours and everything below the table.
"""

from __future__ import annotations

import dataclasses

import numpy as np
from matplotlib.colors import LinearSegmentedColormap, PowerNorm
from matplotlib.patches import FancyBboxPatch

from . import fig3_sim as S
from .fig3_scan import CARD_H, CARD_W, GREY, SLATE_D, W, add, card, panel_ancestry, panel_dosage, txt
from .fig3_scan_chapter_cards import panel_matrix_h
from .style import ANC_COLORS, ANC_TEXT, new_figure, save, stamp_mockup

A_TOP = 1.08  # matrix top (inches): room for the heading above the rotated chapter names
R2_TOP_REF, H_REF = 3.98, 9.26  # layout 1 geometry (row 2 top, figure height)
P_SPECIFIC = 0.30  # share of pooled signals that are driven by a single group (placeholder)
SPECIFIC_LOCI = {2: "AFR"}  # APOL1 G1/G2 is AFR-driven
Z_GW = S._z_from_logp(S.GW_LOGP)

# cell fill by number of significant associations: tan (none) -> green (many), continuous on a square-root scale;
# counts of CNT_MAX or more share the darkest green. White = nothing testable.
CNT_MAX = 8
TAN_GREEN = LinearSegmentedColormap.from_list("tan_green", ["#EFE5CC", "#CFD39F", "#97B672", "#538F62", "#1F5A3E"])
NORM = PowerNorm(0.5, vmin=0, vmax=CNT_MAX, clip=True)
NOT_TESTABLE = "#FFFFFF"


def _lum(rgb) -> float:
    lin = [x / 12.92 if x <= 0.04045 else ((x + 0.055) / 1.055) ** 2.4 for x in rgb[:3]]
    return 0.2126 * lin[0] + 0.7152 * lin[1] + 0.0722 * lin[2]


def sig_style(cov: S.Coverage, sig: dict):
    """Cell style for panel_matrix_h: colour = significant-association count; white = nothing testable."""
    def style(ch, a, k):
        n = sig.get((ch, a), 0)
        assert k > 0 or n == 0, (ch, a, k, n)
        if k == 0:
            return NOT_TESTABLE, "#B0B0B0", "0"
        fc = TAN_GREEN(NORM(n))
        return fc, ("#FFFFFF" if _lum(fc) < 0.30 else "#555555" if n == 0 else "#111111"), f"{k}"
    return style


def _legend_strip(fig, x_right, y, lead, w=1.7, h=0.075):
    """Lead text, then a gradient strip with the end annotations either side (position = association count).

    ``y`` is a shared baseline (inches from the top): the three texts sit on it and the strip's bottom edge rests on
    it, the strip being as tall as the capitals of the title (0.075 in). The whole element ends at ``x_right``
    (inches), so it can be right-aligned to the table above.
    """
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    if not hasattr(fig.canvas, "get_renderer"):
        FigureCanvasAgg(fig)
    rd = fig.canvas.get_renderer()

    def width(s):
        t = txt(fig, 0, 0, s, fontsize=7)
        wd = t.get_window_extent(rd).width / fig.dpi
        t.remove()
        return wd
    parts = [(lead, width(lead)), ("0", width("0")), (None, w), (f"{CNT_MAX}+", width(f"{CNT_MAX}+"))]
    pads = [0.12, 0.06, 0.06]
    x = x_right - sum(wd for _, wd in parts) - sum(pads)
    for i, (s, wd) in enumerate(parts):
        if s is None:
            ax = add(fig, x, y - h, w, h)
            ax.imshow(NORM(np.linspace(0, CNT_MAX, 256))[None, :], cmap=TAN_GREEN, vmin=0, vmax=1, aspect="auto")
            ax.axis("off")
        else:
            txt(fig, x, y, s, fontsize=7, color="#555555", va="baseline", ha="left")
        x += wd + (pads[i] if i < len(pads) else 0)


# ---- panel B: calibration on null phenotypes, one curve per case-count stratum (all placeholders)
STRATA = [("50\u201399", "#C5D6E0"), ("100\u2013199", "#93AEBF"), ("200\u2013499", "#557A91"), ("\u2265500", SLATE_D)]
N_NULL = 3_000_000  # tests per stratum in the mock (null phenotypes x SVs)


def sim_null_qq(seed_offset: int = 6):  # seed chosen so no placeholder point lands on the 5e-8 line
    """Per stratum: sorted null p-values and the median-based lambda_GC. Placeholder: every stratum is calibrated."""
    rng = np.random.default_rng(S.SEED + seed_offset)
    out = []
    for name, col in STRATA:
        z2 = rng.chisquare(1, N_NULL)  # null: no inflation built in, so lambda_GC differs only by sampling noise
        p = np.sort(S.chi2_sf_1df(z2))
        out.append((name, col, p, float(np.median(z2) / S.CHI2_MEDIAN_1DF)))
    return out


def panel_null_qq(ax):
    for name, col, p, lam in sim_null_qq():
        n = p.size
        exp = -np.log10(np.arange(1, n + 1) / (n + 1))
        obs = -np.log10(p)
        keep = np.unique(np.concatenate([np.arange(0, 400), np.unique(np.geomspace(400, n, 320).astype(int)) - 1]))
        keep = keep[keep < n]
        ax.scatter(exp[keep], obs[keep], s=2.6, color=col, lw=0, label=f"{name}  {lam:.2f}")
    ax.plot([0, 8], [0, 8], color="#444444", lw=0.6, zorder=1)
    ax.axhline(S.GW_LOGP, color=GREY, lw=0.6, ls=(0, (3, 2)), zorder=1)
    ax.set_xlim(0, 7.2)
    ax.set_ylim(0, 8.4)
    ax.set_xlabel("Expected $-\\log_{10}P$", labelpad=1)
    ax.set_ylabel("Observed $-\\log_{10}P$", labelpad=1)
    ax.legend(loc="lower right", fontsize=7, handletextpad=0.1, markerscale=2.4, labelspacing=0.25, borderaxespad=0.1,
              title="Cases  $\\lambda_{GC}$", alignment="left", title_fontsize=7)


# ---- panel E: the loci shown (full-set card number -> number shown); C and D mark only these
SHOWN = {1: 1, 2: 2, 4: 3}


def sig_counts(cov: S.Coverage, hits: list[S.Hit]) -> dict:
    """Simulated number of significant associations per (chapter, row). All = pooled hits (as in layout 1's hit map).

    A group row counts a hit when the phecode is testable in the group (>= 50 cases) and the within-group z clears
    5e-8. Shared effects dilute to z_pooled * sqrt(group share); an effect confined to one group reaches
    z_pooled / sqrt(group share) there and is null elsewhere.
    """
    rng = np.random.default_rng(S.SEED + 11)
    share = {g: cov.totals[g] / cov.n_total for g in cov.keys}
    out = {(c, k): 0 for c in cov.chapters for k in cov.keys + ["ALL"]}
    for h in hits:
        out[(h.chapter, "ALL")] += 1
        z_pool = S._z_from_logp(h.logp)
        t_all = max(cov.testable[h.chapter]["ALL"], 1)
        p_t = {g: cov.testable[h.chapter][g] / t_all for g in cov.keys}
        spec = SPECIFIC_LOCI.get(h.num)
        if spec is None and rng.random() < P_SPECIFIC:
            w = np.array([share[g] * p_t[g] for g in cov.keys])
            spec = cov.keys[int(rng.choice(len(w), p=w / w.sum()))] if w.sum() > 0 else None
        for g in cov.keys:
            if spec == g:
                z_g = z_pool / np.sqrt(share[g])
            elif spec is not None:
                z_g = rng.normal(0, 1)
            elif rng.random() < p_t[g]:
                z_g = z_pool * np.sqrt(share[g]) + rng.normal(0, 0.8)
            else:
                continue
            out[(h.chapter, g)] += int(abs(z_g) > Z_GW)
    return out


HEAD_DX = 0.28  # letter-to-title offset, the same for every panel
A_HEAD_Y = 0.15  # baseline of panel A's title and of the colour legend (in from the top)
LETTER_DY = 0.010  # in: the 9 pt bold letter sits this far below the title baseline so its cap-height centre meets the title's


def _head(fig, ch, text, x, y, tags=()):
    """Panel letter, title and italic tags on one line (y = title baseline); letter centred on the title cap height."""
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    if not hasattr(fig.canvas, "get_renderer"):
        FigureCanvasAgg(fig)
    rd = fig.canvas.get_renderer()
    txt(fig, x, y + LETTER_DY, ch, fontsize=9, fontweight="bold", ha="left", va="baseline")
    t = txt(fig, x + HEAD_DX, y, text, fontsize=7.5, ha="left", va="baseline")
    x_end = x + HEAD_DX + t.get_window_extent(rd).width / fig.dpi
    for s, real in tags:
        tt = txt(fig, x_end + 0.09, y, s, fontsize=7, style="italic", color="#1E6E8C" if real else "#8A8A8A",
                 ha="left", va="baseline")
        x_end += 0.09 + tt.get_window_extent(rd).width / fig.dpi


def key_card_annotated(fig, left, top):
    """The card key, drawn as an example card (same geometry as card()) with in-place labels and leader lines.

    Placeholders stand where a real card has text ("n", "Locus", "SV class", "Phecode", effect unit, P_het). Labels
    carry a leader line only where there is room. The three cards shown all have an SV lead, so no SNV-lead marker.
    """
    pl, pw = left + 0.27, CARD_W - 0.31
    R = CARD_W - 0.03  # right text margin, inside the dashed border
    ov = add(fig, left, top, CARD_W, CARD_H, zorder=10)  # overlay in inches from the card's top-left (y downward)
    ov.set_xlim(0, CARD_W)
    ov.set_ylim(CARD_H, 0)
    ov.axis("off")
    ov.add_patch(FancyBboxPatch((-0.03, -0.04), CARD_W + 0.03, CARD_H + 0.06, boxstyle="round,pad=0,rounding_size=0.04",
                                fc="none", ec="#C4C4C4", lw=0.5, ls=(0, (3, 2)), clip_on=False))
    ov.text(0.0, 0.02, "n", fontsize=7.5, fontweight="bold", va="top", ha="left")
    ov.text(0.12, 0.02, "Locus", fontsize=7.5, style="italic", va="top", ha="left")
    ov.text(R, 0.03, "SV class", fontsize=7, color="#666666", va="top", ha="right")
    ov.text(0.12, 0.165, "Phecode", fontsize=7, color="#333333", va="top", ha="left")
    # the word LEGEND, vertical in the gutter where a real card has its y tick labels: tracked, all caps, grey,
    # one letter at a time so the spacing is exact; centred on the regional plot and strip (y 0.31-1.30 in)
    for i, ch in enumerate("LEGEND"):
        ov.text(0.12, 1.05 - 0.10 * i, ch, rotation=90, fontsize=7, fontweight="bold", color="#8A8A8A", ha="center",
                va="center")

    # ---- regional plot (schematic): axes 0.40 in tall; data x in [-0.52, 0.52], y in [0, 1]
    rng = np.random.default_rng(7)
    x = np.concatenate([rng.uniform(-0.5, 0.5, 70), rng.normal(0, 0.09, 45)])
    x = x[np.abs(x) < 0.5]
    y = 0.03 + 0.07 * rng.random(x.size) + 0.58 * np.exp(-0.5 * (x / 0.10) ** 2) * (0.55 + 0.45 * rng.random(x.size))
    x = np.append(x, 0.19)  # the SNV the "SNVs" leader points to
    y = np.append(y, 0.27)
    ax = add(fig, pl, top + 0.31, pw, 0.40)
    ax.scatter(x, y, s=3.2, color="#B9BEC2", lw=0, zorder=2)
    ax.axhline(0.40, color=GREY, lw=0.5, ls=(0, (3, 2)), zorder=1)
    ax.scatter([0], [0.86], marker="D", s=26, fc=SLATE_D, ec="white", lw=0.4, zorder=5, clip_on=False)
    ax.set_xlim(-0.52, 0.52)
    ax.set_ylim(0, 1)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.spines["bottom"].set_linewidth(0.5)

    def rx(v):  # regional data x -> inches from the card's left
        return 0.27 + (v + 0.52) / 1.04 * pw

    def ry(v):  # regional data y -> inches from the card's top
        return 0.31 + (1 - v) * 0.40
    lead = dict(color="#666666", lw=0.5, solid_capstyle="butt")
    ov.text(rx(-0.14), ry(0.86), "Lead SV", fontsize=7, ha="right", va="center")
    ov.plot([rx(-0.125), rx(-0.045)], [ry(0.86)] * 2, **lead)
    ov.text(rx(-0.50), ry(0.40) - 0.012, r"5$\times$10$^{-8}$", fontsize=7, ha="left", va="bottom")
    ov.text(0.27 + 0.97 * pw, 0.31 + 0.03 * 0.40, "lead P", fontsize=7, ha="right", va="top")
    ov.text(rx(0.50), ry(0.25), "SNVs", fontsize=7, ha="right", va="center")
    ov.plot([rx(0.50) - 0.30, rx(0.19) + 0.03], [ry(0.25), ry(0.27)], **lead)

    # ---- effect strip (schematic): shared effect, then the seven groups in fixed order; OTH drawn as <50 cases
    akeys = S.ancestry_keys()
    xs = [0.0] + [1.7 + i for i in range(len(akeys))]
    vals = {"ALL": (0.48, 0.12), "AFR": (0.30, 0.12), "AMR": (0.55, 0.16), "EAS": (0.22, 0.10), "EUR": (0.42, 0.12),
            "MID": (0.40, 0.14), "SAS": (0.50, 0.18), "OTH": None}
    lo, hi, st_top, st_h = -0.95, 1.40, 0.75, 0.50
    ax2 = add(fig, pl, top + st_top, pw, st_h)
    for key, xx in zip(["ALL"] + akeys, xs):
        colr = SLATE_D if key == "ALL" else ANC_COLORS[key]
        if vals[key] is None:
            ax2.scatter([xx], [0], s=9, fc="white", ec="#9A9A9A", lw=0.6, zorder=3)
            continue
        b, ci = vals[key]
        ax2.plot([xx, xx], [0, b], color=colr, lw=1.0, zorder=2, solid_capstyle="butt")
        ax2.plot([xx, xx], [b - ci, b + ci], color=colr, lw=0.7, alpha=0.65, zorder=2)
        ax2.scatter([xx], [b], s=11, color=colr, lw=0, zorder=4)
    ax2.axhline(0, color="#333333", lw=0.5, zorder=1)
    ax2.plot([0.85, 0.85], [lo, 0.80], color="#D6D6D6", lw=0.5, zorder=0)  # All | groups, below the label band
    ax2.set_xlim(-0.6, xs[-1] + 0.6)
    ax2.set_ylim(lo, hi)
    ax2.set_xticks([])
    ax2.set_yticks([])
    ax2.spines["bottom"].set_visible(False)

    def sx(v):  # strip data x -> inches from the card's left
        return 0.27 + (v + 0.6) / (xs[-1] + 1.2) * pw

    def sy(v):  # strip data y -> inches from the card's top
        return st_top + (hi - v) / (hi - lo) * st_h
    band = sy(1.15)
    ov.text(0.30, band, "All = shared", fontsize=7, ha="left", va="center")
    ov.plot([sx(0), sx(0)], [band + 0.05, sy(vals["ALL"][0] + vals["ALL"][1]) - 0.02], **lead)
    ov.text(R - 0.03, band, "<50 cases", fontsize=7, ha="right", va="center")
    ov.plot([sx(xs[-1])] * 2, [band + 0.05, sy(0) - 0.045], **lead)
    for key, xx in zip(["ALL"] + akeys, xs):  # group names under their lollipops, in the card's colours
        ov.text(sx(xx), sy(0) + 0.045, "All" if key == "ALL" else key, rotation=90, fontsize=7, ha="center", va="top",
                color=SLATE_D if key == "ALL" else ANC_TEXT[key], fontweight="bold" if key == "ALL" else "normal")
    ov.text(0.12, 1.33, "Effect, 95% CI", fontsize=7, color="#666666", va="top", ha="left")
    ov.text(R, 1.33, r"$P_{het}$", fontsize=7, color="#666666", va="top", ha="right")


def render():
    cov = S.build_coverage("testable")  # chapters ordered by the Pooled cell, largest first (ties: alphabetical)
    pooled = [cov.testable[ch]["ALL"] for ch in cov.chapters]
    assert all(a >= b for a, b in zip(pooled, pooled[1:])), "chapter columns must be sorted by the Pooled cell"
    assert cov.mode == "group", "significance colours are defined for the 7-group view"
    hits = S.build_hits(cov)
    sig = sig_counts(cov, hits)
    rng = np.random.default_rng(S.SEED + 9)

    # figure height: matrix rows (+ the two gaps) set where B-D start; everything below follows layout 1
    th_expect = 0.150 * (len(cov.keys) + 2) + 0.04  # one gap, under the chapter-count row; Pooled is an ordinary row
    r2_top, r2_h = A_TOP + th_expect + 0.62, 1.35
    fig = new_figure(H_REF - (R2_TOP_REF - r2_top) - (CARD_H + 0.10))  # layout 1's bottom margin, one card row instead of two

    # ---- A: table coloured by significance
    _head(fig, "A", "Testable phecodes by ancestry group and PhecodeX chapter*", 0.02, A_HEAD_Y)  # * = not real data yet
    th = panel_matrix_h(fig, cov, left=1.30, top=A_TOP, style=sig_style(cov, sig), pooled_label="Pooled", cw=0.294,
                       x_pop=-1.22, pooled_gap=0.0,
                       count_labels=("Participants with\nPacBio data", "In analysis"))
    assert abs(th - th_expect) < 1e-9, (th, th_expect)
    ax_tab = fig.axes[0]  # the matrix axes: its right edge is the table's right edge
    x_right = ax_tab.get_position().x1 * fig.get_size_inches()[0]
    _legend_strip(fig, x_right, A_HEAD_Y, "Significant associations")  # on the title's baseline, right-aligned to the table

    # ---- B-D: calibration and method comparisons (layout 1's C-E)
    for ch, ttl, fn, x_title, w in (  # each y-axis sits at the left edge of its panel title
        ("B", "Calibration on null phenotypes", panel_null_qq, 0.30, 2.08),
        ("C", "Dosage vs biallelic coding", lambda ax: panel_dosage(ax, SHOWN), 2.74, 1.93),
        ("D", "Ancestry-aware vs pooled", lambda ax: panel_ancestry(ax, SHOWN), 5.04, 1.93),
    ):
        ax = add(fig, x_title, r2_top, w, r2_h)
        fn(ax)
        _head(fig, ch, ttl + "*", x_title - HEAD_DX, r2_top - 0.27)

    # ---- E: loci
    f_top = r2_top + r2_h + 0.72
    _head(fig, "E", "Representative loci: regional association and effect in each ancestry*", 0.02, f_top - 0.24)
    gap = (W - 0.10 - 4 * CARD_W) / 3
    by_num = {l.num: l for l in S.LOCI}
    for i, (orig, new) in enumerate(SHOWN.items()):
        card(fig, 0.05 + i * (CARD_W + gap), f_top, dataclasses.replace(by_num[orig], num=new), rng)
    key_card_annotated(fig, 0.05 + len(SHOWN) * (CARD_W + gap), f_top)
    stamp_mockup(fig, "*Not real data yet: panel A shows observed case counts with simulated colours; panels B\u2013E are simulated.")
    return save(fig, "fig3_scan_merged")
