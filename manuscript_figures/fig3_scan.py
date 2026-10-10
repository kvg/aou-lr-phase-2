"""Figure 3 mock-up, layout 1: sufficiency matrix + hit map on shared chapter rows,
method panels, and locus cards.

A  Testable phecodes per PheCodeX chapter x ancestry group (>= MIN_CASES cases).
B  Genome-wide significant SV-phecode associations; rows are the same chapters as A.
C  Calibration (QQ, lambda_GC) for the three test encodings.
D  Repeat-length dosage vs biallelic coding at repeat loci.
E  Local-ancestry-aware joint test vs pooled test.
F  Locus cards: regional plot plus per-ancestry effects.

Real: PheCodeX structure, phenotype-table denominators. Everything else is
simulated; see fig3_sim.py for the replacement data contract.
"""

from __future__ import annotations

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, PowerNorm
from matplotlib.patches import Rectangle

from . import fig3_sim as S
from .style import (
    ANC_COLORS, ANC_ORDER, ANC_TEXT, FIG_W, HG38_MB, chrom_offsets,
    new_figure, save, stamp_mockup,
)

W, H = FIG_W, 9.26
LETTERS = True  # lettered plate for design review; the repo convention is letters only on real-data plates

SLATE_D = "#2F4858"  # focal: repeat loci / what the new features recover
SLATE_L = "#A9BCC8"  # other SVs
GREY = "#9A9A9A"
LIGHT = "#E4E7EA"
CMAP = LinearSegmentedColormap.from_list("slate", ["#F3F5F7", "#B9C9D3", "#2F4858"])

STAMP_REAL = {
    "group": "Panel A: observed case counts. All other panels are simulated.",
    "local": ("Panel A: observed case counts; local-ancestry columns are case haplotypes expected from chr1 global ancestry. "
              "All other panels are simulated."),
}
SUBTITLE = {
    "group": "At least {c} cases in the ancestry group (All: pooled)",
    "local": "At least {c} cases (All) or {h} case haplotypes of a local ancestry",
}
STAMP_SIM = "Case counts, statistics and loci are placeholders; denominators from covariates.v8."
ROWH = 0.145  # inches per chapter row in A/B
MAT_TOP = 0.60


def add(fig, left, top, w, h, **kw):
    """Axes at (left, top) with size (w, h), all in inches from the figure's top-left."""
    fw, fh = fig.get_size_inches()
    return fig.add_axes([left / fw, 1 - (top + h) / fh, w / fw, h / fh], **kw)


def txt(fig, x, y, s, **kw):
    fw, fh = fig.get_size_inches()
    return fig.text(x / fw, 1 - y / fh, s, **kw)


def letter(fig, ch, x, y):
    if LETTERS:
        txt(fig, x, y, ch, fontsize=9, fontweight="bold", ha="left", va="top")


def title(fig, s, x, y, tag=None, real=False):
    """Panel title; ``tag`` marks the panel's data as observed (real=True) or simulated."""
    t = txt(fig, x, y, s, fontsize=7.5, ha="left", va="top")
    if tag:
        from matplotlib.backends.backend_agg import FigureCanvasAgg
        if not hasattr(fig.canvas, "get_renderer"):
            FigureCanvasAgg(fig)
        w_in = t.get_window_extent(fig.canvas.get_renderer()).width / fig.dpi
        txt(fig, x + w_in + 0.09, y + 0.012, tag, fontsize=7, style="italic",
            color="#1E6E8C" if real else "#8A8A8A", ha="left", va="top")


def fmt_p(logp: float) -> str:
    """P = 10**-logp as mantissa x power of ten, to one significant figure."""
    e = int(np.ceil(logp))
    m = round(10 ** (e - logp))
    if m >= 10:
        m, e = 1, e - 1
    return rf"10$^{{-{e}}}$" if m == 1 else rf"{m}$\times$10$^{{-{e}}}$"


# ------------------------------------------------------------------ A
def panel_matrix(fig, cov: S.Coverage, left=0.98):
    """Chapters x ancestry columns + All. Cell = testable phecodes; colour = share of the chapter."""
    keys = cov.keys
    group = cov.mode == "group"
    cols = keys + ["ALL"]
    gap, all_w, tot_w = 0.06, 0.38, 0.40
    c_w = 0.31 if group else 0.40
    x0 = {a: 0.04 + c_w * i for i, a in enumerate(keys)}
    x0["ALL"] = 0.04 + c_w * len(keys) + gap
    x0["TOTAL"] = x0["ALL"] + all_w + gap
    total_w = x0["TOTAL"] + tot_w
    nrow = len(cov.chapters)
    fmt_total = (lambda v: f"{v:,.0f}") if group else (lambda v: f"{v / 1000:.1f}k")
    ax = add(fig, left, MAT_TOP, total_w, ROWH * nrow)
    ax.set_xlim(0, total_w)
    ax.set_ylim(nrow, 0)
    ax.axis("off")
    norm = PowerNorm(0.5, vmin=0, vmax=0.55)
    for r, ch in enumerate(cov.chapters):
        ax.text(-0.04, r + 0.5, ch, ha="right", va="center", fontsize=7, clip_on=False)
        ax.text(x0["TOTAL"] + tot_w / 2, r + 0.5, f"{cov.n_in_chapter[ch]}", ha="center", va="center", fontsize=7,
                color="#333333")
        for a in cols:
            k = cov.testable[ch][a]
            frac = k / cov.n_in_chapter[ch]
            w = all_w if a == "ALL" else c_w
            ax.add_patch(Rectangle((x0[a], r), w, 1, fc=CMAP(norm(frac)), ec="white", lw=0.8))
            dark = norm(frac) > 0.5
            ax.text(x0[a] + w / 2, r + 0.5, f"{k}", ha="center", va="center", fontsize=7,
                    color="white" if dark else ("#B0B0B0" if k == 0 else "#111111"))
    ax.add_patch(Rectangle((x0["ALL"], 0), all_w, nrow, fc="none", ec="#111111", lw=0.7))
    # column headers: ancestry code, then participants (groups) or haplotypes (local ancestries)
    ax.text(x0["TOTAL"] + tot_w / 2, -1.55, "Total", ha="center", va="center", fontsize=7, fontweight="bold", clip_on=False)
    ax.plot([x0["TOTAL"] - gap / 2] * 2, [0, nrow], color="#C8C8C8", lw=0.5)
    ax.text(0.0, -0.5, "Participants" if group else "Haplotypes", ha="right", va="center", fontsize=7, color="#555555", clip_on=False)
    for a in keys:
        ax.text(x0[a] + c_w / 2, -1.55, a, ha="center", va="center", fontsize=7, fontweight="bold",
                color=ANC_TEXT[a], clip_on=False)
        ax.text(x0[a] + c_w / 2, -0.5, fmt_total(cov.totals[a]), ha="center", va="center", fontsize=7,
                color="#555555", clip_on=False)
    ax.text(x0["ALL"] + all_w / 2, -1.55, "All", ha="center", va="center", fontsize=7, fontweight="bold", clip_on=False)
    ax.text(x0["ALL"] + all_w / 2, -0.5, fmt_total(cov.n_total if group else 2 * cov.n_total), ha="center", va="center",
            fontsize=7, color="#555555", clip_on=False)
    # colour key
    cb_top = MAT_TOP + ROWH * nrow + 0.20
    cax = add(fig, left + 0.04, cb_top, 1.30, 0.075)
    grad = np.linspace(0, 0.55, 256)[None, :]
    cax.imshow(norm(grad), cmap=CMAP, aspect="auto", extent=[0, 0.55, 0, 1], vmin=0, vmax=1)
    cax.set_yticks([])
    cax.set_xticks([0, 0.25, 0.5])
    cax.set_xticklabels(["0", "25%", "50%"])
    cax.tick_params(length=2, pad=1.5)
    for sp in cax.spines.values():
        sp.set_visible(False)
    txt(fig, left + 0.0, cb_top + 0.04, "Share of chapter", fontsize=7, va="center", ha="right")
    return total_w


# ------------------------------------------------------------------ B
def panel_hits(fig, cov: S.Coverage, hits: list[S.Hit], left: float, w: float):
    nrow = len(cov.chapters)
    ax = add(fig, left, MAT_TOP, w, ROWH * nrow)
    off = chrom_offsets()
    xmax = off[22] + HG38_MB[22] + 4
    row = {c: i for i, c in enumerate(cov.chapters)}
    for c in HG38_MB:
        if c % 2 == 1:
            ax.axvspan(off[c] - 4, off[c] + HG38_MB[c] + 4, color="#F1F3F5", lw=0, zorder=0)
    for r in range(nrow + 1):
        ax.axhline(r, color="#E1E4E7", lw=0.4, zorder=0.5)
    xs = np.array([off[h.chrom] + h.pos_mb for h in hits])
    ys = np.array([row[h.chapter] + 0.5 for h in hits])
    lp = np.array([min(h.logp, 18) for h in hits])
    size = 5 + (lp - S.GW_LOGP) * 3.3
    rep = np.array([h.repeat for h in hits])
    for flag, fc, ec, z in ((False, SLATE_L, "#6F8A99", 2), (True, SLATE_D, "#16242D", 3)):
        m = rep == flag
        ax.scatter(xs[m], ys[m], s=size[m], fc=fc, ec=ec, lw=0.3, zorder=z)
    u_in = xmax / w  # data units per inch
    for h, x, y, sz in zip(hits, xs, ys, size):
        if h.num:
            edge = x > xmax - 60
            dx = np.sqrt(sz / np.pi) / 72 * u_in + 0.03 * u_in
            ax.text(x + (-dx if edge else dx), y - 0.03, str(h.num), fontsize=7, fontweight="bold", va="center",
                    ha="right" if edge else "left", zorder=5,
                    bbox=dict(boxstyle="square,pad=0.05", fc="white", ec="none", alpha=0.75))
    ax.set_xlim(-8, xmax + 6)
    ax.set_ylim(nrow, 0)
    ax.set_yticks([])
    ticks = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 14, 16, 18, 22]
    ax.set_xticks([off[c] + HG38_MB[c] / 2 for c in ticks])
    ax.set_xticklabels([str(c) for c in ticks], fontsize=7)
    ax.tick_params(axis="x", length=0, pad=2)
    for s in ("left", "top", "right"):
        ax.spines[s].set_visible(False)
    ax.set_xlabel("Chromosome", fontsize=7, labelpad=1)
    # legend in the empty bottom rows (offsets in inches, converted to data units)
    lx, ly = xmax * 0.38, nrow - 3.0
    for i, (lab, fc, ec) in enumerate((("Repeat locus, tested as dosage", SLATE_D, "#16242D"),
                                      ("Other SV, tested as carrier", SLATE_L, "#6F8A99"))):
        ax.scatter([lx], [ly + 0.5 + i], s=14, fc=fc, ec=ec, lw=0.3, clip_on=False, zorder=6)
        ax.text(lx + 0.09 * u_in, ly + 0.5 + i, lab, fontsize=7, va="center", ha="left", zorder=6)
    ax.text(lx - 0.03 * u_in, ly + 2.5, "$-\\log_{10}P$", fontsize=7, va="center", ha="left", zorder=6)
    for j, v in enumerate((8, 12, 16)):
        cx = lx + (0.62 + 0.36 * j) * u_in
        ax.scatter([cx], [ly + 2.5], s=5 + (v - S.GW_LOGP) * 3.3, fc="#D5DEE4", ec="#6F8A99", lw=0.3, zorder=6)
        ax.text(cx + 0.08 * u_in, ly + 2.5, str(v), fontsize=7, va="center", ha="left", zorder=6)
    return ax


# ------------------------------------------------------------------ C
def panel_qq(ax):
    arms = S.sim_qq()
    cols = {"Pooled SAIGE": "#8A8F94", "FELIX, biallelic": "#7C9BAD", "FELIX, repeat dosage": SLATE_D}
    for name, p, lam in arms:
        n = p.size
        exp = -np.log10((np.arange(1, n + 1)) / (n + 1))
        obs = -np.log10(p)
        keep = np.unique(np.concatenate([np.arange(0, 400), np.unique(np.geomspace(400, n, 260).astype(int)) - 1]))
        keep = keep[keep < n]
        ax.scatter(exp[keep], obs[keep], s=2.6, color=cols[name], lw=0, zorder=3 if "dosage" in name else 2,
                   label=f"{name.replace('repeat dosage', 'dosage')}  {lam:.2f}")
    ax.plot([0, 8], [0, 8], color="#444444", lw=0.6, zorder=1)
    ax.axhline(S.GW_LOGP, color=GREY, lw=0.6, ls=(0, (3, 2)), zorder=1)
    ax.set_xlim(0, 5.2)
    ax.set_ylim(0, 14)
    ax.set_xlabel("Expected $-\\log_{10}P$", labelpad=1)
    ax.set_ylabel("Observed $-\\log_{10}P$", labelpad=1)
    ax.legend(loc="upper left", fontsize=7, handletextpad=0.1, markerscale=2.4, labelspacing=0.25,
              borderaxespad=0.1, title=r"$\lambda_{GC}$", alignment="left", title_fontsize=7)


# ------------------------------------------------------------------ D / E
def _scatter_frame(ax, lim=17):
    ax.plot([0, lim], [0, lim], color="#444444", lw=0.6, zorder=1)
    ax.axhline(S.GW_LOGP, color=GREY, lw=0.6, ls=(0, (3, 2)), zorder=1)
    ax.axvline(S.GW_LOGP, color=GREY, lw=0.6, ls=(0, (3, 2)), zorder=1)
    ax.set_xlim(0, lim)
    ax.set_ylim(0, lim)


def _number(ax, x, y, n, dx=0.5, dy=0.0):
    ax.text(x + dx, y + dy, str(n), fontsize=7, fontweight="bold", ha="left", va="center", zorder=6,
            bbox=dict(boxstyle="square,pad=0.05", fc="white", ec="none", alpha=0.75))


def panel_dosage(ax, numbers=None):
    x, y = S.sim_dosage_vs_biallelic()
    x, y = np.minimum(x, 16.5), np.minimum(y, 16.5)
    named = {1: (6.5, 14.8), 4: (4.1, 9.9)}  # the repeat-locus cards
    if numbers is not None:  # {card number in the full set: number shown}; loci without a card are not marked
        named = {numbers[n]: xy for n, xy in named.items() if n in numbers}
    gain = (y > S.GW_LOGP) & (x <= S.GW_LOGP)
    ax.scatter(x[~gain], y[~gain], s=4, color="#B8BDC1", lw=0, zorder=2)
    ax.scatter(x[gain], y[gain], s=9, color=SLATE_D, lw=0, zorder=3)
    for n, (nx, ny) in named.items():
        ax.scatter([nx], [ny], s=14, color=SLATE_D, ec="white", lw=0.4, zorder=4)
        _number(ax, nx, ny, n)
    _scatter_frame(ax)
    ax.set_xlabel("Biallelic coding, $-\\log_{10}P$", labelpad=1)
    ax.set_ylabel("Repeat dosage, $-\\log_{10}P$", labelpad=1)
    ax.scatter([8.0], [1.6], s=9, color=SLATE_D, lw=0, clip_on=False)
    ax.text(8.7, 1.6, "Significant only\nwith dosage", fontsize=7, color=SLATE_D, ha="left", va="center", linespacing=0.95)


def panel_ancestry(ax, numbers=None):
    x, y, het = S.sim_ancestry_vs_pooled()
    x, y = np.minimum(x, 16.5), np.minimum(y, 16.5)
    sig = het > S.HET_LOGP
    ax.scatter(x[~sig], y[~sig], s=2.6, color="#C3C7CA", lw=0, zorder=2, rasterized=True)
    ax.scatter(x[sig], y[sig], s=8, color=SLATE_D, lw=0, zorder=3)
    named = {2: (6.4, 11.6), 3: (12.5, 13.1), 4: (7.9, 9.9), 6: (6.8, 8.1)}
    if numbers is not None:
        named = {numbers[n]: xy for n, xy in named.items() if n in numbers}
    for n, (nx, ny) in named.items():
        ax.scatter([nx], [ny], s=14, color=SLATE_D, ec="white", lw=0.4, zorder=4)
        _number(ax, nx, ny, n)
    _scatter_frame(ax)
    ax.set_xlabel("Pooled SAIGE, $-\\log_{10}P$", labelpad=1)
    ax.set_ylabel("FELIX local ancestry, $-\\log_{10}P$", labelpad=1)
    ax.scatter([8.0], [1.6], s=9, color=SLATE_D, lw=0, clip_on=False)
    ax.text(8.7, 1.6, "Heterogeneous\nacross ancestries", fontsize=7, color=SLATE_D, ha="left", va="center", linespacing=0.95)


# ------------------------------------------------------------------ F
def _nice(x: float) -> float:
    """Largest 1-2-5 value not above x."""
    e = 10 ** np.floor(np.log10(x))
    return next(m * e for m in (5, 2, 1) if m * e <= x)


KIND_LABEL = {"repeat": "repeat", "del": "deletion", "dup": "copy number"}
CARD_W, CARD_H = 1.66, 1.44


def card(fig, left, top, loc: S.Locus, rng):
    eff, n_est, p_het = S.card_effects(loc)
    txt(fig, left, top + 0.02, f"{loc.num}", fontsize=7.5, fontweight="bold", va="top", ha="left")
    txt(fig, left + 0.12, top + 0.02, loc.gene, fontsize=7.5, style="italic", va="top", ha="left")
    txt(fig, left + CARD_W, top + 0.03, KIND_LABEL[loc.kind], fontsize=7, color="#666666", va="top", ha="right")
    txt(fig, left + 0.12, top + 0.165, loc.trait, fontsize=7, va="top", ha="left", color="#333333")
    pl, pw = left + 0.27, CARD_W - 0.31
    # regional plot
    ax = add(fig, pl, top + 0.36, pw, 0.50)
    x, y, sv = S.regional(loc, rng)
    ax.scatter(x, y, s=3.2, color="#B9BEC2", lw=0, zorder=2)
    ymax = int(np.ceil((max(y.max(), sv, S.GW_LOGP) + 1.0) / 5) * 5)
    ax.axhline(S.GW_LOGP, color=GREY, lw=0.5, ls=(0, (3, 2)), zorder=1)
    if loc.lead == "SV":
        ax.scatter([0], [sv], marker="D", s=26, fc=SLATE_D, ec="white", lw=0.4, zorder=5, clip_on=False)
    else:
        ax.scatter([0], [sv], marker="D", s=26, fc="white", ec=SLATE_D, lw=0.9, zorder=5, clip_on=False)
    ax.set_xlim(-0.52, 0.52)
    ax.set_ylim(0, ymax)
    ax.set_yticks([0, ymax])
    ax.set_xticks([])
    ax.tick_params(axis="y", length=2, pad=1)
    ax.spines["bottom"].set_linewidth(0.5)
    ax.text(0.97, 0.97, fmt_p(loc.logp), transform=ax.transAxes, fontsize=7, ha="right", va="top")
    if loc.lead == "SNV":
        ax.text(0.03, 0.97, "SNV lead", transform=ax.transAxes, fontsize=7, ha="left", va="top", color="#555555")
    # effect strip: shared effect first, then each local ancestry (hollow = not estimable)
    ax2 = add(fig, pl, top + 0.93, pw, 0.36)
    akeys = S.ancestry_keys()
    keys = ["ALL"] + akeys
    xs = [0.0] + [1.7 + i for i in range(len(akeys))]
    m = max(abs(b) + 1.96 * se for b, se in eff.values() if np.isfinite(b)) * 1.08
    for key, xx in zip(keys, xs):
        b, se = eff[key]
        colr = SLATE_D if key == "ALL" else ANC_COLORS[key]
        if not np.isfinite(b):
            ax2.scatter([xx], [0], s=9, fc="white", ec="#9A9A9A", lw=0.6, zorder=3)
            continue
        ax2.plot([xx, xx], [0, b], color=colr, lw=1.0, zorder=2, solid_capstyle="butt")
        ax2.plot([xx, xx], [b - 1.96 * se, b + 1.96 * se], color=colr, lw=0.7, alpha=0.65, zorder=2)
        ax2.scatter([xx], [b], s=11, color=colr, lw=0, zorder=4)
    ax2.axhline(0, color="#333333", lw=0.5, zorder=1)
    ax2.axvline(0.85, color="#D6D6D6", lw=0.5, zorder=0)
    ax2.set_xlim(-0.6, xs[-1] + 0.6)
    ax2.set_ylim(-m, m)
    ax2.set_xticks([])
    t = _nice(m * 0.95)
    ax2.set_yticks([-t, 0, t])
    ax2.set_yticklabels([f"{-t:g}".replace("-", "\u2212"), "", f"{t:g}"])
    ax2.tick_params(axis="y", length=2, pad=1)
    ax2.spines["bottom"].set_visible(False)
    if p_het is None:
        het = "n/a"
    elif p_het >= 0.05:
        het = "n.s."
    else:
        het = fmt_p(-np.log10(max(p_het, 1e-300)))
    txt(fig, left + 0.12, top + 1.33, loc.unit, fontsize=7, color="#666666", va="top", ha="left")
    txt(fig, left + CARD_W, top + 1.33, rf"$P_{{het}}$ {het}", fontsize=7, color="#666666", va="top", ha="right")


def key_card(fig, left, top):
    ax = add(fig, left, top, CARD_W, CARD_H)
    ax.set_xlim(0, CARD_W)
    ax.set_ylim(CARD_H, 0)
    ax.axis("off")
    ax.text(0.0, 0.08, "How to read a card", fontsize=7.5, va="center")
    xt = 0.42
    y = 0.31
    ax.scatter([0.08], [y], marker="D", s=20, fc=SLATE_D, ec="white", lw=0.4)
    ax.scatter([0.25], [y], marker="D", s=20, fc="white", ec=SLATE_D, lw=0.9)
    ax.text(xt, y, "Lead SV; hollow = SNV lead", fontsize=7, va="center")
    y += 0.165
    ax.scatter([0.06, 0.13, 0.20, 0.27], [y + 0.03, y - 0.03, y + 0.01, y + 0.04], s=4, color="#B9BEC2", lw=0)
    ax.text(xt, y, "SNVs, 1 Mb, \u2212log$_{10}$P", fontsize=7, va="center")
    y += 0.165
    ax.plot([0.03, 0.32], [y, y], color=GREY, lw=0.6, ls=(0, (3, 2)))
    ax.text(xt, y, "Genome-wide significance", fontsize=7, va="center")
    y += 0.22
    for xx, colr, b_ in ((0.06, SLATE_D, 0.05), (0.19, ANC_COLORS["AFR"], -0.05), (0.28, ANC_COLORS["EUR"], 0.07)):
        ax.plot([xx, xx], [y, y - b_ * 2.0], color=colr, lw=1.0)
        ax.plot([xx, xx], [y - b_ * 2.0 - 0.05, y - b_ * 2.0 + 0.05], color=colr, lw=0.6, alpha=0.7)
        ax.scatter([xx], [y - b_ * 2.0], s=8, color=colr, lw=0)
    ax.text(xt, y, "Effect, 95% CI;\nAll = shared", fontsize=7, va="center", linespacing=0.95)
    y += 0.24
    ax.scatter([0.17], [y], s=9, fc="white", ec="#9A9A9A", lw=0.6)
    ax.text(xt, y, f"<{S.MIN_CASES} cases" if S.ANC_MODE == "group" else f"<{S.MIN_HAP} case haplotypes", fontsize=7, va="center")
    y += 0.20
    labs = ([("All", SLATE_D)] if S.ANC_MODE != "group" else []) + [(a_, ANC_TEXT[a_]) for a_ in S.ancestry_keys()]
    pitch = 0.225 if len(labs) >= 7 else 0.26
    for i, (lab, colr) in enumerate(labs):
        ax.text(0.10 + pitch * i, y, lab, fontsize=7, ha="center", va="center", color=colr, fontweight="bold" if lab == "All" else "normal")


# ------------------------------------------------------------------ assemble
def render():
    cov = S.build_coverage("testable")
    hits = S.build_hits(cov)
    fig = new_figure(H)
    rng = np.random.default_rng(S.SEED + 9)

    # Row 1
    letter(fig, "A", 0.02, 0.0)
    a_real = not cov.simulated
    title(fig, "Testable phecodes", 0.30, 0.02, tag="observed counts" if a_real else "simulated", real=a_real)
    txt(fig, 0.30, 0.20, SUBTITLE[cov.mode].format(c=S.MIN_CASES, h=S.MIN_HAP) + "; Total: phecodes in the chapter",
        fontsize=7, color="#555555", va="top", ha="left")
    total_w = panel_matrix(fig, cov)
    bl = 0.98 + total_w + 0.22
    letter(fig, "B", bl - 0.16, 0.0)
    title(fig, "Significant SV\u2013phecode associations", bl + 0.04, 0.02, tag="simulated")
    panel_hits(fig, cov, hits, bl, W - bl - 0.06)

    # Row 2
    r2_top, r2_h = 3.98, 1.35
    for i, (ch, ttl, fn, x0, w) in enumerate((
        ("C", "Calibration", panel_qq, 0.52, 2.0),
        ("D", "Dosage vs biallelic coding", panel_dosage, 2.98, 1.95),
        ("E", "Ancestry-aware vs pooled", panel_ancestry, 5.28, 1.84),
    )):
        ax = add(fig, x0, r2_top, w, r2_h)
        fn(ax)
        letter(fig, ch, x0 - 0.50 if i == 0 else x0 - 0.52, r2_top - 0.37)
        title(fig, ttl, x0 - 0.30, r2_top - 0.36, tag="simulated")

    # Row 3
    f_top = r2_top + r2_h + 0.72
    letter(fig, "F", 0.02, f_top - 0.34)
    title(fig, "Representative loci: regional association and effect in each ancestry", 0.30, f_top - 0.32, tag="simulated")
    gap = (W - 0.10 - 4 * CARD_W) / 3
    for i in range(8):
        r, c = divmod(i, 4)
        left, top = 0.05 + c * (CARD_W + gap), f_top + r * (CARD_H + 0.10)
        if i < len(S.LOCI):
            card(fig, left, top, S.LOCI[i], rng)
        else:
            key_card(fig, left, top)
    stamp_mockup(fig, STAMP_REAL[cov.mode] if a_real else STAMP_SIM)
    return save(fig, "fig3_scan")
