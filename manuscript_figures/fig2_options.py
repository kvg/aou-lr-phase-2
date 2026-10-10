"""Figure 2 layout options: from the catalog to the full association model.

Exploration pass. Length counts come from the merged Terra bins when present
(currently the integrated phased callset, not yet the GLnexus remake). Global
ancestry is the real FLARE chr1 global estimate. Everything else (context
fractions, Phase 1 spectrum, discovery curves, phasing, allelic series,
repeat expansions, haplotype window) is synthetic and anchored to Table 2
where Table 2 has a number.
"""

from __future__ import annotations

import csv
import gzip
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.gridspec import GridSpec
from matplotlib.lines import Line2D
from matplotlib.patches import Patch, Rectangle

from .fig2_full_model import (
    _accumulate_abs_counts,
    _length_bins_path,
    _load_length_bins_tsv,
    _paint_hap,
)
from .style import (
    ANC_COLORS,
    ANC_ORDER,
    ANC_TEXT as ANC_TEXT_COLORS,
    HG38_MB,
    MIN_FONT,
    OUTDIR,
    P1_N,
    P2_N,
    ROOT,
    SV_COLORS,
    covariates_path,
    despine,
    label_panel,
    new_figure,
    save,
    stamp_mockup,
)

FS = MIN_FONT
LAI = ["AFR", "AMR", "EUR", "EAS", "SAS"]  # FLARE 5-way paint order
# Discovery order: previously well-catalogued non-AFR first, then more
# admixed / underrepresented non-AFR, AFR last (Ebert-style novelty axis).
DISCOVERY_BLOCKS = ("EUR", "EAS", "MID", "SAS", "AMR", "OTH", "AFR")
# MID is ~3% of the cohort; push its neighbours' labels outward so the row stays legible.
STRIP_LABEL_NUDGE = {"EAS": -0.012, "SAS": 0.012}
# FLARE global ancestry. Prefer the genome-wide estimate (hg38 length-weighted
# mean over the 22 autosomes from the production run); fall back to the earlier
# chr1-only file. Single-chromosome noise depresses each participant's maximum
# component, so chr1 overstates how admixed the cohort is.
GLOBAL_ANC_CANDIDATES = (
    ROOT / "aou_lr_phase2_v1.genomewide.global.anc.gz",
    ROOT / "aou_lr_phase2_v1.chr1.global.anc.gz",
)
GLOBAL_ANC = next(
    (p for p in GLOBAL_ANC_CANDIDATES if p.is_file()), GLOBAL_ANC_CANDIDATES[-1]
)

CTX = ("US", "RM", "SD", "SR")
CTX_COLORS = {"US": "#D4D4D4", "RM": "#7F9DC1", "SD": "#D49A55", "SR": "#4E9A80"}
CTX_LABELS = {
    "US": "US (unique)",
    "RM": "RM (RepeatMasker)",
    "SD": "SD (segmental dup.)",
    "SR": "SR (simple repeat)",
}

P2_FILL = "#A9B4CF"
P2_UL_FILL = "#D7B9AE"
P1_FILL = "#3E4C75"
BND_FILL = "#8C8C8C"
P1_LINE = "#222222"
FOLD_COLORS = {"INS": "#222222", "DEL": "#9A9A9A"}
P2_LINE = "#3C5A99"
HPRC_FILL = "#BDBDBD"
NOVEL_FILL = "#3C5A99"

N_SNV_FALLBACK = 1.07e8
N_BND = 130_740

def _log_bp_edges(lo: float, hi: float, frac: float = 0.01) -> np.ndarray:
    """Half-integer edges ~frac of length wide, never narrower than 1 bp."""
    n = int(np.ceil(np.log(hi / lo) / np.log1p(frac)))
    return np.unique(np.round(np.geomspace(lo, hi, n + 1) - 0.5)) + 0.5


SPECS = (
    {
        "title": "≤100 bp",
        "xlim": (0.5, 100),
        "xscale": None,
        "xticks": [1, 20, 50, 100],
        "labels": ["1", "20", "50", "100"],
        "layer": "integrated",
        "width": 1.7,
        "edges": np.arange(0.5, 100.5 + 1e-9, 1.0),
        "snv": False,
    },
    {
        "title": "100 bp–10 kb",
        "xlim": (100, 1e4),
        "xscale": "log",
        "xticks": [100, 335, 1_000, 6_000, 10_000],
        "labels": ["", "335", "1k", "6k", ""],
        "layer": "integrated",
        "width": 2.0,
        "edges": _log_bp_edges(100.5, 1e4 + 0.5),
        "landmarks": ((335, "Alu"), (2_750, "SVA"), (6_000, "L1")),
    },
    {
        "title": "≥10 kb",
        "xlim": (1e4, 3e8),
        "xscale": "log",
        "xticks": [1e4, 1e5, 1e6, 1e7, 1e8],
        "labels": ["10k", "", "1M", "", "100M"],
        "layer": "ultralong",
        "width": 1.0,
        "edges": _log_bp_edges(1e4 + 0.5, 3e8),
    },
)


# ---------------------------------------------------------------------------
# Synthetic per-length annotations (Table 2 anchors where available)
# ---------------------------------------------------------------------------

_CTX_ANCHORS = {
    # log10 |L| → US, RM, SD, SR (%). 20–49 bp and ≥50 bp averages match Table 2.
    "INS": (
        (0.0, (30, 24, 3, 43)),
        (0.7, (18, 16, 2, 64)),
        (1.28, (10, 15, 2, 73)),
        (1.54, (8.1, 17.2, 1.8, 72.9)),
        (2.0, (7, 8, 1.2, 83.8)),
        (3.0, (7, 8, 1.2, 83.8)),
        (3.999, (9, 10, 1.5, 79.5)),
        (4.0, (10, 22, 20, 48)),
        (6.3, (12, 25, 35, 28)),
    ),
    "DEL": (
        (0.0, (32, 25, 3, 40)),
        (0.7, (20, 18, 2, 60)),
        (1.28, (13, 13, 1.6, 72.4)),
        (1.54, (11.5, 12.8, 1.6, 74.1)),
        (2.0, (10, 15, 2, 73)),
        (3.0, (10, 20, 2.5, 67.5)),
        (3.999, (11, 28, 6, 55)),
        (4.0, (14, 35, 22, 29)),
        (6.3, (12, 30, 38, 20)),
    ),
}
_SNV_CTX = np.array([47.0, 45.0, 5.0, 3.0])


def _bump(lg: np.ndarray, center_bp: float, sd: float) -> np.ndarray:
    return np.exp(-0.5 * ((lg - np.log10(center_bp)) / sd) ** 2)


def ctx_fractions(length: np.ndarray, direction: str) -> np.ndarray:
    """(n, 4) US/RM/SD/SR fractions for absolute lengths ≥1.

    MEI-sized insertions pull away from SR toward the genome background
    (US/RM); MEI-sized deletions are reference elements, so they move to RM.
    """
    lg = np.log10(np.clip(np.asarray(length, dtype=float), 1, None))
    anchors = _CTX_ANCHORS[direction]
    xs = np.array([a[0] for a in anchors])
    ys = np.array([a[1] for a in anchors], dtype=float)
    out = np.column_stack([np.interp(lg, xs, ys[:, j]) for j in range(4)])
    for center, sd, amp_ins, amp_del in ((310, 0.05, 0.50, 0.55), (2_000, 0.06, 0.20, 0.15), (6_000, 0.04, 0.40, 0.50)):
        g = _bump(lg, center, sd)
        amp = amp_ins if direction == "INS" else amp_del
        moved = out[:, 3] * amp * g
        out[:, 3] -= moved
        if direction == "INS":
            out[:, 0] += 0.5 * moved
            out[:, 1] += 0.5 * moved
        else:
            out[:, 1] += moved
    return out / out.sum(axis=1, keepdims=True)


def p1_ratio(length: np.ndarray, direction: str) -> np.ndarray:
    """Phase 1 / Phase 2 site ratio from Table 2; NaN where Phase 1 has no catalog."""
    L = np.asarray(length, dtype=float)
    r = np.full(L.shape, np.nan)
    small = (L >= 1) & (L < 20)
    r[small] = 3_932_966 / 13_218_208 if direction == "INS" else 4_734_797 / 17_707_931
    sv = (L >= 50) & (L < 1e4)
    base = 498_090 / 1_875_567 if direction == "INS" else 165_717 / 599_220
    # Lower depth in Phase 1 costs more sensitivity for longer events.
    r[sv] = base * (1.15 - 0.15 * np.log10(L[sv] / 50) / 2.3)
    return r


def multiallelic_pct(length: np.ndarray) -> np.ndarray:
    lg = np.log10(np.clip(np.asarray(length, dtype=float), 1, None))
    xs = (0.0, 0.7, 1.28, 1.7, 2.5, 3.0, 3.8, 4.0, 6.3)
    ys = (12, 25, 33, 38, 30, 35, 25, 18, 8)
    out = np.interp(lg, xs, ys)
    out -= 15 * _bump(lg, 310, 0.05) + 12 * _bump(lg, 6_000, 0.04)
    return out


def phased_pct(length: np.ndarray) -> np.ndarray:
    lg = np.log10(np.clip(np.asarray(length, dtype=float), 1, None))
    xs = (0.0, 1.28, 1.7, 3.0, 3.999, 4.0, 6.3)
    ys = (99.0, 98.2, 97.0, 95.5, 93.0, 85.0, 70.0)
    return np.interp(lg, xs, ys)


# ---------------------------------------------------------------------------
# Length spectrum data
# ---------------------------------------------------------------------------


def _bin_centers(edges: np.ndarray, log: bool) -> np.ndarray:
    if log:
        return np.sqrt(edges[:-1] * edges[1:])
    return 0.5 * (edges[:-1] + edges[1:])


from .style import MANUSCRIPT_ROOT

REGION_BINS_CANDIDATES = (
    MANUSCRIPT_ROOT / "data" / "fig2_length_spectrum.region.bins.tsv",
    MANUSCRIPT_ROOT / "scratch" / "chr6.bins.tsv",
)
P1_BINS_PATH = MANUSCRIPT_ROOT / "data" / "fig2_length_spectrum.phase1.bins.tsv"
# Table 2 genome-wide GLnexus small-variant total (SNV + INS<20 + DEL<20).
GENOME_SMALL_TOTAL = 245_289_070 + 13_218_208 + 17_707_931


def load_region_bins(path: Path) -> dict:
    """VariantLengthSpectrum output: chrom partition signed_len region_class n_sites.

    Returns {(partition, signed_len): counts in CTX order} plus the chroms seen.
    """
    from collections import defaultdict

    counts: dict[tuple[str, int], np.ndarray] = defaultdict(lambda: np.zeros(len(CTX)))
    chroms: set[str] = set()
    with path.open(encoding="utf-8") as fh:
        header = fh.readline()
        if not header.startswith("chrom\tpartition\tsigned_len\tregion_class"):
            raise ValueError(f"unexpected region bins header in {path}: {header!r}")
        for line in fh:
            chrom, part, L, cls, n = line.rstrip("\n").split("\t")
            counts[(part, int(L))][CTX.index(cls)] += float(n)
            chroms.add(chrom)
    return {"counts": dict(counts), "chroms": sorted(chroms), "path": path}


def _region_bins_path() -> Path | None:
    return next((p for p in REGION_BINS_CANDIDATES if p.is_file()), None)


def _bin_region_counts(region: dict, edges: np.ndarray, sign: int, partitions) -> np.ndarray:
    """(n_bin, 4) context counts on absolute-length edges for one direction."""
    out = np.zeros((len(edges) - 1, len(CTX)))
    for (part, L), c in region["counts"].items():
        if part not in partitions or L == 0 or np.sign(L) != sign:
            continue
        i = int(np.searchsorted(edges, abs(L), side="right")) - 1
        if 0 <= i < out.shape[0]:
            out[i] += c
    return out


MIN_CELL = 20


def _merge_sparse_bins(edges, ins, dele, ins_frac, del_frac, log: bool):
    """Merge adjacent bins until every drawn INS/DEL bar is 0 or ≥ MIN_CELL sites."""
    ins_c, del_c = ins_frac * ins[:, None], del_frac * dele[:, None]
    ok = lambda v: v == 0 or v >= MIN_CELL
    groups, start = [], 0
    for j in range(len(ins)):
        if ok(ins[start : j + 1].sum()) and ok(dele[start : j + 1].sum()):
            groups.append((start, j + 1))
            start = j + 1
    if start < len(ins):
        if groups:
            groups[-1] = (groups[-1][0], len(ins))
        else:
            groups.append((0, len(ins)))
    new_edges = np.array([edges[a] for a, _b in groups] + [edges[groups[-1][1]]])
    gi = np.array([ins[a:b].sum() for a, b in groups])
    gd = np.array([dele[a:b].sum() for a, b in groups])
    fi = np.array([ins_c[a:b].sum(axis=0) for a, b in groups])
    fd = np.array([del_c[a:b].sum(axis=0) for a, b in groups])
    fi = np.where(gi[:, None] > 0, fi / np.maximum(gi, 1)[:, None], ins_frac[0])
    fd = np.where(gd[:, None] > 0, fd / np.maximum(gd, 1)[:, None], del_frac[0])
    return new_edges, _bin_centers(new_edges, log), gi, gd, fi, fd


def spectrum_data(rng: np.random.Generator, *, region_path: Path | None = None) -> dict:
    """Length-spectrum bars. Region-class bins (real context) are used wherever
    they exist; lengths they do not cover fall back to the legacy integrated
    bins, rescaled to the region file's share of the genome."""
    path = _length_bins_path()
    real = _load_length_bins_tsv(path) if path else None
    region = load_region_bins(region_path) if region_path else None
    region_parts = {p for p, _L in region["counts"]} if region else set()
    p1_region = load_region_bins(P1_BINS_PATH) if P1_BINS_PATH.is_file() else None
    has_sv = bool(region_parts & {"sv", "ultralong"})
    scale = 1.0
    if region and not has_sv:
        small_total = sum(c.sum() for (p, _L), c in region["counts"].items() if p == "small")
        scale = small_total / GENOME_SMALL_TOTAL
    facets = []
    n_snv = 0.0
    snv_frac = _SNV_CTX / _SNV_CTX.sum()
    for spec in SPECS:
        edges = spec["edges"]
        centers = _bin_centers(edges, spec["xscale"] == "log")
        if real is not None:
            ins_d, del_d, snv = _accumulate_abs_counts(real, spec["layer"], edges)
            ins = sum(ins_d.values()) * scale
            dele = sum(del_d.values()) * scale
            n_snv += snv * scale
        else:
            base = 1.2e6 if spec["layer"] == "integrated" else 2e3
            ins = base * centers ** -1.1 * rng.uniform(0.8, 1.2, centers.size)
            dele = 1.3 * ins
        ins_frac = ctx_fractions(centers, "INS")
        del_frac = ctx_fractions(centers, "DEL")
        real_mask = np.zeros(centers.size, dtype=bool)
        if region:
            parts = ("small", "sv") if spec["layer"] == "integrated" else ("ultralong",)
            for sign, tgt, frac in ((1, "ins", ins_frac), (-1, "del", del_frac)):
                rc = _bin_region_counts(region, edges, sign, parts)
                tot = rc.sum(axis=1)
                hit = tot > 0
                vals = ins if tgt == "ins" else dele
                vals[hit] = tot[hit]
                if has_sv:
                    vals[~hit] = 0.0
                frac[hit] = rc[hit] / tot[hit, None]
                real_mask |= hit
            if has_sv:
                edges, centers, ins, dele, ins_frac, del_frac = _merge_sparse_bins(edges, ins, dele, ins_frac, del_frac, spec["xscale"] == "log")
                real_mask = np.ones(centers.size, dtype=bool)
        facet = {"spec": spec, "edges": edges, "centers": centers, "ins": ins, "del": dele, "ins_frac": ins_frac, "del_frac": del_frac, "real_ctx": real_mask}
        if p1_region:
            parts = ("small", "sv") if spec["layer"] == "integrated" else ("ultralong",)
            for sign, key in ((1, "p1_ins"), (-1, "p1_del")):
                v = _bin_region_counts(p1_region, edges, sign, parts).sum(axis=1)
                # Bins under MIN_CELL are blanked (disclosure); empty bins leave a gap.
                facet[key] = np.where(v >= MIN_CELL, v, np.nan)
        facets.append(facet)
    if region and ("small", 0) in region["counts"]:
        c = region["counts"][("small", 0)]
        n_snv = float(c.sum())
        snv_frac = c / c.sum()
    return {
        "facets": facets,
        "n_snv": n_snv or N_SNV_FALLBACK,
        "snv_frac": snv_frac,
        "real": real is not None,
        "region": region,
        "p1_real": p1_region is not None,
        "p1_n_snv": float(p1_region["counts"][("small", 0)].sum()) if p1_region and ("small", 0) in p1_region["counts"] else None,
        "scale": scale,
        "n_bnd": float(region["counts"][("bnd", 0)].sum()) if region and ("bnd", 0) in region["counts"] else N_BND * scale,
        "bnd_frac": region["counts"][("bnd", 0)] / region["counts"][("bnd", 0)].sum() if region and ("bnd", 0) in region["counts"] else None,
    }


# ---------------------------------------------------------------------------
# Spectrum with aligned tracks
# ---------------------------------------------------------------------------

TRACK_HEIGHT = {"counts": 3.2, "binw": 1.0, "ctx_ins": 0.62, "ctx_del": 0.62, "fold": 1.0, "phased": 0.9, "multi": 0.9}
TRACK_LABEL = {
    "ctx_ins": "INS",
    "ctx_del": "DEL",
    "fold": "Phase 2 /\nPhase 1",
    "phased": "Phased\n(% het)",
    "multi": "Multi-\nallelic (%)",
}


def _style_facet_x(ax, spec, *, show_labels: bool, first: bool) -> None:
    if spec["xscale"]:
        ax.set_xscale(spec["xscale"])
    ax.set_xlim(*spec["xlim"])
    ax.set_xticks(spec["xticks"])
    labels = list(spec["labels"])
    if not first and labels and labels[0] and spec["xscale"] is None:
        labels[0] = ""
    if show_labels:
        ax.set_xticklabels(labels, fontsize=FS)
    else:
        ax.set_xticklabels([])
    ax.xaxis.set_minor_locator(plt.NullLocator())


def draw_spectrum(
    fig: plt.Figure,
    cell,
    data: dict,
    *,
    tracks: tuple[str, ...] = ("ctx_ins", "ctx_del"),
    p1: bool = True,
    legend: bool = True,
    stack_ctx: bool = False,
    bnd_right: bool = False,
    scope: bool = False,
    scope_brackets: bool = True,
    per_bp: bool = False,
    bin_notes: bool = False,
    yscale: str = "symlog",
) -> plt.Axes:
    """Faceted INS-up / DEL-down spectrum with aligned per-length tracks below.

    With ``stack_ctx`` each bar is split into US / RM / SD / SR in proportion
    to its share of the bar's drawn (symlog) height; total height stays the
    site count.

    Returns the leftmost counts axis (for the panel letter).
    """
    rows = ("counts",) + tuple(tracks)
    widths = [s["width"] for s in SPECS] + [0.045, 0.30] if bnd_right else [0.30, 0.05] + [s["width"] for s in SPECS]
    sub = cell.subgridspec(
        len(rows),
        len(widths),
        width_ratios=widths,
        height_ratios=[TRACK_HEIGHT[r] for r in rows],
        wspace=0.07,
        hspace=0.18,
    )
    y_max = 0.0
    y_min = np.inf
    for f in data["facets"]:
        div = np.diff(f["edges"]) if per_bp else 1.0
        for v in (f["ins"] / div, f["del"] / div):
            y_max = max(y_max, float(np.max(v)))
            if np.any(v > 0):
                y_min = min(y_min, float(np.min(v[v > 0])))
    if not per_bp:
        y_max = max(y_max, data.get("n_bnd", N_BND))
    y_max *= 1.6
    linthresh = 10 ** np.floor(np.log10(y_min)) if per_bp else 1.0

    first_ax = None
    nf = len(SPECS)
    for r, row in enumerate(rows):
        if bnd_right:
            facet_axes = [fig.add_subplot(sub[r, i]) for i in range(nf)]
            spacer = fig.add_subplot(sub[r, nf])
            lead = fig.add_subplot(sub[r, nf + 1])
        else:
            lead = fig.add_subplot(sub[r, 0])
            spacer = fig.add_subplot(sub[r, 1])
            facet_axes = [fig.add_subplot(sub[r, i + 2]) for i in range(nf)]
        last_row = r == len(rows) - 1
        if row == "counts":
            first_ax = facet_axes[0] if bnd_right else lead
            _draw_counts_row(lead, facet_axes, data, y_max=y_max, p1=p1, show_x=last_row, stack_ctx=stack_ctx, bnd_right=bnd_right, per_bp=per_bp, linthresh=linthresh, bin_notes=bin_notes, yscale=yscale)
            if scope and bnd_right:
                _draw_scope(fig, facet_axes, lead, spacer, brackets=last_row and scope_brackets, axis_label=last_row)
            else:
                spacer.set_axis_off()
                if scope:
                    _draw_scope(fig, facet_axes, lead, None, brackets=last_row and scope_brackets, axis_label=last_row)
        elif row == "binw":
            _draw_binw_row(lead, facet_axes, data, show_x=last_row, bnd_right=bnd_right)
            if scope and bnd_right:
                _draw_scope(fig, facet_axes, lead, spacer, brackets=last_row and scope_brackets, axis_label=last_row)
            else:
                spacer.set_axis_off()
        elif row.startswith("ctx_"):
            spacer.set_axis_off()
            _draw_ctx_row(lead, facet_axes, data, row[4:].upper(), show_x=last_row)
        else:
            _draw_line_row(lead, facet_axes, data, row, show_x=last_row, bnd_right=bnd_right)
            if scope and bnd_right:
                _draw_scope(fig, facet_axes, lead, spacer, brackets=last_row and scope_brackets, axis_label=last_row)
            else:
                spacer.set_axis_off()

    if legend:
        if stack_ctx:
            handles = [Patch(facecolor=CTX_COLORS[k], label=k) for k in CTX_STACK[::-1]]
        else:
            handles = [Patch(facecolor=P2_FILL, label="Phase 2"), Patch(facecolor=P2_UL_FILL, label="Phase 2 ultralong")]
        if p1:
            handles.append(Line2D([], [], color=P1_LINE, lw=0.8, label="Phase 1"))
        if not stack_ctx and any(t.startswith("ctx") for t in tracks):
            handles += [Patch(facecolor=CTX_COLORS[k], label=k) for k in CTX]
        if not stack_ctx and any(t in ("fold", "phased", "multi") for t in tracks):
            handles += [
                Line2D([], [], color=SV_COLORS["INS"], lw=1.0, label="INS"),
                Line2D([], [], color=SV_COLORS["DEL"], lw=1.0, label="DEL"),
            ]
        bbox = first_ax.get_position()
        fig.legend(
            handles=handles,
            loc="lower left",
            bbox_to_anchor=(bbox.x0 - 0.01, bbox.y1 + 0.022),
            ncol=len(handles),
            frameon=False,
            fontsize=FS,
            handlelength=0.9,
            handletextpad=0.35,
            columnspacing=0.9,
        )
    return first_ax


CTX_STACK = ("SR", "SD", "RM", "US")  # axis outward


def _scale_frac(transform, values, lo: float, hi: float) -> np.ndarray:
    """Map data values to 0–1 through an axis scale, without needing a drawn canvas."""
    v = np.asarray(transform.transform(np.atleast_1d(np.asarray(values, dtype=float))), dtype=float)
    ends = np.asarray(transform.transform(np.array([lo, hi], dtype=float)), dtype=float)
    return (v - ends[0]) / (ends[1] - ends[0])


# Rows cover the symlog range at print resolution; columns are fine enough that a
# 1%-of-length bin in the ≥10 kb facet is more than a pixel, so nothing is left
# to blend with the background.
_CTX_PX_X = 2400
_CTX_PX_Y = 1200
_CTX_RGB = np.array(
    [[int(CTX_COLORS[k][i : i + 2], 16) for i in (1, 3, 5)] for k in CTX_STACK], dtype=np.float32
) / 255.0


def _paint_ctx(ax, parts) -> None:
    """Draw context stacks as one opaque image in axes coordinates.

    Vector bars in the ≥10 kb facet are a fraction of a pixel wide. PDF readers
    rasterize each rectangle on its own and blend it with the background, so
    the thin bars look paler than the wide ones. Compositing into one image
    makes every pixel an opaque mix of the bars that cover it.
    """
    nx, ny = _CTX_PX_X, _CTX_PX_Y
    rgb = np.zeros((ny, nx, 3), dtype=np.float32)
    wgt = np.zeros((ny, nx), dtype=np.float32)
    idx = [CTX.index(k) for k in CTX_STACK]
    x_lo, x_hi = ax.get_xlim()
    y_lo, y_hi = ax.get_ylim()
    y0 = float(_scale_frac(ax.yaxis.get_transform(), [0.0], y_lo, y_hi)[0])

    def splat(x0f, x1f, y0f, y1f, color) -> None:
        xa, xb = (x0f, x1f) if x1f >= x0f else (x1f, x0f)
        ya, yb = (y0f, y1f) if y1f >= y0f else (y1f, y0f)
        xa, xb = xa * nx, xb * nx
        ya, yb = ya * ny, yb * ny
        c0, c1 = max(0, int(np.floor(xa))), min(nx, int(np.ceil(xb)))
        r0, r1 = max(0, int(np.floor(ya))), min(ny, int(np.ceil(yb)))
        if c1 <= c0 or r1 <= r0:
            return
        cols = np.arange(c0, c1)
        rows = np.arange(r0, r1)
        xw = np.clip(np.minimum(xb, cols + 1) - np.maximum(xa, cols), 0, None)
        yw = np.clip(np.minimum(yb, rows + 1) - np.maximum(ya, rows), 0, None)
        w = yw[:, None] * xw[None, :]
        rgb[r0:r1, c0:c1] += w[..., None] * color
        wgt[r0:r1, c0:c1] += w

    for x0, widths, totals, frac, sign in parts:
        totals = np.asarray(totals, dtype=float)
        frac = np.asarray(frac, dtype=float)
        if frac.ndim == 1:
            frac = frac[None, :]
        left = np.asarray(x0, dtype=float)
        xf0 = _scale_frac(ax.xaxis.get_transform(), left, x_lo, x_hi)
        xf1 = _scale_frac(ax.xaxis.get_transform(), left + np.asarray(widths, dtype=float), x_lo, x_hi)
        cum = np.cumsum(frac[:, idx], axis=1)
        bounds = np.column_stack([np.zeros(len(totals)), cum])
        y_top = _scale_frac(ax.yaxis.get_transform(), sign * totals, y_lo, y_hi)
        for i, total in enumerate(totals):
            if total <= 0:
                continue
            xf_i = (xf0[i], xf1[i])
            for j in range(len(CTX_STACK)):
                if bounds[i, j + 1] <= bounds[i, j]:
                    continue
                span = float(y_top[i] - y0)
                splat(float(xf_i[0]), float(xf_i[1]), y0 + span * bounds[i, j], y0 + span * bounds[i, j + 1], _CTX_RGB[j])

    if not np.any(wgt > 0.5):
        return
    img = np.zeros((ny, nx, 4), dtype=np.float32)
    m = wgt > 0.5
    img[m, :3] = rgb[m] / wgt[m, None]
    img[m, 3] = 1.0
    ax.imshow(
        img,
        origin="lower",
        extent=(0, 1, 0, 1),
        transform=ax.transAxes,
        aspect="auto",
        interpolation="bilinear",
        zorder=2,
        clip_on=True,
    )


SCOPE_SHADE = "#F1F1F1"
SCOPE_COLOR = "#333333"


def _draw_scope(fig, facet_axes, bnd_ax, spacer=None, *, brackets: bool = True, axis_label: bool = False) -> None:
    """Shade ≥10 kb + BND. ``brackets`` adds the study/companion labels under the x axis."""
    for ax in (facet_axes[-1], bnd_ax):
        ax.set_facecolor(SCOPE_SHADE)
    if spacer is not None:
        spacer.set_xlim(0, 1)
        spacer.set_ylim(0, 1)
        spacer.set_xticks([])
        spacer.set_yticks([])
        for side in ("top", "bottom", "left", "right"):
            spacer.spines[side].set_visible(False)
        spacer.set_facecolor("white")
        spacer.patch.set_facecolor("white")
        pos = spacer.get_position()
        fw, fh = fig.get_size_inches()
        w_in, h_in = pos.width * fw, pos.height * fh
        rise = h_in / w_in  # axes-x run for a 45° line spanning the full height
        step = 0.045 / w_in
        for t in np.arange(-rise, 1.0 + step, step):
            spacer.plot([t, t + rise], [0, 1], transform=spacer.transAxes, color="#8A8A8A", lw=0.45, solid_capstyle="butt", clip_on=True, zorder=1)
    if not brackets and not axis_label:
        return
    fig.canvas.draw()
    tick_labels = [t for a in (*facet_axes, bnd_ax) for t in a.get_xticklabels() if t.get_text()]
    tick_bottom = fig.transFigure.inverted().transform((0, min(t.get_window_extent().y0 for t in tick_labels)))[1]
    xl0, xl1 = facet_axes[0].get_position().x0, facet_axes[-1].get_position().x1
    fig.text((xl0 + xl1) / 2, tick_bottom - 0.006, "Variant length (bp)", ha="center", va="top", fontsize=FS)
    if not brackets:
        return
    groups = (
        (facet_axes[0], facet_axes[-2], "Association analyses (this study)"),
        (facet_axes[-1], bnd_ax, "Companion manuscript"),
    )
    for a0, a1, label in groups:
        p0, p1 = a0.get_position(), a1.get_position()
        y = tick_bottom - 0.046
        x0, x1 = p0.x0 + 0.002, p1.x1 - 0.002
        tick = 0.008
        fig.add_artist(Line2D([x0, x0, x1, x1], [y + tick, y, y, y + tick], transform=fig.transFigure, color=SCOPE_COLOR, lw=0.6))
        fig.text((x0 + x1) / 2, y - 0.004, label, ha="center", va="top", fontsize=FS, color=SCOPE_COLOR)


def _bin_note(spec: dict, drawn_edges: np.ndarray) -> str:
    """Bin width for a facet title: fixed bp for linear facets, % of length for log ones."""
    base = spec["edges"]
    if spec["xscale"] is None:
        return f"{base[1] - base[0]:g} bp bins"
    return "log bins"


def _lin_ticklabels(ax) -> None:
    lo, hi = ax.get_ylim()
    step = 10 ** np.floor(np.log10(max(hi, 1) / 2.5))
    ticks = np.arange(np.ceil(lo / step) * step, hi, step)
    ticks = ticks[np.abs(ticks) > step / 2]

    def fmt(v: float) -> str:
        a = abs(v)
        if a >= 1e6:
            return f"{v / 1e6:g}M"
        if a >= 1e3:
            return f"{v / 1e3:g}k"
        return f"{v:g}"

    ax.set_yticks(ticks)
    ax.set_yticklabels([fmt(t) for t in ticks], fontsize=FS)


def _draw_counts_row(lead, axes, data, *, y_max: float, p1: bool, show_x: bool, stack_ctx: bool = False, bnd_right: bool = False, per_bp: bool = False, linthresh: float = 1.0, bin_notes: bool = False, yscale: str = "symlog") -> None:
    for i, (ax, f) in enumerate(zip(axes, data["facets"])):
        spec = f["spec"]
        edges = f["edges"]
        fill = P2_FILL if spec["layer"] == "integrated" else P2_UL_FILL
        w = np.diff(edges)
        div = w if per_bp else 1.0
        f_ins, f_del = f["ins"] / div, f["del"] / div
        if yscale == "linear_facet":
            local = max(float(np.max(f["ins"])), float(np.max(f["del"])), 1.0) * 1.25
            ax.set_yscale("linear")
            ax.set_ylim(-local, local)
        elif yscale == "linear":
            ax.set_yscale("linear")
            ax.set_ylim(-y_max, y_max)
        else:
            ax.set_yscale("symlog", linthresh=linthresh, linscale=0.3, base=10)
            ax.set_ylim(-y_max, y_max)
        if spec["xscale"]:
            ax.set_xscale(spec["xscale"])
        ax.set_xlim(*spec["xlim"])
        if stack_ctx:
            parts = [
                (edges[:-1], w, f_ins, f.get("ins_frac", ctx_fractions(f["centers"], "INS")), 1.0),
                (edges[:-1], w, f_del, f.get("del_frac", ctx_fractions(f["centers"], "DEL")), -1.0),
            ]
            if spec.get("snv"):
                snv_frac = np.asarray(data.get("snv_frac", _SNV_CTX / _SNV_CTX.sum()))[None, :]
                for sign in (1.0, -1.0):
                    parts.append((np.array([-0.5]), np.array([1.0]), [data["n_snv"]], snv_frac, sign))
            _paint_ctx(ax, parts)
        else:
            ax.bar(edges[:-1], f_ins, width=w, align="edge", color=fill, lw=0, zorder=2)
            ax.bar(edges[:-1], -f_del, width=w, align="edge", color=fill, lw=0, zorder=2)
        if p1:
            for d, sign in (("INS", 1.0), ("DEL", -1.0)):
                if "p1_ins" in f:
                    y = sign * f["p1_ins" if d == "INS" else "p1_del"] / div
                else:
                    vals = f_ins if d == "INS" else f_del
                    y = sign * vals * p1_ratio(f["centers"], d)
                ax.step(edges, np.r_[y, y[-1]], where="post", color=P1_LINE, lw=0.4, alpha=0.7, zorder=4)
        if spec.get("snv"):
            if not stack_ctx:
                ax.bar([0], [data["n_snv"]], width=1.0, color=fill, lw=0, zorder=2)
                ax.bar([0], [-data["n_snv"]], width=1.0, color=fill, lw=0, zorder=2)
            if p1:
                p1_snv = data.get("p1_n_snv") or data["n_snv"] * 54_790_068 / 245_289_070
                for sign in (1.0, -1.0):
                    ax.plot([-0.5, 0.5], [sign * p1_snv] * 2, color=P1_LINE, lw=0.55, zorder=4)
            ax.annotate(
                "SNV",
                xy=(0.5, data["n_snv"]),
                xytext=(5, -1),
                textcoords="offset points",
                fontsize=FS,
                va="center",
                color="#444444",
            )
        for x, lab in spec.get("landmarks", ()):
            ax.axvline(x, color="#B8B8B8", lw=0.45, ls="--", zorder=1)
            ax.text(x, 0.97, lab, transform=ax.get_xaxis_transform(), ha="center", va="top", fontsize=FS, color="#555555", clip_on=True)
        ax.axhline(0, color="#999999", lw=0.4, zorder=4)
        _style_facet_x(ax, spec, show_labels=show_x, first=i == 0)
        # Common title row across facets (avoid set_title pad differences).
        title = spec["title"]
        if bin_notes:
            title = f"{title}, {_bin_note(spec, edges)}"
        ax.text(0.5, 1.04, title, transform=ax.transAxes, ha="center", va="bottom", fontsize=FS)
        despine(ax, left=False)
        if yscale == "linear_facet":
            ax.tick_params(left=True, labelleft=True)
            _lin_ticklabels(ax)
        else:
            ax.tick_params(left=False, labelleft=False)
        if not show_x:
            ax.spines["bottom"].set_visible(False)
            ax.tick_params(bottom=False)
    axes[0].text(0.98, 0.97, "INS", transform=axes[0].transAxes, ha="right", va="top", fontsize=FS, color="#444444", fontweight="bold")
    axes[0].text(0.98, 0.03, "DEL", transform=axes[0].transAxes, ha="right", va="bottom", fontsize=FS, color="#444444", fontweight="bold")

    n_bnd = data.get("n_bnd", N_BND)
    lead.set_xlim(-0.6, 0.6)
    if yscale == "linear_facet":
        lead.set_yscale("linear")
        lead.set_ylim(-n_bnd * 1.25, n_bnd * 1.25)
    elif yscale == "linear":
        lead.set_yscale("linear")
        lead.set_ylim(-y_max, y_max)
    else:
        lead.set_yscale("symlog", linthresh=1, linscale=0.3, base=10)
    if per_bp:
        # Breakends have no length, so their bar is a count on its own scale.
        lead.set_ylim(-n_bnd * 40, n_bnd * 40)
        lead.text(0, n_bnd * 1.6, f"{n_bnd:,.0f}", ha="center", va="bottom", fontsize=FS, color="#333333")
    elif yscale != "linear_facet":
        lead.set_ylim(-y_max, y_max)
    # The context split is done in scale coordinates, so the y scale must be set first.
    if stack_ctx and data.get("bnd_frac") is not None:
        _paint_ctx(lead, [(np.array([-0.3]), np.array([0.6]), [n_bnd], np.asarray(data["bnd_frac"])[None, :], 1.0)])
    else:
        lead.bar([0], [n_bnd], width=0.6, color=BND_FILL, lw=0, zorder=2)
    lead.axhline(0, color="#999999", lw=0.4)
    lead.set_xticks([0])
    lead.set_xticklabels(["BND"] if show_x else [], fontsize=FS)
    lead.text(0.5, 1.04, "no length", transform=lead.transAxes, ha="center", va="bottom", fontsize=FS)
    yax = axes[0] if bnd_right else lead
    if yscale == "symlog":
        ticks = [10.0 ** k for k in range(int(np.ceil(np.log10(linthresh))) + 2, int(np.log10(y_max)) + 1) if k % 4 == 0] if per_bp else [1e3, 1e6]
        yax.tick_params(left=True, labelleft=True)
        yax.spines["left"].set_visible(True)
        yax.set_yticks([-t for t in reversed(ticks)] + [0] + ticks)
        yax.set_yticklabels(
            [rf"$10^{{{int(np.log10(t))}}}$" for t in reversed(ticks)] + ["0"] + [rf"$10^{{{int(np.log10(t))}}}$" for t in ticks],
            fontsize=FS,
        )
        yax.yaxis.set_minor_locator(plt.NullLocator())
    elif yscale == "linear":
        yax.tick_params(left=True, labelleft=True)
        yax.spines["left"].set_visible(True)
        _lin_ticklabels(yax)
    else:
        lead.tick_params(right=True, labelright=True)
        lead.yaxis.tick_right()
        _lin_ticklabels(lead)
    yax.set_ylabel("Sites per bp" if per_bp else "Sites", fontsize=FS, labelpad=1)
    despine(lead, left=not bnd_right)
    if not show_x:
        lead.spines["bottom"].set_visible(False)
        lead.tick_params(bottom=False)
    if bnd_right:
        lead.spines["left"].set_visible(False)
        lead.tick_params(left=False, labelleft=False)


def _ctx_draw(ax, x0, widths, frac) -> None:
    bottom = np.zeros(len(x0))
    for j, k in enumerate(CTX):
        ax.bar(x0, frac[:, j], width=widths, bottom=bottom, align="edge", color=CTX_COLORS[k], lw=0)
        bottom += frac[:, j]


def _draw_ctx_row(lead, axes, data, direction: str, *, show_x: bool) -> None:
    for i, (ax, f) in enumerate(zip(axes, data["facets"])):
        spec = f["spec"]
        edges = f["edges"]
        frac = ctx_fractions(f["centers"], direction)
        _ctx_draw(ax, edges[:-1], np.diff(edges), frac)
        if spec.get("snv"):
            _ctx_draw(ax, np.array([-0.5]), np.array([1.0]), (_SNV_CTX / _SNV_CTX.sum())[None, :])
        ax.set_ylim(0, 1)
        _style_facet_x(ax, spec, show_labels=show_x, first=i == 0)
        ax.set_yticks([])
        for side in ("left", "top", "right"):
            ax.spines[side].set_visible(False)
        if not show_x:
            ax.spines["bottom"].set_visible(False)
            ax.tick_params(bottom=False)
    lead.set_axis_off()
    lead.text(1.0, 0.5, TRACK_LABEL[f"ctx_{direction.lower()}"], transform=lead.transAxes, ha="right", va="center", fontsize=FS)


def _draw_binw_row(lead, axes, data, *, show_x: bool, bnd_right: bool) -> None:
    """Bin width (bp) of the drawn bars as a function of variant length."""
    widths = [np.diff(f["edges"]) for f in data["facets"]]
    hi = 10 ** np.ceil(np.log10(max(w.max() for w in widths)))
    for i, (ax, f, w) in enumerate(zip(axes, data["facets"], widths)):
        spec = f["spec"]
        ax.step(f["edges"], np.r_[w, w[-1]], where="post", color="#333333", lw=0.7)
        ax.set_yscale("log")
        ax.set_ylim(0.5, hi)
        _style_facet_x(ax, spec, show_labels=show_x, first=i == 0)
        despine(ax, left=False)
        ax.tick_params(left=False, labelleft=False)
        ax.yaxis.set_minor_locator(plt.NullLocator())
    yax = axes[0] if bnd_right else lead
    ticks, labels = [1, 1e3, 1e6], ["1", "1k", "1M"]
    yax.set_yticks(ticks)
    yax.set_yticklabels(labels, fontsize=FS)
    yax.tick_params(left=True, labelleft=True)
    yax.spines["left"].set_visible(True)
    yax.set_ylabel("Bin (bp)", fontsize=FS, labelpad=1)
    lead.set_xlim(-0.6, 0.6)
    lead.set_yticks([])
    lead.set_xticks([0])
    lead.set_xticklabels(["BND"] if show_x else [], fontsize=FS)
    lead.text(0, 0.5, "—", transform=lead.get_xaxis_transform(), ha="center", va="center", fontsize=FS, color="#777777")
    despine(lead, left=True)


def _draw_line_row(lead, axes, data, row: str, *, show_x: bool, bnd_right: bool = False) -> None:
    ylims = {"fold": (1, 30), "phased": (60, 100), "multi": (0, 45)}
    for i, (ax, f) in enumerate(zip(axes, data["facets"])):
        spec = f["spec"]
        c = f["centers"]
        if row == "fold":
            for d in ("INS", "DEL"):
                if "p1_ins" in f:
                    p2v, p1v = (f["ins"], f["p1_ins"]) if d == "INS" else (f["del"], f["p1_del"])
                    with np.errstate(divide="ignore", invalid="ignore"):
                        y = np.where(p1v > 0, p2v / p1v, np.nan)
                    ax.step(f["edges"], np.r_[y, y[-1]], where="post", color=FOLD_COLORS[d], lw=0.6)
                    if i == 1:
                        ax.text(2_500, 19 if d == "DEL" else 2.1, d, ha="center", va="center", fontsize=FS, color=FOLD_COLORS[d], fontweight="bold")
                else:
                    ax.plot(c, 1.0 / p1_ratio(c, d), color=SV_COLORS[d], lw=0.9)
            if spec.get("snv"):
                ax.plot([0], [245_289_070 / 54_790_068], "o", ms=2.5, color="#333333")
            ax.axhline(P2_N / P1_N, color="#777777", lw=0.5, ls=":")
            ax.set_yscale("log")
            ax.set_yticks([1, 3, 10, 30])
            if i == 0:
                ax.text(spec["xlim"][1], P2_N / P1_N * 1.08, f"{P2_N / P1_N:.0f}× participants", ha="right", va="bottom", fontsize=FS, color="#555555")
                if "p1_ins" not in f:
                    ax.axvspan(19.5, 49.5, color="#F0F0F0", lw=0, zorder=0)
            if i == 2:
                ax.text(0.5, 0.5, "not in Phase 1", transform=ax.transAxes, ha="center", va="center", fontsize=FS, color="#777777")
        elif row == "phased":
            ax.plot(c, phased_pct(c), color="#333333", lw=0.9)
            if spec.get("snv"):
                ax.plot([0], [99.4], "o", ms=2.5, color="#333333")
            ax.set_yticks([60, 80, 100])
        else:
            ax.plot(c, multiallelic_pct(c), color="#333333", lw=0.9)
            if spec.get("snv"):
                ax.plot([0], [2.0], "o", ms=2.5, color="#333333")
            ax.set_yticks([0, 20, 40])
        ax.set_ylim(*ylims[row])
        ax.yaxis.set_minor_locator(plt.NullLocator())
        _style_facet_x(ax, spec, show_labels=show_x, first=i == 0)
        ax.grid(axis="y", color="#EEEEEE", lw=0.4)
        despine(ax)
        if i > 0:
            ax.tick_params(labelleft=False)
        else:
            ax.tick_params(axis="y", labelsize=FS)
            if row == "fold":
                ax.set_yticklabels(["1", "3", "10", "30"])
        if not show_x:
            ax.tick_params(bottom=False)
    if bnd_right:
        axes[0].set_ylabel(TRACK_LABEL[row], fontsize=FS)
        lead.set_xlim(-0.6, 0.6)
        lead.set_yticks([])
        lead.set_xticks([0])
        lead.set_xticklabels(["BND"] if show_x else [], fontsize=FS)
        lead.tick_params(bottom=show_x)
        despine(lead, left=True)
        return
    lead.set_axis_off()
    lead.text(1.0, 0.5, TRACK_LABEL[row], transform=lead.transAxes, ha="right", va="center", fontsize=FS)


def spectrum_xlabel(fig, ax_ref: plt.Axes) -> None:
    pos = ax_ref.get_position()
    fig.text((pos.x0 + 0.99) / 2, pos.y0 - 0.032, "|Length| (bp)", ha="center", va="top", fontsize=FS)


# ---------------------------------------------------------------------------
# Other panels
# ---------------------------------------------------------------------------


def panel_discovery(ax, *, freq: bool = False) -> None:
    """Cumulative SVs ≥50 bp vs genomes; split into HPRC/HGSVC3-known and new."""
    k = np.unique(np.geomspace(1, P2_N, 400).astype(int))
    total = 25_000 * k ** 0.488
    known_frac = np.clip(0.95 * k ** -0.22, 0, 1)
    known = total * known_frac
    novel = total - known
    if freq:
        sing = novel * np.clip(0.15 + 0.35 * np.log10(k) / 4.1, 0, 0.5)
        ax.stackplot(k, known, novel - sing, sing, colors=[HPRC_FILL, NOVEL_FILL, "#9DB3DD"], lw=0)
    else:
        ax.stackplot(k, known, novel, colors=[HPRC_FILL, NOVEL_FILL], lw=0)
    k1 = k[k <= P1_N]
    p1 = 11_000 * k1 ** 0.55
    ax.plot(k1, p1, color=P1_LINE, lw=0.9, ls="--")
    ax.plot([P1_N], [p1[-1]], "o", ms=2.5, color=P1_LINE)
    ax.annotate(
        "Phase 1",
        xy=(P1_N, p1[-1]),
        xytext=(-6, 14),
        textcoords="offset points",
        ha="right",
        fontsize=FS,
        color=P1_LINE,
        arrowprops=dict(arrowstyle="-", lw=0.4, color=P1_LINE),
    )
    ax.axvline(291, color="#666666", lw=0.45, ls=":")
    ax.text(291 * 0.85, total[-1] * 0.98, "291 HPRC +\nHGSVC3 genomes", ha="right", va="top", fontsize=FS, color="#444444")
    ax.text(3_500, 25_000 * 3_500 ** 0.488 * 0.55, "new", ha="center", va="center", fontsize=FS, color="white")
    ax.text(3_500, 25_000 * 3_500 ** 0.488 * 0.95 * 3_500 ** -0.22 * 0.4, "known", ha="center", va="center", fontsize=FS, color="#333333")
    ax.set_xscale("log")
    ax.set_xlim(1, P2_N)
    ax.set_ylim(0, total[-1] * 1.04)
    ax.set_xticks([1, 10, 100, 1_000, 10_000])
    ax.set_xticklabels(["1", "10", "100", "1k", "10k"])
    ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda v, _p: f"{v / 1e6:.1f}"))
    ax.set_xlabel("Genomes")
    ax.set_ylabel("SVs ≥50 bp absent from\nHPRC/HGSVC3 (M)")
    despine(ax)


EBERT_CLASSES = (
    # figS8 stratum key, label, colour (HGSVC3 Fig. 1f palette)
    ("shared", "Fixed", "#B5403A"),
    ("major", "Major", "#7B2D8E"),
    ("polymorphic", "Minor", "#0B1F78"),
    ("singleton", "Singleton", "#79A6AC"),
)
# Panel B counts SVs with no carrier among the HPRC/HGSVC3 controls in the joint
# callset; prefix-recomputed carrier classes, axis outward (dark → light).
DISCOVERY_CLASSES = (
    ("recurrent", "≥2 carriers", "#4F6878"),
    ("singleton", "Singleton", "#CFDAE0"),
)
# Real panel B bands (fig2_01_panel_b_discovery), stacked outward from the axis.
# Carrier class is recomputed at each prefix, as in HGSVC3 Fig. 1f.
DISCOVERY_CLASSES_P1 = (
    ("in_controls", "In HPRC/HGSVC3", "#8E8E8E"),
    ("in_phase1", "In Phase 1", "#D4D4D4"),
    ("new_recurrent", "New, ≥2 carriers", "#4F6878"),
    ("new_singleton", "New, singleton", "#CFDAE0"),
)
PANEL_B_PATH = MANUSCRIPT_ROOT / "data" / "fig2_panel_b_discovery.tsv"


def load_panel_b(n: int) -> tuple[dict, dict] | None:
    """Cumulative counts per (svtype, class) interpolated onto participants 1..n, plus control-carried totals."""
    if not PANEL_B_PATH.is_file():
        return None
    grid: dict[tuple[str, str], list[tuple[int, int]]] = {}
    with PANEL_B_PATH.open(encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            grid.setdefault((r["svtype"], r["class"]), []).append((int(r["n_participants"]), int(r["cumulative"])))
    summary_path = PANEL_B_PATH.with_suffix(".summary.json")
    summary = json.loads(summary_path.read_text()) if summary_path.is_file() else {}
    if summary.get("n_participants", n) != n:
        return None
    x = np.arange(1, n + 1)
    curves = {}
    for (sv, cls), pts in grid.items():
        k, v = np.array(sorted(pts)).T
        curves.setdefault(sv, {})[cls] = np.interp(x, np.r_[0, k], np.r_[0, v])
    counts = summary.get("counts", {})
    ctrl = {sv: counts.get(f"in_controls_{sv}") for sv in ("INS", "DEL")}
    return curves, {"ctrl": ctrl, "summary": summary}
N_CONTROLS = 293
# Mock: SVs ≥50 bp carried by ≥1 control (the assembly-panel catalog within this callset).
CTRL_INS_50 = 600_000
CTRL_DEL_50 = 250_000
HGSVC3_SV = 188_000
P1_SV_50 = 165_717 + 498_090
P1_DEL_50 = 165_717
P1_INS_50 = 498_090
# Table 2 Phase 2 ≥50 bp totals (for mocking the butterfly split).
P2_DEL_50 = 599_220
P2_INS_50 = 1_875_567


def panel_discovery_ebert(ax, rng: np.random.Generator) -> None:
    """Ebert-style cumulative SVs ≥50 bp, samples added non-AFR first, then AFR.

    Frequency classes are recomputed after each added sample (figS8 model).
    Reference levels: HGSVC3 catalog and the Phase 1 catalog.
    """
    from .figS8_ebert_discovery import _incremental, _p2_blocks, _strata
    from .style import SV_N

    blocks = _p2_blocks()
    afr_start = next(s for a, s, _ in blocks if a == "AFR")
    new = _incremental(P2_N, 26_000, float(SV_N["SV_50"]), afr_start, 1.70, rng)
    st = _strata(new, shared_inf=1_200, major_plat=28_000, major_tau=80, sing_tau=1600, sing_max=0.50, afr_start=afr_start)
    x = np.arange(1, P2_N + 1)
    ax.stackplot(x, *[st[k] for k, _l, _c in EBERT_CLASSES], colors=[c for *_x, c in EBERT_CLASSES], labels=[l for _k, l, _c in EBERT_CLASSES], lw=0, rasterized=True)
    ax.axvline(afr_start + 0.5, color="#555555", ls="--", lw=0.6)
    ax.text(afr_start * 0.5, 1.01, "Non-AFR", transform=ax.get_xaxis_transform(), ha="center", va="bottom", fontsize=FS)
    ax.text(afr_start + (P2_N - afr_start) * 0.5, 1.01, "AFR", transform=ax.get_xaxis_transform(), ha="center", va="bottom", fontsize=FS)
    for y, lab in ((P1_SV_50, "Phase 1 catalog"), (HGSVC3_SV, "HGSVC3 catalog")):
        ax.axhline(y, color="white", lw=0.6, ls=(0, (2, 1.5)))
        ax.text(afr_start * 0.5, y, lab, ha="center", va="bottom", fontsize=FS, color="white")
    ax.set_xlim(1, P2_N)
    ax.set_ylim(0, st["total"][-1] * 1.04)
    ax.set_xticks([1, 5_000, 10_000])
    ax.set_xticklabels(["1", "5k", "10k"])
    ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda v, _p: f"{v / 1e6:.1f}"))
    ax.set_ylabel("SVs ≥50 bp absent from\nHPRC/HGSVC3 (M)")
    ax.legend(loc="upper left", fontsize=FS, handlelength=0.8, handletextpad=0.3, borderaxespad=0.1, labelspacing=0.2, reverse=True)
    despine(ax)
    bar = ax.inset_axes([0, -0.055, 1, 0.04])
    for anc, start, end in blocks:
        bar.add_patch(Rectangle((start + 1, 0), end - start, 1, facecolor=ANC_COLORS[anc], edgecolor="none"))
    bar.set_xlim(1, P2_N)
    bar.set_ylim(0, 1)
    bar.set_axis_off()
    ax.tick_params(axis="x", pad=9, length=0)
    ax.set_xlabel("Discovery sample (ancestry-ordered)")


def load_discovery_cohort(*, with_ids: bool = False):
    """Join Fig 1 computed-ancestry labels with FLARE chr1 global proportions.

    Returns (props n×5 in LAI order, labels in ANC_ORDER codes), plus the
    research ids when ``with_ids``; None if inputs are missing.
    """
    path = covariates_path()
    if path is None or not GLOBAL_ANC.is_file():
        return None

    def _truthy(v) -> bool:
        return str(v or "").strip().lower() in {"1", "true", "t", "yes", "y"}

    flare: dict[str, np.ndarray] = {}
    with gzip.open(GLOBAL_ANC, "rt") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            flare[str(r["SAMPLE"]).strip()] = np.array([float(r[a.lower()]) for a in LAI])
    labels, props, ids = [], [], []
    with gzip.open(path, "rt") as fh:
        for row in csv.DictReader(fh):
            if not _truthy(row.get("final_releasable_v9", "")):
                continue
            if row.get("technology") != "PacBio":
                continue
            rid = str(row.get("research_id") or "").strip()
            if rid not in flare:
                continue
            anc = (row.get("ancestry_pred_other") or "oth").upper()
            if anc not in ANC_ORDER:
                anc = "OTH"
            labels.append(anc)
            props.append(flare[rid])
            ids.append(rid)
    if not labels:
        return None
    if with_ids:
        return np.vstack(props), np.array(labels), np.array(ids)
    return np.vstack(props), np.array(labels)


def write_discovery_order(out: Path) -> Path:
    """Panel B participant order (one research_id per line, with block) for the Workbench discovery run."""
    loaded = load_discovery_cohort(with_ids=True)
    if loaded is None:
        raise SystemExit("covariates or FLARE global ancestry missing; cannot build discovery order")
    q, labels, ids = loaded
    q = q / q.sum(axis=1, keepdims=True)
    order, blocks = discovery_order(q, labels)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("wt", encoding="utf-8") as fh:
        fh.write("research_id\tblock\n")
        for i, b in zip(order, blocks):
            fh.write(f"{ids[i]}\t{DISCOVERY_BLOCKS[b]}\n")
    return out


def discovery_order(q: np.ndarray, labels: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Computed-ancestry blocks (Non-AFR first, AFR last), then rising AFR fraction."""
    block = np.array([DISCOVERY_BLOCKS.index(a) for a in labels])
    afr = q[:, LAI.index("AFR")]
    order = np.lexsort((afr, block))
    return order, block[order]


def panel_discovery_structure(fig, cell, rng: np.random.Generator) -> plt.Axes:
    """Butterfly of cumulative SVs absent from HPRC/HGSVC3 (INS +y, DEL −y) over the ancestry strip.

    Grouping matches Figure 1 (ancestry_pred_other); the paint is real FLARE
    chr1 global ancestry. The discovery curves and control catalog are mocked.
    """
    loaded = load_discovery_cohort()
    if loaded is not None:
        q, labels = loaded
    else:
        q = rng.dirichlet(np.ones(len(LAI)), P2_N)
        labels = np.array([DISCOVERY_BLOCKS[i % len(DISCOVERY_BLOCKS)] for i in range(P2_N)])
    q = q / q.sum(axis=1, keepdims=True)
    order, blocks = discovery_order(q, labels)
    q = q[order]
    n = q.shape[0]
    afr = q[:, LAI.index("AFR")]

    # Per-sample sites absent from the controls; INS share ~ Table 2 ≥50 bp mix.
    weight = np.arange(1, n + 1, dtype=float) ** -0.55 * (1.0 + 0.8 * afr) * rng.uniform(0.93, 1.07, n)
    ins_frac = 0.76 + 0.04 * afr
    k = np.arange(1, n + 1, dtype=float)
    sing_frac = 0.58 + 0.42 * np.exp(-(k - 1.0) / 900.0)
    real_b = load_panel_b(n) if loaded is not None else None
    classes = DISCOVERY_CLASSES_P1 if real_b else DISCOVERY_CLASSES
    ctrl_ins, ctrl_del = CTRL_INS_50, CTRL_DEL_50
    if real_b:
        strata, _meta = real_b
        for key in strata:
            for c, _l, _col in classes:
                strata[key].setdefault(c, np.zeros(n))
            strata[key]["total"] = sum(strata[key][c] for c, _l, _col in classes)
    else:
        strata = {}
        for key, frac, total in (("INS", ins_frac, P2_INS_50 - CTRL_INS_50), ("DEL", 1.0 - ins_frac, P2_DEL_50 - CTRL_DEL_50)):
            new = weight * frac
            cum = np.cumsum(new * total / new.sum())
            strata[key] = {"singleton": cum * sing_frac, "recurrent": cum * (1.0 - sing_frac), "total": cum}

    sub = cell.subgridspec(2, 1, height_ratios=[1.0, 0.30], hspace=0.06)
    ax = fig.add_subplot(sub[0, 0])
    ax_q = fig.add_subplot(sub[1, 0], sharex=ax)
    n_gaps = len(np.unique(blocks)) - 1
    width_pt = ax.get_position().width * fig.get_figwidth() * 72
    gap = n / (width_pt / 0.7 - n_gaps)  # 0.7 pt per gap on the page
    segments = []
    for b in range(len(DISCOVERY_BLOCKS)):
        idx = np.where(blocks == b)[0]
        if idx.size:
            segments.append((b, idx.min(), idx.max() + 1, gap * len(segments)))
    keys = [k for k, _l, _c in classes]
    colors = [c for *_x, c in classes]
    for _b, i0, i1, off in segments:
        xseg = np.arange(i0, i1) + 1 + off
        ax.stackplot(xseg, *[strata["INS"][k][i0:i1] for k in keys], colors=colors, lw=0, rasterized=True)
        ax.stackplot(xseg, *[-strata["DEL"][k][i0:i1] for k in keys], colors=colors, lw=0, rasterized=True)
        ax.plot([xseg[0], xseg[-1]], [0, 0], color="#999999", lw=0.4, zorder=4, solid_capstyle="butt")
    if not real_b:
        x_lab = n * 0.012
        for y in (ctrl_ins, -ctrl_del):
            ax.axhline(y, color="#111111", lw=0.6, ls=(0, (2.5, 1.5)), zorder=4)
        ax.text(x_lab, ctrl_ins, "HPRC/HGSVC3 catalog", ha="left", va="bottom", fontsize=FS, color="#111111", zorder=5)
    else:
        ctrl_ins = ctrl_del = 0
    ins_top = strata["INS"]["total"][-1]
    del_bot = strata["DEL"]["total"][-1]
    ax.set_ylim(-max(del_bot, ctrl_del) * 1.32, max(ins_top, ctrl_ins) * 1.12)
    step = 0.5e6 if max(ins_top, del_bot) > 1.6e6 else 0.25e6
    yt = np.arange(-np.floor(max(del_bot, ctrl_del) * 1.3 / step) * step, max(ins_top, ctrl_ins) * 1.1, step)
    ax.set_yticks(yt)
    ax.set_yticklabels([f"{abs(t) / 1e6:g}" for t in yt])
    ax.set_ylabel("Phase 2 SVs ≥50 bp (M)" if real_b else "SVs ≥50 bp absent from\nHPRC/HGSVC3 (M)")
    ax.text(0.99, 0.98, "INS", transform=ax.transAxes, ha="right", va="top", fontsize=FS, color="#444444", fontweight="bold")
    ax.text(0.99, 0.02, "DEL", transform=ax.transAxes, ha="right", va="bottom", fontsize=FS, color="#444444", fontweight="bold")
    handles = [Patch(facecolor=c, label=l) for _k, l, c in classes[::-1]]
    ax.legend(handles=handles, loc="upper left", bbox_to_anchor=(0.0, 1.0), fontsize=FS, handlelength=0.9, handleheight=0.9, handletextpad=0.35, borderaxespad=0.2, labelspacing=0.25)
    ax.tick_params(labelbottom=False, bottom=False)
    despine(ax)
    ax.spines["bottom"].set_color("black")

    for _b, i0, i1, off in segments:
        bottom = np.zeros(i1 - i0)
        xs = np.arange(i0, i1 + 1) + 0.5 + off
        for j, a in enumerate(LAI):
            top = bottom + q[i0:i1, j]
            ax_q.fill_between(xs, np.r_[bottom, bottom[-1]], np.r_[top, top[-1]], step="post", color=ANC_COLORS[a], lw=0, rasterized=True)
            bottom = top
    ax_q.set_ylim(0, 1)
    ax_q.set_yticks([0, 1])
    ax_q.set_ylabel("Global\nancestry")
    ax_q.set_xticks([])
    for b, i0, i1, off in segments:
        name = DISCOVERY_BLOCKS[b]
        cx = (i0 + i1) / 2 + 0.5 + off + STRIP_LABEL_NUDGE.get(name, 0.0) * n
        ax_q.text(cx, -0.08, name, transform=ax_q.get_xaxis_transform(), ha="center", va="top", fontsize=FS - 0.5, color=ANC_TEXT_COLORS.get(name, "#111111"))
    ax_q.set_xlim(0.5, n + 0.5 + segments[-1][3])
    ax_q.set_xlabel(f"Participants in discovery order (n = {n:,})", labelpad=12)
    despine(ax_q)
    ax_q.spines["bottom"].set_color("black")
    return ax


def panel_allelic_series(ax) -> None:
    """Loci with ≥k distinct alleles (CCDF), Phase 1 vs Phase 2; tails <20 cut."""
    k = np.geomspace(2, 600, 200)
    p2 = 185_000 * (k / 2) ** -1.3
    p1 = 42_000 * (k / 2) ** -1.8
    p2_tr = 1_900_000 * (k / 2) ** -1.15
    for y, color, ls, lab in (
        (p2_tr, P2_LINE, ":", "Phase 2 TR"),
        (p2, P2_LINE, "-", "Phase 2 SV"),
        (p1, P1_LINE, "--", "Phase 1 SV"),
    ):
        m = y >= 20
        ax.plot(k[m], y[m], color=color, ls=ls, lw=1.0, label=lab)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(2, 600)
    ax.set_ylim(20, 4e6)
    ax.set_xticks([2, 10, 100])
    ax.set_xticklabels(["2", "10", "100"])
    ax.set_xlabel("Distinct alleles per locus (≥k)")
    ax.set_ylabel("Loci")
    ax.legend(loc="upper right", fontsize=FS, handlelength=1.6)
    despine(ax)


def panel_phasing(ax) -> None:
    classes = ["SNV", "Indel", "SV", "TR"]
    strata = (
        ("Phase 1", "#BBBBBB", "o", (96.0, 93.5, 85.0, 80.0)),
        ("Mid-pass", "#A9BBC9", "o", (99.2, 98.4, 95.1, 93.0)),
        ("High-pass", "#334E68", "o", (99.7, 99.3, 97.8, 96.5)),
    )
    y = np.arange(len(classes))[::-1]
    for i, cls in enumerate(classes):
        vals = [s[3][i] for s in strata]
        ax.plot([min(vals), max(vals)], [y[i]] * 2, color="#CCCCCC", lw=1.2, zorder=1)
    for name, color, marker, vals in strata:
        ax.scatter(vals, y, s=16, color=color, marker=marker, edgecolor="#333333" if name == "Phase 1" else "none", lw=0.4, zorder=3, label=name)
    ax.set_yticks(y)
    ax.set_yticklabels(classes)
    ax.set_xlim(78, 100.5)
    ax.set_ylim(-0.7, len(classes) - 0.3)
    ax.set_xlabel("Heterozygous calls phased (%)")
    ax.legend(loc="lower center", bbox_to_anchor=(0.45, 1.0), ncol=3, fontsize=FS, handletextpad=0.1, columnspacing=0.8, borderaxespad=0.1)
    despine(ax)


def load_global_ancestry() -> tuple[np.ndarray, list[str]] | None:
    if not GLOBAL_ANC.is_file():
        return None
    rows = []
    with gzip.open(GLOBAL_ANC, "rt") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for r in reader:
            rows.append([float(r[a.lower()]) for a in LAI])
    return np.array(rows), LAI


def panel_global_ancestry(ax, *, phase1_bar: bool = False) -> None:
    """Real FLARE chr1 global ancestry, sorted by dominant component."""
    loaded = load_global_ancestry()
    if loaded is None:
        ax.text(0.5, 0.5, "global ancestry file missing", ha="center", va="center", transform=ax.transAxes)
        return
    q, names = loaded
    q = q / q.sum(axis=1, keepdims=True)
    dom = q.argmax(axis=1)
    order = np.lexsort((-q.max(axis=1), dom))
    q = q[order]
    dom = dom[order]
    n = q.shape[0]
    x = np.arange(n + 1)
    bottom = np.zeros(n)
    for j, a in enumerate(names):
        top = bottom + q[:, j]
        ax.fill_between(x, np.r_[bottom, bottom[-1]], np.r_[top, top[-1]], step="post", color=ANC_COLORS[a], lw=0, rasterized=True)
        bottom = top
    for j, a in enumerate(names):
        idx = np.where(dom == j)[0]
        width_in = ax.get_position().width * ax.figure.get_figwidth()
        if idx.size / n * width_in < 0.16:
            continue
        ax.text(idx.mean(), 1.03, a, ha="center", va="bottom", fontsize=FS, color=ANC_COLORS[a])
        ax.axvline(idx.max() + 1, color="white", lw=0.6)
    admixed = float(np.mean(q.max(axis=1) < 0.8))
    ax.set_xlim(0, n)
    ax.set_ylim(0, 1)
    ax.set_yticks([0, 0.5, 1])
    ax.set_xticks([])
    ax.set_ylabel("Global ancestry")
    ax.set_xlabel(f"Participants (n = {n:,}); {100 * admixed:.0f}% with max component < 0.8")
    despine(ax, bottom=False)


def panel_karyogram(ax, rng: np.random.Generator, *, chroms: int = 22) -> None:
    mix = np.array([0.22, 0.38, 0.36, 0.02, 0.02])
    h = 0.34
    for c in range(1, chroms + 1):
        y = chroms - c
        L = HG38_MB[c]
        for k, off in enumerate((0.04, -0.04 - h)):
            tracts = _paint_hap(L, mix, mean_mb=18.0, flicker=0.05, rng=rng)
            for x0, x1, j in tracts:
                ax.add_patch(Rectangle((x0, y + off), x1 - x0, h, facecolor=ANC_COLORS[LAI[j]], edgecolor="none"))
        if c in (1, 5, 10, 15, 20):
            ax.text(-4, y, str(c), ha="right", va="center", fontsize=FS, color="#444444")
    ax.set_xlim(-12, 252)
    ax.set_ylim(-0.6, chroms - 0.3)
    ax.set_yticks([])
    ax.set_xticks([0, 100, 200])
    ax.set_xlabel("Mb (one admixed participant, two haplotypes)")
    despine(ax, left=False)


def panel_tr_by_lai(ax, rng: np.random.Generator) -> None:
    """VNTR allele frequency by repeat length, one line per local ancestry at the locus.

    Frequencies backed by fewer than 20 haplotypes are not drawn.
    """
    alleles = np.array([3, 4, 5, 6, 7, 8, 9, 10, 12, 14, 16, 20, 24, 30])
    freqs = {
        "AFR": [2, 18, 6, 22, 3, 14, 2, 9, 7, 5, 4, 4, 2, 2],
        "AMR": [1, 30, 3, 34, 1, 16, 1, 6, 3, 2, 1, 1, 0.5, 0.5],
        "EUR": [1, 36, 2, 38, 1, 14, 0.5, 5, 1, 0.5, 0.5, 0.3, 0.1, 0.1],
        "EAS": [0.5, 18, 1, 58, 0.5, 15, 0.5, 4, 1, 0.5, 0.5, 0.3, 0.1, 0.1],
        "SAS": [1, 28, 2, 42, 1, 16, 0.5, 6, 2, 0.6, 0.5, 0.3, 0.1, 0.1],
    }
    lai_haps = {"AFR": 6_200, "AMR": 4_100, "EUR": 9_300, "EAS": 2_500, "SAS": 2_400}
    x = np.arange(alleles.size)
    for a in LAI:
        p = np.array(freqs[a], dtype=float)
        p /= p.sum()
        counts = rng.multinomial(lai_haps[a], p).astype(float)
        af = np.where(counts >= 20, counts / lai_haps[a], np.nan)
        ax.plot(x, af, color=ANC_COLORS[a], lw=0.9, marker="o", ms=2.2, label=a)
    ax.set_yscale("log")
    ax.set_ylim(1e-3, 1)
    ax.set_xticks(x[::2])
    ax.set_xticklabels(alleles[::2])
    ax.set_xlabel("Repeat units (VNTR allele)")
    ax.set_ylabel("Allele frequency in\nlocal-ancestry haplotypes")
    ax.legend(loc="lower center", bbox_to_anchor=(0.5, 1.0), ncol=5, fontsize=FS, handlelength=0.9, handletextpad=0.25, columnspacing=0.6, borderaxespad=0.1)
    despine(ax)


EXPANSION_LOCI = (
    # locus, motif, normal median RU, log sd, pathogenic threshold RU, expanded alleles
    ("FMR1", "CGG", 30, 0.10, 55, 120),
    ("HTT", "CAG", 18, 0.13, 36, 31),
    ("C9orf72", "GGGGCC", 4, 0.35, 30, 9),
    ("FXN", "GAA", 10, 0.35, 66, 210),
    ("DMPK", "CTG", 11, 0.30, 50, 14),
    ("ATXN2", "CAG", 22, 0.04, 33, 26),
)


def panel_expansions(ax, rng: np.random.Generator) -> None:
    """Allele length relative to the pathogenic threshold at clinical TR loci."""
    bins = np.geomspace(0.05, 20, 90)
    n_alleles = 2 * P2_N
    for i, (locus, motif, med, sd, thr, n_exp) in enumerate(EXPANSION_LOCI):
        normal = rng.lognormal(np.log(med), sd, n_alleles - n_exp)
        expanded = thr * rng.lognormal(np.log(1.6), 0.45, n_exp)
        ratio = np.r_[normal, expanded] / thr
        counts, _ = np.histogram(ratio, bins=bins)
        y = len(EXPANSION_LOCI) - 1 - i
        dens = np.log10(counts + 1)
        dens = dens / dens.max()
        for j in np.nonzero(counts)[0]:
            ax.add_patch(
                Rectangle((bins[j], y - 0.32), bins[j + 1] - bins[j], 0.64, facecolor=plt.cm.Greys(0.25 + 0.7 * dens[j]), edgecolor="none")
            )
        beyond = int(np.sum(ratio >= 1))
        ax.text(21, y, f"{beyond:,}" if beyond >= 20 else "<20", ha="left", va="center", fontsize=FS, color="#B03A1A")
    ax.axvline(1, color="#B03A1A", lw=0.7)
    ax.text(1.0, 1.01, "pathogenic threshold", transform=ax.get_xaxis_transform(), ha="center", va="bottom", fontsize=FS, color="#B03A1A")
    ax.set_xscale("log")
    ax.set_xlim(0.05, 20)
    n = len(EXPANSION_LOCI)
    ax.set_ylim(-0.6, n - 0.4)
    ax.set_yticks([n - 1 - i for i in range(n)])
    ax.set_yticklabels([loc for loc, *_ in EXPANSION_LOCI])
    ax.tick_params(axis="y", length=0)
    ax.set_xticks([0.1, 1, 10])
    ax.set_xticklabels(["0.1×", "1×", "10×"])
    ax.xaxis.set_minor_locator(plt.NullLocator())
    ax.set_xlabel("Length / pathogenic threshold")
    despine(ax, left=False)


P1_N_PARTICIPANTS = 1_027
# Haplotype ancestry mix per phase: Phase 1 was African American/Black participants;
# Phase 2 uses the cohort-mean FLARE proportions when available.
P1_HAP_ANC_P = {"AFR": 0.80, "AMR": 0.02, "EUR": 0.17, "EAS": 0.005, "SAS": 0.005}
P2_HAP_ANC_P = {"AFR": 0.30, "AMR": 0.12, "EUR": 0.35, "EAS": 0.12, "SAS": 0.11}
ANC_SD_MULT = {"AFR": 1.25, "AMR": 1.05, "EUR": 1.0, "EAS": 0.95, "SAS": 1.0}
THRESHOLD_COLOR = "#B03A1A"
# Mirrored allelic series: Phase 1 above each locus axis, Phase 2 below.
P1_FILL, P1_TICK = "#CFCFCF", "#7A7A7A"
P2_FILL, P2_TICK_SEEN, P2_TICK_NEW = "#BCCAD3", "#8FA5B2", "#22323F"


def _locus_model(rng, med, sd):
    """Locus-level truth shared by both phases: per-ancestry mode shifts and a
    few ancestry-private alleles (units, within-ancestry frequency)."""
    shift = {a: rng.normal(0, 0.12) for a in LAI}
    second = {a: (rng.choice([-1, 1]) * rng.uniform(0.25, 0.45), rng.uniform(0.0, 0.45)) for a in LAI}
    private = {a: [(int(round(med * np.exp(rng.normal(0, 1.8 * sd)))), rng.uniform(0.004, 0.02)) for _ in range(rng.integers(1, 4))] for a in LAI}
    return shift, second, private


def _locus_haplotypes(rng, n_hap, med, sd, tail_center, exp_rate, model, anc_p):
    shift, second, private = model
    anc = rng.choice(LAI, size=n_hap, p=[anc_p[a] for a in LAI])
    mult = np.array([ANC_SD_MULT[a] for a in anc])
    loc = np.log(med) + np.array([shift[a] for a in anc])
    off, w = np.array([second[a] for a in anc]).T
    loc = loc + np.where(rng.random(n_hap) < w, off, 0.0)
    units = rng.lognormal(loc, sd * mult)
    n_exp = rng.binomial(n_hap, exp_rate)
    if n_exp:
        units[:n_exp] = tail_center * rng.lognormal(0.0, 0.45, n_exp)
    for a, alleles in private.items():
        idx = np.nonzero(anc == a)[0]
        for u, f in alleles:
            units[idx[rng.random(idx.size) < f]] = u
    return np.maximum(1, np.round(units)).astype(int), anc


def _log_density(units: np.ndarray, ref: float, grid: np.ndarray, bw: float = 0.035) -> np.ndarray:
    """Haplotype density over log10(length / ref), Gaussian-smoothed."""
    x = np.log10(units / ref)
    step = grid[1] - grid[0]
    hist, _ = np.histogram(x, bins=np.r_[grid - step / 2, grid[-1] + step / 2])
    k = np.arange(-4 * bw, 4 * bw + step, step)
    kern = np.exp(-0.5 * (k / bw) ** 2)
    return np.convolve(hist / x.size, kern / kern.sum(), mode="same")


# Allelic-series rows: (label, is_gene, median RU, log sd, pathogenic threshold RU or None,
# long-tail haplotype rate, long-tail centre as a multiple of the median).
# Clinical anchors, then Danzi et al. PLVI candidates, then top-PLVI placeholders.
TR_SERIES_LOCI = (
    ("FMR1", True, 30, 0.10, 200, 0.004, 2.2),
    ("HTT", True, 18, 0.13, 40, 0.002, 2.0),
    ("C9orf72", True, 5, 0.40, 30, 0.002, 4.0),
    ("EP400", True, 14, 0.16, None, 0.004, 2.4),
    ("FAM193B", True, 9, 0.22, None, 0.004, 2.8),
    ("PLVI top 1", False, 12, 0.22, None, 0.010, 2.4),
    ("PLVI top 2", False, 24, 0.18, None, 0.012, 2.2),
    ("PLVI top 3", False, 7, 0.32, None, 0.008, 3.0),
)
SERIES_LO, SERIES_HI = 0.05, 5.0
OUTLIER_Q = 99.9
OUTLIER_LINE = "#333333"


def panel_allelic_series_loci(fig, cell, rng: np.random.Generator) -> plt.Axes:
    """Per-locus mirrored allelic series on a length / threshold axis.

    The threshold is the pathogenic threshold where one is known, otherwise the
    Phase 2 99.9th percentile (Danzi et al. outlier boundary), so 1× lines up
    across rows. Above each locus axis: Phase 1
    haplotype density and one tick per distinct allele. Below: the same for
    Phase 2, with alleles absent from Phase 1 drawn dark. Known pathogenic
    thresholds are marked on their own row.
    """
    sub = cell.subgridspec(1, 2, width_ratios=[1.0, 0.20], wspace=0.12)
    ax = fig.add_subplot(sub[0, 0])
    ax_bar = fig.add_subplot(sub[0, 1], sharey=ax)

    loaded = load_discovery_cohort()
    if loaded is not None:
        mean_q = loaded[0].mean(axis=0)
        p2_anc = dict(zip(LAI, mean_q / mean_q.sum()))
    else:
        p2_anc = P2_HAP_ANC_P
    n = len(TR_SERIES_LOCI)
    amp, tick = 0.38, 0.14
    lo, hi = SERIES_LO, SERIES_HI
    grid = np.linspace(np.log10(lo), np.log10(hi), 400)
    xg = 10 ** grid
    counts = []
    for i, (_label, _gene, med, sd, thr, rate, fold) in enumerate(TR_SERIES_LOCI):
        model = _locus_model(rng, med, sd)
        u2, _a2 = _locus_haplotypes(rng, 2 * P2_N, med, sd, med * fold, rate, model, p2_anc)
        u1, _a1 = _locus_haplotypes(rng, 2 * P1_N_PARTICIPANTS, med, sd, med * fold, rate, model, P1_HAP_ANC_P)
        ref = float(thr) if thr is not None else float(np.percentile(u2, OUTLIER_Q))
        p1, p2 = np.unique(u1), np.unique(u2)
        seen = np.isin(p2, p1)
        counts.append((int(p1.size), int(p2.size), int(seen.sum())))
        base = n - 1 - i
        d1, d2 = _log_density(u1, ref, grid), _log_density(u2, ref, grid)
        scale = amp / max(d1.max(), d2.max())
        ax.fill_between(xg, base, base + d1 * scale, color=P1_FILL, lw=0, zorder=1)
        ax.fill_between(xg, base, base - d2 * scale, color=P2_FILL, lw=0, zorder=1)
        in_view = lambda u: u[(u / ref >= lo) & (u / ref <= hi)]
        ax.vlines(in_view(p1) / ref, base, base + tick, color=P1_TICK, lw=0.45, zorder=2)
        ax.vlines(in_view(p2[seen]) / ref, base - tick, base, color=P2_TICK_SEEN, lw=0.45, zorder=2)
        ax.vlines(in_view(p2[~seen]) / ref, base - tick, base, color=P2_TICK_NEW, lw=0.6, zorder=3)
        ax.plot([lo, hi], [base, base], color="#333333", lw=0.5, zorder=4, solid_capstyle="butt")
        if thr is not None:
            ax.plot([1, 1], [base - 0.4, base + 0.4], color=THRESHOLD_COLOR, lw=0.9, zorder=5, solid_capstyle="butt")
        else:
            ax.plot([1, 1], [base - 0.4, base + 0.4], color=OUTLIER_LINE, lw=0.7, ls=(0, (2.5, 1.5)), zorder=5)
        if i == 0:
            ax.text(hi, base + 0.08, "Phase 1", ha="right", va="bottom", fontsize=FS, color="#555555")
            ax.text(hi, base - 0.08, "Phase 2", ha="right", va="top", fontsize=FS, color="#111111")
    ax.set_xscale("log")
    ax.set_xlim(lo / 1.25, hi * 1.15)
    ax.set_ylim(-0.6, n - 0.4)
    ax.set_yticks([n - 1 - i for i in range(n)])
    ax.set_yticklabels([lab for lab, *_ in TR_SERIES_LOCI])
    for t, (_lab, gene, *_r) in zip(ax.get_yticklabels(), TR_SERIES_LOCI):
        t.set_style("italic" if gene else "normal")
    ax.tick_params(axis="y", length=0)
    xt = [0.05, 0.1, 0.2, 0.5, 1, 2, 5]
    ax.set_xticks(xt)
    ax.set_xticklabels([f"{t:g}×" for t in xt])
    ax.xaxis.set_minor_locator(plt.NullLocator())
    ax.set_xlabel("Repeat length / threshold")
    x0, x1 = np.log10(ax.get_xlim())
    ax.xaxis.label.set_x((np.log10(np.sqrt(lo * hi)) - x0) / (x1 - x0))
    despine(ax, left=False)
    ax.spines["left"].set_visible(False)
    ax.spines["bottom"].set_bounds(lo, hi)

    xmax = max(c for row in counts for c in row)
    for i, (n1, n2, n2_seen) in enumerate(counts):
        base = n - 1 - i
        y1, y2 = base + 0.17, base - 0.17
        ax_bar.barh(y1, n1, height=0.22, color=P1_TICK, lw=0)
        ax_bar.barh(y2, n2_seen, height=0.22, color=P2_TICK_SEEN, lw=0)
        ax_bar.barh(y2, n2 - n2_seen, height=0.22, left=n2_seen, color=P2_TICK_NEW, lw=0)
        ax_bar.text(n1 + 0.04 * xmax, y1, f"{n1:,}", ha="left", va="center", fontsize=FS, color="#555555")
        ax_bar.text(n2 + 0.04 * xmax, y2, f"{n2:,}", ha="left", va="center", fontsize=FS, color="#111111")
    ax_bar.set_xlim(0, xmax * 1.55)
    ax_bar.set_xticks([])
    ax_bar.tick_params(labelleft=False, left=False)
    for side in ("top", "right", "left", "bottom"):
        ax_bar.spines[side].set_visible(False)
    ax_bar.text(0.0, n - 0.4, "Alleles", ha="left", va="bottom", fontsize=FS)

    tick_kw = dict(marker="|", ls="none", ms=6, mew=1.0)
    handles = [
        Patch(facecolor=P1_FILL, label="Phase 1 lengths"),
        Patch(facecolor=P2_FILL, label="Phase 2 lengths"),
        Line2D([], [], color=THRESHOLD_COLOR, label="Pathogenic threshold", **tick_kw),
        Line2D([], [], color=P2_TICK_SEEN, label="Allele seen in Phase 1", **tick_kw),
        Line2D([], [], color=P2_TICK_NEW, label="New allele in Phase 2", **tick_kw),
        Line2D([], [], color=OUTLIER_LINE, lw=0.6, ls=(0, (2.5, 1.5)), label=f"Phase 2 {OUTLIER_Q:g}th percentile"),
    ]
    ax.legend(handles=handles, bbox_to_anchor=(0.0, 1.02), ncol=2, frameon=False, fontsize=FS, handlelength=0.8, handleheight=0.7, handletextpad=0.35, columnspacing=0.9, labelspacing=0.25, borderaxespad=0, loc="lower left")
    return ax


def panel_context_p1p2(ax) -> None:
    """Table 2 sequence-context mix for SVs ≥50 bp, Phase 1 vs Phase 2 (real)."""
    rows = (
        ("DEL  Phase 1", (5.0, 12.8, 2.0, 80.2)),
        ("DEL  Phase 2", (10.1, 23.1, 2.7, 64.0)),
        ("INS  Phase 1", (6.7, 9.1, 1.2, 83.0)),
        ("INS  Phase 2", (7.4, 8.4, 1.2, 83.0)),
    )
    ys = [3.2, 2.4, 1.0, 0.2]
    for y, (lab, vals) in zip(ys, rows):
        left = 0.0
        for k, v in zip(CTX, vals):
            ax.barh(y, v, left=left, height=0.62, color=CTX_COLORS[k], lw=0)
            if v >= 8:
                ax.text(left + v / 2, y, f"{v:.0f}", ha="center", va="center", fontsize=FS, color="white" if k != "US" else "#333333")
            left += v
    ax.set_yticks(ys)
    ax.set_yticklabels([r[0] for r in rows])
    ax.set_xlim(0, 100)
    ax.set_xlabel("SVs ≥50 bp (%)")
    despine(ax, left=False)
    ax.tick_params(left=False)


def panel_haplotype_window(fig, cell, rng: np.random.Generator) -> plt.Axes:
    """Phased haplotypes in one 200 kb window: SNVs, an SV, a VNTR allele."""
    sub = cell.subgridspec(1, 3, width_ratios=[0.035, 1.0, 0.32], wspace=0.04)
    ax_lai = fig.add_subplot(sub[0, 0])
    ax = fig.add_subplot(sub[0, 1])
    ax_tr = fig.add_subplot(sub[0, 2])
    n_hap = 48
    weights = np.array([0.30, 0.20, 0.32, 0.09, 0.09])
    anc = np.sort(rng.choice(len(LAI), size=n_hap, p=weights))
    snv_pos = np.sort(rng.uniform(2, 198, 70))
    base_af = rng.beta(0.6, 1.6, snv_pos.size)
    shift = rng.normal(0, 0.25, (len(LAI), snv_pos.size))
    del_af = np.array([0.45, 0.25, 0.10, 0.05, 0.15])
    alu_af = np.array([0.10, 0.30, 0.35, 0.02, 0.20])
    tr_mode = np.array([12, 6, 5, 6, 6])
    tr_len = np.maximum(3, rng.poisson(tr_mode[anc]) + (rng.random(n_hap) < 0.08) * rng.integers(8, 20, n_hap))
    has_del = rng.random(n_hap) < del_af[anc]
    has_alu = rng.random(n_hap) < alu_af[anc]
    order = np.lexsort((tr_len, has_del, anc))
    anc, tr_len, has_del, has_alu = anc[order], tr_len[order], has_del[order], has_alu[order]
    del_span = (84.0, 96.5)
    alu_x = 142.0
    tr_x = 170.0
    for i in range(n_hap):
        y = n_hap - 1 - i
        a = anc[i]
        p = np.clip(base_af + shift[a], 0.01, 0.99)
        alt = rng.random(snv_pos.size) < p
        ax.plot([0, 200], [y, y], color="#E4E4E4", lw=0.5, zorder=1)
        xs = snv_pos[alt]
        if has_del[i]:
            xs = xs[(xs < del_span[0]) | (xs > del_span[1])]
            ax.add_patch(Rectangle((del_span[0], y - 0.35), del_span[1] - del_span[0], 0.7, facecolor="#BFD7EA", edgecolor=SV_COLORS["DEL"], lw=0.4, zorder=2))
        ax.scatter(xs, np.full(xs.size, y), marker="|", s=9, lw=0.6, color="#333333", zorder=3)
        if has_alu[i]:
            ax.scatter([alu_x], [y], marker="v", s=9, color=SV_COLORS["INS"], lw=0, zorder=4)
        ax_lai.add_patch(Rectangle((0, y - 0.5), 1, 1, facecolor=ANC_COLORS[LAI[a]], edgecolor="none"))
    ax.axvline(tr_x, color="#B8B8B8", lw=0.4, ls="--")
    for x, lab in (((del_span[0] + del_span[1]) / 2, "DEL 12 kb"), (alu_x, "Alu INS"), (tr_x, "VNTR")):
        ax.text(x, n_hap + 0.3, lab, ha="center", va="bottom", fontsize=FS, color="#444444")
    ax.set_xlim(0, 200)
    ax.set_ylim(-0.8, n_hap - 0.2)
    ax.set_yticks([])
    ax.set_xticks([0, 50, 100, 150, 200])
    ax.set_xlabel("Position in window (kb)")
    despine(ax, left=False)
    ax_lai.set_xlim(0, 1)
    ax_lai.set_ylim(-0.8, n_hap - 0.2)
    ax_lai.set_axis_off()
    ax_lai.text(-0.6, n_hap / 2, f"{n_hap} phased haplotypes, local ancestry", rotation=90, ha="right", va="center", fontsize=FS, transform=ax_lai.transData)
    ys = n_hap - 1 - np.arange(n_hap)
    ax_tr.barh(ys, tr_len, height=0.8, color=[ANC_COLORS[LAI[a]] for a in anc], lw=0)
    ax_tr.set_ylim(-0.8, n_hap - 0.2)
    ax_tr.set_yticks([])
    ax_tr.set_xlabel("VNTR (RU)")
    ax_tr.set_title("allele length", fontsize=FS, pad=2)
    despine(ax_tr, left=False)
    return ax_lai


def anc_legend(fig, x: float, y: float) -> None:
    handles = [Patch(facecolor=ANC_COLORS[a], label=a) for a in LAI]
    fig.legend(handles=handles, loc="lower left", bbox_to_anchor=(x, y), ncol=5, frameon=False, fontsize=FS, handlelength=0.8, handletextpad=0.3, columnspacing=0.8)


# ---------------------------------------------------------------------------
# Plates
# ---------------------------------------------------------------------------

STAMP = "Length counts: integrated phased callset; context, Phase 1, phasing, alleles, expansions mocked; global ancestry real (chr1 FLARE)"


def letter(fig, ax, s: str, *, dx: float = 0.055, dy: float = 0.012, note: str | None = None) -> None:
    pos = ax.get_position()
    t = fig.text(pos.x0 - dx, pos.y1 + dy, s, fontsize=9, fontweight="bold", ha="left", va="bottom")
    if note:
        ax.annotate(note, xy=(1, 0), xycoords=t, xytext=(3, 0), textcoords="offset points",
                    fontsize=FS, ha="left", va="bottom", annotation_clip=False)


def _spectrum_letter(fig, ax_lead) -> None:
    pos = ax_lead.get_position()
    fig.text(pos.x0 - 0.055, pos.y1 + 0.024, "A", fontsize=9, fontweight="bold", ha="left", va="bottom")


def option_a(rng) -> tuple[Path, Path]:
    """A · Catalog → haplotypes → alleles (three tiers)."""
    data = spectrum_data(rng)
    fig = new_figure(8.9)
    gs = GridSpec(3, 1, figure=fig, left=0.075, right=0.975, top=0.935, bottom=0.05, height_ratios=[1.05, 1.0, 1.0], hspace=0.48)
    ax_a = draw_spectrum(fig, gs[0], data, tracks=(), p1=True, stack_ctx=True)
    _spectrum_letter(fig, ax_a)
    row2 = gs[1].subgridspec(1, 3, wspace=0.42, width_ratios=[1.1, 1.0, 1.0])
    ax = fig.add_subplot(row2[0, 0])
    panel_discovery_ebert(ax, rng)
    letter(fig, ax, "B")
    ax = fig.add_subplot(row2[0, 1])
    panel_allelic_series(ax)
    letter(fig, ax, "C")
    ax = fig.add_subplot(row2[0, 2])
    panel_phasing(ax)
    letter(fig, ax, "D")
    row3 = gs[2].subgridspec(1, 3, wspace=0.42, width_ratios=[1.1, 1.0, 1.0])
    ax = fig.add_subplot(row3[0, 0])
    panel_global_ancestry(ax)
    letter(fig, ax, "E", dy=0.02)
    ax = fig.add_subplot(row3[0, 1])
    panel_tr_by_lai(ax, rng)
    letter(fig, ax, "F", dy=0.02)
    ax = fig.add_subplot(row3[0, 2])
    panel_expansions(ax, rng)
    letter(fig, ax, "G", dx=0.1)
    stamp_mockup(fig, extra=STAMP)
    return save(fig, "fig2_option_a_tiers")


def _match_xlabel_baseline(fig, source: str, target: str) -> None:
    """Drop the xlabel containing ``target`` onto the baseline of the one containing ``source``."""
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    found = {}
    for ax in fig.axes:
        for key in (source, target):
            if key in ax.get_xlabel():
                found[key] = ax
    if len(found) < 2:
        return

    def baseline(ax):
        extent = ax.xaxis.label.get_window_extent(renderer)
        return extent.transformed(fig.transFigure.inverted()).y0

    gap_pt = (baseline(found[target]) - baseline(found[source])) * fig.get_figheight() * 72
    if gap_pt > 0.3:
        found[target].xaxis.labelpad += gap_pt


def option_b(rng, examples=None, stem: str = "fig2_option_b_length_spine", ratios: tuple[float, float] = (0.55, 1.50), summary: str = "new", row_scales: str = "shared", note: str | None = "(mock-up)") -> tuple[Path, Path]:
    """B · Context-stacked spectrum over discovery and the repeat allelic series.

    ``examples`` is the panel C locus list. None keeps two known loci and two
    high-PLVI stand-ins. The summary above them is Phase 2 lengths by carrier
    count, split into lengths Phase 1 already had and lengths it did not.
    """
    from .fig2_panel_c_options import panel_rarefaction_examples

    data = spectrum_data(rng, region_path=_region_bins_path())
    fig = new_figure(6.2)
    # The scope brackets used to sit in this gap. Without them the axis title
    # only needs a short pad, and the height goes back to panels B and C.
    gs = GridSpec(2, 1, figure=fig, left=0.09, right=0.975, top=0.922, bottom=0.12, height_ratios=[0.552, 1.2], hspace=0.34)
    ax_a = draw_spectrum(fig, gs[0], data, tracks=(), p1=True, stack_ctx=True, bnd_right=True, scope=True, scope_brackets=False, bin_notes=True)
    _spectrum_letter(fig, ax_a)
    # The symlog rows put their haplotype scale and gene names in the gutter, so C starts further left.
    row2 = gs[1].subgridspec(1, 2, wspace=0.30 if row_scales != "symlog" else 0.24, width_ratios=[1.0, 1.05 if row_scales != "symlog" else 1.12])
    ax = panel_discovery_structure(fig, row2[0, 0], rng)
    letter(fig, ax, "B", dx=0.075, dy=0.02)
    ax, ax_rows = panel_rarefaction_examples(fig, row2[0, 1], rng, examples=examples, ratios=ratios, summary=summary, row_scales=row_scales)
    letter(fig, ax, "C", dx=0.075, dy=0.02, note=note)
    if summary != "none":
        # D shares C's left edge and sits in the gap above the allelic series.
        fig.text(ax.get_position().x0 - 0.075, ax_rows.get_position().y1 + 0.012, "D",
                 fontsize=9, fontweight="bold", ha="left", va="bottom")
    _match_xlabel_baseline(fig, "Repeat units from reference", "Participants in discovery order")
    return save(fig, stem)


def option_c(rng) -> tuple[Path, Path]:
    """C · Phase 1 → Phase 2 in every panel."""
    data = spectrum_data(rng)
    fig = new_figure(8.6)
    gs = GridSpec(3, 1, figure=fig, left=0.085, right=0.975, top=0.935, bottom=0.05, height_ratios=[1.0, 1.0, 1.0], hspace=0.5)
    ax_a = draw_spectrum(fig, gs[0], data, tracks=("fold",), p1=True, stack_ctx=True)
    _spectrum_letter(fig, ax_a)
    row2 = gs[1].subgridspec(1, 3, wspace=0.45, width_ratios=[1.1, 1.0, 1.0])
    ax = fig.add_subplot(row2[0, 0])
    panel_discovery_ebert(ax, rng)
    letter(fig, ax, "B")
    ax = fig.add_subplot(row2[0, 1])
    panel_allelic_series(ax)
    letter(fig, ax, "C")
    ax = fig.add_subplot(row2[0, 2])
    panel_context_p1p2(ax)
    handles = [Patch(facecolor=CTX_COLORS[k], label=k) for k in CTX]
    ax.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.4, 1.0), ncol=4, fontsize=FS, handlelength=0.8, columnspacing=0.7)
    letter(fig, ax, "D", dx=0.1, dy=0.02)
    row3 = gs[2].subgridspec(1, 3, wspace=0.45, width_ratios=[1.1, 1.0, 1.0])
    ax = fig.add_subplot(row3[0, 0])
    panel_global_ancestry(ax)
    letter(fig, ax, "E", dy=0.02)
    ax = fig.add_subplot(row3[0, 1])
    panel_phasing(ax)
    letter(fig, ax, "F", dy=0.02)
    ax = fig.add_subplot(row3[0, 2])
    panel_expansions(ax, rng)
    letter(fig, ax, "G", dx=0.1)
    stamp_mockup(fig, extra=STAMP)
    return save(fig, "fig2_option_c_phase1_pairs")


def option_d(rng) -> tuple[Path, Path]:
    """D · Haplotype window as the centerpiece."""
    data = spectrum_data(rng)
    fig = new_figure(8.8)
    gs = GridSpec(3, 1, figure=fig, left=0.075, right=0.975, top=0.935, bottom=0.05, height_ratios=[0.85, 1.35, 0.9], hspace=0.45)
    ax_a = draw_spectrum(fig, gs[0], data, tracks=(), p1=False, stack_ctx=True)
    _spectrum_letter(fig, ax_a)
    row2 = gs[1].subgridspec(1, 2, wspace=0.3, width_ratios=[1.7, 1.0])
    ax_lai = panel_haplotype_window(fig, row2[0, 0], rng)
    letter(fig, ax_lai, "B", dx=0.04, dy=0.02)
    sub = row2[0, 1].subgridspec(2, 1, height_ratios=[0.75, 1.25], hspace=0.55)
    ax = fig.add_subplot(sub[0, 0])
    panel_global_ancestry(ax)
    ax.set_xlabel("Participants")
    letter(fig, ax, "C", dy=0.02)
    ax = fig.add_subplot(sub[1, 0])
    panel_tr_by_lai(ax, rng)
    ax.get_legend().remove()
    letter(fig, ax, "D")
    anc_legend(fig, ax_lai.get_position().x0 + 0.03, ax_lai.get_position().y1 + 0.03)
    row3 = gs[2].subgridspec(1, 3, wspace=0.45, width_ratios=[1.1, 1.0, 1.0])
    ax = fig.add_subplot(row3[0, 0])
    panel_discovery_ebert(ax, rng)
    letter(fig, ax, "E")
    ax = fig.add_subplot(row3[0, 1])
    panel_allelic_series(ax)
    letter(fig, ax, "F")
    ax = fig.add_subplot(row3[0, 2])
    panel_expansions(ax, rng)
    letter(fig, ax, "G", dx=0.1)
    stamp_mockup(fig, extra=STAMP)
    return save(fig, "fig2_option_d_haplotypes")


def contact_sheet(pngs: list[Path]) -> Path:
    titles = ("A · Three tiers", "B · Length spine", "C · Phase 1 → Phase 2")
    fig, axes = plt.subplots(1, len(pngs), figsize=(4 * len(pngs), 5.6))
    fig.subplots_adjust(left=0.01, right=0.99, top=0.93, bottom=0.01, wspace=0.03)
    for ax, png, title in zip(axes, pngs, titles):
        ax.imshow(plt.imread(png))
        ax.set_title(title, loc="left", fontsize=10)
        ax.axis("off")
    out = OUTDIR / "fig2_options_contact_sheet.png"
    fig.savefig(out, dpi=200)
    plt.close(fig)
    return out


def spectrum_panel(rng) -> tuple[Path, Path]:
    """Standalone panel A: context-stacked counts with the Phase 1 outline."""
    data = spectrum_data(rng)
    fig = new_figure(2.2)
    gs = GridSpec(1, 1, figure=fig, left=0.075, right=0.975, top=0.78, bottom=0.12)
    draw_spectrum(fig, gs[0], data, tracks=(), p1=True, stack_ctx=True)
    stamp_mockup(fig, extra="context fractions and Phase 1 mocked")
    return save(fig, "fig2_length_spectrum_ctx")


def render_main() -> tuple[Path, Path]:
    """Figure 2 for the manuscript: the symlog plate copied to figures/fig2_catalog.{pdf,png}."""
    import shutil

    from .fig2_panel_c_options import load_chosen_loci

    dest = MANUSCRIPT_ROOT / "figures"
    dest.mkdir(exist_ok=True)
    out = []
    for src in option_b(
        np.random.default_rng(7),
        examples=load_chosen_loci(),
        stem="fig2_option_b_symlog",
        ratios=(0.48, 1.62),
        summary="new",
        row_scales="symlog",
        note=None,
    ):
        dst = dest / f"fig2_catalog{Path(src).suffix}"
        shutil.copyfile(src, dst)
        out.append(dst)
    return out[0], out[1]


def render() -> list[tuple[Path, Path]]:
    out = [spectrum_panel(np.random.default_rng(7))]
    for fn in (option_a, option_b, option_c):
        out.append(fn(np.random.default_rng(7)))
    contact_sheet([png for _pdf, png in out[1:]])
    return out


if __name__ == "__main__":
    for pdf, png in render():
        print(png)
