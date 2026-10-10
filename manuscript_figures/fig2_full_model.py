"""Figure 2 full-model mock-ups: FELIXla → FELIXassoc narrative.

Exploration pass only. Every panel and plate is stamped MOCK-UP and uses
deterministic synthetic values for visual grammar. Do not treat numbers as
results. The existing fig2_associations mock-up is left untouched.
"""

from __future__ import annotations

from pathlib import Path
from collections import Counter

import numpy as np
from matplotlib.gridspec import GridSpec
from matplotlib.patches import FancyBboxPatch, Rectangle
import matplotlib.pyplot as plt

from .style import (
    ANC_COLORS,
    HG38_MB,
    OKABE,
    OUTDIR,
    P2_N,
    SV_COLORS,
    SV_N,
    ax_panel,
    despine,
    label_panel,
    new_figure,
    save,
    stamp_mockup,
)

# Real-data bins from variant_length/ WDL + merge (optional).
LENGTH_BINS_CANDIDATES = [
    OUTDIR / "data" / "fig2_length_spectrum.bins.tsv",
]

# Association FLARE panels (num_ancs = 5). No MID in the LAI model.
LAI = ["AFR", "AMR", "EUR", "EAS", "SAS"]

# Mock class counts — order-of-magnitude placeholders, all ≥20.
CLASS_COUNTS = {
    "SNV": 78_400_000,
    "indel": 9_600_000,
    "biallelic SV": SV_N["SV_50"],
    "multi-allelic SV": 185_000,
}

CLASS_COLORS = {
    "SNV": OKABE["blue"],
    "indel": OKABE["sky"],
    "biallelic SV": OKABE["orange"],
    "multi-allelic SV": OKABE["pink"],
    "DEL": SV_COLORS["DEL"],
    "INS": SV_COLORS["INS"],
    "other": SV_COLORS["INV"],
}

LENGTH_COLORS = {
    "DeepVariant": OKABE["sky"],
    "SV callset": OKABE["blue"],
    "Ultralong": OKABE["vermillion"],
}


def _fmt_count(n: float) -> str:
    if n >= 1_000_000:
        v = n / 1_000_000
        return f"{v:.1f}M" if v < 10 else f"{v:.0f}M"
    if n >= 1_000:
        return f"{n / 1_000:.0f}k"
    return f"{n:.0f}"


def _round_box(ax, x, y, w, h, text, color, fontsize=6.2, alpha=0.16):
    ax.add_patch(
        FancyBboxPatch(
            (x, y),
            w,
            h,
            boxstyle="round,pad=0.03,rounding_size=0.08",
            facecolor=color,
            edgecolor="none",
            alpha=alpha,
            mutation_aspect=0.55,
        )
    )
    ax.add_patch(
        FancyBboxPatch(
            (x, y),
            w,
            h,
            boxstyle="round,pad=0.03,rounding_size=0.08",
            facecolor="none",
            edgecolor=color,
            linewidth=0.9,
        )
    )
    ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", fontsize=fontsize, color="#222222")


def _arrow(ax, x0, x1, y, color="#444444"):
    ax.annotate(
        "",
        xy=(x1, y),
        xytext=(x0, y),
        arrowprops=dict(arrowstyle="-|>", color=color, lw=0.9, mutation_scale=8),
    )


def _paint_hap(length: float, mix: np.ndarray, mean_mb: float, flicker: float, rng: np.random.Generator):
    tracts = []
    x = 0.0
    mix = np.asarray(mix, dtype=float)
    mix = mix / mix.sum()
    while x < length - 1e-6:
        anc = int(rng.choice(len(LAI), p=mix))
        w = float(rng.exponential(mean_mb))
        if rng.random() < flicker:
            w = min(w, rng.uniform(0.04, 0.35))
        w = max(0.08, min(w, length - x))
        tracts.append((x, x + w, anc))
        x += w
    return tracts


def _draw_hap(ax, y, h, tracts):
    for x0, x1, j in tracts:
        ax.add_patch(
            Rectangle((x0, y), x1 - x0, h, facecolor=ANC_COLORS[LAI[j]], edgecolor="none", zorder=2)
        )


# ---------------------------------------------------------------------------
# Panel alternatives — variant classes
# ---------------------------------------------------------------------------


# Placeholder catalog sizes for the length spectrum.
LENGTH_CATALOG = {
    "snv": 95_000_000,
    "integrated": SV_N["SV_20"] + 12_900_000,  # phased SNV/indel/SV ≤10 kb
    "ultralong": 225_000,
    "bnd": 85_000,  # callout only — not on the length axis
}

# Muted fills for thin-bar overlays (Science-like, not slide-bright).
LENGTH_FILL = {
    "Integrated": "#5B6B9A",
    "Ultralong": "#B07A6A",
}

# Plot layers from encoding-based partitions. The integrated phased VCF is
# sequence REF/ALT alleles (no symbolic/SVTYPE structural sites in chr1–22
# merge); companion ultralong is the only separate length callset.
STACK_ORDER = ("integrated", "ultralong")
STACK_FILL = {
    "integrated": "#5B6B9A",
    "ultralong": "#B07A6A",
}
STACK_LABEL = {
    "integrated": "Integrated phased callset",
    "ultralong": "Ultralong (companion)",
}
BND_FILL = "#8C8C8C"
BND_LABEL = "Breakends (companion)"
# Partitions that belong to the integrated phased callset (≤10 kb).
_INTEGRATED_PARTITIONS = (
    "deepvariant",
    "main_sv",
    "sequence",
    "structural",
    "integrated",
)


def _simulate_length_spectrum(rng: np.random.Generator) -> dict[str, np.ndarray]:
    """Synthetic unsigned lengths for layout; peaks are visual cues only."""
    n = 120_000

    short = np.clip(rng.geometric(0.34, int(n * 0.55)) + 1, 1, 35).astype(float)
    mid = np.clip(rng.lognormal(np.log(48), 0.62, int(n * 0.20)), 36, 999)
    alu = rng.lognormal(np.log(310), 0.11, int(n * 0.12))
    l1 = rng.lognormal(np.log(6_000), 0.12, int(n * 0.08))
    other = np.clip(rng.lognormal(np.log(1_400), 0.35, int(n * 0.05)), 1_000, 9_999)
    integrated = np.concatenate([short, mid, alu, l1, other])
    integrated = np.clip(integrated, 1, 9_999)

    weights = rng.random(n)
    ultralong = np.empty(n)
    components = [
        (weights < 0.58, 24_000, 0.55),
        ((weights >= 0.58) & (weights < 0.88), 150_000, 0.75),
        (weights >= 0.88, 1_100_000, 0.85),
    ]
    for mask, median, sigma in components:
        ultralong[mask] = rng.lognormal(np.log(median), sigma, mask.sum())
    ultralong = np.clip(ultralong, 10_000, 2_000_000)

    def _split(values: np.ndarray, p_ins: float) -> tuple[np.ndarray, np.ndarray]:
        mask = rng.random(values.size) < p_ins
        return values[mask], values[~mask]

    int_ins, int_del = _split(integrated, 0.55)
    ul_ins, ul_del = _split(ultralong, 0.22)
    return {
        "int_ins": int_ins,
        "int_del": int_del,
        "ul_ins": ul_ins,
        "ul_del": ul_del,
    }


def _signed_series(data: dict[str, np.ndarray]) -> dict[str, np.ndarray]:
    """Deletions negative, insertions positive (reference-figure convention)."""
    return {
        "Integrated": np.concatenate([-data["int_del"], data["int_ins"]]),
        "Ultralong": np.concatenate([-data["ul_del"], data["ul_ins"]]),
    }


def _catalog_n(key: str) -> float:
    if key == "Integrated":
        return LENGTH_CATALOG["integrated"]
    return LENGTH_CATALOG["ultralong"]


def _thin_hist(
    ax,
    values: np.ndarray,
    n_true: float,
    n_draw: int,
    bins: np.ndarray,
    color: str,
) -> None:
    """Fine-binned step-filled histogram; scale uses full draw size, not the facet subset."""
    if values.size == 0:
        return
    counts, edges = np.histogram(values, bins=bins)
    scale = n_true / max(n_draw, 1)
    heights = np.clip(counts * scale, 0, None)
    widths = np.diff(edges)
    # Needle bars: draw as zero-edge rectangles so fine bins read as hair, not chunks.
    ax.bar(
        edges[:-1],
        heights,
        width=widths,
        align="edge",
        color=color,
        alpha=0.55,
        linewidth=0.0,
        edgecolor="none",
        zorder=2,
    )


def _length_facet_specs() -> list[dict]:
    """Absolute-length facets; INS on +y, DEL on −y. Callset split only at 10 kb."""
    return [
        {
            "title": r"$|$length$|$ $\leq$ 100",
            "xlim": (0, 100),
            "xticks": [0, 20, 50, 100],
            "xticklabels": ["0", "20", "50", "100"],
            "series": "Integrated",
            "layer": "integrated",
            "snv": True,
            "width": 1.6,
        },
        {
            "title": r"100 $<$ $|$length$|$ $<$ 10,000",
            "xlim": (100, 10_000),
            "xscale": "log",
            "xticks": [100, 310, 1_000, 6_000, 10_000],
            "xticklabels": ["100", "310", "1000", "6000", "10000"],
            "series": "Integrated",
            "layer": "integrated",
            "landmarks": ((310, "Alu"), (6_000, "L1")),
            "width": 2.0,
        },
        {
            "title": r"$|$length$|$ $\geq$ 10,000",
            "xlim": (10_000, 2.0e6),
            "xscale": "log",
            "xticks": [10_000, 1e6, 2e6],
            "xticklabels": [r"$10^{4}$", r"$10^{6}$", r"$2{\times}10^{6}$"],
            "series": "Ultralong",
            "layer": "ultralong",
            "width": 1.0,
        },
    ]


def _load_length_bins_tsv(path: Path) -> dict[str, object]:
    """Load merged length counts.

    New schema: ``chrom partition signed_len n_sites`` (1 bp).
    Legacy schema: ``chrom facet partition bin_index bin_lo bin_hi n_sites``.
    """
    from collections import defaultdict

    with path.open(encoding="utf-8") as fh:
        header = fh.readline()
        if header.startswith("chrom\tpartition\tsigned_len"):
            counts: dict[tuple[str, int], float] = Counter()
            for line in fh:
                parts = line.rstrip("\n").split("\t")
                if len(parts) != 4:
                    continue
                _chrom, partition, signed_len, n_sites = parts
                counts[(partition, int(signed_len))] += float(n_sites)
            return {"mode": "1bp", "counts": counts}

        if not header.startswith("chrom\tfacet"):
            raise ValueError(f"unexpected bins header in {path}: {header!r}")
        rows: dict[tuple[str, str], list[tuple[int, float, float, int]]] = defaultdict(list)
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) != 7:
                continue
            _chrom, facet, partition, bin_index, bin_lo, bin_hi, n_sites = parts
            rows[(facet, partition)].append((int(bin_index), float(bin_lo), float(bin_hi), int(n_sites)))
    out = {}
    for key, items in rows.items():
        items.sort()
        out[key] = (
            np.array([x[0] for x in items], dtype=int),
            np.array([x[1] for x in items], dtype=float),
            np.array([x[2] for x in items], dtype=float),
            np.array([x[3] for x in items], dtype=float),
        )
    return {"mode": "legacy", "rows": out}


def _length_bins_path() -> Path | None:
    for path in LENGTH_BINS_CANDIDATES:
        if path.suffix == ".tsv" and path.is_file():
            return path
    return None


def _iter_legacy_bars(rows, layer: str):
    """Yield (partition, bin_lo, bin_hi, n) from pre-faceted Terra bins."""
    facets = (
        "del_1k_10k",
        "del_indel_1k",
        "snv_indel",
        "ins_indel_1k",
        "ins_1k_10k",
        "le1k",
        "del_20_50",
        "ins_20_50",
        "le20",
        "del_50_1k",
        "ins_50_1k",
        "del_50_10k",
        "ins_50_10k",
        "le10k",
        "del_ge10k",
        "ins_ge10k",
    )
    parts = _INTEGRATED_PARTITIONS if layer == "integrated" else (layer,)
    for facet in facets:
        for part in parts:
            packed = rows.get((facet, part))
            if packed is None:
                continue
            _idx, bin_lo, bin_hi, n_sites = packed
            for a, b, n in zip(bin_lo, bin_hi, n_sites):
                if n > 0:
                    yield part, float(a), float(b), float(n)


def _facet_abs_edges(spec: dict) -> np.ndarray:
    lo, hi = spec["xlim"]
    if hi <= 100:
        start = 0.5 if lo <= 0 else lo + 0.5
        return np.arange(start, hi + 0.5 + 1e-9, 1.0)
    if hi <= 10_000:
        return np.geomspace(max(lo, 100), hi, 160)
    return np.geomspace(max(lo, 10_000), hi, 80)


def _accumulate_abs_counts(
    real: dict[str, object],
    layer: str,
    edges: np.ndarray,
) -> tuple[dict[str, np.ndarray], dict[str, np.ndarray], float]:
    """Bin onto absolute-length display edges, split by INS/DEL."""
    n_bin = len(edges) - 1
    ins = {p: np.zeros(n_bin, dtype=float) for p in STACK_ORDER}
    dele = {p: np.zeros(n_bin, dtype=float) for p in STACK_ORDER}
    n_snv = 0.0
    abs_lo, abs_hi = float(edges[0]), float(edges[-1])

    def _add_interval(part: str, lo_abs: float, hi_abs: float, n: float, target: dict[str, np.ndarray]) -> None:
        if part in _INTEGRATED_PARTITIONS:
            stack = "integrated"
        elif part == "ultralong":
            stack = "ultralong"
        else:
            return
        if hi_abs <= abs_lo or lo_abs >= abs_hi:
            return
        i0 = max(0, int(np.searchsorted(edges, lo_abs, side="right") - 1))
        i1 = min(n_bin - 1, int(np.searchsorted(edges, hi_abs, side="left") - 1))
        if i1 < i0:
            return
        overlaps = []
        for i in range(i0, i1 + 1):
            left = max(lo_abs, edges[i])
            right = min(hi_abs, edges[i + 1])
            overlaps.append(max(0.0, right - left))
        total = sum(overlaps)
        if total <= 0:
            return
        for i, w in zip(range(i0, i1 + 1), overlaps):
            if w > 0:
                target[stack][i] += n * (w / total)

    if real["mode"] == "1bp":
        counts: dict[tuple[str, int], float] = real["counts"]  # type: ignore[assignment]
        for (part, slen), n in counts.items():
            if layer == "integrated":
                if part not in _INTEGRATED_PARTITIONS:
                    continue
            elif part != layer:
                continue
            if slen == 0:
                n_snv += n
                continue
            if slen < 0:
                _add_interval(part, float(-slen) - 0.5, float(-slen) + 0.5, n, dele)
            else:
                _add_interval(part, float(slen) - 0.5, float(slen) + 0.5, n, ins)
    else:
        for part, a, b, n in _iter_legacy_bars(real["rows"], layer):
            if a < 0 <= b or a <= 0 < b:
                n_snv += n
                continue
            if b <= 0:
                _add_interval(part, -b, -a, n, dele)
            else:
                _add_interval(part, a, b, n, ins)
    return ins, dele, n_snv


def _draw_stack(ax, edges: np.ndarray, layers: dict[str, np.ndarray], *, sign: float) -> None:
    bottoms = np.zeros(len(edges) - 1, dtype=float)
    widths = np.diff(edges)
    for part in STACK_ORDER:
        h = layers.get(part)
        if h is None or not np.any(h):
            continue
        ax.bar(
            edges[:-1],
            sign * h,
            width=widths,
            bottom=sign * bottoms,
            align="edge",
            color=STACK_FILL[part],
            alpha=0.78,
            linewidth=0.0,
            edgecolor="none",
            zorder=2,
        )
        bottoms += h



def _draw_length_butterfly(fig: plt.Figure, rng: np.random.Generator, *, rect=(0.07, 0.16, 0.915, 0.78)) -> None:
    """Faceted length spectrum: INS on +y, DEL on −y (symlog); BND detached."""
    import json

    real_path = _length_bins_path()
    real = _load_length_bins_tsv(real_path) if real_path else None
    data = None if real else _signed_series(_simulate_length_spectrum(rng))
    specs = _length_facet_specs()
    x0, y0, w, h = rect
    # BND (undefined length) + spacer, then the length facets.
    widths = [0.22, 0.06] + [float(s.get("width", 1.0)) for s in specs]
    gs = GridSpec(
        1,
        len(widths),
        figure=fig,
        left=x0,
        right=x0 + w,
        bottom=y0,
        top=y0 + h,
        wspace=0.08,
        width_ratios=widths,
    )
    ax_bnd = fig.add_subplot(gs[0, 0])
    fig.add_subplot(gs[0, 1]).set_axis_off()
    axes = [fig.add_subplot(gs[0, i + 2]) for i in range(len(specs))]

    n_bnd = LENGTH_CATALOG["bnd"]
    if real_path is not None:
        summary_path = real_path.with_name(real_path.name.replace(".bins.tsv", ".summary.json"))
        if summary_path.is_file():
            summary = json.loads(summary_path.read_text(encoding="utf-8"))
            n_bnd = int(summary.get("n_bnd", n_bnd))

    # Precompute stacks so all facets share one y-scale.
    prepared: list[dict] = []
    y_max = 1.2e7
    for spec in specs:
        edges = _facet_abs_edges(spec)
        entry: dict = {"spec": spec, "edges": edges, "ins": None, "dele": None, "n_snv": 0.0}
        if real is not None:
            ins, dele, n_snv = _accumulate_abs_counts(real, spec["layer"], edges)
            entry.update(ins=ins, dele=dele, n_snv=n_snv)
            for layer in (*ins.values(), *dele.values()):
                if layer.size:
                    y_max = max(y_max, float(layer.max()) * 1.15)
            if n_snv > 0:
                y_max = max(y_max, n_snv * 1.05)
        prepared.append(entry)
    y_max = max(y_max, float(n_bnd) * 1.15)

    yticks = [1e2, 1e4, 1e6]
    for i, (ax, entry) in enumerate(zip(axes, prepared)):
        spec = entry["spec"]
        lo, hi = spec["xlim"]
        edges = entry["edges"]
        if real is not None:
            _draw_stack(ax, edges, entry["ins"], sign=1.0)
            _draw_stack(ax, edges, entry["dele"], sign=-1.0)
            if spec.get("snv") and entry["n_snv"] > 0:
                ax.vlines(
                    0.0,
                    -entry["n_snv"],
                    entry["n_snv"],
                    colors=STACK_FILL["integrated"],
                    lw=0.9,
                    alpha=0.85,
                    zorder=3,
                )
        else:
            key = spec["series"]
            color = LENGTH_FILL[key]
            vals = data[key]
            abs_v = np.abs(vals)
            mask = (abs_v >= lo) & (abs_v <= hi) & (vals != 0)
            for signed, flip in ((vals[mask & (vals > 0)], 1.0), (vals[mask & (vals < 0)], -1.0)):
                if signed.size == 0:
                    continue
                counts, e = np.histogram(np.abs(signed), bins=edges)
                scale = _catalog_n(key) / max(vals.size, 1)
                ax.bar(
                    e[:-1],
                    flip * counts * scale,
                    width=np.diff(e),
                    align="edge",
                    color=color,
                    alpha=0.55,
                    linewidth=0.0,
                    zorder=2,
                )
            if spec.get("snv"):
                ax.vlines(0.0, -LENGTH_CATALOG["snv"], LENGTH_CATALOG["snv"], colors=STACK_FILL["integrated"], lw=0.9, alpha=0.8, zorder=3)

        ax.axhline(0, color="#BBBBBB", lw=0.4, zorder=1)
        for x, lab in spec.get("landmarks", ()):
            ax.axvline(x, color="#B8B8B8", lw=0.45, ls="--", zorder=1)
            ax.text(x, 0.96, lab, transform=ax.get_xaxis_transform(), ha="center", va="top", fontsize=5.2, color="#666666")

        if spec.get("xscale"):
            ax.set_xscale(spec["xscale"])
        ax.set_xlim(lo, hi)
        ax.set_yscale("symlog", linthresh=1, linscale=0.35, base=10)
        ax.set_ylim(-y_max, y_max)
        ax.set_xticks(spec["xticks"])
        # Facet boundaries are shared with the previous panel; label them once.
        labels = list(spec["xticklabels"])
        if i > 0:
            labels[0] = ""
        ax.set_xticklabels(labels, fontsize=5.5)
        ax.set_title(spec["title"], fontsize=6.0, pad=3)
        despine(ax)
        ax.spines["bottom"].set_position(("data", 0))
        ax.tick_params(labelleft=False)
        ax.spines["left"].set_visible(False)
        ax.yaxis.set_ticks_position("none")

    # Breakends: undefined length and neither INS nor DEL, so one upward bar
    # outside the length axis.
    ax_bnd.bar([0.22], [n_bnd], width=0.55, color=BND_FILL, linewidth=0.0, zorder=2)
    ax_bnd.axhline(0, color="#BBBBBB", lw=0.4, zorder=1)
    ax_bnd.set_xlim(-0.6, 0.6)
    ax_bnd.set_yscale("symlog", linthresh=1, linscale=0.35, base=10)
    ax_bnd.set_ylim(-y_max, y_max)
    ax_bnd.set_xticks([0.22])
    ax_bnd.set_xticklabels(["BND"], fontsize=5.5)
    ax_bnd.set_title("no length", fontsize=6.0, pad=3)
    despine(ax_bnd)
    ax_bnd.spines["bottom"].set_position(("data", 0))
    ax_bnd.set_ylabel("Variant count", fontsize=6.0, labelpad=1.5)
    for y_frac, lab in ((0.94, "INS"), (0.06, "DEL")):
        axes[0].text(0.03, y_frac, lab, transform=axes[0].transAxes, ha="left", va="center", fontsize=5.2, color="#555555")
    ax_bnd.set_yticks([-t for t in reversed(yticks)] + [0] + yticks)
    ax_bnd.set_yticklabels(
        [rf"$-10^{{{int(np.log10(t))}}}$" if t >= 10 else str(int(-t)) for t in reversed(yticks)]
        + ["0"]
        + [rf"$10^{{{int(np.log10(t))}}}$" if t >= 10 else str(int(t)) for t in yticks],
        fontsize=5.2,
    )

    handles = [
        Rectangle((0, 0), 1, 1, facecolor=STACK_FILL[p], edgecolor="none", alpha=0.85, label=STACK_LABEL[p])
        for p in STACK_ORDER
    ] + [Rectangle((0, 0), 1, 1, facecolor=BND_FILL, edgecolor="none", label=BND_LABEL)]
    fig.legend(
        handles=handles,
        loc="lower center",
        ncol=3,
        frameon=False,
        fontsize=5.8,
        bbox_to_anchor=(0.5, 0.015),
        handlelength=1.1,
        columnspacing=1.4,
    )


def _panel_variant_ribbon(ax) -> None:
    """Compact four-class icon ribbon."""
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 3.2)
    ax.axis("off")
    ax.set_title("MOCK · Variant classes in the full model", loc="left", pad=2)

    items = [
        ("SNV", "A/G", CLASS_COLORS["SNV"], CLASS_COUNTS["SNV"]),
        ("indel", "ins/del\n<50 bp", CLASS_COLORS["indel"], CLASS_COUNTS["indel"]),
        ("biallelic SV", "DEL / INS\n≥50 bp", CLASS_COLORS["biallelic SV"], CLASS_COUNTS["biallelic SV"]),
        ("multi-allelic SV", "repeat /\ncopy-number", CLASS_COLORS["multi-allelic SV"], CLASS_COUNTS["multi-allelic SV"]),
    ]
    for i, (name, glyph, color, n) in enumerate(items):
        x = 0.35 + i * 2.4
        _round_box(ax, x, 1.05, 2.1, 1.55, f"{name}\n{glyph}\n{_fmt_count(n)} sites", color, fontsize=6.0)
        if i < 3:
            _arrow(ax, x + 2.15, x + 2.35, 1.82)
    ax.text(
        0.35,
        0.35,
        f"Phase 2 PacBio association cohort  ·  n = {P2_N:,}  ·  counts are mock placeholders",
        fontsize=5.8,
        color="#555555",
    )


def _panel_variant_counts(ax, rng: np.random.Generator) -> None:
    """Data-rich class counts + SV length densities."""
    fig = ax.figure
    ax.set_axis_off()
    pos = ax.get_position()
    # Left: counts. Right: length densities.
    ax_c = fig.add_axes([pos.x0, pos.y0, pos.width * 0.48, pos.height])
    ax_l = fig.add_axes([pos.x0 + pos.width * 0.56, pos.y0, pos.width * 0.44, pos.height])

    labels = list(CLASS_COUNTS)
    vals = np.array([CLASS_COUNTS[k] for k in labels], dtype=float)
    colors = [CLASS_COLORS[k] for k in labels]
    y = np.arange(len(labels))
    ax_c.barh(y, vals, color=colors, height=0.62, edgecolor="none")
    ax_c.set_yticks(y)
    ax_c.set_yticklabels(labels, fontsize=6.0)
    ax_c.set_xscale("log")
    ax_c.set_xlabel("Sites (mock)")
    ax_c.set_title("MOCK · Class counts", loc="left", pad=2)
    ax_c.invert_yaxis()
    despine(ax_c)
    for yi, v in zip(y, vals):
        ax_c.text(v * 1.08, yi, _fmt_count(v), va="center", fontsize=5.5, color="#333333")

    # Length densities for biallelic SV types (log10 bp).
    for typ, mu, color in [
        ("DEL", 2.6, CLASS_COLORS["DEL"]),
        ("INS", 2.4, CLASS_COLORS["INS"]),
        ("other", 3.1, CLASS_COLORS["other"]),
    ]:
        x = rng.normal(mu, 0.45, 4000)
        x = np.clip(x, 1.7, 5.5)
        dens, edges = np.histogram(x, bins=40, density=True)
        centers = 0.5 * (edges[:-1] + edges[1:])
        ax_l.plot(centers, dens, color=color, lw=1.1, label=typ)
        ax_l.fill_between(centers, dens, color=color, alpha=0.12, linewidth=0)
    ax_l.set_xlabel(r"$\log_{10}$ SV length (bp)")
    ax_l.set_ylabel("Density")
    ax_l.set_title("MOCK · SV length", loc="left", pad=2)
    ax_l.legend(loc="upper right", fontsize=5.4)
    despine(ax_l)


def _panel_variant_counts_simple(ax) -> None:
    """Horizontal bar counts only (fits tight layouts better than nested axes)."""
    labels = list(CLASS_COUNTS)
    vals = np.array([CLASS_COUNTS[k] for k in labels], dtype=float)
    colors = [CLASS_COLORS[k] for k in labels]
    y = np.arange(len(labels))
    ax.barh(y, vals, color=colors, height=0.62, edgecolor="none")
    ax.set_yticks(y)
    ax.set_yticklabels(labels)
    ax.set_xscale("log")
    ax.set_xlabel("Sites (mock)")
    ax.set_title("MOCK · Variant class counts", loc="left", pad=2)
    ax.invert_yaxis()
    despine(ax)
    for yi, v in zip(y, vals):
        ax.text(v * 1.08, yi, _fmt_count(v), va="center", fontsize=5.5, color="#333333")


def _panel_af_spectrum(ax, rng: np.random.Generator) -> None:
    """Optional AF spectrum strip per class."""
    bins = np.array([1e-4, 1e-3, 0.01, 0.05, 0.5, 1.0])
    labels = ["sing.", "rare", "LF", "common", "MA"]
    x = np.arange(len(labels))
    width = 0.18
    for i, (cls, shape) in enumerate(
        [
            ("SNV", (2.2, 8.0)),
            ("indel", (1.8, 6.0)),
            ("biallelic SV", (1.4, 9.0)),
            ("multi-allelic SV", (2.8, 4.0)),
        ]
    ):
        draws = rng.beta(shape[0], shape[1], 20_000)
        hist, _ = np.histogram(draws, bins=bins)
        hist = hist / hist.sum()
        ax.bar(
            x + (i - 1.5) * width,
            hist,
            width=width * 0.92,
            color=CLASS_COLORS[cls],
            label=cls,
            edgecolor="none",
        )
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.set_ylabel("Fraction of sites")
    ax.set_title("MOCK · Allele-frequency spectrum by class", loc="left", pad=2)
    despine(ax)
    ax.legend(loc="upper right", fontsize=5.2, ncol=2)


# ---------------------------------------------------------------------------
# Panel alternatives — local ancestry
# ---------------------------------------------------------------------------


def _panel_karyograms(ax, rng: np.random.Generator) -> None:
    people = [
        ("Person A · AMR-admixed", np.array([0.12, 0.48, 0.32, 0.04, 0.04])),
        ("Person B · AFR–EUR", np.array([0.55, 0.05, 0.38, 0.01, 0.01])),
        ("Person C · SAS–EUR", np.array([0.05, 0.08, 0.28, 0.04, 0.55])),
    ]
    length = HG38_MB[1]
    # Example variant overlays (Mb on chr1).
    markers = [
        (22.5, "SNV", CLASS_COLORS["SNV"]),
        (48.0, "indel", CLASS_COLORS["indel"]),
        (96.0, "DEL", CLASS_COLORS["DEL"]),
        (155.0, "INS", CLASS_COLORS["INS"]),
        (210.0, "mSV", CLASS_COLORS["multi-allelic SV"]),
    ]
    ax.set_xlim(-18, length + 4)
    ax.set_ylim(-0.6, 11.8)
    ax.axis("off")
    ax.set_title("MOCK · Local-ancestry painting with example variants (chr1)", loc="left", pad=2)
    y = 10.9
    for name, mix in people:
        mix = mix / mix.sum()
        ax.text(0, y + 0.18, name, fontsize=6.0, color="#222222", va="bottom")
        for hap_lab in ("h1", "h2"):
            ax.text(-1.5, y - 0.35, hap_lab, ha="right", va="center", fontsize=5.3, color="#555555")
            y -= 0.72
            tracts = _paint_hap(length, mix, mean_mb=22.0, flicker=0.04, rng=rng)
            _draw_hap(ax, y, 0.55, tracts)
        # Variant ticks on both haplotypes.
        for xpos, lab, color in markers:
            ax.plot([xpos, xpos], [y, y + 1.55], color=color, lw=0.7, zorder=4)
            ax.scatter([xpos], [y + 1.65], s=10, color=color, zorder=5, edgecolors="white", linewidths=0.25)
        y -= 0.55
    # Axis
    ax.plot([0, length], [0.35, 0.35], color="#222222", lw=0.5)
    for tick, lab in [(0, "0"), (50, "50"), (100, "100"), (150, "150"), (200, "200"), (248, "Mb")]:
        ax.plot([tick, tick], [0.25, 0.35], color="#222222", lw=0.5)
        ax.text(tick, -0.05, lab, ha="center", va="top", fontsize=5.3, color="#444444")
    # Legend for markers
    lx = 0.5
    for lab, color in [("SNV", CLASS_COLORS["SNV"]), ("indel", CLASS_COLORS["indel"]), ("DEL/INS", CLASS_COLORS["DEL"]), ("mSV", CLASS_COLORS["multi-allelic SV"])]:
        ax.scatter([lx], [-0.35], s=12, color=color, zorder=5)
        ax.text(lx + 3.5, -0.35, lab, fontsize=5.2, va="center", color="#333333")
        lx += 28
    # Switch-spanning SV callout (Methods, not central).
    ax.text(
        length - 2,
        10.6,
        "SVs spanning an ancestry\nswitch in a carrier: ignored\n(Methods audit)",
        ha="right",
        va="top",
        fontsize=5.2,
        color="#666666",
        style="italic",
    )


def _panel_dosage_zoom(ax) -> None:
    """Zoomed locus: two haplotypes → five ancestry dosage arms."""
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 6.2)
    ax.axis("off")
    ax.set_title("MOCK · FELIXla: haplotype × local ancestry → dosage arms", loc="left", pad=2)

    # Haplotype tracks; single focal variant at POS = 5.0.
    segs_h1 = [(0.4, 3.2, "AFR"), (3.2, 9.4, "EUR")]
    segs_h2 = [(0.4, 4.8, "AMR"), (4.8, 9.4, "EUR")]
    alleles = [("hap 1", 4.7, "ALT"), ("hap 2", 3.7, "ALT")]
    for lab, y, allele in alleles:
        segs = segs_h1 if lab.endswith("1") else segs_h2
        ax.text(0.15, y + 0.28, lab, fontsize=5.6, ha="left", va="center", color="#444444")
        for x0, x1, anc in segs:
            ax.add_patch(
                Rectangle((x0, y), x1 - x0, 0.55, facecolor=ANC_COLORS[anc], edgecolor="none", alpha=0.9)
            )
            ax.text((x0 + x1) / 2, y + 0.27, anc, ha="center", va="center", fontsize=5.4, color="white", fontweight="bold")
        color = OKABE["vermillion"] if allele == "ALT" else "#555555"
        ax.scatter([5.0], [y + 0.27], s=28, color=color, zorder=4, edgecolors="white", linewidths=0.3)
        ax.text(5.0, y + 0.72, allele, ha="center", va="bottom", fontsize=5.2, color=color)

    ax.annotate(
        "",
        xy=(5.0, 2.85),
        xytext=(5.0, 3.55),
        arrowprops=dict(arrowstyle="-|>", color="#333333", lw=0.9, mutation_scale=8),
    )
    ax.text(5.25, 3.15, "variant at POS", fontsize=5.4, color="#333333", va="center")

    # Dosage columns.
    dosages = [("DS_AFR", 0, OKABE["blue"]), ("DS_AMR", 0, ANC_COLORS["AMR"]), ("DS_EUR", 2, ANC_COLORS["EUR"]), ("DS_EAS", 0, ANC_COLORS["EAS"]), ("DS_SAS", 0, ANC_COLORS["SAS"]), ("DSALL", 2, "#333333")]
    for i, (name, val, color) in enumerate(dosages):
        x = 0.5 + i * 1.55
        _round_box(ax, x, 0.55, 1.35, 1.7, f"{name}\n\n{val}", color, fontsize=6.0, alpha=0.14)
    ax.text(
        0.4,
        0.2,
        "K = 5 FLARE panels  ·  ALT on EUR hap1 + EUR hap2 → DS_EUR = 2, DSALL = 2",
        fontsize=5.5,
        color="#555555",
    )


def _panel_admixed_summary(ax, rng: np.random.Generator) -> None:
    """Cohort local-ancestry proportions + retained-beyond-clustering concept."""
    # Ternary-ish scatter of max global prob vs second component — simplified:
    # histogram of max ancestry proportion, with a threshold line.
    max_p = np.concatenate(
        [
            rng.beta(12, 2, 7_000),  # mostly continental
            rng.beta(3.5, 3.5, 5_261),  # admixed-ish
        ]
    )
    max_p = np.clip(max_p, 0.35, 0.999)
    ax.hist(max_p, bins=36, color=OKABE["blue"], alpha=0.75, edgecolor="none")
    ymax = ax.get_ylim()[1]
    for thr, ls in [(0.8, "--"), (0.9, ":"), (0.95, "-.")]:
        ax.axvline(thr, color=OKABE["vermillion"], ls=ls, lw=0.9)
        ax.text(
            thr - 0.012,
            ymax * 0.96,
            f"max p < {thr:g}",
            fontsize=5.1,
            color=OKABE["vermillion"],
            rotation=90,
            va="top",
            ha="right",
        )
    # Mock retained fraction at 0.8 — large enough to avoid small-count issues.
    ax.text(
        0.38,
        ymax * 0.72,
        "MOCK retained beyond\nglobal clustering\n(max p < 0.8): ~18%\n(FELIX preprint: 12.1%)",
        fontsize=5.6,
        color="#333333",
        bbox=dict(boxstyle="round,pad=0.25", facecolor="#F7F7F7", edgecolor="#CCCCCC", linewidth=0.5),
    )
    ax.set_xlim(0.35, 1.0)
    ax.set_xlabel("Max global ancestry probability")
    ax.set_ylabel("Participants")
    ax.set_title("MOCK · Admixed participants retained by local-ancestry model", loc="left", pad=2)
    despine(ax)


# ---------------------------------------------------------------------------
# Panel alternatives — multi-allelic placeholder
# ---------------------------------------------------------------------------


def _panel_multiallelic(ax, rng: np.random.Generator) -> None:
    loci = [
        ("Locus 1 · VNTR", [0.42, 0.28, 0.18, 0.08, 0.04]),
        ("Locus 2 · STR", [0.55, 0.22, 0.12, 0.07, 0.04]),
        ("Locus 3 · CNV", [0.35, 0.25, 0.20, 0.12, 0.08]),
    ]
    ax.set_xlim(0, 11.5)
    ax.set_ylim(0, 4.6)
    ax.axis("off")
    ax.set_title("MOCK · Multi-allelic SV allele distributions (encoding TBD)", loc="left", pad=2)

    for i, (name, fracs) in enumerate(loci):
        x0 = 0.4 + i * 2.8
        ax.text(x0 + 1.0, 4.15, name, ha="center", fontsize=5.8, color="#222222")
        bottom = 0.9
        alleles = list(range(len(fracs)))
        colors = [plt.cm.PuOr(0.15 + 0.7 * j / max(1, len(fracs) - 1)) for j in alleles]
        # Vertical stacked bar of allele fractions (cohort).
        for j, (f, c) in enumerate(zip(fracs, colors)):
            h = f * 2.8
            ax.add_patch(Rectangle((x0 + 0.35, bottom), 1.3, h, facecolor=c, edgecolor="white", linewidth=0.4))
            if f >= 0.08:
                ax.text(x0 + 1.0, bottom + h / 2, f"a{j}", ha="center", va="center", fontsize=5.2, color="#222222")
            bottom += h
        # Tiny ancestry-split sparkline (mock, all cells ≥20 conceptually via %).
        ay = 0.35
        for k, anc in enumerate(["AFR", "EUR", "AMR"]):
            ax.add_patch(
                Rectangle((x0 + 0.2 + k * 0.7, ay), 0.55, 0.35, facecolor=ANC_COLORS[anc], edgecolor="none", alpha=0.85)
            )
        ax.text(x0 + 1.0, 0.12, "by local anc. (%)", ha="center", fontsize=4.8, color="#666666")

    # Placeholder box.
    _round_box(
        ax,
        8.7,
        1.0,
        2.5,
        2.8,
        "→ FELIX\ndosage\ncolumns\n\nTBD",
        OKABE["grey"],
        fontsize=6.4,
        alpha=0.10,
    )
    # Dashed emphasis.
    ax.add_patch(
        FancyBboxPatch(
            (8.7, 1.0),
            2.5,
            2.8,
            boxstyle="round,pad=0.03,rounding_size=0.08",
            facecolor="none",
            edgecolor="#666666",
            linewidth=1.0,
            linestyle=(0, (2, 2)),
        )
    )
    ax.text(9.95, 0.55, "encoding not invented", ha="center", fontsize=5.2, color="#666666", style="italic")


# ---------------------------------------------------------------------------
# Panel alternatives — FELIX association tests
# ---------------------------------------------------------------------------


def _panel_felix_tests(ax) -> None:
    """Shared vs ancestry-specific effects + CCT combination."""
    ax.set_xlim(0, 12)
    ax.set_ylim(0, 6.5)
    ax.axis("off")
    ax.set_title("MOCK · FELIXassoc: shared, ancestry-specific, and combined tests", loc="left", pad=2)

    # Two example loci side by side.
    examples = [
        (
            "Shared-effect locus",
            [("AFR", 0.42), ("AMR", 0.38), ("EUR", 0.40), ("EAS", 0.35), ("SAS", 0.37)],
            r"$P_{\mathrm{hom}} \ll P_{\mathrm{het}}$",
            OKABE["green"],
        ),
        (
            "Ancestry-specific locus",
            [("AFR", 0.62), ("AMR", 0.10), ("EUR", -0.05), ("EAS", 0.08), ("SAS", 0.12)],
            r"$P_{\mathrm{het}} \ll P_{\mathrm{hom}}$",
            OKABE["vermillion"],
        ),
    ]
    for i, (title, betas, note, accent) in enumerate(examples):
        x0 = 0.4 + i * 5.8
        ax.text(x0 + 2.4, 6.15, title, ha="center", fontsize=6.2, fontweight="bold", color=accent)
        # Forest-like points.
        ax.plot([x0 + 2.4, x0 + 2.4], [5.5, 5.5 - 4 * 0.55], color="#CCCCCC", lw=0.6, zorder=1)
        for j, (anc, beta) in enumerate(betas):
            y = 5.4 - j * 0.55
            ax.plot([x0 + 1.2, x0 + 3.6], [y, y], color="#EEEEEE", lw=0.6, zorder=1)
            se = 0.08
            ax.plot([x0 + 2.4 + beta - se, x0 + 2.4 + beta + se], [y, y], color=ANC_COLORS[anc], lw=1.2, zorder=2)
            ax.scatter([x0 + 2.4 + beta], [y], s=16, color=ANC_COLORS[anc], zorder=3, edgecolors="white", linewidths=0.3)
            ax.text(x0 + 0.9, y, anc, ha="right", va="center", fontsize=5.3, color="#333333")
        ax.text(x0 + 2.4, 2.55, note, ha="center", fontsize=5.5, color=accent)
        _round_box(ax, x0 + 1.3, 1.55, 2.2, 0.75, "P_hom   P_het", accent, fontsize=5.8, alpha=0.12)

    # CCT funnel.
    ax.annotate(
        "",
        xy=(6.0, 0.95),
        xytext=(3.2, 1.55),
        arrowprops=dict(arrowstyle="-|>", color="#444444", lw=0.85, mutation_scale=8),
    )
    ax.annotate(
        "",
        xy=(6.0, 0.95),
        xytext=(8.8, 1.55),
        arrowprops=dict(arrowstyle="-|>", color="#444444", lw=0.85, mutation_scale=8),
    )
    _round_box(ax, 4.4, 0.25, 3.2, 0.85, "P_cct_admixed_c\n(Cauchy combination)", OKABE["blue"], fontsize=6.0)
    ax.text(
        11.7,
        0.4,
        "No categorical\n“model selected”\nfield in FELIX",
        ha="right",
        va="bottom",
        fontsize=5.2,
        color="#666666",
        style="italic",
    )


def _panel_felix_compact(ax) -> None:
    """Compact one-row version for hybrid top/bottom layouts."""
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 3.4)
    ax.axis("off")
    ax.set_title("MOCK · Adaptive shared vs ancestry-specific tests", loc="left", pad=2)

    # Shared betas (flat).
    ax.text(1.6, 3.05, "Shared", ha="center", fontsize=6.0, color=OKABE["green"], fontweight="bold")
    for i, (anc, b) in enumerate([("AFR", 0.4), ("AMR", 0.38), ("EUR", 0.41), ("EAS", 0.36), ("SAS", 0.39)]):
        ax.scatter([0.6 + i * 0.5], [2.3 + b * 0.8], s=14, color=ANC_COLORS[anc], zorder=3)
        ax.plot([0.6 + i * 0.5, 0.6 + i * 0.5], [2.3, 2.3 + b * 0.8], color=ANC_COLORS[anc], lw=1.0)
    ax.axhline(2.3, color="#CCCCCC", lw=0.5)
    ax.text(1.6, 1.85, r"$P_{\mathrm{hom}}$", ha="center", fontsize=6.0, color=OKABE["green"])

    # Heterogeneous.
    ax.text(5.2, 3.05, "Ancestry-specific", ha="center", fontsize=6.0, color=OKABE["vermillion"], fontweight="bold")
    for i, (anc, b) in enumerate([("AFR", 0.7), ("AMR", 0.15), ("EUR", -0.1), ("EAS", 0.05), ("SAS", 0.12)]):
        ax.scatter([4.2 + i * 0.5], [2.3 + b * 0.8], s=14, color=ANC_COLORS[anc], zorder=3)
        ax.plot([4.2 + i * 0.5, 4.2 + i * 0.5], [2.3, 2.3 + b * 0.8], color=ANC_COLORS[anc], lw=1.0)
    ax.axhline(2.3, color="#CCCCCC", lw=0.5)
    ax.text(5.2, 1.85, r"$P_{\mathrm{het}}$", ha="center", fontsize=6.0, color=OKABE["vermillion"])

    _arrow(ax, 3.3, 3.7, 2.5)
    _arrow(ax, 6.9, 7.3, 2.5)
    _round_box(ax, 7.4, 1.7, 2.3, 1.4, r"$P_{\mathrm{cct}}$" + "\nCauchy", OKABE["blue"], fontsize=6.4)
    ax.text(
        0.3,
        0.45,
        "Both tests always computed; CCT is the primary association statistic.",
        fontsize=5.5,
        color="#555555",
    )
    # Ancestry legend.
    for i, anc in enumerate(LAI):
        ax.add_patch(Rectangle((0.3 + i * 0.85, 0.05), 0.25, 0.22, facecolor=ANC_COLORS[anc], edgecolor="none"))
        ax.text(0.58 + i * 0.85, 0.16, anc, fontsize=5.0, va="center", color="#333333")


# ---------------------------------------------------------------------------
# Pipeline schematic pieces
# ---------------------------------------------------------------------------


def _panel_pipeline_flow(ax) -> None:
    ax.set_xlim(0, 12)
    ax.set_ylim(0, 3.6)
    ax.axis("off")
    ax.set_title("MOCK · Full association model", loc="left", pad=2)
    boxes = [
        (0.2, 1.35, 2.5, 1.5, "Variant classes\nSNV · indel\nbiallelic · multi-allelic", CLASS_COLORS["biallelic SV"]),
        (3.5, 1.35, 2.7, 1.5, "FELIXla\nlocal-ancestry\ndosage (K = 5)", OKABE["blue"]),
        (7.0, 1.35, 2.5, 1.5, "FELIXassoc\n$P_{hom}$ · $P_{het}$\n→ $P_{cct}$", OKABE["green"]),
        (10.0, 1.35, 1.8, 1.5, "Hits\n+ showcase\nloci", OKABE["pink"]),
    ]
    for x, y, w, h, text, color in boxes:
        _round_box(ax, x, y, w, h, text, color, fontsize=6.0)
    for x0, x1 in [(2.75, 3.45), (6.25, 6.95), (9.55, 9.95)]:
        _arrow(ax, x0, x1, 2.1)
    # Multi-allelic callout under first arrow.
    ax.add_patch(
        FancyBboxPatch(
            (3.6, 0.25),
            2.5,
            0.85,
            boxstyle="round,pad=0.02,rounding_size=0.06",
            facecolor="#FAFAFA",
            edgecolor="#888888",
            linewidth=0.8,
            linestyle=(0, (2, 2)),
        )
    )
    ax.text(4.85, 0.67, "multi-allelic dosage\nencoding — TBD", ha="center", va="center", fontsize=5.4, color="#555555")
    ax.text(
        0.2,
        0.35,
        f"Null covariates (Fig. 1): age, sex, PCs, coverage, genome center  ·  n = {P2_N:,}",
        fontsize=5.5,
        color="#555555",
    )


# ---------------------------------------------------------------------------
# Full-figure versions
# ---------------------------------------------------------------------------


def _render_version1_pipeline() -> tuple[Path, Path]:
    rng = np.random.default_rng(11)
    fig = new_figure(8.6)
    gs = GridSpec(
        4,
        2,
        figure=fig,
        left=0.07,
        right=0.985,
        top=0.97,
        bottom=0.04,
        wspace=0.28,
        hspace=0.42,
        height_ratios=[0.85, 1.25, 1.15, 1.35],
    )
    ax_a = fig.add_subplot(gs[0, :])
    ax_b = fig.add_subplot(gs[1, 0])
    ax_c = fig.add_subplot(gs[1, 1])
    ax_d = fig.add_subplot(gs[2, :])
    ax_e = fig.add_subplot(gs[3, :])

    _panel_pipeline_flow(ax_a)
    label_panel(ax_a, "A", x=-0.01, y=1.08)
    _panel_variant_ribbon(ax_b)
    label_panel(ax_b, "B", x=-0.04, y=1.06)
    _panel_dosage_zoom(ax_c)
    label_panel(ax_c, "C", x=-0.04, y=1.06)
    _panel_multiallelic(ax_d, rng)
    label_panel(ax_d, "D", x=-0.01, y=1.06)
    _panel_felix_tests(ax_e)
    label_panel(ax_e, "E", x=-0.01, y=1.04)

    stamp_mockup(fig, extra="Version 1 — pipeline schematic")
    return save(fig, "fig2_full_model_v1_pipeline")


def _render_version2_evidence() -> tuple[Path, Path]:
    rng = np.random.default_rng(11)
    fig = new_figure(9.0)
    gs = GridSpec(
        3,
        2,
        figure=fig,
        left=0.08,
        right=0.985,
        top=0.97,
        bottom=0.04,
        wspace=0.30,
        hspace=0.38,
        height_ratios=[1.05, 1.25, 1.35],
    )
    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[0, 1])
    ax_c = fig.add_subplot(gs[1, 0])
    ax_d = fig.add_subplot(gs[1, 1])
    ax_e = fig.add_subplot(gs[2, :])

    _panel_variant_counts_simple(ax_a)
    label_panel(ax_a, "A", x=-0.08, y=1.05)
    _panel_af_spectrum(ax_b, rng)
    label_panel(ax_b, "B", x=-0.08, y=1.05)
    _panel_karyograms(ax_c, rng)
    label_panel(ax_c, "C", x=-0.04, y=1.04)
    _panel_admixed_summary(ax_d, rng)
    label_panel(ax_d, "D", x=-0.08, y=1.05)
    # Bottom: multi-allelic + FELIX side by side via nested.
    # Redraw as split bottom row instead.
    ax_e.remove()
    gs_b = gs[2, :].subgridspec(1, 2, wspace=0.28)
    ax_e1 = fig.add_subplot(gs_b[0, 0])
    ax_e2 = fig.add_subplot(gs_b[0, 1])
    _panel_multiallelic(ax_e1, rng)
    label_panel(ax_e1, "E", x=-0.04, y=1.05)
    _panel_felix_compact(ax_e2)
    label_panel(ax_e2, "F", x=-0.04, y=1.05)

    stamp_mockup(fig, extra="Version 2 — evidence-first")
    return save(fig, "fig2_full_model_v2_evidence")


def _render_length_butterfly_panel(rng: np.random.Generator) -> tuple[Path, Path]:
    """Primary full-width length-spectrum panel (faceted signed butterfly)."""
    from .style import OUTDIR, PANELDIR

    fig = new_figure(1.375)
    _draw_length_butterfly(fig, rng, rect=(0.07, 0.24, 0.915, 0.66))
    # Real Terra bins → no MOCK stamp; synthetic fallback keeps the stamp.
    if _length_bins_path() is None:
        stamp_mockup(fig)
    OUTDIR.mkdir(parents=True, exist_ok=True)
    PANELDIR.mkdir(parents=True, exist_ok=True)
    pdf = OUTDIR / "fig2_length_spectrum.pdf"
    png = OUTDIR / "fig2_length_spectrum.png"
    fig.savefig(pdf, dpi=600, bbox_inches="tight", pad_inches=0.04)
    fig.savefig(png, dpi=400, bbox_inches="tight", pad_inches=0.04)
    # Mirror into panels/ for remixing.
    fig.savefig(PANELDIR / "fig2_MOCK_variant_length_butterfly.pdf", dpi=600, bbox_inches="tight", pad_inches=0.05)
    fig.savefig(PANELDIR / "fig2_MOCK_variant_length_butterfly.png", dpi=400, bbox_inches="tight", pad_inches=0.05)
    plt.close(fig)
    return pdf, png


def _render_comparison_sheet(paths: list[Path]) -> Path:
    apply = __import__("matplotlib.pyplot", fromlist=["plt"])
    from .style import apply_style, OUTDIR

    apply_style()
    fig, axes = plt.subplots(1, 2, figsize=(8.4, 5.2))
    fig.subplots_adjust(left=0.02, right=0.98, top=0.90, bottom=0.04, wspace=0.04)
    titles = [
        "V1  Pipeline schematic",
        "V2  Evidence-first",
    ]
    for ax, path, title in zip(axes, paths, titles):
        img = plt.imread(path)
        ax.imshow(img)
        ax.set_title(title, loc="left", fontsize=8, pad=3)
        ax.axis("off")
    fig.suptitle("MOCK · Figure 2 full-model layout options (V3 dropped)", fontsize=10, fontweight="bold", y=0.98)
    stamp_mockup(fig, extra="length spectrum iterated separately as butterfly")
    OUTDIR.mkdir(parents=True, exist_ok=True)
    out = OUTDIR / "fig2_full_model_comparison.png"
    fig.savefig(out, dpi=220, bbox_inches="tight", pad_inches=0.04)
    fig.savefig(OUTDIR / "fig2_full_model_comparison.pdf", dpi=400, bbox_inches="tight", pad_inches=0.04)
    plt.close(fig)
    return out


def _write_requirements() -> Path:
    from .style import OUTDIR

    OUTDIR.mkdir(parents=True, exist_ok=True)
    path = OUTDIR / "fig2_full_model_requirements.md"
    path.write_text(
        """# Figure 2 full-model mock-up — data requirements

Exploration artifact. Numbers in the mock panels are **not** results.
V3 plate and ridgeline length spectrum were dropped; length panel is the butterfly.

Disclosure rule: any count in \\{1 … 19\\} must appear as `<20` (or be aggregated).

## Visible elements → required fields

| Element | Mock panel(s) | Fields / statistic needed | Existing producer | New work? | Disclosure-sensitive? | Droppable? |
|---|---|---|---|---|---|---|
| Variant length spectrum + analysis scope | `fig2_length_spectrum` butterfly | Length histogram for DeepVariant SNVs/indels (≤35 bp mocked), main SVs (20 bp–<10 kb), ultralong (≥10 kb); INS/DEL split; BND count as footnote | SNV/indel VCF + `AnnotateSvCallset` site table + ultralong/BND callsets | One harmonized binned summary | No (histograms and total counts) | No — primary length panel |
| Class counts (SNV, indel, biallelic SV, multi-allelic SV) | ribbon, counts bars (V1/V2) | Site counts by class; biallelic SV split DEL/INS/other | `sv_annotation` site table + `BcftoolsGlnexusStats` merge | Thin summarizer script | No if totals ≫20 | Likely subsumed by length butterfly |
| AF spectrum by class | V2 panel B | Per-site AF (or AC/AN) by class | Site table `af`; SNV/indel from bcftools or VCF INFO | Summarizer only | No if binned | **Yes** |
| Karyogram painting (2–3 people) | V2 C | FLARE `AN1`/`AN2` tracts for allow-listed admixed IDs; example variant POS by class | `FlareByPopulation` `*.anc.vcf.gz` | Small LAI extract for ≤3 samples | Sample IDs must not be released; use anonymous labels | Zoom dosage panel can replace |
| Dosage-arm schematic | V1 C | Conceptual: genotype × ancestry → `DS1…DS5`, `DSALL` | None (schematic) | No data | No | Keep as schematic even with real data |
| Admixed retained fraction | V2 D | Join cohort to `anc_df` `prob1–prob6`; fraction with `max(prob)<{0.8,0.9,0.95}` | `covariates.v8` + `anc_df.csv.gz` | Local script | Report % / private counts ≥20 | Threshold choice is manuscript decision |
| Multi-allelic allele distributions | V1 D, V2 E | 2–3 loci with several common alleles; allele counts; optional split by local ancestry | Repeat / RU-annotated SV VCF; FLARE on site | Example extractor | **Yes** for per-ancestry carrier counts | Encoding box stays TBD |
| Shared vs ancestry-specific forest | V1 E, V2 F | Per-ancestry `BETA_c_anc{j}`; `P_hom_admixed_c`, `P_het_admixed_c`, `P_cct_admixed_c` | `FelixGenome` / `summarize_felix_results.py` | Thin classifier optional | Effect plots OK; n-case cells may be sensitive | Do **not** invent a “model selected” chip |

## Real-data swap-in for the butterfly

Produce one binned TSV (or three partition files) with:
`partition` ∈ {deepvariant, main_sv, ultralong}, `direction` ∈ {INS, DEL, SNV},
and either per-site `svlen_bp` or pre-binned `bin_lo, bin_hi, n_sites`.
Exclude BNDs from the length histogram; report BND site count separately.
Confirm DeepVariant endpoint (mocked ≤35 bp) before final rendering.

## Mock files written

| File | Role |
|---|---|
| `fig2_length_spectrum.*` | Primary butterfly length panel |
| `panels/fig2_MOCK_variant_length_butterfly.*` | Same panel under panels/ |
| `panels/fig2_MOCK_*.pdf/png` | Other remixable panels |
| `fig2_full_model_v1_pipeline.*` | Version 1 plate |
| `fig2_full_model_v2_evidence.*` | Version 2 plate |
| `fig2_full_model_comparison.*` | V1 vs V2 only |
| `fig2_full_model_requirements.md` | This file |
""",
        encoding="utf-8",
    )
    return path


def _export_individual_panels(rng: np.random.Generator) -> list[tuple[str, Path, Path]]:
    out = []
    specs = [
        ("fig2_MOCK_variant_ribbon", 1.55, _panel_variant_ribbon, dict(left=0.03, right=0.99, top=0.86, bottom=0.08)),
        ("fig2_MOCK_variant_counts", 2.4, _panel_variant_counts_simple, dict(left=0.22, right=0.96, top=0.88, bottom=0.16)),
        ("fig2_MOCK_af_spectrum", 2.5, lambda ax: _panel_af_spectrum(ax, rng), dict(left=0.12, right=0.98, top=0.88, bottom=0.16)),
        ("fig2_MOCK_karyograms", 3.2, lambda ax: _panel_karyograms(ax, rng), dict(left=0.05, right=0.99, top=0.90, bottom=0.10)),
        ("fig2_MOCK_dosage_zoom", 2.8, _panel_dosage_zoom, dict(left=0.04, right=0.99, top=0.88, bottom=0.06)),
        ("fig2_MOCK_admixed_summary", 2.6, lambda ax: _panel_admixed_summary(ax, rng), dict(left=0.14, right=0.96, top=0.88, bottom=0.16)),
        ("fig2_MOCK_multiallelic", 2.7, lambda ax: _panel_multiallelic(ax, rng), dict(left=0.03, right=0.99, top=0.88, bottom=0.06)),
        ("fig2_MOCK_felix_tests", 3.3, _panel_felix_tests, dict(left=0.03, right=0.99, top=0.90, bottom=0.05)),
        ("fig2_MOCK_felix_compact", 1.9, _panel_felix_compact, dict(left=0.03, right=0.99, top=0.86, bottom=0.08)),
        ("fig2_MOCK_pipeline_flow", 1.7, _panel_pipeline_flow, dict(left=0.03, right=0.99, top=0.86, bottom=0.06)),
    ]
    for stem, height, drawer, kw in specs:
        pdf, png = ax_panel(stem, height, drawer, **kw)
        out.append((stem, pdf, png))
        print(f"wrote {png}")
    # Nested counts+length panel as its own figure (uses add_axes).
    fig = new_figure(2.8)
    ax = fig.add_axes([0.0, 0.0, 1.0, 1.0])
    ax.set_axis_off()
    _panel_variant_counts(ax, rng)
    stamp_mockup(fig)
    from .style import save_panel

    pdf, png = save_panel(fig, "fig2_MOCK_variant_counts_length", stamp=False)
    out.append(("fig2_MOCK_variant_counts_length", pdf, png))
    print(f"wrote {png}")
    return out


def render() -> tuple[Path, Path]:
    """Render butterfly length panel, remaining panel alternatives, V1/V2 plates."""
    rng = np.random.default_rng(11)

    pdf_len, png_len = _render_length_butterfly_panel(rng)
    print(f"wrote {png_len}")

    _export_individual_panels(rng)

    pdf1, png1 = _render_version1_pipeline()
    print(f"wrote {png1}")
    pdf2, png2 = _render_version2_evidence()
    print(f"wrote {png2}")

    sheet = _render_comparison_sheet([png1, png2])
    print(f"wrote {sheet}")

    req = _write_requirements()
    print(f"wrote {req}")

    # Primary return: the iterated length-spectrum butterfly.
    return pdf_len, png_len
