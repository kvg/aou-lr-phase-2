"""Figure 1: haplotype-resolved cohort, local ancestry, short-read vs HiFi tracts."""

from __future__ import annotations

import csv
import gzip
from pathlib import Path

import numpy as np
from matplotlib.colors import to_rgb
from matplotlib.gridspec import GridSpec
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
from matplotlib.transforms import blended_transform_factory

from .figS16_zip_map import GEO_DIR, _load_people, _rings
from .style import (
    ANC_COLORS,
    ANC_ORDER,
    ANC_TEXT,
    COL_IN,
    FIG_W,
    MIN_FONT,
    HG38_MB,
    OKABE,
    P2_N,
    ROOT,
    ax_panel,
    chrom_offsets,
    covariates_path,
    despine,
    manhattan_x,
    new_figure,
    save,
)

# FLARE 5-way local ancestry (no MID in the Phase 2 FLARE model).
LAI = ["AFR", "AMR", "EUR", "EAS", "SAS"]
LAI_RGB = {a: np.array(to_rgb(ANC_COLORS[a])) for a in LAI}

# Typical five-way mix by global label (placeholder until FLARE global.anc is wired).
_MIX = {
    "AFR": np.array([0.84, 0.07, 0.06, 0.01, 0.02]),
    "AMR": np.array([0.12, 0.48, 0.32, 0.04, 0.04]),
    "EAS": np.array([0.02, 0.03, 0.04, 0.88, 0.03]),
    "EUR": np.array([0.03, 0.05, 0.88, 0.02, 0.02]),
    "MID": np.array([0.08, 0.10, 0.62, 0.04, 0.16]),
    "SAS": np.array([0.05, 0.08, 0.12, 0.06, 0.69]),
    "OTH": np.array([0.22, 0.18, 0.28, 0.14, 0.18]),
}


def _flare_path() -> Path | None:
    for path in (
        ROOT / "aou_lr_phase2_v1.chr1.global.anc.gz",
        ROOT / "data" / "aou_lr_phase2_v1.chr1.global.anc.gz",
    ):
        if path.is_file():
            return path
    return None


def _load_flare() -> dict[str, np.ndarray]:
    """chr1 FLARE global.anc: five-way mean local ancestry (no MID)."""
    path = _flare_path()
    if path is None:
        return {}
    out: dict[str, np.ndarray] = {}
    with gzip.open(path, "rt") as fh:
        header = next(fh).rstrip("\n").split("\t")
        idx = {h.lower(): i for i, h in enumerate(header)}
        cols = ["afr", "amr", "eur", "eas", "sas"]  # LAI order
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            sid = parts[idx["sample"]]
            vec = np.array([float(parts[idx[c]]) for c in cols], dtype=float)
            s = float(vec.sum())
            if s > 0:
                vec /= s
            out[sid] = vec
    return out


def _truthy(v) -> bool:
    return str(v).strip().lower() in {"true", "1", "yes"}


def _float_or_nan(v) -> float:
    try:
        if v in (None, ""):
            return np.nan
        return float(v)
    except (TypeError, ValueError):
        return np.nan


def _load_cohort() -> dict | None:
    path = covariates_path()
    if path is None:
        return None
    rows = []
    with gzip.open(path, "rt") as fh:
        for row in csv.DictReader(fh):
            if not _truthy(row.get("final_releasable_v9", "")):
                continue
            if row.get("technology") != "PacBio":
                continue
            rid = str(row.get("research_id") or "").strip()
            anc = (row.get("ancestry_pred_other") or "oth").upper()
            if anc not in ANC_ORDER:
                anc = "OTH"
            try:
                cov = float(row["coverage"]) if row.get("coverage") not in ("", None) else np.nan
            except (TypeError, ValueError):
                cov = np.nan
            try:
                ont = float(row["ont_coverage"]) if row.get("ont_coverage") not in ("", None) else np.nan
            except (TypeError, ValueError):
                ont = np.nan
            rows.append(
                (
                    rid,
                    anc,
                    cov,
                    ont,
                    _truthy(row.get("has_ehr_data")),
                    _truthy(row.get("has_ehr_data")) and _truthy(row.get("in_phenotype_table")),
                    _truthy(row.get("has_rna")),
                    _truthy(row.get("has_proteomics")),
                    _truthy(row.get("has_methylation")),
                    _truthy(row.get("is_AIAN")),
                    _truthy(row.get("in_pedigree")),
                    _truthy(row.get("has_ONT")),
                    str(row.get("zip3") or "").strip().split(".")[0],
                    _float_or_nan(row.get("lr_PC1")),
                    _float_or_nan(row.get("lr_PC2")),
                )
            )
    if not rows:
        return None
    return {
        "rid": np.array([r[0] for r in rows]),
        "anc": np.array([r[1] for r in rows]),
        "coverage": np.array([r[2] for r in rows], dtype=float),
        "ont": np.array([r[3] for r in rows], dtype=float),
        "ehr": np.array([r[4] for r in rows], dtype=bool),
        "pheno": np.array([r[5] for r in rows], dtype=bool),
        "rna": np.array([r[6] for r in rows], dtype=bool),
        "olink": np.array([r[7] for r in rows], dtype=bool),
        "meth": np.array([r[8] for r in rows], dtype=bool),
        "aian": np.array([r[9] for r in rows], dtype=bool),
        "ped": np.array([r[10] for r in rows], dtype=bool),
        "has_ont": np.array([r[11] for r in rows], dtype=bool),
        "zip3": np.array([r[12].zfill(3)[:3] if r[12] and r[12].lower() not in {"nan", "none", "na"} else "" for r in rows]),
        "pc1": np.array([r[13] for r in rows], dtype=float),
        "pc2": np.array([r[14] for r in rows], dtype=float),
        "real": True,
    }


def _synthetic_cohort(rng: np.random.Generator) -> dict:
    counts = [2816, 1966, 1309, 2177, 377, 1197, 2419]
    labels = np.concatenate([np.repeat(a, c) for a, c in zip(ANC_ORDER, counts)])
    n = labels.size
    high = rng.random(n) < (1133 / n)
    cov = np.empty(n)
    cov[~high] = rng.normal(16.3, 3.6, int((~high).sum()))
    cov[high] = rng.normal(33.2, 4.3, int(high.sum()))
    ont = np.full(n, np.nan)
    ont_ix = rng.choice(np.flatnonzero(high), size=min(876, int(high.sum())), replace=False)
    ont[ont_ix] = rng.normal(34.5, 5.1, ont_ix.size)
    ehr = rng.random(n) < 0.865
    return {
        "rid": np.array([f"sim{i}" for i in range(n)]),
        "anc": labels,
        "coverage": cov,
        "ont": ont,
        "ehr": ehr,
        "pheno": ehr & (rng.random(n) < 0.893),
        "rna": rng.random(n) < 0.679,
        "olink": rng.random(n) < 0.754,
        "meth": np.ones(n, dtype=bool),
        "aian": rng.random(n) < (800 / n),
        "ped": rng.random(n) < (134 / n),
        "has_ont": ~np.isnan(ont),
        "zip3": np.array([""] * n),
        "pc1": rng.normal(size=n),
        "pc2": rng.normal(size=n),
        "real": False,
    }


def _admixture(data: dict, rng: np.random.Generator) -> tuple[np.ndarray, bool]:
    """chr1 FLARE global.anc when the cohort is real; Dirichlet only for the mock cohort."""
    flare = _load_flare()
    n = data["anc"].size
    rids = data.get("rid")
    if data.get("real") and flare:
        missing = [str(rids[i]) for i in range(n) if str(rids[i]) not in flare]
        if missing:
            raise RuntimeError(
                f"{len(missing)} Phase 2 samples are missing from chr1 global.anc "
                f"(e.g. {missing[:3]})"
            )
        props = np.vstack([flare[str(rid)] for rid in rids])
        return props, True
    props = np.zeros((n, len(LAI)))
    for i, a in enumerate(data["anc"]):
        base = _MIX.get(a, _MIX["OTH"])
        props[i] = rng.dirichlet(np.clip(base * 28, 0.4, None))
    return props, False


def _sort_cohort(data: dict, props: np.ndarray) -> tuple[dict, np.ndarray]:
    """Ancestry groups, then mid-pass before high-pass, then decreasing African fraction."""
    afr = props[:, 0]
    order_key = np.array([ANC_ORDER.index(a) for a in data["anc"]])
    high = np.isfinite(data["coverage"]) & (data["coverage"] >= PB_HIGH_AT)
    ix = np.lexsort((-afr, high.astype(int), order_key))
    out = {k: (v[ix] if isinstance(v, np.ndarray) else v) for k, v in data.items()}
    return out, props[ix]


def _mosaic_image(props: np.ndarray, height: int = 64) -> np.ndarray:
    n, k = props.shape
    cum = np.cumsum(props, axis=1)
    fr = (np.arange(height)[:, None] + 0.5) / height
    colors = np.stack([LAI_RGB[a] for a in LAI])
    layer = np.zeros((height, n), dtype=int)
    for j in range(k):
        prev = 0.0 if j == 0 else cum[:, j - 1]
        layer[(fr >= prev) & (fr < cum[:, j])] = j
    return colors[layer]


def _paint_hap(length: float, mix: np.ndarray, mean_mb: float, flicker: float, rng: np.random.Generator):
    tracts = []
    x = 0.0
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


def _draw_hap(ax, y, h, tracts, length):
    for x0, x1, j in tracts:
        ax.add_patch(
            Rectangle((x0, y), x1 - x0, h, facecolor=ANC_COLORS[LAI[j]], edgecolor="none", zorder=2)
        )


# Empty columns between blocks, in participant units. The ancestry-group gap
# is wide enough to survive downsampling of ~12k columns onto a 183 mm plate.
# The tier gap is smaller so mid-pass and high-pass read as one group.
_GROUP_GAP = 220
_TIER_GAP = 36

# PacBio points are colored by coverage, not by the dual-technology cohort rule.
# 25× sits in the valley between the mid-pass and high-pass modes.
PB_HIGH_AT = 25.0
PB_MID_COLOR = "#F6D4E8"
PB_MID_LINE = "#C45A8A"
PB_HIGH_COLOR = "#DF1995"


def _tier_spans(anc: np.ndarray, high: np.ndarray) -> list[tuple[str, str, int, int]]:
    """Contiguous (ancestry, tier, start, end) blocks. Mid-pass precedes high-pass."""
    spans = []
    i = 0
    for a in ANC_ORDER:
        for tier, flag in (("mid", False), ("high", True)):
            c = int(((anc == a) & (high == flag)).sum())
            if not c:
                continue
            if not np.all(anc[i : i + c] == a) or not np.all(high[i : i + c] == flag):
                raise RuntimeError(f"{a} {tier} is not a contiguous block after sorting")
            spans.append((a, tier, i, i + c))
            i += c
    if i != anc.size:
        raise RuntimeError("ancestry groups are not contiguous after sorting")
    return spans


def _span_gaps(spans) -> list[int]:
    gaps = []
    for left, right in zip(spans, spans[1:]):
        gaps.append(_TIER_GAP if left[0] == right[0] else _GROUP_GAP)
    return gaps


def _gapped(img: np.ndarray, spans, gaps: list[int], fill: float) -> np.ndarray:
    pieces = []
    for k, (_a, _tier, i0, i1) in enumerate(spans):
        pieces.append(img[:, i0:i1])
        if k < len(spans) - 1:
            pieces.append(np.full((img.shape[0], gaps[k]) + img.shape[2:], fill, dtype=float))
    return np.concatenate(pieces, axis=1)


def _layout(spans, gaps: list[int]):
    """Person x positions, per-block ranges, and ancestry labels centered on each pair."""
    parts = []
    ranges = []
    cursor = 0.0
    for k, (a, tier, i0, i1) in enumerate(spans):
        c = i1 - i0
        parts.append(np.arange(c) + cursor + 0.5)
        ranges.append((a, tier, cursor, cursor + c, c))
        cursor += c + (gaps[k] if k < len(spans) - 1 else 0)
    labels = []
    i = 0
    while i < len(ranges):
        a = ranges[i][0]
        j = i + 1
        while j < len(ranges) and ranges[j][0] == a:
            j += 1
        x0 = ranges[i][2]
        x1 = ranges[j - 1][3]
        n = sum(ranges[k][4] for k in range(i, j))
        labels.append((a, (x0 + x1) / 2.0, n, x0, x1))
        i = j
    return np.concatenate(parts), ranges, labels, cursor


def _coverage_groups(data: dict) -> tuple[np.ndarray, np.ndarray, float]:
    """Masks for PacBio coverage below and at or above PB_HIGH_AT, plus mean ONT coverage."""
    cov = data["coverage"]
    ont = data["ont"]
    finite = np.isfinite(cov)
    high = finite & (cov >= PB_HIGH_AT)
    mid = finite & ~high
    return mid, high, float(ont[np.isfinite(ont)].mean())


# Assay combinations for the upset plot. Methylation is universal, so it is
# omitted here and stated in the caption. Combinations smaller than this are
# dropped; they are unlike one another and add up to a few percent of the cohort.
_UPSET_SETS = (
    ("EHR", "ehr"),
    ("Has phenotypes", "pheno"),
    ("Olink", "olink"),
    ("RNA-seq", "rna"),
)
_UPSET_MIN = 300


def _upset_intersections(data: dict) -> tuple[list[tuple[int, list[int]]], int]:
    bits = np.zeros(data["ehr"].size, dtype=np.int16)
    for i, (_, key) in enumerate(_UPSET_SETS):
        bits |= np.asarray(data[key], dtype=bool).astype(np.int16) << i
    counts: dict[int, int] = {}
    for bit in bits.tolist():
        counts[bit] = counts.get(bit, 0) + 1
    shown = []
    omitted = 0
    for bit, n in sorted(counts.items(), key=lambda item: -item[1]):
        if n < _UPSET_MIN:
            omitted += n
            continue
        members = [i for i in range(len(_UPSET_SETS)) if bit & (1 << i)]
        shown.append((n, members))
    return shown, omitted


def _panel_upset(
    fig,
    spec,
    data: dict,
    combos: list[list[int]],
    mid_pct: list[float],
    high_pct: list[float],
) -> None:
    """Each combination is two bars: percent of the mid-pass tier and of the high-pass tier."""
    inner = spec.subgridspec(2, 1, height_ratios=[1.15, 1.35], hspace=0.06)
    ax_b = fig.add_subplot(inner[0, 0])
    ax_d = fig.add_subplot(inner[1, 0], sharex=ax_b)
    xs = np.arange(len(combos))
    ink = "#3E4C59"
    w = 0.36
    ax_b.bar(
        xs - w / 2, mid_pct, width=w, color=PB_MID_COLOR, edgecolor=PB_MID_LINE,
        linewidth=0.6, zorder=2,
    )
    ax_b.bar(xs + w / 2, high_pct, width=w, color=PB_HIGH_COLOR, zorder=2)
    ax_b.set_xlim(-0.6, max(len(combos) - 0.4, 0.6))
    ymax = max(max(mid_pct, default=0), max(high_pct, default=0)) * 1.28
    ax_b.set_ylim(0, ymax)
    ax_b.set_ylabel("Percent")
    ax_b.tick_params(axis="x", bottom=False, labelbottom=False)
    for x, pct in zip(xs - w / 2, mid_pct):
        if pct >= 8:
            ax_b.text(x, pct, f"{pct:.0f}", ha="center", va="bottom", fontsize=MIN_FONT, color="black")
    for x, pct in zip(xs + w / 2, high_pct):
        if pct >= 8:
            ax_b.text(x, pct, f"{pct:.0f}", ha="center", va="bottom", fontsize=MIN_FONT, color="black")
    ax_b.legend(
        handles=[
            Line2D([0], [0], marker="s", color="none", markerfacecolor=PB_MID_COLOR, markeredgecolor=PB_MID_LINE, markersize=6, label="Mid-pass"),
            Line2D([0], [0], marker="s", color="none", markerfacecolor=PB_HIGH_COLOR, markeredgecolor=PB_HIGH_COLOR, markersize=6, label="High-pass"),
        ],
        loc="upper right",
        frameon=False,
        fontsize=MIN_FONT,
        handletextpad=0.3,
        borderaxespad=0.1,
    )
    despine(ax_b)
    n_sets = len(_UPSET_SETS)
    for x, members in enumerate(combos):
        if len(members) >= 2:
            ax_d.plot([x, x], [min(members), max(members)], color=ink, lw=0.9, zorder=1, solid_capstyle="round")
        if members:
            ax_d.scatter(np.full(len(members), x), members, s=16, c=ink, linewidths=0, zorder=2)
    ax_d.set_ylim(n_sets - 0.5, -0.5)
    ax_d.set_yticks(range(n_sets))
    ax_d.set_yticklabels(
        [f"{name}  {100.0 * np.mean(data[key]):.0f}%" for name, key in _UPSET_SETS]
    )
    ax_d.tick_params(axis="y", length=0, pad=2)
    ax_d.tick_params(axis="x", bottom=False, labelbottom=False)
    for sp in ax_d.spines.values():
        sp.set_visible(False)


def _subset(data: dict, props: np.ndarray, mask: np.ndarray) -> tuple[dict, np.ndarray]:
    n = mask.shape[0]
    out = {}
    for key, value in data.items():
        if isinstance(value, np.ndarray) and value.shape[:1] == (n,):
            out[key] = value[mask]
        else:
            out[key] = value
    return out, props[mask]


def _exact_percents(data: dict, combos: list[list[int]]) -> list[float]:
    n = data["ehr"].size
    masks = [np.asarray(data[key], dtype=bool) for _, key in _UPSET_SETS]
    out = []
    for members in combos:
        keep = set(members)
        match = np.ones(n, dtype=bool)
        for i, mask in enumerate(masks):
            match &= mask if i in keep else ~mask
        out.append(100.0 * float(match.sum()) / n if n else 0.0)
    return out


def _panel_cohort(fig, spec, data: dict, props: np.ndarray) -> None:
    """One strip: within each ancestry group, mid-pass sits immediately left of high-pass."""
    rows = (
        "mosaic", "assaygap",
        "olink", "trackgap1", "rna", "trackgap2", "pheno",
        "covgap", "pb", "gap", "ont",
    )
    height_of = {
        "mosaic": 0.640625,
        "assaygap": 0.12,
        "olink": 0.11,
        "trackgap1": 0.04,
        "rna": 0.11,
        "trackgap2": 0.04,
        "pheno": 0.11,
        "covgap": 0.14,
        "pb": 0.38,
        "gap": 0.08,
        "ont": 0.34,
    }
    inner = spec.subgridspec(len(rows), 1, height_ratios=[height_of[k] for k in rows], hspace=0.04)
    n = int(data["anc"].size)
    high = np.isfinite(data["coverage"]) & (data["coverage"] >= PB_HIGH_AT)
    spans = _tier_spans(data["anc"], high)
    gaps = _span_gaps(spans)
    x, _ranges, label_pos, width = _layout(spans, gaps)
    ax_m = fig.add_subplot(inner[0, 0])
    img = _gapped(_mosaic_image(props), spans, gaps, 1.0)
    ax_m.imshow(img, aspect="auto", interpolation="nearest", origin="upper", extent=(0, width, 0, 1))
    ax_m.set_yticks([])
    ax_m.set_xticks([])
    ax_m.set_ylabel(
        "Global ancestry", rotation=0, ha="right", va="center",
        fontsize=MIN_FONT, labelpad=4, color="black",
    )
    for a, xc, c, _x0, _x1 in label_pos:
        if c / n < 0.025:
            continue
        ax_m.annotate(
            f"{c:,}", xy=(xc, 1.02), xytext=(0, 1), textcoords="offset points",
            ha="center", va="bottom", fontsize=MIN_FONT, color="#444444",
            annotation_clip=False,
        )
        ax_m.annotate(
            a, xy=(xc, 1.02), xytext=(0, 12), textcoords="offset points",
            ha="center", va="bottom", fontsize=MIN_FONT, color=ANC_TEXT[a],
            fontweight="bold", annotation_clip=False,
        )
    ax_m.set_ylim(0, 1.02)
    ax_m.set_xlim(0, width)
    for sp in ax_m.spines.values():
        sp.set_visible(False)
    # Right-hand column. Same axes-fraction x on every row, right-aligned.
    right_x = 1.09
    ax_m.annotate(
        "Total", xy=(right_x, 1.02), xycoords=blended_transform_factory(ax_m.transAxes, ax_m.transData),
        xytext=(0, 12), textcoords="offset points",
        ha="right", va="bottom", fontsize=MIN_FONT, color="black",
        fontweight="bold", annotation_clip=False,
    )
    ax_m.annotate(
        f"{n:,}", xy=(right_x, 1.02), xycoords=blended_transform_factory(ax_m.transAxes, ax_m.transData),
        xytext=(0, 1), textcoords="offset points",
        ha="right", va="bottom", fontsize=MIN_FONT, color="#444444",
        annotation_clip=False,
    )

    ont_color = "#007FA8"
    point_s = 1.0
    cov = data["coverage"]
    finite = np.isfinite(cov)
    ax_pb = fig.add_subplot(inner[rows.index("pb"), 0], sharex=ax_m)
    colors = np.where(high[finite], PB_HIGH_COLOR, PB_MID_COLOR)
    ax_pb.scatter(x[finite], cov[finite], s=point_s, c=colors, alpha=0.9, linewidths=0, rasterized=True, zorder=2)
    ax_pb.set_ylim(0, 60)
    ax_pb.set_yticks([0, 20, 40, 60])
    ax_pb.set_ylabel(
        "PacBio ×", rotation=0, ha="right", va="center",
        fontsize=MIN_FONT, labelpad=10, color="black",
    )
    ax_pb.set_xticks([])
    despine(ax_pb)
    ax_pb.legend(
        handles=[
            Line2D(
                [0], [0], marker="o", ls="none", markerfacecolor=PB_MID_COLOR,
                markeredgecolor=PB_MID_LINE, markeredgewidth=0.7, markersize=5.5, label="Mid-pass",
            ),
            Line2D(
                [0], [0], marker="o", ls="none", markerfacecolor=PB_HIGH_COLOR,
                markeredgecolor=PB_HIGH_COLOR, markersize=5.5, label="High-pass",
            ),
        ],
        loc="upper left",
        ncol=2,
        frameon=False,
        fontsize=MIN_FONT,
        borderaxespad=0.15,
        handletextpad=0.25,
        columnspacing=0.8,
    )
    for mask, color in ((~high & finite, PB_MID_LINE), (high, PB_HIGH_COLOR)):
        if not mask.any():
            continue
        y = round(float(cov[mask].mean()), 1)
        ax_pb.axhline(y, color=color, ls=":", lw=0.9, zorder=3)
        ax_pb.annotate(
            f"{y:.1f}×",
            xy=(right_x, y),
            xycoords=blended_transform_factory(ax_pb.transAxes, ax_pb.transData),
            ha="right",
            va="center",
            fontsize=MIN_FONT,
            color="black",
            annotation_clip=False,
        )

    ax_ont = fig.add_subplot(inner[rows.index("ont"), 0], sharex=ax_m)
    ont_ok = np.isfinite(data["ont"])
    ax_ont.scatter(
        x[ont_ok], data["ont"][ont_ok], s=point_s, c=ont_color, alpha=0.9,
        linewidths=0, rasterized=True,
    )
    ax_ont.set_ylim(0, 60)
    ax_ont.set_yticks([0, 20, 40, 60])
    ax_ont.set_ylabel(
        "ONT ×", rotation=0, ha="right", va="center",
        fontsize=MIN_FONT, labelpad=10, color="black",
    )
    if ont_ok.any():
        ont_y = round(float(data["ont"][ont_ok].mean()), 1)
        ax_ont.axhline(ont_y, color=ont_color, ls=":", lw=0.9, zorder=3)
        ax_ont.annotate(
            f"{ont_y:.1f}×",
            xy=(right_x, ont_y),
            xycoords=blended_transform_factory(ax_ont.transAxes, ax_ont.transData),
            ha="right",
            va="center",
            fontsize=MIN_FONT,
            color="black",
            annotation_clip=False,
        )
    despine(ax_ont)
    ax_ont.set_xlabel("Participants")

    present = np.array(to_rgb("#3E4C59"))
    absent = np.array(to_rgb("#E4E8EC"))
    for key, name in (("olink", "Olink"), ("rna", "RNA-seq"), ("pheno", "Phenotypes")):
        flags = np.asarray(data[key], dtype=bool)
        ax = fig.add_subplot(inner[rows.index(key), 0], sharex=ax_m)
        rgb = np.where(flags[:, None], present, absent)
        ax.imshow(
            _gapped(rgb[None, :, :], spans, gaps, 1.0),
            aspect="auto", interpolation="nearest", origin="upper",
            extent=(0, width, 0, 1),
        )
        ax.set_yticks([])
        ax.set_xticks([])
        ax.set_ylabel(
            name, rotation=0, ha="right", va="center",
            fontsize=MIN_FONT, labelpad=10, color="black",
        )
        ax.annotate(
            f"{100.0 * float(flags.mean()):.1f}%",
            xy=(right_x, 0.5),
            xycoords=blended_transform_factory(ax.transAxes, ax.transData),
            ha="right",
            va="center",
            fontsize=MIN_FONT,
            color="black",
            annotation_clip=False,
        )
        for sp in ax.spines.values():
            sp.set_visible(False)


def _panel_karyograms(ax, rng: np.random.Generator) -> None:
    people = [
        ("AMR, highly admixed", _MIX["AMR"]),
        ("AFR–EUR, two-way", np.array([0.55, 0.05, 0.38, 0.01, 0.01])),
        ("SAS, residual EUR", _MIX["SAS"]),
    ]
    length = HG38_MB[1]
    ax.set_xlim(-12, length)
    ax.set_ylim(-0.45, 14.2)
    ax.axis("off")
    ax.set_title("Short-read vs HiFi local ancestry  ·  chr1 (placeholder tracts)", loc="left", pad=2)
    y = 13.3
    for name, mix in people:
        mix = mix / mix.sum()
        ax.add_patch(Rectangle((-11, y - 3.55), length + 11, 3.85, facecolor="#FAFAFA", edgecolor="none", zorder=0))
        ax.text(0, y + 0.12, name, fontsize=6.1, color="#222222", va="bottom")
        for mean_mb, flicker, lab in [(4.2, 0.34, "SR"), (28.0, 0.018, "LR")]:
            ax.text(-1.2, y - 0.85, lab, ha="right", va="center", fontsize=5.5, color="#444444")
            for _hap in range(2):
                y -= 0.68
                tracts = _paint_hap(length, mix, mean_mb, flicker, rng)
                _draw_hap(ax, y, 0.52, tracts, length)
            y -= 0.20
        y -= 0.42
    ax.plot([0, length], [0.22, 0.22], color="#222222", lw=0.5)
    for tick, lab in [(0, "0"), (50, "50"), (100, "100"), (150, "150"), (200, "200"), (248, "chr1 Mb")]:
        ax.plot([tick, tick], [0.12, 0.22], color="#222222", lw=0.5)
        ax.text(tick, -0.12, lab, ha="center", va="top", fontsize=5.4, color="#444444")


def _tract_lengths(rng: np.random.Generator, n, mean, flicker):
    w = rng.exponential(mean, n)
    flick = rng.random(n) < flicker
    w[flick] = rng.uniform(0.05, 0.4, flick.sum())
    return np.clip(w, 0.05, 80)


def _panel_tract_length(ax, rng: np.random.Generator) -> None:
    for arr, color, lab in [
        (_tract_lengths(rng, 4000, 4.2, 0.34), "#888888", "Short-read LAI"),
        (_tract_lengths(rng, 1800, 28.0, 0.018), OKABE["blue"], "HiFi FLARE"),
    ]:
        xs = np.sort(arr)
        ax.plot(xs, np.linspace(0, 1, xs.size), color=color, lw=1.3, label=lab)
    ax.set_xscale("log")
    ax.set_xlim(0.08, 90)
    ax.set_xlabel("Ancestry tract length (Mb)")
    ax.set_ylabel("ECDF")
    ax.set_title("Tract length", loc="left", pad=2)
    despine(ax)
    ax.legend(loc="lower right", fontsize=5.5)


def _panel_switch_rate(ax, rng: np.random.Generator) -> None:
    cats = ["Unique", "TR", "SD", "MHC", "Cen"]
    x = np.arange(len(cats))
    w = 0.36
    ax.bar(x - w / 2, [0.42, 0.91, 1.35, 1.82, 2.40], width=w, color="#C8C8C8", label="SR", zorder=2)
    ax.bar(x + w / 2, [0.11, 0.22, 0.38, 0.44, 0.70], width=w, color=OKABE["blue"], label="HiFi", zorder=2)
    ax.set_xticks(x)
    ax.set_xticklabels(cats)
    ax.set_ylabel("Switches / Mb")
    ax.set_title("Switch rate by context", loc="left", pad=2)
    despine(ax)
    ax.legend(loc="upper left", fontsize=5.5)


def _panel_genome(ax, rng: np.random.Generator) -> None:
    off = chrom_offsets()
    xs, stack = [], []
    for c, length in HG38_MB.items():
        pos = np.arange(0.4, length, 1.6)
        base = np.array([0.28, 0.16, 0.26, 0.14, 0.16])
        p = np.tile(base, (pos.size, 1)) + rng.normal(0, 0.012, (pos.size, 5))
        if c == 6:
            bump = np.exp(-0.5 * ((pos - 31.0) / 2.2) ** 2)
            p[:, 0] += 0.12 * bump
            p[:, 2] += 0.08 * bump
            p[:, 1] -= 0.10 * bump
        if c == 2:
            bump = np.exp(-0.5 * ((pos - 136.0) / 1.8) ** 2)
            p[:, 2] += 0.16 * bump
            p[:, 0] -= 0.08 * bump
        if c == 16:
            bump = np.exp(-0.5 * ((pos - 0.4) / 0.9) ** 2)
            p[:, 0] += 0.18 * bump
            p[:, 2] -= 0.08 * bump
        p = np.clip(p, 0.01, None)
        p = p / p.sum(axis=1, keepdims=True)
        xs.append(off[c] + pos)
        stack.append(p)
    x = np.concatenate(xs)
    p = np.vstack(stack)
    bottom = np.zeros(x.size)
    for j in [0, 2, 1, 3, 4]:
        ax.fill_between(x, bottom, bottom + p[:, j], color=ANC_COLORS[LAI[j]], linewidth=0, alpha=0.92, label=LAI[j])
        bottom += p[:, j]
    ax.set_ylim(0, 1)
    ax.set_xlim(0, x.max())
    ax.set_ylabel("Mean local ancestry")
    ax.set_xlabel("Chromosome")
    ax.set_title("Cohort-mean local ancestry along the genome (placeholder FLARE)", loc="left", pad=2)
    despine(ax)
    centers = {c: off[c] + HG38_MB[c] / 2 for c in HG38_MB}
    ax.set_xticks([centers[c] for c in range(1, 23)])
    ax.set_xticklabels([str(c) if c not in (11, 13, 15, 17, 19, 21) else "" for c in range(1, 23)], fontsize=5.5)
    for chrom, pos_mb, text in [(6, 32.0, "MHC"), (2, 136.0, "LCT"), (16, 0.4, "HBA")]:
        xx = float(manhattan_x(np.array([chrom]), np.array([pos_mb]))[0])
        ax.annotate(text, xy=(xx, 1.0), xytext=(0, 3), textcoords="offset points", fontsize=5.6, ha="center", va="bottom")
    ax.legend(loc="upper right", ncol=5, bbox_to_anchor=(1.0, 1.18), fontsize=5.6)


def _sort_ancestry(data: dict, props: np.ndarray) -> tuple[dict, np.ndarray]:
    """Ancestry groups, then decreasing African fraction within each group."""
    afr = props[:, 0]
    order_key = np.array([ANC_ORDER.index(a) for a in data["anc"]])
    ix = np.lexsort((-afr, order_key))
    out = {k: (v[ix] if isinstance(v, np.ndarray) else v) for k, v in data.items()}
    return out, props[ix]


def _ancestry_spans(anc: np.ndarray) -> list[tuple[str, str, int, int]]:
    spans = []
    i = 0
    for a in ANC_ORDER:
        c = int((anc == a).sum())
        if not c:
            continue
        spans.append((a, "all", i, i + c))
        i += c
    if i != anc.size:
        raise RuntimeError("ancestry groups are not contiguous after sorting")
    return spans


def _panel_admixture(ax, data: dict, props: np.ndarray) -> None:
    """Global-ancestry stacks, one block per computed ancestry group."""
    n = int(data["anc"].size)
    spans = _ancestry_spans(data["anc"])
    gaps = _span_gaps(spans)
    _x, _ranges, label_pos, width = _layout(spans, gaps)
    img = _gapped(_mosaic_image(props), spans, gaps, 1.0)
    ax.imshow(img, aspect="auto", interpolation="nearest", origin="upper", extent=(0, width, 0, 1))
    ax.set_yticks([])
    ax.set_xticks([])
    ax.set_ylabel(
        "Global ancestry", rotation=0, ha="right", va="center",
        fontsize=MIN_FONT, labelpad=6, color="black",
    )
    for a, xc, c, _x0, _x1 in label_pos:
        if c / n < 0.025:
            continue
        ax.annotate(
            f"{c:,}", xy=(xc, 1.02), xytext=(0, 1), textcoords="offset points",
            ha="center", va="bottom", fontsize=MIN_FONT, color="#444444",
            annotation_clip=False,
        )
        ax.annotate(
            a, xy=(xc, 1.02), xytext=(0, 12), textcoords="offset points",
            ha="center", va="bottom", fontsize=MIN_FONT, color=ANC_TEXT[a],
            fontweight="bold", annotation_clip=False,
        )
    right_x = 1.09
    ax.annotate(
        "Total", xy=(right_x, 1.02), xycoords=blended_transform_factory(ax.transAxes, ax.transData),
        xytext=(0, 12), textcoords="offset points",
        ha="right", va="bottom", fontsize=MIN_FONT, color="black",
        fontweight="bold", annotation_clip=False,
    )
    ax.annotate(
        f"{n:,}", xy=(right_x, 1.02), xycoords=blended_transform_factory(ax.transAxes, ax.transData),
        xytext=(0, 1), textcoords="offset points",
        ha="right", va="bottom", fontsize=MIN_FONT, color="#444444",
        annotation_clip=False,
    )
    ax.set_ylim(0, 1.02)
    ax.set_xlim(0, width)
    for sp in ax.spines.values():
        sp.set_visible(False)


# State of the 3-digit ZIP centroid. A tile is colored only when the most
# common ancestry in that state itself has at least this many participants.
_MIN_STATE_N = 20
_STATE_POSTAL = {
    "Alabama": "AL", "Alaska": "AK", "Arizona": "AZ", "Arkansas": "AR",
    "California": "CA", "Colorado": "CO", "Connecticut": "CT", "Delaware": "DE",
    "District of Columbia": "DC", "Florida": "FL", "Georgia": "GA", "Hawaii": "HI",
    "Idaho": "ID", "Illinois": "IL", "Indiana": "IN", "Iowa": "IA",
    "Kansas": "KS", "Kentucky": "KY", "Louisiana": "LA", "Maine": "ME",
    "Maryland": "MD", "Massachusetts": "MA", "Michigan": "MI", "Minnesota": "MN",
    "Mississippi": "MS", "Missouri": "MO", "Montana": "MT", "Nebraska": "NE",
    "Nevada": "NV", "New Hampshire": "NH", "New Jersey": "NJ", "New Mexico": "NM",
    "New York": "NY", "North Carolina": "NC", "North Dakota": "ND", "Ohio": "OH",
    "Oklahoma": "OK", "Oregon": "OR", "Pennsylvania": "PA", "Puerto Rico": "PR", "Rhode Island": "RI",
    "South Carolina": "SC", "South Dakota": "SD", "Tennessee": "TN", "Texas": "TX",
    "Utah": "UT", "Vermont": "VT", "Virginia": "VA", "Washington": "WA",
    "West Virginia": "WV", "Wisconsin": "WI", "Wyoming": "WY",
}


_CENSUS_DIVISION = {
    "CT": "New England", "ME": "New England", "MA": "New England",
    "NH": "New England", "RI": "New England", "VT": "New England",
    "NJ": "Mid-Atlantic", "NY": "Mid-Atlantic", "PA": "Mid-Atlantic",
    "IL": "E. North Central", "IN": "E. North Central", "MI": "E. North Central",
    "OH": "E. North Central", "WI": "E. North Central",
    "IA": "W. North Central", "KS": "W. North Central", "MN": "W. North Central",
    "MO": "W. North Central", "NE": "W. North Central", "ND": "W. North Central",
    "SD": "W. North Central",
    "DE": "South Atlantic", "DC": "South Atlantic", "FL": "South Atlantic",
    "GA": "South Atlantic", "MD": "South Atlantic", "NC": "South Atlantic",
    "SC": "South Atlantic", "VA": "South Atlantic", "WV": "South Atlantic",
    "AL": "E. South Central", "KY": "E. South Central", "MS": "E. South Central",
    "TN": "E. South Central",
    "AR": "W. South Central", "LA": "W. South Central", "OK": "W. South Central",
    "TX": "W. South Central",
    "AZ": "Mountain", "CO": "Mountain", "ID": "Mountain", "MT": "Mountain",
    "NV": "Mountain", "NM": "Mountain", "UT": "Mountain", "WY": "Mountain",
    "AK": "Pacific", "CA": "Pacific", "HI": "Pacific", "OR": "Pacific", "WA": "Pacific",
}
_DIVISION_REGION = {
    "New England": "Northeast", "Mid-Atlantic": "Northeast",
    "E. North Central": "Midwest", "W. North Central": "Midwest",
    "South Atlantic": "South", "E. South Central": "South", "W. South Central": "South",
    "Mountain": "West", "Pacific": "West",
}
_REGION_ORDER = ("Northeast", "Midwest", "South", "West")
_DIVISION_ORDER = (
    "New England", "Mid-Atlantic",
    "E. North Central", "W. North Central",
    "South Atlantic", "E. South Central", "W. South Central",
    "Mountain", "Pacific",
)


def _label_ink(color: str) -> str:
    r, g, b = to_rgb(color)
    return "#222222" if 0.2126 * r + 0.7152 * g + 0.0722 * b > 0.62 else "white"


def _assign_states(lon: np.ndarray, lat: np.ndarray) -> np.ndarray:
    """State name for each ZIP3 centroid. None when the centroid is not in a state."""
    import json
    from matplotlib.path import Path as MPath

    features = json.loads((GEO_DIR / "us-states.json").read_text())["features"]
    states = []
    for feat in features:
        name = feat["properties"]["name"]
        if name not in _STATE_POSTAL:
            continue
        paths, cents = [], []
        for ring in _rings(feat["geometry"]):
            arr = np.asarray(ring, dtype=float)
            paths.append(MPath(arr))
            cents.append(arr.mean(axis=0))
        states.append((name, paths, np.mean(cents, axis=0)))

    def one(x: float, y: float) -> str | None:
        pt = (x, y)
        for name, paths, _c in states:
            if any(path.contains_point(pt) for path in paths):
                return name
        best, dist = None, 1e9
        for name, _paths, center in states:
            d = (center[0] - x) ** 2 + (center[1] - y) ** 2
            if d < dist:
                best, dist = name, d
        return best if dist < 4.0 else None

    return np.array([one(float(x), float(y)) for x, y in zip(lon, lat)], dtype=object)


def _panel_state_hex(ax, postal_counts: dict[str, dict[str, int]]) -> None:
    """US state hex cartogram. Color is the plurality ancestry, or blank."""
    import json

    features = json.loads((GEO_DIR / "us_states_hexgrid.geojson.json").read_text())["features"]
    blank = "#E6E6E6"
    neutral = "#C5C9CE"
    for feat in features:
        postal = feat["properties"]["iso3166_2"]
        ring = np.asarray(next(_rings(feat["geometry"])), dtype=float)
        counts = postal_counts.get(postal, {})
        n = sum(counts.values())
        top = max(counts, key=counts.get) if counts else None
        top_n = counts.get(top, 0) if top else 0
        if top is not None and top_n >= _MIN_STATE_N:
            face = ANC_COLORS[top]
            ink = _label_ink(face)
            label = f"{postal}\n{n:,}"
        elif n >= _MIN_STATE_N:
            face, ink, label = neutral, "#333333", f"{postal}\n{n:,}"
        else:
            face, ink, label = blank, "#A0A0A0", postal
        ax.fill(ring[:, 0], ring[:, 1], facecolor=face, edgecolor="white", lw=1.1, zorder=1, closed=True)
        ax.text(
            ring[:, 0].mean(), ring[:, 1].mean(), label,
            ha="center", va="center", fontsize=MIN_FONT, color=ink, zorder=2, linespacing=0.95,
        )
    ax.set_aspect("equal")
    ax.axis("off")
    ax.text(
        0.5, -0.02, "State of 3-digit ZIP code  ·  blank tiles have fewer than 20 participants",
        transform=ax.transAxes, ha="center", va="top",
        fontsize=MIN_FONT, color="black", clip_on=False,
    )


def _hex_aspect() -> float:
    import json

    features = json.loads((GEO_DIR / "us_states_hexgrid.geojson.json").read_text())["features"]
    xs, ys = [], []
    for feat in features:
        ring = np.asarray(next(_rings(feat["geometry"])), dtype=float)
        xs.append(ring[:, 0])
        ys.append(ring[:, 1])
    x = np.concatenate(xs)
    y = np.concatenate(ys)
    return float((x.max() - x.min()) / (y.max() - y.min()))


def _panel_mix(ax, rows: list[tuple[str, dict[str, int]]], *, show_legend: bool, show_xlabel: bool) -> None:
    """Horizontal stacked bars of computed-ancestry percent. Count sits at the right."""
    y = np.arange(len(rows))[::-1]
    left = np.zeros(len(rows))
    for anc in ANC_ORDER:
        vals = []
        for _name, counts in rows:
            total = sum(counts.values())
            vals.append(100.0 * counts.get(anc, 0) / total if total else 0.0)
        ax.barh(y, vals, left=left, height=0.72, color=ANC_COLORS[anc], linewidth=0, label=anc)
        left = left + np.asarray(vals)
    ax.set_yticks(y)
    ax.set_yticklabels([name for name, _counts in rows])
    ax.set_xlim(0, 100)
    ax.set_xticks([0, 25, 50, 75, 100])
    if show_xlabel:
        ax.set_xlabel("Percent of participants with a 3-digit ZIP code")
    else:
        ax.tick_params(labelbottom=False)
    for yi, (_name, counts) in zip(y, rows):
        ax.annotate(
            f"{sum(counts.values()):,}",
            xy=(100, yi), xytext=(4, 0), textcoords="offset points",
            ha="left", va="center", fontsize=MIN_FONT, color="#444444",
            annotation_clip=False,
        )
    if show_legend:
        ax.legend(
            loc="lower center", bbox_to_anchor=(0.5, 1.16), ncol=7, frameon=False,
            fontsize=MIN_FONT, handlelength=0.8, handletextpad=0.3, columnspacing=0.8,
            borderaxespad=0,
        )
    despine(ax)


def _zip3_postal() -> dict[str, str]:
    """Map a 3-digit ZIP code to the state that contains its centroid."""
    import csv

    cents_lon, cents_lat, keys = [], [], []
    with (GEO_DIR / "zip3_centroids.tsv").open() as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            keys.append(row["zip3"].zfill(3)[:3])
            cents_lat.append(float(row["lat"]))
            cents_lon.append(float(row["lon"]))
    names = _assign_states(np.asarray(cents_lon), np.asarray(cents_lat))
    out = {}
    for key, name in zip(keys, names):
        postal = _STATE_POSTAL.get(name)
        if postal:
            out[key] = postal
    return out


# Square cartogram. Row 0 is the top. Matches the usual state-tile arrangement:
# Alaska and Maine on the top corners, Hawaii and Puerto Rico on the bottom.
_STATE_TILES = {
    "AK": (0, 0), "ME": (10, 0),
    "WI": (5, 1), "VT": (9, 1), "NH": (10, 1),
    "WA": (0, 2), "ID": (1, 2), "MT": (2, 2), "ND": (3, 2), "MN": (4, 2),
    "IL": (5, 2), "MI": (6, 2), "NY": (8, 2), "MA": (9, 2),
    "OR": (0, 3), "NV": (1, 3), "WY": (2, 3), "SD": (3, 3), "IA": (4, 3),
    "IN": (5, 3), "OH": (6, 3), "PA": (7, 3), "NJ": (8, 3), "CT": (9, 3), "RI": (10, 3),
    "CA": (0, 4), "UT": (1, 4), "CO": (2, 4), "NE": (3, 4), "MO": (4, 4),
    "WV": (5, 4), "VA": (6, 4), "MD": (7, 4), "DC": (8, 4), "DE": (9, 4),
    "AZ": (1, 5), "NM": (2, 5), "KS": (3, 5), "AR": (4, 5), "KY": (5, 5),
    "TN": (6, 5), "NC": (7, 5), "SC": (8, 5),
    "OK": (3, 6), "LA": (4, 6), "MS": (5, 6), "AL": (6, 6), "GA": (7, 6),
    "HI": (0, 7), "TX": (3, 7), "FL": (8, 7), "PR": (10, 7),
}
_TILE_NCOL = 11
_TILE_NROW = 8
_TILE_MIN_N = 20


def _panel_state_tiles(ax, postal: np.ndarray, props: np.ndarray) -> None:
    """One ancestry stack per state. Tiles below the count threshold stay gray."""
    import matplotlib.patheffects as pe

    gap = 0.14
    stroke = [pe.withStroke(linewidth=0.8, foreground="white")]
    for code, (col, row) in _STATE_TILES.items():
        x0 = col * (1 + gap)
        y0 = (_TILE_NROW - 1 - row) * (1 + gap)
        mask = postal == code
        n = int(mask.sum())
        if n < _TILE_MIN_N:
            ax.add_patch(Rectangle((x0, y0), 1, 1, facecolor="#D5D5D5", edgecolor="white", lw=0.6, zorder=1))
        else:
            sub = props[mask]
            sub = sub[np.argsort(-sub[:, 0])]
            img = _mosaic_image(sub, height=40)
            ax.imshow(
                img, extent=(x0, x0 + 1, y0, y0 + 1), origin="upper",
                interpolation="nearest", aspect="auto", zorder=1,
            )
        ax.text(
            x0 + 0.06, y0 + 0.94, code,
            ha="left", va="top", fontsize=MIN_FONT, color="#222222", zorder=3,
            path_effects=stroke,
        )
    width = _TILE_NCOL + (_TILE_NCOL - 1) * gap
    height = _TILE_NROW + (_TILE_NROW - 1) * gap
    ax.set_xlim(-0.08, width - 1 + 1.08)
    ax.set_ylim(-0.08, height - 1 + 1.08)
    ax.set_aspect("equal")
    ax.axis("off")
    ax.text(
        0.5, -0.02,
        "State of 3-digit ZIP code  ·  gray tiles have fewer than 20 participants",
        transform=ax.transAxes, ha="center", va="top",
        fontsize=MIN_FONT, color="black", clip_on=False,
    )


def _tile_aspect() -> float:
    gap = 0.14
    width = _TILE_NCOL + (_TILE_NCOL - 1) * gap
    height = _TILE_NROW + (_TILE_NROW - 1) * gap
    return width / height


def _panel_pca(ax, data: dict) -> None:
    """Long-read global PC1 vs PC2. Colors match the ancestry strip."""
    from matplotlib.ticker import MultipleLocator

    x = np.asarray(data["pc1"], dtype=float)
    y = np.asarray(data["pc2"], dtype=float)
    anc = data["anc"]
    # Larger groups first so smaller clusters stay visible.
    order = sorted(ANC_ORDER, key=lambda a: -int((anc == a).sum()))
    for name in order:
        mask = (anc == name) & np.isfinite(x) & np.isfinite(y)
        ax.scatter(
            x[mask], y[mask],
            s=4.0, c=ANC_COLORS[name], alpha=0.55,
            linewidths=0, rasterized=True, zorder=2,
        )
    finite = np.isfinite(x) & np.isfinite(y)
    span = max(float(np.ptp(x[finite])), float(np.ptp(y[finite]))) * 1.06
    half = span / 2
    xc = float(np.min(x[finite]) + np.max(x[finite])) / 2
    yc = float(np.min(y[finite]) + np.max(y[finite])) / 2
    ax.set_xlim(xc - half, xc + half)
    ax.set_ylim(yc - half, yc + half)
    ax.set_aspect("equal", adjustable="box")
    ax.xaxis.set_major_locator(MultipleLocator(0.2))
    ax.yaxis.set_major_locator(MultipleLocator(0.2))
    ax.set_xlabel("PC1", fontsize=MIN_FONT)
    ax.set_ylabel("PC2", rotation=0, ha="right", va="center", fontsize=MIN_FONT, labelpad=8)
    ax.tick_params(length=2.5, width=0.6, labelsize=MIN_FONT)
    for sp in ax.spines.values():
        sp.set_linewidth(0.6)


def _prepared(rng: np.random.Generator) -> tuple[dict, np.ndarray, bool]:
    loaded = _load_cohort()
    data = loaded if loaded is not None else _synthetic_cohort(rng)
    props, flare_real = _admixture(data, rng)
    data["flare_real"] = flare_real
    data, props = _sort_ancestry(data, props)
    return data, props, flare_real


def render(*, placeholders: bool = True) -> tuple[Path, Path]:
    rng = np.random.default_rng(21)
    data, props, flare_real = _prepared(rng)

    if placeholders and not flare_real:
        ax_panel("fig1_karyograms", 3.15, lambda ax: _panel_karyograms(ax, rng), left=0.06, right=0.98, top=0.90, bottom=0.10)
        ax_panel("fig1_tract_length", 2.85, lambda ax: _panel_tract_length(ax, rng), width=COL_IN * 1.05, left=0.18, right=0.96, top=0.88, bottom=0.16)
        ax_panel("fig1_switch_rate", 2.85, lambda ax: _panel_switch_rate(ax, rng), width=COL_IN * 1.05, left=0.16, right=0.96, top=0.88, bottom=0.16)
        ax_panel("fig1_genome_lai", 2.35, lambda ax: _panel_genome(ax, rng), left=0.10, right=0.98, top=0.82, bottom=0.16)

    left = 1.18
    width = FIG_W - left - 0.62
    bar_h = 0.58
    top_room = 0.42
    bottom = 0.06
    fig_h = bottom + bar_h + top_room
    fig = new_figure(fig_h)
    fig_w, fig_h_in = fig.get_size_inches()
    ax_m = fig.add_axes([
        left / fig_w,
        bottom / fig_h_in,
        width / fig_w,
        bar_h / fig_h_in,
    ])
    _panel_admixture(ax_m, data, props)
    return save(fig, "fig1_cohort")


def render_pca() -> tuple[Path, Path]:
    """Long-read PCA, kept as a supplemental figure."""
    rng = np.random.default_rng(21)
    data, _props, _flare_real = _prepared(rng)
    side = 4.60
    left = (FIG_W - side) / 2 + 0.20
    bottom = 0.42
    top = 0.16
    fig_h = bottom + side + top
    fig = new_figure(fig_h)
    fig_w, fig_h_in = fig.get_size_inches()
    ax = fig.add_axes([
        left / fig_w,
        bottom / fig_h_in,
        side / fig_w,
        side / fig_h_in,
    ])
    _panel_pca(ax, data)
    return save(fig, "figS_pca")


def render_state_tiles() -> tuple[Path, Path]:
    """State-tile ancestry map, kept separate from Figure 1."""
    rng = np.random.default_rng(21)
    data, props, _flare_real = _prepared(rng)
    lookup = _zip3_postal()
    postal = np.array([lookup.get(z, "") for z in data["zip3"]], dtype=object)

    map_left = 0.15
    map_w = FIG_W - 0.30
    map_h = map_w / _tile_aspect()
    top_room = 0.32
    bottom = 0.28
    fig_h = bottom + map_h + top_room
    fig = new_figure(fig_h)
    fig_w, fig_h_in = fig.get_size_inches()
    ax = fig.add_axes([
        map_left / fig_w,
        bottom / fig_h_in,
        map_w / fig_w,
        map_h / fig_h_in,
    ])
    _panel_state_tiles(ax, postal, props)
    for i, name in enumerate(LAI):
        fig.text(
            0.50 + (i - 2) * 0.09,
            (bottom + map_h + 0.06) / fig_h,
            name,
            ha="center", va="bottom", fontsize=MIN_FONT,
            color=ANC_TEXT[name], fontweight="bold",
        )
    return save(fig, "fig_state_tiles")
