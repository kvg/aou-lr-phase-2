"""Figure 1: cohort and association-analysis design.

Age and sex, long-read principal components (with Phase 1 projected in),
PacBio coverage, and sequencing center, plus multi-omic data availability.
Local ancestry is reserved for the genotype figure.
"""

from __future__ import annotations

import csv
import gzip
import importlib.util
from collections import Counter

import numpy as np
from matplotlib import patheffects
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.legend_handler import HandlerTuple
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.ticker import StrMethodFormatter
from matplotlib.transforms import Bbox

from .style import (
    ANC_COLORS,
    ANC_ORDER,
    ANC_TEXT,
    MIN_FONT,
    ROOT,
    covariates_path,
    new_figure,
    save,
)

_CENTERS = (
    ("BI", "Broad"),
    ("HA", "HudsonAlpha"),
    ("BCM", "Baylor"),
    ("UW", "UW"),
    ("JHU", "Johns Hopkins"),
)
_MULTIOMIC_COLOR = "#788585"
_MULTIOMIC_PALE = "#C5CDD1"
_MID_PASS_COLOR = "#A9BBC9"
_HIGH_PASS_COLOR = "#334E68"
_ONT_COLOR = "#007FA8"
_ONT_ONLY_COLOR = "#BFE3EF"
_FEMALE = "#1F6F78"
_FEMALE_PALE = "#B9D7DA"
_MALE = "#C9A227"
_MALE_PALE = "#EEE0B0"
_PHASE1 = "#9C3D1A"
_HIGH_PASS_N = 1_133
# Shared metrics so the 2- and 3-row legends above panels C–F share baselines.
_ROW_LEGEND = dict(
    frameon=False,
    fontsize=MIN_FONT,
    handlelength=1.0,
    handleheight=0.7,
    handletextpad=0.3,
    borderaxespad=0,
    borderpad=0.4,
    labelspacing=0.5,
)
_ROW_LEGEND_Y = 1.01
_ILLUMINA_COUNTS = ROOT / "mockup_figures" / "data" / "illumina_wgs_age_sex_counts.tsv"
_THOUSANDS = StrMethodFormatter("{x:,.0f}")
_N_PCS = 10
_PC_ZLIM = 3
_PC_CMAP = LinearSegmentedColormap.from_list(
    "pc_diverging", ["#2C5F7C", "#9DBFD1", "#F7F7F5", "#E1B48F", "#9C4A1E"]
)


def _pyramid():
    path = ROOT / "mockup_figures" / "plot_age_sex_pyramid.py"
    spec = importlib.util.spec_from_file_location("age_sex_pyramid", path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"could not load {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _load_model_rows() -> list[dict[str, str]]:
    path = covariates_path()
    if path is None:
        raise RuntimeError("covariates were not found")
    rows = []
    with gzip.open(path, "rt", newline="") as handle:
        for row in csv.DictReader(handle):
            if row.get("technology") != "PacBio":
                continue
            if str(row.get("final_releasable_v9", "")).strip() not in {"True", "true"}:
                continue
            if row.get("lr_phase") not in ("phase_2", "phase_1_phase_2"):
                continue
            rows.append(row)
    return rows


def _load_phase1_rows() -> list[dict[str, str]]:
    """Phase 1 PacBio discovery cohort, including the 14 people resequenced in Phase 2."""
    path = covariates_path()
    if path is None:
        raise RuntimeError("covariates were not found")
    rows = []
    with gzip.open(path, "rt", newline="") as handle:
        for row in csv.DictReader(handle):
            if row.get("technology") != "PacBio":
                continue
            if row.get("lr_phase") not in ("phase_1", "phase_1_phase_2"):
                continue
            rows.append(row)
    return rows


def _phase1_coverage(rows: list[dict[str, str]]) -> np.ndarray:
    """Phase 1 sequencing depth.

    Overlap samples keep their Phase 2 depth on the Phase 2 histogram, so the
    Phase 1 distribution uses the original Phase 1 mosdepth value for those 14.
    """
    values = []
    for row in rows:
        key = "coverage" if row.get("lr_phase") == "phase_1" else "phase1_mosdepth_cov"
        raw = row.get(key)
        if raw in ("", None):
            continue
        values.append(float(raw))
    return np.asarray(values, dtype=float)


def _load_ont_only_counts() -> Counter:
    """Releasable Phase 2 participants sequenced only with ONT, by center."""
    path = covariates_path()
    counts: Counter = Counter()
    with gzip.open(path, "rt", newline="") as handle:
        for row in csv.DictReader(handle):
            if row.get("technology") != "ONT":
                continue
            if str(row.get("final_releasable_v9", "")).strip() not in {"True", "true"}:
                continue
            counts[row.get("GC") or ""] += 1
    return counts


def _load_illumina_reference() -> tuple[Counter, Counter, int]:
    """Illumina WGS age×sex counts from the Workbench CDR extract."""
    female: Counter = Counter()
    male: Counter = Counter()
    with _ILLUMINA_COUNTS.open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            side = row["sex"]
            if side not in {"Female", "Male"}:
                continue
            age = min(int(row["age"]), 90)
            n = int(row["n"])
            (female if side == "Female" else male)[age] += n
    return female, male, sum(female.values()) + sum(male.values())


def _truthy(value: str | None) -> bool:
    return str(value or "").strip().lower() in {"true", "1", "yes"}


def _panel_age_sex(ax, rows, illumina_f, illumina_m, illumina_n, pyr, phase1_rows=None):
    """Rotated age–sex comparison, split by phenotype availability.

    ``phase1_rows`` draws the Phase 1 age–sex outline on the same axis.
    """
    ages = np.arange(20, 91)
    layers = {
        (side, has_pheno): Counter()
        for side in ("Female", "Male")
        for has_pheno in (True, False)
    }
    for row in rows:
        side, _source = pyr.placement(row)
        if side not in {"Female", "Male"}:
            continue
        age = min(int(float(row["age"])), 90)
        layers[(side, _truthy(row.get("in_phenotype_table")))][age] += 1

    ref_f = np.array([illumina_f[int(age)] for age in ages], dtype=float)
    ref_m = np.array([illumina_m[int(age)] for age in ages], dtype=float)
    ax_ref = ax.twinx()
    # Year-wide filled steps, so the Illumina counts read as a silhouette
    # around the narrower long-read bars. Dual axes keep both magnitudes
    # readable; a shared linear scale would shrink the long-read bars away.
    for values, sign in ((ref_m, 1.0), (ref_f, -1.0)):
        ax_ref.fill_between(
            ages, 0, sign * values, step="mid",
            facecolor=pyr.POP, edgecolor="none", linewidth=0, zorder=1,
        )

    sex_style = {
        "Female": (_FEMALE, _FEMALE_PALE, -1.0),
        "Male": (_MALE, _MALE_PALE, 1.0),
    }
    totals = {"Female": np.zeros(ages.size), "Male": np.zeros(ages.size)}
    for side in ("Female", "Male"):
        color, pale, sign = sex_style[side]
        bottom = np.zeros(ages.size)
        for has_pheno in (True, False):
            values = np.array(
                [layers[(side, has_pheno)][int(age)] for age in ages],
                dtype=float,
            )
            face = color if has_pheno else pale
            ax.bar(
                ages,
                sign * values,
                bottom=sign * bottom,
                width=0.58,
                facecolor=face,
                edgecolor="none",
                linewidth=0,
                zorder=3,
            )
            bottom += values
        totals[side] = bottom

    if phase1_rows is not None:
        # One outline per sex, not split by the Phase 2 phenotype table.
        # Phase 1 had EHRs but is absent from that table.
        p1_layers = {"Female": Counter(), "Male": Counter()}
        for row in phase1_rows:
            side, _source = pyr.placement(row)
            if side not in {"Female", "Male"}:
                continue
            age = min(int(float(row["age"])), 90)
            p1_layers[side][age] += 1
        halo = [patheffects.withStroke(linewidth=2.6, foreground="white")]
        for side, sign in (("Male", 1.0), ("Female", -1.0)):
            values = np.array([p1_layers[side][int(age)] for age in ages], dtype=float)
            ax.step(
                ages, sign * values, where="mid", color=_PHASE1, lw=1.25,
                zorder=6, path_effects=halo,
            )

    ax.axhline(0, color="#222222", linewidth=0.7, zorder=4)
    ax.set_zorder(ax_ref.get_zorder() + 1)
    ax.patch.set_visible(False)
    focal_lim = pyr.nice_ceil(max(totals["Female"].max(), totals["Male"].max()) * 1.08, 25)
    ref_lim = pyr.nice_ceil(max(ref_f.max(), ref_m.max()) * 1.08, 2_500)
    ax.set_xlim(19.4, 93.0)
    ax.set_ylim(-focal_lim, focal_lim)
    ax_ref.set_ylim(-ref_lim, ref_lim)
    ax.set_xticks(list(range(20, 91, 10)))
    ax.set_xticklabels(["20", "30", "40", "50", "60", "70", "80", "90+"])
    ax.set_xlabel("Age (years)")
    ax.set_ylabel("Long-read participants" if phase1_rows is not None else "Phase 2 participants")
    # The gray right-axis ticks are identified by the legend; omitting a
    # second vertical label keeps this compact panel clear of panel B.
    ax.yaxis.set_major_formatter(pyr.plt.FuncFormatter(pyr.cohort_tick))
    ax_ref.yaxis.set_major_formatter(pyr.plt.FuncFormatter(pyr.cohort_tick))
    ax_ref.tick_params(axis="y", colors=pyr.GRAY)
    ax_ref.spines["top"].set_visible(False)
    ax_ref.spines["left"].set_visible(False)
    ax_ref.spines["right"].set_color(pyr.GRAY)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    for text, y, va, color in (("Male", 0.97, "top", _MALE), ("Female", 0.03, "bottom", _FEMALE)):
        ax.text(0.98, y, text, transform=ax.transAxes, ha="right", va=va,
                fontsize=MIN_FONT, color=color, fontweight="bold")

    n_with = sum(sum(layers[(side, True)].values()) for side in ("Female", "Male"))
    n_without = sum(sum(layers[(side, False)].values()) for side in ("Female", "Male"))

    def swatch(*colors):
        return tuple(Patch(facecolor=c, edgecolor="none") for c in colors)

    # Columns fill first. Without Phase 1, the blank entry puts Illumina
    # alone on the first row. With Phase 1, the long-read cohorts fill the
    # left column (left axis) and Illumina sits alone on the right (right axis).
    if phase1_rows is None:
        handles = [
            Patch(facecolor=pyr.POP, edgecolor="none"),
            swatch(_FEMALE, _MALE),
            Patch(facecolor="none", edgecolor="none"),
            swatch(_FEMALE_PALE, _MALE_PALE),
        ]
        labels = [
            f"Illumina WGS ({illumina_n:,})",
            f"Long-read, with phenotypes ({n_with:,})",
            "",
            f"Long-read, without phenotypes ({n_without:,})",
        ]
    else:
        handles = [
            Line2D([0], [0], color=_PHASE1, lw=1.25),
            swatch(_FEMALE, _MALE),
            swatch(_FEMALE_PALE, _MALE_PALE),
            Patch(facecolor=pyr.POP, edgecolor="none"),
            Patch(facecolor="none", edgecolor="none"),
        ]
        labels = [
            f"Phase 1 ({len(phase1_rows):,})",
            f"Phase 2, with phenotypes ({n_with:,})",
            f"Phase 2, without phenotypes ({n_without:,})",
            f"Illumina WGS ({illumina_n:,})",
            "",
        ]
    return ax_ref, handles, labels


def _pc_data(rows: list[dict[str, str]]) -> tuple[np.ndarray, np.ndarray]:
    """PC1–PC10 for participants with all ten, and their computed ancestry labels."""
    keys = [f"PC{i}" for i in range(1, _N_PCS + 1)]
    kept = [r for r in rows if all(r.get(k) not in ("", None) for k in keys)]
    pcs = np.array([[float(r[k]) for k in keys] for r in kept])
    labels = np.array([
        anc if (anc := (r.get("ancestry_pred_other") or "oth").upper()) in ANC_COLORS else "OTH"
        for r in kept
    ])
    return pcs, labels


def _draw_order(labels: np.ndarray) -> list[str]:
    """Largest groups first so smaller ones stay visible on top."""
    return sorted((a for a in ANC_ORDER if (labels == a).any()), key=lambda a: -(labels == a).sum())


def _scatter_pair(ax, pcs, labels, i, j, size=3.0, square=True) -> None:
    for name in _draw_order(labels):
        mask = labels == name
        ax.scatter(
            pcs[mask, i], pcs[mask, j],
            s=size, c=ANC_COLORS[name], alpha=0.45,
            linewidths=0, rasterized=True, zorder=2,
        )
    if square:
        x, y = pcs[:, i], pcs[:, j]
        half = max(float(np.ptp(x)), float(np.ptp(y))) * 1.08 / 2
        xc, yc = (x.min() + x.max()) / 2, (y.min() + y.max()) / 2
        ax.set_xlim(xc - half, xc + half)
        ax.set_ylim(yc - half, yc + half)
        ax.set_aspect("equal", adjustable="box")


def _anc_handles(markersize=4):
    return [
        Line2D([0], [0], marker="o", ls="none", color=ANC_COLORS[a], markersize=markersize, label=a)
        for a in ANC_ORDER
    ]


def _anc_legend_above(ax) -> None:
    # Legends fill column by column; interleave so the rows read in ANC_ORDER.
    ncol = 4
    handles = _anc_handles()
    handles = [h for c in range(ncol) for h in handles[c::ncol]]
    ax.legend(
        handles=handles,
        loc="lower center",
        bbox_to_anchor=(0.5, 1.08),
        ncol=ncol,
        frameon=False,
        fontsize=MIN_FONT,
        handlelength=0.8,
        handletextpad=0.15,
        columnspacing=0.6,
        borderaxespad=0.0,
    )


def _panel_pcs_scatter(ax, rows) -> tuple:
    pcs, labels = _pc_data(rows)
    _scatter_pair(ax, pcs, labels, 0, 1)
    ax.set_xlabel("PC1")
    ax.set_ylabel("PC2")
    _anc_legend_above(ax)
    return ()


def _panel_pcs_inset(ax, rows) -> tuple:
    """PC1 vs PC2, with PC3 vs PC5 inset in the emptiest corner."""
    pcs, labels = _pc_data(rows)
    _scatter_pair(ax, pcs, labels, 0, 1)
    ax.set_xlabel("PC1")
    ax.set_ylabel("PC2")
    _anc_legend_above(ax)
    xy = ax.transAxes.inverted().transform(ax.transData.transform(pcs[:, :2]))
    side = 0.42
    corners = {
        (0.0, 0.0): (xy[:, 0] < side) & (xy[:, 1] < side),
        (1 - side, 0.0): (xy[:, 0] > 1 - side) & (xy[:, 1] < side),
        (0.0, 1 - side): (xy[:, 0] < side) & (xy[:, 1] > 1 - side),
        (1 - side, 1 - side): (xy[:, 0] > 1 - side) & (xy[:, 1] > 1 - side),
    }
    x0, y0 = min(corners, key=lambda c: corners[c].sum())
    pad = 0.03
    x0 = x0 + pad if x0 == 0 else x0 - pad
    y0 = y0 + pad if y0 == 0 else y0 - pad
    inset = ax.inset_axes([x0, y0, side, side])
    _scatter_pair(inset, pcs, labels, 2, 4, size=1.2)
    inset.set_xticks([])
    inset.set_yticks([])
    for spine in inset.spines.values():
        spine.set_visible(True)
        spine.set_linewidth(0.5)
        spine.set_color("#666666")
    inset.text(0.97, 0.04, "PC3", transform=inset.transAxes, ha="right", va="bottom", fontsize=MIN_FONT)
    inset.text(0.04, 0.04, "PC5", transform=inset.transAxes, ha="left", va="bottom",
               rotation=90, fontsize=MIN_FONT)
    inset.set_facecolor("white")
    inset.set_zorder(5)
    return (inset,)


def _load_reference_pcs(keys) -> tuple[np.ndarray, np.ndarray]:
    """1000 Genomes controls projected into the long-read PC space."""
    refs, labels = [], []
    with gzip.open(covariates_path(), "rt", newline="") as handle:
        for row in csv.DictReader(handle):
            if not _truthy(row.get("is_reference_control")):
                continue
            if any(row.get(k) in ("", None) for k in keys) or row.get("population") == "OTH":
                continue
            refs.append([float(row[k]) for k in keys])
            labels.append(row["population"])
    return np.array(refs), np.array(labels)


def _panel_pcs_clines(ax, rows) -> tuple:
    """Participants spread along clines between tight reference clusters."""
    keys = ("lr_PC1", "lr_PC2")
    kept = [r for r in rows if all(r.get(k) not in ("", None) for k in keys)]
    pcs = np.array([[float(r[k]) for k in keys] for r in kept])
    labels = np.array([
        anc if (anc := (r.get("ancestry_pred_other") or "oth").upper()) in ANC_COLORS else "OTH"
        for r in kept
    ])
    _scatter_pair(ax, pcs, labels, 0, 1, size=2.4)
    ref, ref_labels = _load_reference_pcs(keys)
    ax.scatter(ref[:, 0], ref[:, 1], s=2.2, c="#111111", linewidths=0, zorder=4)
    offsets = {
        "AFR": ((-4, 9), "left", "bottom"),
        "EAS": ((-8, 0), "right", "center"),
        "SAS": ((-14, 0), "right", "center"),
        "EUR": ((-10, -6), "right", "top"),
    }
    halo = [patheffects.withStroke(linewidth=2.2, foreground="white")]
    for name, (offset, ha, va) in offsets.items():
        core = np.median(ref[ref_labels == name], axis=0)
        ax.annotate(
            f"1kGP {name}", xy=core, xytext=offset, textcoords="offset points",
            ha=ha, va=va, fontsize=MIN_FONT, color=ANC_TEXT[name],
            fontweight="bold", zorder=6, path_effects=halo,
        )
    ax.set_xlabel("Long-read PC1")
    ax.set_ylabel("Long-read PC2")
    ax.set_xticks([])
    ax.set_yticks([])

    ncol = 4
    handles = _anc_handles() + [
        Line2D([0], [0], marker="o", ls="none", color="#111111", markersize=2.5, label="1kGP")
    ]
    handles = [h for c in range(ncol) for h in handles[c::ncol]]
    ax.legend(
        handles=handles, loc="lower center", bbox_to_anchor=(0.5, 1.08), ncol=ncol,
        frameon=False, fontsize=MIN_FONT, handlelength=0.8, handletextpad=0.15,
        columnspacing=0.6, borderaxespad=0.0,
    )
    return ()


def _clines_data(rows, keys):
    kept = [r for r in rows if all(r.get(k) not in ("", None) for k in keys)]
    pcs = np.array([[float(r[k]) for k in keys] for r in kept])
    labels = np.array([
        anc if (anc := (r.get("ancestry_pred_other") or "oth").upper()) in ANC_COLORS else "OTH"
        for r in kept
    ])
    ref, ref_labels = _load_reference_pcs(keys)
    return pcs, labels, ref, ref_labels


def _project_phase1_lr_pcs(phase2_rows, phase1_rows) -> tuple[np.ndarray, np.ndarray]:
    """Place Phase 1 in long-read PC1–PC3 using a linear map fit on Phase 2.

    Phase 2 participants with both short-read PC1–PC5 and long-read PCs define
    the map. Those five short-read axes already track the long-read coordinates
    (holdout R² above 0.99). Phase 1 samples that were resequenced in Phase 2
    keep their long-read coordinates and are not projected again.
    """
    sr = [f"PC{i}" for i in range(1, 6)]
    lr = ["lr_PC1", "lr_PC2", "lr_PC3"]

    def vec(row, keys):
        out = []
        for key in keys:
            raw = row.get(key)
            if raw in ("", None):
                return None
            out.append(float(raw))
        return out

    train_x, train_y = [], []
    for row in phase2_rows:
        x, y = vec(row, sr), vec(row, lr)
        if x is not None and y is not None:
            train_x.append(x)
            train_y.append(y)
    design = np.column_stack([np.ones(len(train_x)), np.asarray(train_x, dtype=float)])
    coef, *_ = np.linalg.lstsq(design, np.asarray(train_y, dtype=float), rcond=None)
    proj, labels = [], []
    for row in phase1_rows:
        if row.get("lr_phase") != "phase_1":
            continue
        x = vec(row, sr)
        if x is not None:
            proj.append(x)
            anc = (row.get("ancestry_pred_other") or "oth").upper()
            labels.append(anc if anc in ANC_COLORS else "OTH")
    predictors = np.asarray(proj, dtype=float)
    return np.column_stack([np.ones(len(predictors)), predictors]) @ coef, np.array(labels)


def _nearest_point(points: np.ndarray, target: np.ndarray) -> np.ndarray:
    """Return the row of *points* closest to *target*."""
    delta = points - np.asarray(target, dtype=float)
    return points[np.argmin(np.einsum("ij,ij->i", delta, delta))]


def _tip_toward(points: np.ndarray, axes: tuple[int, int], offset: tuple[float, float]) -> np.ndarray:
    """Pick a real sample near the cluster edge facing the label offset."""
    xy = points[:, list(axes)]
    center = np.median(xy, axis=0)
    direction = np.array([-offset[0], -offset[1]], dtype=float)
    if not np.any(direction):
        return _nearest_point(xy, center)
    direction /= np.linalg.norm(direction)
    score = (xy - center) @ direction
    # Prefer points on the label-facing side, still near the dense core.
    rank = score - 0.35 * np.linalg.norm(xy - center, axis=1)
    return xy[np.argmax(rank)]


def _phase1_tips(phase1, phase1_labels, i, j) -> list[np.ndarray]:
    """Representative AFR and OTH Phase 1 samples to hang leaders on."""
    tips = []
    afr = phase1[phase1_labels == "AFR"]
    oth = phase1[phase1_labels == "OTH"]
    if len(afr):
        tips.append(_tip_toward(afr, (i, j), (0, -12)))
    if len(oth):
        tips.append(_nearest_point(oth[:, [i, j]], np.median(oth[:, [i, j]], axis=0)))
    return tips


def _annotate_phase1(ax, tips, halo, *, lower_panel=False) -> None:
    """One Phase 1 label with leaders forking to AFR and OTH samples."""
    if not tips:
        return
    tips = np.asarray(tips, dtype=float)
    y0, y1 = ax.get_ylim()
    span = y1 - y0
    # The lower PCA has open space beneath its central cloud; use it instead of
    # laying the label over the points.
    junction_y = y0 + 0.28 * span if lower_panel else np.min(tips[:, 1]) - 0.07 * span
    junction = np.array([float(np.mean(tips[:, 0])), float(junction_y)])
    for tip in tips:
        ax.annotate(
            "", xy=tip, xytext=junction, textcoords="data",
            arrowprops=dict(
                arrowstyle="-", color="#111111", lw=0.55, shrinkA=1.5, shrinkB=0,
                path_effects=halo,
            ),
            zorder=6,
        )
    ax.annotate(
        "Phase 1", xy=junction, xytext=(0, -2), textcoords="offset points",
        ha="center", va="top", fontsize=MIN_FONT, color="#111111",
        fontweight="bold", zorder=6, path_effects=halo,
    )


def _clines_view(
    ax, pcs, labels, ref, ref_labels, i, j, notes,
    phase1=None, phase1_labels=None, show_ref: bool = True,
) -> None:
    """One PC projection: participants by ancestry, references in black, labeled.

    ``show_ref`` draws the 1kGP controls and lets them widen the axis limits.
    With it off the view is set by the participants alone.
    """
    for name in _draw_order(labels):
        mask = labels == name
        ax.scatter(
            pcs[mask, i], pcs[mask, j], s=1.8, c=ANC_COLORS[name], alpha=0.45,
            linewidths=0, rasterized=True, zorder=2,
        )
    if phase1 is not None and len(phase1):
        # Larger than the Phase 2 dots, filled with ancestry color, ringed in black.
        for name in _draw_order(phase1_labels):
            mask = phase1_labels == name
            ax.scatter(
                phase1[mask, i], phase1[mask, j], s=4.5, c=ANC_COLORS[name],
                edgecolors="#111111", linewidths=0.3, marker="o",
                rasterized=True, zorder=3,
            )
    if show_ref:
        ax.scatter(ref[:, i], ref[:, j], s=2.0, c="#111111", linewidths=0, zorder=4)
    pad = 0.05
    for axis, k in ((ax.set_xlim, i), (ax.set_ylim, j)):
        cols = [pcs[:, k]] + ([ref[:, k]] if show_ref else [])
        if phase1 is not None and len(phase1):
            cols.append(phase1[:, k])
        lo = min(float(np.min(c)) for c in cols)
        hi = max(float(np.max(c)) for c in cols)
        axis(lo - pad * (hi - lo), hi + pad * (hi - lo))
    halo = [patheffects.withStroke(linewidth=2.2, foreground="white")]
    for text, xy, offset, ha, va, color in notes:
        ax.annotate(
            text, xy=xy, xytext=offset, textcoords="offset points", ha=ha, va=va,
            fontsize=MIN_FONT, color=color, fontweight="bold", zorder=6, path_effects=halo,
            arrowprops=dict(
                arrowstyle="-", color=color, lw=0.55, shrinkA=1, shrinkB=1.5,
                path_effects=halo,
            ),
        )
    if phase1 is not None and len(phase1) and phase1_labels is not None:
        _annotate_phase1(
            ax, _phase1_tips(phase1, phase1_labels, i, j), halo,
            lower_panel=(j == 2),
        )
    ax.set_xticks([])
    ax.set_yticks([])


def _panel_pcs_clines2(
    ax,
    rows,
    phase1_xy: np.ndarray | None = None,
    phase1_labels: np.ndarray | None = None,
    with_reference: bool = False,
) -> tuple:
    """PC1 against PC2 and PC3: continental clines, then the Indigenous American axis.

    ``with_reference`` overlays the 291 1kGP controls in black and labels their
    cluster cores. The participants already carry their similarity group in
    colour, so the default view leaves them out; set it True to show where the
    reference populations sit in this PC space.
    """
    pcs, labels, ref, ref_labels = _clines_data(rows, ("lr_PC1", "lr_PC2", "lr_PC3"))
    ax.axis("off")
    top = ax.inset_axes([0, 0.52, 1, 0.48])
    bottom = ax.inset_axes([0, 0, 1, 0.48])
    top_notes = [
        ("1kGP AFR", (0, 1), (-2, 12), "left", "bottom", ANC_TEXT["AFR"], "AFR"),
        ("1kGP EAS", (0, 1), (-14, 2), "right", "center", ANC_TEXT["EAS"], "EAS"),
        ("1kGP SAS", (0, 1), (-16, -4), "right", "center", ANC_TEXT["SAS"], "SAS"),
        ("1kGP EUR", (0, 1), (-14, -8), "right", "top", ANC_TEXT["EUR"], "EUR"),
    ]
    _clines_view(top, pcs, labels, ref, ref_labels, 0, 1, [
        (text, _tip_toward(ref[ref_labels == name], axes, offset), offset, ha, va, color)
        for text, axes, offset, ha, va, color, name in top_notes
    ] if with_reference else [], phase1=phase1_xy, phase1_labels=phase1_labels,
        show_ref=with_reference)
    amr = ref[ref_labels == "AMR"]
    amr_offset = (-14, 0)
    _clines_view(bottom, pcs, labels, ref, ref_labels, 0, 2, [
        (
            "1kGP AFR",
            _tip_toward(ref[ref_labels == "AFR"], (0, 2), (-2, 12)),
            (-2, 12), "left", "bottom", ANC_TEXT["AFR"],
        ),
        (
            "1kGP AMR",
            _tip_toward(amr, (0, 2), amr_offset),
            amr_offset, "right", "center", ANC_TEXT["AMR"],
        ),
    ] if with_reference else [], phase1=phase1_xy, phase1_labels=phase1_labels,
        show_ref=with_reference)
    top.set_ylabel("PC2")
    bottom.set_ylabel("PC3")
    bottom.set_xlabel("Long-read PC1")
    for sub in (top, bottom):
        sub.yaxis.label.set_rotation(90)
    by_label = {h.get_label(): h for h in _anc_handles()}
    if phase1_xy is not None:
        # Three rows, bottoms aligned with panel A's legend (bbox y=1.08, default
        # labelspacing). Column-major fill for ncol=4:
        #   AFR  AMR  EAS  OTH
        #   MID  SAS  EUR
        #  1kGP  Phase 1     (1kGP only when with_reference)
        blank = Line2D([0], [0], ls="none", marker="", label=" ")
        comparison = {
            "1kGP": Line2D(
                [0], [0], marker="o", ls="none", color="#111111", markersize=2.5, label="1kGP",
            ),
            "Phase 1": Line2D(
                [0], [0], marker="o", ls="none",
                markerfacecolor="white", markeredgecolor="#111111",
                markersize=4.2, markeredgewidth=0.7, label="Phase 1",
            ),
        }
        bottom_row = (
            [comparison["1kGP"], comparison["Phase 1"], blank]
            if with_reference
            else [comparison["Phase 1"], blank, blank]
        )
        handles = [
            by_label["AFR"], by_label["MID"], bottom_row[0],
            by_label["AMR"], by_label["SAS"], bottom_row[1],
            by_label["EAS"], by_label["EUR"], bottom_row[2],
            by_label["OTH"], blank, blank,
        ]
        ax.legend(
            handles=handles, loc="lower center", bbox_to_anchor=(0.5, 1.02), ncol=4,
            frameon=False, fontsize=MIN_FONT, handlelength=0.8, handletextpad=0.15,
            columnspacing=0.6, borderaxespad=0.0, labelspacing=0.5,
        )
    else:
        handles = _anc_handles() + [
            Line2D([0], [0], marker="o", ls="none", color="#111111", markersize=2.5, label="1kGP")
        ]
        handles = [h for c in range(4) for h in handles[c::4]]
        ax.legend(
            handles=handles, loc="lower center", bbox_to_anchor=(0.5, 1.02), ncol=4,
            frameon=False, fontsize=MIN_FONT, handlelength=0.8, handletextpad=0.15,
            columnspacing=0.6, borderaxespad=0.0, labelspacing=0.5,
        )
    return top, bottom


def _panel_pcs_multiples(ax, rows) -> tuple:
    """Five mini-scatters, one per consecutive PC pair, plus a legend tile."""
    pcs, labels = _pc_data(rows)
    ax.axis("off")
    gap_x, gap_y = 0.12, 0.16
    w = (1 - 2 * gap_x) / 3
    h = (1 - gap_y) / 2
    minis = []
    for k in range(6):
        col, row = k % 3, k // 3
        rect = [col * (w + gap_x), (1 - row) * (h + gap_y), w, h]
        sub = ax.inset_axes(rect)
        if k == 5:
            sub.axis("off")
            sub.legend(
                handles=_anc_handles(), loc="center", frameon=False, fontsize=MIN_FONT,
                handlelength=0.8, handletextpad=0.3, labelspacing=0.25, borderaxespad=0.0,
            )
        else:
            i, j = 2 * k, 2 * k + 1
            _scatter_pair(sub, pcs, labels, i, j, size=0.9, square=False)
            sub.set_xticks([])
            sub.set_yticks([])
            sub.set_xlabel(f"PC{i + 1}", labelpad=1)
            sub.set_ylabel(f"PC{j + 1}", labelpad=1)
            for spine in ("top", "right"):
                sub.spines[spine].set_visible(False)
        minis.append(sub)
    return tuple(minis)


def _panel_pcs_strips(ax, rows) -> tuple:
    """Every participant's standardized score on each null-model PC."""
    pcs, labels = _pc_data(rows)
    z = (pcs - pcs.mean(axis=0)) / pcs.std(axis=0)
    rng = np.random.default_rng(7)
    jitter = rng.uniform(-0.34, 0.34, size=z.shape)
    for name in _draw_order(labels):
        mask = labels == name
        for k in range(_N_PCS):
            ax.scatter(
                z[mask, k], np.full(mask.sum(), k) + jitter[mask, k],
                s=0.5, c=ANC_COLORS[name], alpha=0.35,
                linewidths=0, rasterized=True, zorder=2,
            )
    lim = float(np.ceil(np.percentile(np.abs(z), 99.9)))
    ax.set_xlim(-lim, lim)
    ax.set_ylim(_N_PCS - 0.5, -0.5)
    ax.set_yticks(range(_N_PCS))
    ax.set_yticklabels([f"PC{k + 1}" for k in range(_N_PCS)])
    ax.tick_params(axis="y", length=0)
    ax.axvline(0, color="#BBBBBB", linewidth=0.5, zorder=1)
    ax.set_xlabel("Standardized score (SD)")
    _anc_legend_above(ax)
    return ()


def _panel_pcs_heatmap(ax, rows) -> tuple:
    """Mean of each null-model PC by ancestry group, with between-group R² above."""
    pcs, labels = _pc_data(rows)
    z = (pcs - pcs.mean(axis=0)) / pcs.std(axis=0)
    groups = [a for a in ANC_ORDER if (labels == a).any()]
    means = np.array([z[labels == a].mean(axis=0) for a in groups])
    sizes = np.array([(labels == a).sum() for a in groups])
    r2 = (sizes[:, None] * means**2).sum(axis=0) / len(z)

    image = ax.imshow(
        means, cmap=_PC_CMAP, vmin=-_PC_ZLIM, vmax=_PC_ZLIM,
        aspect="auto", interpolation="nearest",
    )
    ax.set_xticks(range(_N_PCS))
    ax.set_xticklabels([str(i) for i in range(1, _N_PCS + 1)])
    ax.set_yticks(range(len(groups)))
    ax.set_yticklabels(groups)
    for tick, name in zip(ax.get_yticklabels(), groups):
        tick.set_color(ANC_TEXT[name])
        tick.set_fontweight("bold")
    ax.tick_params(length=0, pad=2)
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.set_xlabel("Genetic principal component")

    ax_r2 = ax.inset_axes([0, 1.04, 1, 0.24])
    ax_r2.bar(range(_N_PCS), r2, width=0.78, color="#5A6670", linewidth=0)
    ax_r2.set_xlim(-0.5, _N_PCS - 0.5)
    ax_r2.set_ylim(0, 1)
    ax_r2.set_yticks([0, 1])
    ax_r2.set_xticks([])
    ax_r2.spines["bottom"].set_visible(False)
    ax_r2.set_ylabel("Group R²", rotation=0, ha="right", va="center", labelpad=2)

    ax_bar = ax.inset_axes([1.04, 0, 0.045, 1])
    bar = ax.figure.colorbar(image, cax=ax_bar, extend="both")
    bar.set_ticks([-_PC_ZLIM, 0, _PC_ZLIM])
    bar.outline.set_linewidth(0.4)
    bar.ax.tick_params(length=2, pad=1.5)
    bar.set_label("Group mean (SD)", labelpad=2)
    return ax_r2, ax_bar


_PC_PANELS = {
    "scatter": _panel_pcs_scatter,
    "inset": _panel_pcs_inset,
    "clines": _panel_pcs_clines,
    "clines2": _panel_pcs_clines2,
    "multiples": _panel_pcs_multiples,
    "strips": _panel_pcs_strips,
    "heatmap": _panel_pcs_heatmap,
}


def _high_pass_mask(rows: list[dict[str, str]]) -> np.ndarray:
    cov = np.array(
        [float(row["coverage"]) if row.get("coverage") not in ("", None) else np.nan for row in rows],
        dtype=float,
    )
    high = np.array([_truthy(row.get("has_ONT")) for row in rows], dtype=bool)
    n_needed = _HIGH_PASS_N - int(high.sum())
    if n_needed > 0:
        eligible = np.flatnonzero(~high & np.isfinite(cov))
        chosen = eligible[np.argsort(cov[eligible])[-n_needed:]]
        high[chosen] = True
    return high


def _panel_coverage(ax, rows: list[dict[str, str]], phase1_cov: np.ndarray | None = None) -> None:
    cov = np.array(
        [float(row["coverage"]) if row.get("coverage") not in ("", None) else np.nan for row in rows],
        dtype=float,
    )
    finite = np.isfinite(cov)
    high = _high_pass_mask(rows) & finite
    mid = ~high & finite
    bins = np.arange(4 if phase1_cov is not None else 7, 53, 1)
    ax.hist(
        cov[mid], bins=bins, color=_MID_PASS_COLOR,
        alpha=0.9, linewidth=0, label=f"Mid-pass (n={int(mid.sum()):,})",
    )
    ax.hist(
        cov[high], bins=bins, color=_HIGH_PASS_COLOR,
        alpha=0.85, linewidth=0, label=f"High-pass (n={int(high.sum()):,})",
    )
    if phase1_cov is not None and phase1_cov.size:
        ax.hist(
            phase1_cov, bins=bins, color=_PHASE1, alpha=0.55, linewidth=0, zorder=3,
        )
        ax.hist(
            phase1_cov, bins=bins, histtype="step", color=_PHASE1, lw=1.15,
            zorder=4, label=f"Phase 1 (n={int(phase1_cov.size):,})",
        )
        p1_mean = float(phase1_cov.mean())
        p1_counts, _edges = np.histogram(phase1_cov, bins=bins)
        ax.axvline(p1_mean, color=_PHASE1, ls=":", lw=0.9, zorder=4)
        # Just above the Phase 1 mode. The mid-pass and high-pass means stay at the top.
        ax.annotate(
            f"{p1_mean:.1f}×",
            xy=(p1_mean, float(p1_counts.max())),
            xytext=(3, 1),
            textcoords="offset points",
            ha="left",
            va="bottom",
            fontsize=MIN_FONT,
            color=_PHASE1,
            path_effects=[patheffects.withStroke(linewidth=2.0, foreground="white")],
            zorder=5,
        )
    mid_mean = float(cov[mid].mean())
    high_mean = float(cov[high].mean())
    ax.axvline(mid_mean, color="#526C7C", ls=":", lw=0.9)
    ax.axvline(high_mean, color=_HIGH_PASS_COLOR, ls=":", lw=0.9)
    ax.set_xlim(4 if phase1_cov is not None else 7, 52)
    ax.set_ylim(0, ax.get_ylim()[1] * 1.18)
    ax.yaxis.set_major_formatter(_THOUSANDS)
    ax.set_xlabel("PacBio coverage (×)")
    ax.set_ylabel("Participants")
    ax.annotate(
        f"{mid_mean:.1f}×",
        xy=(mid_mean, 1.0),
        xycoords=("data", "axes fraction"),
        xytext=(4, -2),
        textcoords="offset points",
        ha="left",
        va="top",
        fontsize=MIN_FONT,
        color="#526C7C",
    )
    ax.annotate(
        f"{high_mean:.1f}×",
        xy=(high_mean, 0.82),
        xycoords=("data", "axes fraction"),
        xytext=(4, 0),
        textcoords="offset points",
        ha="left",
        va="top",
        fontsize=MIN_FONT,
        color=_HIGH_PASS_COLOR,
    )
    handles, labels = ax.get_legend_handles_labels()
    if phase1_cov is not None and len(handles) == 3:
        handles = [handles[2], handles[0], handles[1]]
        labels = [labels[2], labels[0], labels[1]]
    ax.legend(
        handles,
        labels,
        loc="lower right",
        bbox_to_anchor=(1.0, _ROW_LEGEND_Y),
        **_ROW_LEGEND,
    )


def _panel_center(
    ax,
    rows: list[dict[str, str]],
    ont_only: Counter,
    phase1_n: int | None = None,
    phase1_center: str = "HA",
) -> None:
    counts = {code: 0 for code, _name in _CENTERS}
    for row in rows:
        code = row.get("GC") or ""
        if code in counts:
            counts[code] += 1
    ordered = sorted(_CENTERS, key=lambda item: (-counts[item[0]], -ont_only[item[0]]))
    y = np.arange(len(ordered))
    high_mask = _high_pass_mask(rows)
    high_counts = {
        code: sum(mask and row.get("GC") == code for mask, row in zip(high_mask, rows))
        for code, _name in _CENTERS
    }
    ont_counts = {
        code: sum(_truthy(row.get("has_ONT")) and row.get("GC") == code for row in rows)
        for code, _name in _CENTERS
    }
    # Stack left to right: Phase 1 (HudsonAlpha only), then PacBio-only
    # (mid, then high-only), PacBio+ONT, and ONT only. ONT-overlap PacBio
    # samples are high-pass by construction, so they are not double-counted
    # in the PacBio-only high segment.
    mid_vals = [counts[code] - high_counts[code] for code, _name in ordered]
    high_vals = [high_counts[code] for code, _name in ordered]
    ont_vals = [ont_counts[code] for code, _name in ordered]
    high_only_vals = [high - ont for high, ont in zip(high_vals, ont_vals)]
    ont_only_vals = [ont_only[code] for code, _name in ordered]
    pacbio_vals = [mid + high for mid, high in zip(mid_vals, high_vals)]
    phase1_vals = [
        (phase1_n or 0) if code == phase1_center else 0
        for code, _name in ordered
    ]
    left_high = [p1 + mid for p1, mid in zip(phase1_vals, mid_vals)]
    left_ont = [p1 + mid + high_only for p1, mid, high_only in zip(phase1_vals, mid_vals, high_only_vals)]
    left_ont_only = [p1 + pacbio for p1, pacbio in zip(phase1_vals, pacbio_vals)]
    if phase1_n:
        ax.barh(
            y, phase1_vals, color=_PHASE1, height=0.72, linewidth=0,
            label="Phase 1",
        )
    ax.barh(y, mid_vals, left=phase1_vals, color=_MID_PASS_COLOR, height=0.72, linewidth=0, label="Mid-pass")
    ax.barh(
        y, high_only_vals, left=left_high, color=_HIGH_PASS_COLOR,
        height=0.72, linewidth=0, label="High-pass",
    )
    ax.barh(
        y, ont_vals, left=left_ont, color=_HIGH_PASS_COLOR, edgecolor=_ONT_COLOR,
        height=0.72, linewidth=0.55, hatch="////", label="ONT overlap",
    )
    ax.barh(
        y, ont_only_vals, left=left_ont_only, color=_ONT_ONLY_COLOR,
        height=0.72, linewidth=0, label="ONT only",
    )
    ax.set_yticks(y)
    ax.set_yticklabels([name for _code, name in ordered])
    ax.invert_yaxis()
    ax.set_xlabel("Participants")
    ax.set_xlim(0, max(p1 + v + o for p1, v, o in zip(phase1_vals, pacbio_vals, ont_only_vals)) * 1.18)
    ax.xaxis.set_major_formatter(_THOUSANDS)
    for yi, val, extra, p1 in zip(y, pacbio_vals, ont_only_vals, phase1_vals):
        end = val + extra + p1
        if p1:
            # Phase 1 count in the Phase 1 color, then the Phase 2 count in black,
            # same construction as the PacBio + ONT labels.
            anchor = ax.annotate(
                f"{p1:,}",
                xy=(end, yi),
                xytext=(3, 0),
                textcoords="offset points",
                ha="left",
                va="center",
                fontsize=MIN_FONT,
                color=_PHASE1,
            )
            anchor = ax.annotate(
                f" + {val:,}",
                xy=(1, 0.5),
                xycoords=anchor,
                ha="left",
                va="center",
                fontsize=MIN_FONT,
                color="#111111",
            )
        else:
            anchor = ax.annotate(
                f"{val:,}" if val else "",
                xy=(end, yi),
                xytext=(3, 0),
                textcoords="offset points",
                ha="left",
                va="center",
                fontsize=MIN_FONT,
                color="#111111",
            )
        if extra:
            ax.annotate(
                f" + {extra:,}" if val else f"{extra:,}",
                xy=(1, 0.5),
                xycoords=anchor,
                ha="left",
                va="center",
                fontsize=MIN_FONT,
                color=_ONT_COLOR,
            )
    handles, labels = ax.get_legend_handles_labels()
    # One legend, so its row pitch matches the other bottom-row legends.
    # Columns fill first. The blank keeps Phase 1 alone on the top row, and
    # Mid-pass / High-pass sit on the same two lower rows as panels D and F:
    # Phase 1 / blank, Mid-pass / ONT overlap, High-pass / ONT only.
    if phase1_n:
        blank = Patch(facecolor="none", edgecolor="none", linewidth=0)
        legend_handles = [
            handles[0], handles[1], handles[2],
            blank, handles[3], handles[4],
        ]
        legend_labels = [
            labels[0], labels[1], labels[2],
            " ", labels[3], labels[4],
        ]
    else:
        legend_order = [0, 2, 1, 3]
        legend_handles = [handles[i] for i in legend_order]
        legend_labels = [labels[i] for i in legend_order]
    ax.legend(
        legend_handles,
        legend_labels,
        loc="lower left",
        bbox_to_anchor=(0.0, _ROW_LEGEND_Y),
        ncol=2,
        columnspacing=0.7,
        **_ROW_LEGEND,
    )


def _panel_molecular_data(ax, rows: list[dict[str, str]]) -> None:
    """Horizontal UpSet of assayed modalities, stacked by QTL-analysis membership.

    Matrix rows are intersections of assayed 5mC / Olink / RNA-seq. Bars split
    into participants featured in the multi-omics QTL analyses (``in_pqtl`` for
    Olink; ``in_eqtl`` or ``in_sqtl`` for RNA-seq) versus assayed but not in
    those analysis sets.
    """
    avail: dict[tuple[bool, bool, bool], int] = {}
    featured: dict[tuple[bool, bool, bool], int] = {}
    for row in rows:
        key = (
            _truthy(row.get("has_methylation")),
            _truthy(row.get("has_proteomics")),
            _truthy(row.get("has_rna")),
        )
        if not any(key):
            continue
        avail[key] = avail.get(key, 0) + 1
        olink_ok = (not key[1]) or _truthy(row.get("in_pqtl"))
        rna_ok = (not key[2]) or _truthy(row.get("in_eqtl")) or _truthy(row.get("in_sqtl"))
        if olink_ok and rna_ok:
            featured[key] = featured.get(key, 0) + 1
    combos = sorted(avail, key=lambda key: -avail[key])
    feat_vals = [featured.get(key, 0) for key in combos]
    other_vals = [avail[key] - featured.get(key, 0) for key in combos]
    totals = [avail[key] for key in combos]
    y = np.arange(len(combos))

    ax.set_axis_off()
    ax_mat = ax.inset_axes([0.0, 0.0, 0.34, 1.0])
    ax_bar = ax.inset_axes([0.40, 0.0, 0.60, 1.0], sharey=ax_mat)
    ax_bar.barh(y, feat_vals, height=0.68, color=_MULTIOMIC_COLOR, linewidth=0, label="In QTL analysis")
    ax_bar.barh(
        y, other_vals, left=feat_vals, height=0.68,
        color=_MULTIOMIC_PALE, linewidth=0, label="Assayed only",
    )
    ax_bar.set_xlabel("Participants")
    ax_bar.set_ylim(len(combos) - 0.55, -0.55)
    ax_bar.set_yticks([])
    ax_bar.set_xlim(0, max(totals) * 1.28)
    ax_bar.xaxis.set_major_formatter(_THOUSANDS)
    for yi, total, feat, other in zip(y, totals, feat_vals, other_vals):
        label = f"{feat:,}" if other == 0 else f"{feat:,}+{other:,}"
        ax_bar.annotate(
            label,
            xy=(total, yi),
            xytext=(3, 0),
            textcoords="offset points",
            ha="left",
            va="center",
            fontsize=MIN_FONT,
        )
    ax_bar.spines["top"].set_visible(False)
    ax_bar.spines["right"].set_visible(False)
    ax_bar.spines["left"].set_visible(False)
    ax_bar.legend(
        loc="lower right",
        bbox_to_anchor=(1.0, _ROW_LEGEND_Y),
        ncol=1,
        **_ROW_LEGEND,
    )

    names = ["5mC", "Olink", "RNA-seq"]
    for yi, combo in zip(y, combos):
        present = [i for i, value in enumerate(combo) if value]
        if len(present) > 1:
            ax_mat.plot([min(present), max(present)], [yi, yi], color=_MULTIOMIC_COLOR, lw=1.0)
        for xi, value in enumerate(combo):
            ax_mat.scatter(
                xi,
                yi,
                s=18,
                color=_MULTIOMIC_COLOR if value else "#D9DEE2",
                edgecolors="none",
                zorder=2,
            )
    ax_mat.set_xticks(range(len(names)))
    ax_mat.set_xticklabels(names, rotation=90, ha="center", va="top")
    ax_mat.set_xlim(-0.5, len(names) - 0.5)
    ax_mat.set_yticks([])
    ax_mat.tick_params(axis="x", length=0, pad=2)
    for spine in ax_mat.spines.values():
        spine.set_visible(False)


def _panel_pedigrees(ax, rows: list[dict[str, str]]) -> None:
    """Family-size counts, stacked by whether members are mid-pass, high-pass, or mixed."""
    high = _high_pass_mask(rows)
    families: dict[str, dict[str, object]] = {}
    for row, is_high in zip(rows, high):
        if not _truthy(row.get("in_pedigree")):
            continue
        family_id = row.get("pedigree_family_id") or ""
        family_size = row.get("pedigree_family_size") or ""
        if not family_id or not family_size:
            continue
        entry = families.setdefault(
            family_id, {"size": int(float(family_size)), "n_mid": 0, "n_high": 0},
        )
        if is_high:
            entry["n_high"] = int(entry["n_high"]) + 1
        else:
            entry["n_mid"] = int(entry["n_mid"]) + 1

    sizes = np.array([3, 4, 5, 6])
    strata = ("mid", "mixed", "high")
    colors = {"mid": _MID_PASS_COLOR, "mixed": "#7A8F9E", "high": _HIGH_PASS_COLOR}
    labels = {"mid": "Mid-pass", "mixed": "Mixed", "high": "High-pass"}
    counts = {name: np.zeros(sizes.size, dtype=int) for name in strata}
    for entry in families.values():
        size = int(entry["size"])
        if size not in sizes:
            continue
        n_mid, n_high = int(entry["n_mid"]), int(entry["n_high"])
        if n_high == 0:
            kind = "mid"
        elif n_mid == 0:
            kind = "high"
        else:
            kind = "mixed"
        counts[kind][int(np.where(sizes == size)[0][0])] += 1

    bottom = np.zeros(sizes.size)
    for name in strata:
        ax.bar(
            sizes, counts[name], bottom=bottom, width=0.68,
            color=colors[name], linewidth=0, label=labels[name],
        )
        bottom += counts[name]
    totals = bottom
    ax.set_xticks(sizes)
    ax.set_xlabel("Family size")
    ax.set_ylabel("Families")
    ax.set_ylim(0, max(totals) * 1.28 if totals.max() else 1)
    for size, total in zip(sizes, totals):
        if total == 0:
            continue
        ax.annotate(
            f"{int(total):,}",
            xy=(size, total),
            xytext=(0, 2),
            textcoords="offset points",
            ha="center",
            va="bottom",
            fontsize=MIN_FONT,
        )
    handles, labels = ax.get_legend_handles_labels()
    # Mixed on the top row so Mid-pass and High-pass share the lower two
    # rows with panels C and D. Bar stacking is unchanged.
    by_label = dict(zip(labels, handles))
    ax.legend(
        [by_label["Mixed"], by_label["Mid-pass"], by_label["High-pass"]],
        ["Mixed", "Mid-pass", "High-pass"],
        loc="lower left",
        bbox_to_anchor=(0.0, _ROW_LEGEND_Y),
        ncol=1,
        **_ROW_LEGEND,
    )


# Spacing in inches, measured between the outermost ink of adjacent panels
# (tick labels, axis labels, legends), not between axes frames.
_MARGIN = 0.06
_GUTTER = 0.22
# Panel B's width as a multiple of the row height. It was square (1.0) while
# the ancestry legend read AFR/AMR/...; the longer population descriptors
# ("AFR-like", "Unassigned") push that four-column legend to about 2.28 in,
# which would overhang a 2.00 in panel. The extra width comes out of panel A,
# whose own legends leave roughly 0.9 in of slack between their two blocks.
_PC_ASPECT = 1.18
_ROW_GAP = 0.16
_LETTER_GAP = 0.04
_LETTER_SIZE = 9


class _Panel:
    """One panel in a row.

    ``weight`` shares out the row's flexible width. ``square`` instead fixes
    the width from the row height; ``aspect`` (width / height) widens or
    narrows that fixed box, which is how panel B buys room for a legend
    wider than its plot. Fixed width is taken out before the flexible
    panels are sized, so widening B narrows A by the same amount.
    """

    def __init__(
        self,
        letter: str,
        ax,
        extra=(),
        weight: float = 1.0,
        square: bool = False,
        aspect: float = 1.0,
    ):
        self.letter = letter
        self.ax = ax
        self.axes = [ax, *extra]
        self.weight = weight
        self.square = square
        self.aspect = aspect

    def fixed_width(self, height: float) -> float:
        return height * self.aspect

    def pads(self, renderer, dpi: float, fig_w: float, fig_h: float):
        ink = Bbox.union([a.get_tightbbox(renderer) for a in self.axes])
        pos = self.ax.get_position()
        return (
            pos.x0 * fig_w - ink.x0 / dpi,
            ink.x1 / dpi - pos.x1 * fig_w,
            pos.y0 * fig_h - ink.y0 / dpi,
            ink.y1 / dpi - pos.y1 * fig_h,
        )


def _layout(fig, rows_of_panels, heights) -> list[tuple[str, float, float]]:
    """Tile panels so their ink is separated by fixed gutters; return letter anchors.

    Rows are listed bottom to top. Axes bottoms align within a row, and every
    letter in a row shares one baseline just above that row's tallest ink.
    """
    fig_w, fig_h = fig.get_size_inches()
    dpi = fig.dpi
    letters: list[tuple[str, float, float]] = []
    for _ in range(3):
        fig.canvas.draw()
        renderer = fig.canvas.get_renderer()
        letter_h = _LETTER_SIZE / 72
        floor = _MARGIN
        letters = []
        for panels, height in zip(rows_of_panels, heights):
            pads = [p.pads(renderer, dpi, fig_w, fig_h) for p in panels]
            fixed = sum(p.fixed_width(height) for p in panels if p.square)
            flexible = (
                fig_w - 2 * _MARGIN
                - _GUTTER * (len(panels) - 1)
                - sum(left + right for left, right, _b, _t in pads)
                - fixed
            )
            total_weight = sum(p.weight for p in panels if not p.square)
            y = floor + max(bottom for _l, _r, bottom, _t in pads)
            row_top = y + height + max(top for _l, _r, _b, top in pads)
            cursor = _MARGIN
            for panel, (left, right, _b, _t) in zip(panels, pads):
                width = (
                    panel.fixed_width(height)
                    if panel.square
                    else flexible * panel.weight / total_weight
                )
                panel.ax.set_position([
                    (cursor + left) / fig_w, y / fig_h, width / fig_w, height / fig_h,
                ])
                letters.append((panel.letter, cursor, row_top + _LETTER_GAP))
                cursor += left + width + right + _GUTTER
            floor = row_top + _LETTER_GAP + letter_h + _ROW_GAP
    return letters


def render(
    pc_style: str = "clines2", *, with_phase1: bool = True, with_reference: bool = False,
) -> tuple:
    pyr = _pyramid()
    illumina_f, illumina_m, illumina_n = _load_illumina_reference()
    rows = _load_model_rows()
    phase1_rows = _load_phase1_rows() if with_phase1 else None

    # Model covariates and cohort structure above; sequencing design,
    # molecular follow-up, and pedigrees below.
    fig = new_figure(5.45 if with_phase1 else 5.90)
    ax_pyr = fig.add_axes([0.1, 0.5, 0.4, 0.4])
    ax_ref, handles, labels = _panel_age_sex(
        ax_pyr, rows, illumina_f, illumina_m, illumina_n, pyr, phase1_rows=phase1_rows,
    )
    legend_kwargs = dict(
        frameon=False,
        fontsize=MIN_FONT,
        handlelength=1.6,
        handleheight=0.7,
        borderaxespad=0.0,
        labelspacing=0.5,
        handler_map={tuple: HandlerTuple(ndivide=None, pad=0)},
    )
    if with_phase1:
        # Long-read items keep markers on the left. Illumina sits on the right
        # column with its gray swatch after the text (right y-axis). Blank rows
        # keep Illumina on the top shared baseline with Phase 1.
        blank = Patch(facecolor="none", edgecolor="none")
        left = ax_pyr.legend(
            handles[:3], labels[:3],
            loc="lower left", bbox_to_anchor=(0.0, 1.02), ncol=1,
            **legend_kwargs,
        )
        ax_pyr.add_artist(left)
        ax_pyr.legend(
            [handles[3], blank, blank], [labels[3], " ", " "],
            loc="lower right", bbox_to_anchor=(1.0, 1.02), ncol=1,
            markerfirst=False, **legend_kwargs,
        )
    else:
        ax_pyr.legend(
            handles, labels,
            loc="lower center", bbox_to_anchor=(0.5, 1.02), ncol=2,
            columnspacing=0.9, **legend_kwargs,
        )
    ax_pc, ax_gc, ax_cov, ax_molecular, ax_pedigree = (
        fig.add_axes([0.1 + 0.15 * i, 0.1, 0.1, 0.2]) for i in range(5)
    )
    if with_phase1:
        phase1_xy, phase1_labels = _project_phase1_lr_pcs(rows, phase1_rows)
        pc_extras = _panel_pcs_clines2(
            ax_pc, rows, phase1_xy=phase1_xy, phase1_labels=phase1_labels,
            with_reference=with_reference,
        )
    else:
        pc_extras = _PC_PANELS[pc_style](ax_pc, rows)
    _panel_center(
        ax_gc, rows, _load_ont_only_counts(),
        phase1_n=len(phase1_rows) if phase1_rows is not None else None,
    )
    _panel_coverage(
        ax_cov, rows,
        phase1_cov=_phase1_coverage(phase1_rows) if phase1_rows is not None else None,
    )
    _panel_molecular_data(ax_molecular, rows)
    _panel_pedigrees(ax_pedigree, rows)

    bottom_row = [
        _Panel("C", ax_gc, weight=1.35),
        _Panel("D", ax_cov, weight=1.25),
        _Panel("E", ax_molecular, weight=1.25),
        _Panel("F", ax_pedigree, weight=0.85),
    ]
    top_row = [
        _Panel("A", ax_pyr, extra=(ax_ref,)),
        _Panel("B", ax_pc, extra=pc_extras, square=True, aspect=_PC_ASPECT),
    ]
    fig_w, fig_h = fig.get_size_inches()
    for letter, x, y in _layout(fig, [bottom_row, top_row], heights=[1.15, 2.00]):
        fig.text(
            x / fig_w, y / fig_h, letter,
            fontsize=_LETTER_SIZE, fontweight="bold", ha="left", va="bottom", color="#111111",
        )
    return save(fig, "fig1_cohort")


def render_phase1() -> tuple:
    """Alias for the current Figure 1 (Phase 1 context included)."""
    return render(with_phase1=True)
