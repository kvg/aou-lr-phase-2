"""Phenotype prevalence by computed ancestry.

One row per phecodeX category. The cell is the percent of that ancestry
group with at least one code in the category. The count in parentheses is
the number of parent phecodes in the category. Cells with fewer than 20
cases are left blank.
"""

from __future__ import annotations

import gzip
from pathlib import Path

import numpy as np
from matplotlib.colors import LinearSegmentedColormap

from .fig1_cohort import _load_cohort
from .style import ANC_ORDER, ANC_TEXT, FIG_W, MIN_FONT, REPO_ROOT, new_figure, save

_MIN_CASES = 20
_CAT_ORDER = (
    "ID", "CA", "BI", "EM", "MB", "NS", "SO", "CV", "RE",
    "GI", "GU", "DE", "MS", "CM", "SS", "NB", "PP", "GE",
)
_CAT_NAME = {
    "ID": "Infections",
    "CA": "Neoplasms",
    "BI": "Blood/Immune",
    "EM": "Endocrine",
    "MB": "Mental",
    "NS": "Neurological",
    "SO": "Sense organs",
    "CV": "Cardiovascular",
    "RE": "Respiratory",
    "GI": "Gastrointestinal",
    "GU": "Genitourinary",
    "DE": "Dermatologic",
    "MS": "Musculoskeletal",
    "CM": "Congenital",
    "SS": "Symptoms",
    "NB": "Neonatal",
    "PP": "Pregnancy",
    "GE": "Genetic",
}
_PHENO_CANDIDATES = (
    REPO_ROOT / "tractor_mix" / "resources" / "AoU_Phase2_Phenotype.csv.gz",
    REPO_ROOT.parent / "long-read-pipelines" / "AoU_Phase2_Phenotype.csv.gz",
)


def _pheno_path() -> Path:
    path = next((p for p in _PHENO_CANDIDATES if p.is_file()), None)
    if path is None:
        raise RuntimeError("AoU_Phase2_Phenotype.csv.gz was not found")
    return path


def _category_prevalence() -> tuple[list[str], np.ndarray, np.ndarray, dict[str, int]]:
    """Category labels, percent matrix, case counts, ancestry denominators."""
    cohort = _load_cohort()
    if cohort is None:
        raise RuntimeError("covariates were not found")
    anc = {cohort["rid"][i]: cohort["anc"][i] for i in range(cohort["rid"].size)}
    cats = list(_CAT_ORDER)
    cat_index = {c: i for i, c in enumerate(cats)}
    with gzip.open(_pheno_path(), "rt") as fh:
        header = fh.readline().rstrip("\n").split(",")
        pid = header.index("person_id")
        groups: list[list[int]] = [[] for _ in cats]
        n_parent = [0 for _ in cats]
        for col, name in enumerate(header):
            if "_" not in name:
                continue
            pref, body = name.split("_", 1)
            row = cat_index.get(pref)
            if row is None or not body[:1].isdigit():
                continue
            groups[row].append(col)
            if "." not in body:
                n_parent[row] += 1
        den = {a: 0 for a in ANC_ORDER}
        cases = np.zeros((len(cats), len(ANC_ORDER)))
        for line in fh:
            parts = line.rstrip("\n").split(",")
            group = anc.get(parts[pid])
            if group is None:
                continue
            j = ANC_ORDER.index(group)
            den[group] += 1
            for i, cols in enumerate(groups):
                if any(parts[c] == "TRUE" for c in cols):
                    cases[i, j] += 1
    perc = np.zeros_like(cases)
    for j, group in enumerate(ANC_ORDER):
        perc[:, j] = 100.0 * cases[:, j] / den[group]
    labels = [f"{_CAT_NAME[c]} ({n_parent[i]})" for i, c in enumerate(cats)]
    return labels, perc, cases, den


def render() -> tuple[Path, Path]:
    labels, perc, cases, den = _category_prevalence()
    masked = np.ma.masked_where(cases < _MIN_CASES, perc)
    n = len(labels)
    left = 1.70
    bottom = 0.52
    top = 0.08
    cbar_w = 0.16
    cbar_gap = 0.12
    right = 0.15
    row_h = 0.20
    heat_h = n * row_h
    heat_w = FIG_W - left - cbar_gap - cbar_w - right
    fig_h = bottom + heat_h + top
    fig = new_figure(fig_h)
    fig_w, fig_h_in = fig.get_size_inches()
    ax = fig.add_axes([
        left / fig_w,
        bottom / fig_h_in,
        heat_w / fig_w,
        heat_h / fig_h_in,
    ])
    cmap = LinearSegmentedColormap.from_list(
        "pheno", ["#F4F7FA", "#7FA4C4", "#1E3348"],
    )
    cmap.set_bad("#FFFFFF")
    ax.imshow(
        masked, aspect="auto", cmap=cmap, vmin=0, vmax=70,
        interpolation="nearest", origin="upper",
    )
    ax.set_xticks(range(len(ANC_ORDER)))
    ax.set_xticklabels(list(ANC_ORDER), fontsize=MIN_FONT, fontweight="bold")
    for tick, name in zip(ax.get_xticklabels(), ANC_ORDER):
        tick.set_color(ANC_TEXT[name])
    ax.set_yticks(range(n))
    ax.set_yticklabels(labels, fontsize=MIN_FONT, color="black")
    ax.tick_params(length=0, pad=4)
    for sp in ax.spines.values():
        sp.set_visible(False)
    cax = fig.add_axes([
        (left + heat_w + cbar_gap) / fig_w,
        (bottom + heat_h * 0.20) / fig_h_in,
        cbar_w / fig_w,
        (heat_h * 0.60) / fig_h_in,
    ])
    cbar = fig.colorbar(ax.images[0], cax=cax)
    cbar.set_ticks([0, 35, 70])
    cbar.ax.tick_params(labelsize=MIN_FONT, length=2, width=0.4)
    cbar.outline.set_linewidth(0.4)
    cbar.set_label("Percent", fontsize=MIN_FONT)
    for j, name in enumerate(ANC_ORDER):
        ax.annotate(
            f"{den[name]:,}",
            xy=(j, 0.0), xycoords=("data", "axes fraction"),
            xytext=(0, -15), textcoords="offset points",
            ha="center", va="top", fontsize=MIN_FONT, color="#444444",
            annotation_clip=False,
        )
    return save(fig, "fig_phenotypes")
