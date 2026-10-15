"""Testable phecodes as a function of the case-count cutoff (observed data).

Figure 3 panel A counts a phecode as testable when it has at least
``fig3_sim.MIN_CASES`` cases in the analysis cohort (pooled rule) or in a
participant genetic-ancestry group (group rule); FELIX's per-local-ancestry
tests use the same floor in case haplotypes (``MIN_HAP = 2 x MIN_CASES``).
This figure sweeps those cutoffs instead of fixing them, so the headline count
at the chosen rule can be read against what a stricter or looser rule would
give.

Both panels are observed counts from data/fig3_phecode_case_counts.tsv and
data/phecodeX_info.csv -- the same two files panel A reads. Case haplotypes in
panel B are expectations from chr1 FLARE global ancestry proportions (2 x the
sum of a participant's proportion over cases), not counts of genome-wide local
ancestry calls.
"""

from __future__ import annotations

import shutil
from pathlib import Path

import numpy as np

from .fig3_sim import (
    GROUP_ANC,
    LOCAL_ANC,
    MIN_CASES,
    MIN_HAP,
    load_case_counts,
    load_phecodes,
)
from .style import (
    ANC_COLORS,
    ANC_TEXT,
    DATADIR,
    FIG_W,
    MANUSCRIPT_ROOT,
    MIN_FONT,
    new_figure,
    save,
)

_THR = np.arange(1, 1001)  # candidate cutoffs, in cases (panel A) or haplotypes (panel B)
_XLIM = (10, 1000)
_POOLED_COLOR = "#111111"
_ANY_COLOR = "#777777"


def _counts(values: np.ndarray) -> np.ndarray:
    """Number of phecodes at or above each cutoff in _THR."""
    return (values[None, :] >= _THR[:, None]).sum(axis=1)


def _load() -> dict[str, np.ndarray]:
    """Per-phecode case counts and expected case haplotypes, matched to PheCodeX."""
    real = load_case_counts()
    if real is None:
        raise RuntimeError(
            "data/fig3_phecode_case_counts.tsv is required; this figure has no simulated mode"
        )
    rows = [real[r["phecode"]] for r in load_phecodes() if r["phecode"] in real]
    if not rows:
        raise RuntimeError("no PheCodeX phecode matched the case-count table")
    out = {"pooled": np.array([r["n_cases"] for r in rows], dtype=float)}
    for g in GROUP_ANC:
        out[f"group_{g}"] = np.array([r["cases"][g] for r in rows], dtype=float)
    for k in LOCAL_ANC:
        out[f"hap_{k}"] = np.array([r["hap"][k] for r in rows], dtype=float)
    out["any_group"] = np.max(
        np.vstack([out[f"group_{g}"] for g in GROUP_ANC]), axis=0
    )
    out["n_phecodes"] = np.array([len(rows)], dtype=float)
    return out


def _curves(data: dict[str, np.ndarray]) -> dict[str, np.ndarray]:
    return {key: _counts(val) for key, val in data.items() if key != "n_phecodes"}


def _axis(ax, *, xlabel: str, ylabel: str, cutoff: float, ymax: float) -> None:
    ax.set_xscale("log")
    ax.set_xlim(*_XLIM)
    ax.set_ylim(0, ymax)
    ax.set_xticks([10, 20, 50, 100, 200, 500, 1000])
    ax.set_xticklabels(["10", "20", "50", "100", "200", "500", "1,000"])
    ax.minorticks_off()
    ax.set_xlabel(xlabel, labelpad=2)
    ax.set_ylabel(ylabel, labelpad=2)
    ax.axvline(cutoff, color="#999999", lw=0.6, ls=(0, (2.2, 2.2)), zorder=1)


def _mark(ax, x: float, y: float, text: str, color: str, dy: float = 13.0) -> None:
    ax.plot([x], [y], marker="o", ms=3.0, mfc="white", mec=color, mew=0.9, zorder=6)
    ax.annotate(
        text,
        (x, y),
        xytext=(5.0, dy),
        textcoords="offset points",
        fontsize=MIN_FONT,
        color=color,
        ha="left",
        va="center",
        zorder=7,
    )


def _panel_groups(ax, cur: dict[str, np.ndarray]) -> None:
    ymax = 1750.0
    _axis(
        ax,
        xlabel="Minimum cases in the stratum",
        ylabel="Testable phecodes",
        cutoff=MIN_CASES,
        ymax=ymax,
    )
    for g in GROUP_ANC:
        ax.step(_THR, cur[f"group_{g}"], where="post", color=ANC_COLORS[g], lw=0.9, zorder=3)
    ax.step(_THR, cur["any_group"], where="post", color=_ANY_COLOR, lw=1.0, ls=(0, (3, 1.6)), zorder=4)
    ax.step(_THR, cur["pooled"], where="post", color=_POOLED_COLOR, lw=1.5, zorder=5)
    i = int(MIN_CASES) - 1
    _mark(ax, MIN_CASES, cur["pooled"][i], f"{cur['pooled'][i]:,} joint", _POOLED_COLOR)
    _mark(ax, MIN_CASES, cur["any_group"][i], f"{cur['any_group'][i]:,} in any one group", _ANY_COLOR)
    ax.text(
        MIN_CASES,
        ymax * 0.995,
        f"  {int(MIN_CASES)}-case rule",
        fontsize=MIN_FONT,
        color="#666666",
        ha="left",
        va="top",
    )
    handles = [
        ("Joint", _POOLED_COLOR, "-", 1.5),
        ("Any one group", _ANY_COLOR, (0, (3, 1.6)), 1.0),
    ] + [(g, ANC_COLORS[g], "-", 0.9) for g in GROUP_ANC]
    _legend(ax, handles, ncol=2, loc=(0.53, 0.99))


def _panel_local(ax, cur: dict[str, np.ndarray]) -> None:
    ymax = 1750.0
    _axis(
        ax,
        xlabel="Minimum case haplotypes in the local ancestry",
        ylabel="Testable phecodes",
        cutoff=MIN_HAP,
        ymax=ymax,
    )
    for k in LOCAL_ANC:
        ax.step(_THR, cur[f"hap_{k}"], where="post", color=ANC_COLORS[k], lw=1.0, zorder=3)
    i = int(MIN_HAP) - 1
    best = max(LOCAL_ANC, key=lambda k: cur[f"hap_{k}"][i])
    _mark(ax, MIN_HAP, cur[f"hap_{best}"][i], f"{cur[f'hap_{best}'][i]:,} ({best})", ANC_TEXT[best])
    ax.text(
        MIN_HAP,
        ymax * 0.995,
        f"  {int(MIN_HAP)}-haplotype rule",
        fontsize=MIN_FONT,
        color="#666666",
        ha="left",
        va="top",
    )
    _legend(ax, [(k, ANC_COLORS[k], "-", 1.0) for k in LOCAL_ANC], ncol=1, loc=(0.68, 0.92))


def _legend(ax, entries, *, ncol: int, loc) -> None:
    from matplotlib.lines import Line2D

    handles = [
        Line2D([], [], color=color, ls=ls, lw=lw, label=label) for label, color, ls, lw in entries
    ]
    ax.legend(
        handles=handles,
        loc="upper left",
        bbox_to_anchor=loc,
        ncol=ncol,
        fontsize=MIN_FONT,
        frameon=False,
        handlelength=1.4,
        handletextpad=0.4,
        columnspacing=0.9,
        labelspacing=0.28,
        borderaxespad=0.0,
    )


def write_table(cur: dict[str, np.ndarray], data: dict[str, np.ndarray]) -> Path:
    """Counts at every integer cutoff, for the caption numbers and the text."""
    DATADIR.mkdir(parents=True, exist_ok=True)
    path = DATADIR / "figS_phecode_cutoff.tsv"
    cols = ["pooled", "any_group"] + [f"group_{g}" for g in GROUP_ANC] + [f"hap_{k}" for k in LOCAL_ANC]
    with open(path, "w") as fh:
        fh.write("cutoff\t" + "\t".join(cols) + "\n")
        for j, t in enumerate(_THR):
            fh.write(f"{t}\t" + "\t".join(str(int(cur[c][j])) for c in cols) + "\n")
    return path


def render() -> tuple[Path, Path]:
    data = _load()
    cur = _curves(data)
    write_table(cur, data)
    fig = new_figure(3.05, FIG_W)
    ax_a = fig.add_axes([0.070, 0.165, 0.395, 0.760])
    ax_b = fig.add_axes([0.580, 0.165, 0.395, 0.760])
    _panel_groups(ax_a, cur)
    _panel_local(ax_b, cur)
    for ax, letter in ((ax_a, "A"), (ax_b, "B")):
        ax.text(
            -0.115,
            1.035,
            letter,
            transform=ax.transAxes,
            fontsize=9,
            fontweight="bold",
            ha="left",
            va="bottom",
        )
    pdf, png = save(fig, "figS_phecode_cutoff")
    # Observed data: keep the tracked figures/ pair in step with every re-render.
    dest = MANUSCRIPT_ROOT / "figures"
    dest.mkdir(parents=True, exist_ok=True)
    out = [shutil.copyfile(src, dest / src.name) for src in (pdf, png)]
    return Path(out[0]), Path(out[1])


if __name__ == "__main__":  # pragma: no cover
    print(render())
