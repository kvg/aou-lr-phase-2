"""Figure 1 + Figure 2 merged onto one plate.

Reviewer feedback was that two display items for "what we sequenced" and
"what we found" cost more page than they return. This module composes the
surviving panels of :mod:`fig1_null` and :mod:`fig2_options` into a single
Science-width plate.

Dropped from the merge (both move to the supplement):

* Fig. 1B, the long-read PC scatter. The Fig. 2B global-ancestry strip
  already shows cohort structure, and PC1-10 as null-model covariates is
  conventional enough to relegate.
* Fig. 1F, the pedigree size histogram. No pedigree analysis is carried in
  the main text.

``center`` selects how sequencing center is displayed. Center is a null-model
covariate, so it has a place in the figure, but the counts-ordered bar reads
as a league table of how much each center contributed:

``named``
    The Fig. 1C bar unchanged: centers ordered by participant count.
``omit``
    No center panel. Center stays a covariate, declared in Methods.
``shares``
    One equal-length bar per center, alphabetical, showing the composition
    of what that center produced rather than how much. Counts live in the
    supplement.
``coverage``
    No separate center panel. The coverage panel grows a rail underneath
    carrying each center's median and interquartile coverage, which is the
    batch structure that justifies the covariate in the first place.

Layouts:

``band``
    One cohort row across the top, then the length spectrum, then discovery
    beside the repeat panels.
``tall``
    Same order, but the cohort row is given more height and the allelic
    series is cut to three loci to pay for it.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
from matplotlib import patheffects
from matplotlib.legend_handler import HandlerTuple
from matplotlib.lines import Line2D
from matplotlib.patches import Patch, PathPatch, Rectangle
from matplotlib.path import Path as MplPath
from matplotlib.ticker import MultipleLocator

from . import fig1_null as f1
from . import fig2_options as f2
from .style import MIN_FONT, despine, new_figure, save

_CAT_GAP = 0.50  # default inches above the spectrum axes for its legend and facet titles
_CAT_FOOT = 0.76  # inches below the catalog axes for x labels
_CENTER_ALPHA = ("BCM", "BI", "HA", "JHU", "UW")
_CENTER_NAME = dict(f1._CENTERS)


# ---------------------------------------------------------------------------
# Row placement (fig1_null._layout, but for one row at an arbitrary floor)
# ---------------------------------------------------------------------------


def _place_row(fig, panels, height: float, floor: float, inline_letters: bool = False):
    """Tile one row of :class:`fig1_null._Panel` so ink is separated by gutters.

    ``inline_letters`` drops the letters alongside the top of each panel's ink
    instead of on a band above it. Every legend in the cohort row is anchored
    to its axes, so the letter column at the row's left margin stays clear of
    them and the row reclaims a letter's height.
    """
    fig_w, fig_h = fig.get_size_inches()
    dpi = fig.dpi
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    pads = [p.pads(renderer, dpi, fig_w, fig_h) for p in panels]
    fixed = sum(p.fixed_width(height) for p in panels if p.square)
    flexible = (
        fig_w - 2 * f1._MARGIN
        - f1._GUTTER * (len(panels) - 1)
        - sum(left + right for left, right, _b, _t in pads)
        - fixed
    )
    total_weight = sum(p.weight for p in panels if not p.square) or 1.0
    y = floor + max(bottom for _l, _r, bottom, _t in pads)
    row_top = y + height + max(top for _l, _r, _b, top in pads)
    cursor = f1._MARGIN
    letters = []
    for panel, (left, right, _b, _t) in zip(panels, pads):
        width = (
            panel.fixed_width(height)
            if panel.square
            else flexible * panel.weight / total_weight
        )
        panel.ax.set_position([
            (cursor + left) / fig_w, y / fig_h, width / fig_w, height / fig_h,
        ])
        letter_y = (
            row_top - f1._LETTER_SIZE / 72 if inline_letters
            else row_top + f1._LETTER_GAP
        )
        letters.append((panel.letter, cursor, letter_y))
        cursor += left + width + right + f1._GUTTER
    return letters, row_top


def _settle_row(fig, panels, height: float, top_in: float, inline_letters: bool = False):
    """Place a row so its ink (letters included) ends ``top_in`` below the top edge.

    Returns the letter anchors and the row's floor in inches.
    """
    fig_h = fig.get_size_inches()[1]
    letter_h = f1._LETTER_SIZE / 72
    block = height + 1.0
    letters: list = []
    floor = 0.0
    for _ in range(4):
        floor = fig_h - top_in - block
        letters, row_top = _place_row(fig, panels, height, floor, inline_letters)
        block = (row_top if inline_letters else row_top + f1._LETTER_GAP + letter_h) - floor
    return letters, floor


# ---------------------------------------------------------------------------
# Sequencing-center treatments
# ---------------------------------------------------------------------------


def _center_segments(rows, ont_only):
    """Per center: Phase 1, mid-pass, PacBio-only high-pass, PacBio+ONT, ONT-only."""
    high_mask = f1._high_pass_mask(rows)
    seg = {}
    for code in _CENTER_ALPHA:
        n = sum(row.get("GC") == code for row in rows)
        high = sum(m and row.get("GC") == code for m, row in zip(high_mask, rows))
        ont = sum(f1._truthy(row.get("has_ONT")) and row.get("GC") == code for row in rows)
        seg[code] = dict(
            phase1=0, mid=n - high, high_only=high - ont, ont=ont,
            ont_only=int(ont_only[code]),
        )
    return seg


def _panel_center_shares(ax, rows, ont_only, phase1_n=None, phase1_center="HA") -> None:
    """Composition of each center's output, every bar the same length.

    Removes the counts-ordered ranking while keeping the reason center is a
    covariate visible: the centers ran different parts of the design.
    """
    seg = _center_segments(rows, ont_only)
    if phase1_n:
        seg[phase1_center]["phase1"] = phase1_n
    order = list(_CENTER_ALPHA)
    y = np.arange(len(order))
    keys = ("phase1", "mid", "high_only", "ont", "ont_only")
    labels = ("Phase 1", "Mid-pass", "High-pass", "ONT overlap", "ONT only")
    colors = (f1._PHASE1, f1._MID_PASS_COLOR, f1._HIGH_PASS_COLOR, f1._HIGH_PASS_COLOR, f1._ONT_ONLY_COLOR)
    totals = np.array([sum(seg[c][k] for k in keys) for c in order], dtype=float)
    left = np.zeros(len(order))
    for key, label, color in zip(keys, labels, colors):
        vals = np.array([seg[c][key] for c in order], dtype=float) / totals * 100
        kw = dict(edgecolor=f1._ONT_COLOR, linewidth=0.55, hatch="////") if key == "ont" else dict(linewidth=0)
        ax.barh(y, vals, left=left, color=color, height=0.72, label=label, **kw)
        left += vals
    ax.set_yticks(y)
    ax.set_yticklabels([_CENTER_NAME[c] for c in order])
    ax.invert_yaxis()
    ax.set_xlim(0, 100)
    ax.set_xticks([0, 50, 100])
    ax.set_xlabel("Share of center's participants (%)")
    despine(ax)
    handles, legend_labels = ax.get_legend_handles_labels()
    blank = Patch(facecolor="none", edgecolor="none", linewidth=0)
    ax.legend(
        [handles[0], handles[1], handles[2], blank, handles[3], handles[4]],
        [legend_labels[0], legend_labels[1], legend_labels[2], " ", legend_labels[3], legend_labels[4]],
        loc="lower left", bbox_to_anchor=(0.0, f1._ROW_LEGEND_Y), ncol=2,
        columnspacing=0.7, **f1._ROW_LEGEND,
    )


def _align_letter_to_legend(fig, text, label: str) -> None:
    """Drop ``text`` onto the baseline of the figure legend entry named ``label``."""
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    for legend in fig.legends:
        entries = legend.get_texts()
        if not any(t.get_text() == label for t in entries):
            continue
        y0 = entries[0].get_window_extent(renderer).y0
        text.set_position((text.get_position()[0], y0 / (fig.dpi * fig.get_figheight())))
        return


LAST_GEOMETRY: dict[str, float] = {}
# Genetic-similarity descriptors, following the NASEM population-descriptor
# report: a group is named for the reference panel its members resemble, not
# as a race or ethnicity. OTH is the residual, so it is named as such.
ANC_LIKE = {
    "AFR": "AFR-like", "AMR": "AMR-like", "EAS": "EAS-like", "EUR": "EUR-like",
    "MID": "MID-like", "SAS": "SAS-like", "OTH": "Unassigned",
}
_MIN_DISCLOSED = 20


def _panel_molecular_ont(ax, rows, ont_only_n: int) -> None:
    """UpSet of assayed modalities with ONT as a fourth column.

    ONT is the platform used for validation in this paper rather than an
    assay, so its dots are drawn in the ONT color. Participants sequenced
    only with ONT are not in the PacBio covariate table, so they enter as
    their own row; the panel therefore covers more participants than the
    12,261 PacBio cohort.

    ``In QTL analysis`` counts only participants who have Olink or RNA-seq
    and are in the matching QTL set. Fig. 1E instead treats a missing assay
    as satisfied, which marks 5mC-only participants as in-analysis.
    """
    t = f1._truthy
    avail: dict[tuple[bool, ...], int] = {}
    featured: dict[tuple[bool, ...], int] = {}
    for row in rows:
        key = (t(row.get("has_methylation")), t(row.get("has_proteomics")),
               t(row.get("has_rna")), t(row.get("has_ONT")))
        if not any(key):
            continue
        avail[key] = avail.get(key, 0) + 1
        if key[1] or key[2]:
            pqtl_ok = (not key[1]) or t(row.get("in_pqtl"))
            rna_ok = (not key[2]) or t(row.get("in_eqtl")) or t(row.get("in_sqtl"))
            if pqtl_ok and rna_ok:
                featured[key] = featured.get(key, 0) + 1
    ont_only_key = (False, False, False, True)
    avail[ont_only_key] = avail.get(ont_only_key, 0) + ont_only_n
    combos = sorted(avail, key=lambda key: -avail[key])
    hidden = [k for k in combos if avail[k] < _MIN_DISCLOSED]
    combos = [k for k in combos if avail[k] >= _MIN_DISCLOSED]
    feat_vals = [featured.get(key, 0) for key in combos]
    other_vals = [avail[key] - featured.get(key, 0) for key in combos]
    totals = [avail[key] for key in combos]
    y = np.arange(len(combos))

    ax.set_axis_off()
    ax_mat = ax.inset_axes([0.0, 0.0, 0.32, 1.0])
    ax_bar = ax.inset_axes([0.38, 0.0, 0.62, 1.0], sharey=ax_mat)
    ax_bar.barh(y, feat_vals, height=0.66, color=f1._MULTIOMIC_COLOR, linewidth=0,
                label="In QTL analysis")
    ax_bar.barh(y, other_vals, left=feat_vals, height=0.66, color=f1._MULTIOMIC_PALE,
                linewidth=0, label="Assayed only")
    ax_bar.set_xlabel("Participants")
    ax_bar.set_ylim(len(combos) - 0.55, -0.55)
    ax_bar.set_yticks([])
    ax_bar.set_xlim(0, max(totals) * 1.30)
    ax_bar.xaxis.set_major_formatter(f1._THOUSANDS)
    ax_bar.set_xticks([0, 5000])
    for yi, total, feat, other in zip(y, totals, feat_vals, other_vals):
        label = f"{feat:,}+{other:,}" if feat and other else f"{total:,}"
        ax_bar.annotate(label, xy=(total, yi), xytext=(3, 0), textcoords="offset points",
                        ha="left", va="center", fontsize=MIN_FONT)
    for side in ("top", "right", "left"):
        ax_bar.spines[side].set_visible(False)
    ax_bar.legend(loc="lower right", bbox_to_anchor=(1.0, f1._ROW_LEGEND_Y), ncol=1,
                  **f1._ROW_LEGEND)
    if hidden:
        ax_bar.annotate(
            f"{len(hidden)} further combination{'s' if len(hidden) > 1 else ''} <20",
            xy=(1.0, -0.33), xycoords="axes fraction",
            ha="right", va="top", fontsize=MIN_FONT, color="#8A8A8A",
        )

    names = ["5mC", "Olink", "RNA-seq", "ONT"]
    on_color = [f1._MULTIOMIC_COLOR] * 3 + [f1._ONT_COLOR]
    for yi, combo in zip(y, combos):
        present = [i for i, value in enumerate(combo) if value]
        if len(present) > 1:
            ax_mat.plot([min(present), max(present)], [yi, yi], color=f1._MULTIOMIC_COLOR, lw=1.0)
        for xi, value in enumerate(combo):
            ax_mat.scatter(xi, yi, s=15, color=on_color[xi] if value else "#D9DEE2",
                           edgecolors="none", zorder=2)
    ax_mat.set_xticks(range(len(names)))
    ax_mat.set_xticklabels(names, rotation=90, ha="center", va="top")
    for tick, color in zip(ax_mat.get_xticklabels(), on_color):
        tick.set_color(color if color == f1._ONT_COLOR else "#111111")
    ax_mat.set_xlim(-0.5, len(names) - 0.5)
    ax_mat.set_yticks([])
    ax_mat.tick_params(axis="x", length=0, pad=2)
    for spine in ax_mat.spines.values():
        spine.set_visible(False)


_FREEZE_LABELS = ("In companion analyses", "Release only")
_QTL_LABELS = ("In QTL analysis", "Assayed only")


def _load_ont_only_rows() -> list[dict[str, str]]:
    """Releasable Phase 2 participants sequenced with ONT and not PacBio.

    They are absent from ``_load_model_rows``, so every panel built on that
    covers 12,261 of the 13,530 releasable Phase 2 participants.
    """
    import csv
    import gzip

    from .style import covariates_path

    out = []
    with gzip.open(covariates_path(), "rt", newline="") as handle:
        for row in csv.DictReader(handle):
            if row.get("technology") != "ONT":
                continue
            if str(row.get("final_releasable_v9", "")).strip() in {"True", "true"}:
                out.append(row)
    return out


def _add_ont_only_outline(ax, pyr, ont_rows) -> tuple:
    """Draw the ONT-only cohort as a second outline on the age-sex pyramid.

    Same treatment the panel already gives Phase 1, so it adds a cohort to
    the panel that carries cohorts rather than new furniture to a panel that
    does not.
    """
    from collections import Counter

    ages = np.arange(20, 91)
    layers = {"Female": Counter(), "Male": Counter()}
    for row in ont_rows:
        side, _source = pyr.placement(row)
        if side not in {"Female", "Male"}:
            continue
        layers[side][min(int(float(row["age"])), 90)] += 1
    halo = [patheffects.withStroke(linewidth=2.6, foreground="white")]
    for side, sign in (("Male", 1.0), ("Female", -1.0)):
        values = np.array([layers[side][int(a)] for a in ages], dtype=float)
        ax.step(ages, sign * values, where="mid", color=f1._ONT_COLOR, lw=1.1,
                zorder=5, path_effects=halo)
    return Line2D([0], [0], color=f1._ONT_COLOR, lw=1.1)


def _ont_counts(rows) -> tuple[int, int]:
    """ONT participants, split by whether they also have PacBio.

    The overlap pair carry ``has_ONT`` on their PacBio row and so are already
    counted in the UpSet rows; the ONT-only pair are a separate technology in
    the covariate table and are in no UpSet row at all.
    """
    import csv
    import gzip

    from .style import covariates_path

    with_pacbio = sum(f1._truthy(r.get("has_ONT")) for r in rows)
    ont_only = 0
    with gzip.open(covariates_path(), "rt", newline="") as handle:
        for row in csv.DictReader(handle):
            if row.get("technology") != "ONT":
                continue
            if str(row.get("final_releasable_v9", "")).strip() in {"True", "true"}:
                ont_only += 1
    return with_pacbio, ont_only


def _draw_ont_row(ax_mat, ax_bar, rows, y: float) -> None:
    """One bar for ONT, outside the UpSet's partition of the PacBio cohort.

    Deliberately not a fourth matrix column: the 876 participants with both
    platforms are already inside the rows above, so an ONT column would
    double-count them, and splitting every row by ONT is the seven-row
    version that was too dense to read. What this says instead is that ONT
    data exists, that some of it sits on PacBio participants and some does
    not.
    """
    with_pacbio, ont_only = _ont_counts(rows)
    ax_bar.barh(y, with_pacbio, height=0.68, color=f1._ONT_COLOR, linewidth=0)
    ax_bar.barh(y, ont_only, left=with_pacbio, height=0.68,
                color=f1._ONT_ONLY_COLOR, linewidth=0)
    # Two-tone count, as the Fig. 1 center bars label PacBio + ONT.
    anchor = ax_bar.annotate(
        f"{with_pacbio:,}", xy=(with_pacbio + ont_only, y), xytext=(3, 0),
        textcoords="offset points", ha="left", va="center",
        fontsize=MIN_FONT, color=f1._ONT_COLOR,
    )
    ax_bar.annotate(
        f" + {ont_only:,}", xy=(1, 0.5), xycoords=anchor, ha="left", va="center",
        fontsize=MIN_FONT, color="#5C7A88",
    )
    ax_mat.annotate(
        "ONT", xy=(1.0, y), xycoords=("axes fraction", "data"), xytext=(-1, 0),
        textcoords="offset points", ha="right", va="center",
        fontsize=MIN_FONT, color=f1._ONT_COLOR,
    )


def _panel_molecular_any_qtl(ax, rows, split: bool = True, rule: str = "freeze",
                             ont_row: bool = False) -> None:
    """Assayed modalities, as an UpSet, split by what the companion analyses used.

    The release and the analysis freeze are not the same cohort: data that
    arrived after a companion study's freeze is in the release but not in
    that study. ``rule="freeze"`` shades a participant dark when every assay
    they carry was in its analysis — the rule
    ``fig1_null._panel_molecular_data`` uses. For a participant without an
    assay the corresponding test is vacuously true, which is correct here
    (nothing of theirs was left out) but makes 5mC-only participants dark on
    the assumption that methylation has no freeze gap; there is no meQTL
    sample list to check that against.

    ``rule="any"`` instead asks only whether the participant appears in at
    least one of the pQTL, eQTL and sQTL lists, which reads as QTL membership
    rather than as freeze coverage. Without ``split`` the bars are
    availability only.
    """
    t = f1._truthy
    avail: dict[tuple[bool, bool, bool], int] = {}
    featured: dict[tuple[bool, bool, bool], int] = {}
    for row in rows:
        key = (t(row.get("has_methylation")), t(row.get("has_proteomics")), t(row.get("has_rna")))
        if not any(key):
            continue
        avail[key] = avail.get(key, 0) + 1
        rna_ok = (not key[2]) or t(row.get("in_eqtl")) or t(row.get("in_sqtl"))
        if rule == "freeze":
            used = ((not key[1]) or t(row.get("in_pqtl"))) and rna_ok
        else:
            used = t(row.get("in_pqtl")) or t(row.get("in_eqtl")) or t(row.get("in_sqtl"))
        if used:
            featured[key] = featured.get(key, 0) + 1
    combos = sorted(avail, key=lambda key: -avail[key])
    feat_vals = [featured.get(key, 0) for key in combos]
    other_vals = [avail[key] - featured.get(key, 0) for key in combos]
    totals = [avail[key] for key in combos]
    y = np.arange(len(combos))

    ax.set_axis_off()
    ax_mat = ax.inset_axes([0.0, 0.0, 0.34, 1.0])
    ax_bar = ax.inset_axes([0.40, 0.0, 0.60, 1.0], sharey=ax_mat)
    used_label, rest_label = _FREEZE_LABELS if rule == "freeze" else _QTL_LABELS
    if split:
        ax_bar.barh(y, feat_vals, height=0.68, color=f1._MULTIOMIC_COLOR, linewidth=0,
                    label=used_label)
        ax_bar.barh(y, other_vals, left=feat_vals, height=0.68, color=f1._MULTIOMIC_PALE,
                    linewidth=0, label=rest_label)
    else:
        ax_bar.barh(y, totals, height=0.68, color=f1._MULTIOMIC_COLOR, linewidth=0)
    ax_bar.set_xlabel("Participants")
    ont_y = len(combos) - 1 + 1.35
    ax_bar.set_ylim((ont_y + 0.55) if ont_row else (len(combos) - 0.55), -0.55)
    ax_bar.set_yticks([])
    ax_bar.set_xlim(0, max(totals) * 1.28)
    ax_bar.xaxis.set_major_formatter(f1._THOUSANDS)
    for yi, total, feat, other in zip(y, totals, feat_vals, other_vals):
        label = f"{feat:,}+{other:,}" if split and feat and other else f"{total:,}"
        ax_bar.annotate(label, xy=(total, yi), xytext=(3, 0), textcoords="offset points",
                        ha="left", va="center", fontsize=MIN_FONT)
    for side in ("top", "right", "left"):
        ax_bar.spines[side].set_visible(False)
    if split:
        ax_bar.legend(loc="lower right", bbox_to_anchor=(1.0, f1._ROW_LEGEND_Y), ncol=1,
                      **f1._ROW_LEGEND)

    names = ["5mC", "Olink", "RNA-seq"]
    for yi, combo in zip(y, combos):
        present = [i for i, value in enumerate(combo) if value]
        if len(present) > 1:
            ax_mat.plot([min(present), max(present)], [yi, yi], color=f1._MULTIOMIC_COLOR, lw=1.0)
        for xi, value in enumerate(combo):
            ax_mat.scatter(xi, yi, s=18, color=f1._MULTIOMIC_COLOR if value else "#D9DEE2",
                           edgecolors="none", zorder=2)
    ax_mat.set_xticks(range(len(names)))
    ax_mat.set_xticklabels(names, rotation=90, ha="center", va="top")
    ax_mat.set_xlim(-0.5, len(names) - 0.5)
    ax_mat.set_yticks([])
    ax_mat.tick_params(axis="x", length=0, pad=2)
    for spine in ax_mat.spines.values():
        spine.set_visible(False)
    if ont_row:
        # A rule marks where the partition of the PacBio cohort stops.
        ax_bar.axhline(ont_y - 0.68, color="#D0D5D8", lw=0.6, zorder=0)
        _draw_ont_row(ax_mat, ax_bar, rows, ont_y)


def _panel_cohort_bracket(ax, rows) -> None:
    """Cohort bar, bracket, then the assay UpSet for the bracketed part.

    Reading order is the panel's argument: every Phase 2 participant is on
    the top bar, split by platform; the bracket picks out the ones with
    PacBio; the UpSet below breaks that subset down by assay. The bar sits
    above the UpSet's legend because it is the premise, not one more series
    in it.
    """
    t = f1._truthy
    avail: dict[tuple[bool, bool, bool], int] = {}
    featured: dict[tuple[bool, bool, bool], int] = {}
    for row in rows:
        key = (t(row.get("has_methylation")), t(row.get("has_proteomics")), t(row.get("has_rna")))
        if not any(key):
            continue
        avail[key] = avail.get(key, 0) + 1
        rna_ok = (not key[2]) or t(row.get("in_eqtl")) or t(row.get("in_sqtl"))
        if ((not key[1]) or t(row.get("in_pqtl"))) and rna_ok:
            featured[key] = featured.get(key, 0) + 1
    combos = sorted(avail, key=lambda k: -avail[k])
    feat = [featured.get(k, 0) for k in combos]
    rest = [avail[k] - featured.get(k, 0) for k in combos]
    totals = [avail[k] for k in combos]
    y = np.arange(len(combos))

    with_pb, ont_only = _ont_counts(rows)
    pb_only = len(rows) - with_pb
    segs = [("PacBio only", pb_only, f1._MID_PASS_COLOR, "#6F8697"),
            ("PacBio + ONT", with_pb, f1._HIGH_PASS_COLOR, f1._HIGH_PASS_COLOR),
            ("ONT only", ont_only, f1._ONT_ONLY_COLOR, "#5C8FA8")]
    cohort = pb_only + with_pb + ont_only

    ax.set_axis_off()
    cohort = pb_only + with_pb + ont_only
    pb_end, ont_start = pb_only + with_pb, pb_only
    # The cohort bar sits in the UpSet's bar column and shares its x scale, so
    # the two read against one ruler. Its axis runs past the cohort to a round
    # 15,000; the UpSet's spine below is truncated to the range its bars use.
    xmax = 15_600
    ax_top = ax.inset_axes([0.40, 0.60, 0.60, 0.40])
    ax_top.set_xlim(0, xmax)
    ax_top.set_ylim(-5.6, 4.2)
    for x0, x1, color in ((0, pb_only, f1._MID_PASS_COLOR),
                          (pb_only, pb_end, f1._HIGH_PASS_COLOR),
                          (pb_end, cohort, f1._ONT_ONLY_COLOR)):
        ax_top.barh(0, x1 - x0, left=x0, height=1.5, color=color, linewidth=0)

    def brace(x0, x1, y, tip, text, va):
        ax_top.plot([x0, x0, x1, x1], [y + tip, y, y, y + tip], color="#5A6B78",
                    lw=0.7, solid_joinstyle="miter")
        ax_top.annotate(text, xy=((x0 + x1) / 2, y), xytext=(0, 2 if va == "bottom" else -2),
                        textcoords="offset points", ha="center", va=va,
                        fontsize=MIN_FONT, color="#44525C")

    brace(ont_start, cohort, 1.5, -0.6, "ONT", "bottom")
    brace(0, pb_end, -1.5, 0.6, "PacBio", "top")
    ax_top.set_yticks([])
    ax_top.set_xticks([0, 5000, 10000, 15000])
    ax_top.xaxis.set_major_formatter(f1._THOUSANDS)
    ax_top.tick_params(axis="x", length=2, pad=1.5)
    for side in ("top", "right", "left"):
        ax_top.spines[side].set_visible(False)
    ax_top.spines["bottom"].set_position(("data", -3.8))

    ax_mat = ax.inset_axes([0.0, 0.0, 0.34, 0.46])
    ax_bar = ax.inset_axes([0.40, 0.0, 0.60, 0.46], sharey=ax_mat)
    ax_bar.barh(y, feat, height=0.68, color=f1._MULTIOMIC_COLOR, linewidth=0,
)
    ax_bar.barh(y, rest, left=feat, height=0.68, color=f1._MULTIOMIC_PALE,
                linewidth=0, label=_FREEZE_LABELS[1])
    for yi, total, f_, r_ in zip(y, totals, feat, rest):
        label = f"{f_:,}+{r_:,}" if f_ and r_ else f"{total:,}"
        ax_bar.annotate(label, xy=(total, yi), xytext=(3, 0), textcoords="offset points",
                        ha="left", va="center", fontsize=MIN_FONT)
    ax_bar.set_xlabel("Participants")
    ax_bar.set_ylim(len(combos) - 0.55, -0.55)
    ax_bar.set_yticks([])
    ax_bar.set_xlim(0, xmax)
    ax_bar.xaxis.set_major_formatter(f1._THOUSANDS)
    ax_bar.set_xticks([0, 5000])
    for side in ("top", "right", "left"):
        ax_bar.spines[side].set_visible(False)
    ax_bar.spines["bottom"].set_bounds(0, 5000)
    ax_bar.xaxis.set_label_coords(2500, -0.42, transform=ax_bar.get_xaxis_transform())

    names = ["5mC", "Olink", "RNA-seq"]
    for yi, combo in zip(y, combos):
        present = [i for i, v in enumerate(combo) if v]
        if len(present) > 1:
            ax_mat.plot([min(present), max(present)], [yi, yi], color=f1._MULTIOMIC_COLOR, lw=1.0)
        for xi, v in enumerate(combo):
            ax_mat.scatter(xi, yi, s=18, color=f1._MULTIOMIC_COLOR if v else "#D9DEE2",
                           edgecolors="none", zorder=2)
    ax_mat.set_xticks(range(len(names)))
    ax_mat.set_xticklabels(names, rotation=90, ha="center", va="top")
    ax_mat.set_xlim(-0.5, len(names) - 0.5)
    ax_mat.set_yticks([])
    ax_mat.tick_params(axis="x", length=0, pad=2)
    for spine in ax_mat.spines.values():
        spine.set_visible(False)


_RULE_COLOR = "#5A6B78"


def _square_bracket(ax, x0, x1, y_line, tip, label, *, color=_RULE_COLOR):
    """A square bracket spanning ``x0``-``x1`` with tips of length ``tip``.

    ``tip`` is signed in data units, so on the inverted assay axis a positive
    value points the tips down toward the bar being annotated.
    """
    ax.plot([x0, x0, x1, x1], [y_line + tip, y_line, y_line, y_line + tip],
            color=color, lw=0.7, solid_joinstyle="miter", clip_on=False, zorder=3)
    return 0.5 * (x0 + x1), label


def _curly_brace(ax, x0, x1, y_flat, y_tip, *, color=_RULE_COLOR):
    """A curly brace: flat ends at ``y_flat``, meeting in a tip at ``y_tip``.

    Two cubic arms whose control points sit half way between the two levels,
    which gives the usual brace silhouette -- a short turn at each end, a long
    flat run, a turn into the point -- at any aspect ratio.
    """
    xm = 0.5 * (x0 + x1)
    depth = y_tip - y_flat
    r = 0.10 * (xm - x0)
    verts: list[tuple[float, float]] = []
    codes: list[int] = []
    for sign, xe in ((1.0, x0), (-1.0, x1)):
        # Up-turned end, rounded corner, straight run, rounded turn into the
        # point. Quadratic segments with the corner as the control point.
        verts += [(xe, y_flat - 0.40 * depth), (xe, y_flat), (xe + sign * r, y_flat),
                  (xm - sign * r, y_flat), (xm, y_flat), (xm, y_tip)]
        codes += [MplPath.MOVETO, MplPath.CURVE3, MplPath.CURVE3,
                  MplPath.LINETO, MplPath.CURVE3, MplPath.CURVE3]
    ax.add_patch(PathPatch(MplPath(verts, codes), fill=False, lw=0.7,
                           edgecolor=color, clip_on=False, zorder=3))


def _panel_bar_over_upset(ax, rows, *, brackets: bool = True,
                          connector: str = "brace", counts: str = "under_bar",
                          legend_lift: int = 1, tip_label: str | None = None) -> None:
    """A platform composition bar over the assay UpSet for the PacBio cohort.

    The top bar is the whole releasable Phase 2 cohort split by platform;
    the UpSet below it is the assay detail for the PacBio part only, so the
    two elements answer separate questions on one Participants axis. The
    876 PacBio+ONT and 1,269 ONT-only segments are a few per cent of 13,530
    and far too narrow to label in place, so their counts run as one
    colour-coded line and the caption keys the colours.

    ``connector`` is how the PacBio span of the bar is tied to the UpSet that
    partitions it -- ``brace`` (curly), ``bracket`` (square, with a stem into
    the UpSet), ``tint`` (a pale band behind the gap) or ``none``.
    ``counts`` is where the segment counts go: ``under_bar``,
    ``under_connector``, ``in_brackets`` (folded into two-line platform
    labels, which drops the 876 overlap to the caption) or ``none``.

    Whatever sits directly under the bar is what the connector appears to
    attach to, so the two placements are not interchangeable.
    """
    if connector not in {"brace", "bracket", "tint", "none"}:
        raise ValueError(f"unknown connector {connector!r}")
    if counts not in {"under_bar", "under_connector", "in_brackets",
                      "study", "none"}:
        raise ValueError(f"unknown counts placement {counts!r}")
    t = f1._truthy
    avail: dict[tuple[bool, bool, bool], int] = {}
    featured: dict[tuple[bool, bool, bool], int] = {}
    for row in rows:
        key = (t(row.get("has_methylation")), t(row.get("has_proteomics")), t(row.get("has_rna")))
        if not any(key):
            continue
        avail[key] = avail.get(key, 0) + 1
        rna_ok = (not key[2]) or t(row.get("in_eqtl")) or t(row.get("in_sqtl"))
        if ((not key[1]) or t(row.get("in_pqtl"))) and rna_ok:
            featured[key] = featured.get(key, 0) + 1
    combos = sorted(avail, key=lambda k: -avail[k])
    feat = [featured.get(k, 0) for k in combos]
    rest = [avail[k] - featured.get(k, 0) for k in combos]
    totals = [avail[k] for k in combos]
    y = np.arange(len(combos))

    with_pb, ont_only = _ont_counts(rows)
    segs = [("PacBio only", len(rows) - with_pb, f1._MID_PASS_COLOR),
            ("PacBio + ONT", with_pb, f1._HIGH_PASS_COLOR),
            ("ONT only", ont_only, f1._ONT_ONLY_COLOR)]
    ont_start, pb_end = segs[0][1], segs[0][1] + segs[1][1]
    total = pb_end + segs[2][1]

    # Stack the top element upward from the first UpSet row. TXT is one line
    # of MIN_FONT in row units, which is what every gap in here is made of.
    TXT, GAP, CONN_H = 0.78, 0.12, 0.66
    BAR_H, TIP, DROP = 0.62, 0.26, 0.78
    BAR_GAP = DROP - TIP - BAR_H / 2   # bar edge to bracket tip, above and below
    y_counts = None
    if not brackets:
        connector, counts = "none", "under_bar"
        y_bar, y_counts, y_top = -1.75, -1.13, -2.20
    elif counts == "study":
        # Three statements rather than a platform taxonomy: how big the
        # long-read release is, how much of it is dual-platform, and which
        # part this paper analyses. Everything hangs off the bar, placed so
        # the study bracket, its label and its stem clear the UpSet.
        y_bar = (-0.34 - 0.25 - 0.63 - 0.06 - TXT - 0.18 - TIP - BAR_GAP
                 - BAR_H / 2)
        y_lo = y_bar + BAR_H / 2 + BAR_GAP + TIP
        y_hi = y_bar - BAR_H / 2 - BAR_GAP - TIP
        y_study = y_lo + 0.18
        y_both = y_hi - 0.45
        y_head = y_both - TXT - 0.10
        y_top = y_head - (0.0 if legend_lift else TXT) - 0.10
    else:
        cur = -0.34 - (0.62 if tip_label else 0.52)
        if counts == "under_connector":
            cur -= TXT
            y_counts = cur
            cur -= GAP
        if connector in {"brace", "bracket"}:
            y_tip = cur
            cur -= CONN_H
            y_flat = cur
            # Match the clearance the platform brackets leave above the bar.
            # The brace's ends turn up before the run, so its gap has to be
            # deeper by that curl to look the same.
            cur -= BAR_GAP + (0.40 * CONN_H if connector == "brace" else 0.0)
        if counts == "under_bar":
            cur -= TXT
            y_counts = cur
            cur -= 0.04
        y_bar = cur - 0.31
        y_pb = y_bar - DROP
        # The overlap count has no home in this variant: every spot in the gap
        # between the bracket levels is crossed by one of the two brackets'
        # tips, and masking them reads as an erasure. It goes in the caption.
        y_ont = y_pb - 0.50
        y_label = y_ont - 0.18
        # The top label line is allowed to sit above the axes, in the band the
        # lifted legend vacated; text is not clipped and panel A's four-row
        # legend still sets the row's top, so the plate does not grow.
        lines = 2 if counts == "in_brackets" else 1
        y_top = y_label - (lines - (1 if legend_lift else 0)) * TXT - 0.10

    ax.set_axis_off()
    ax_mat = ax.inset_axes([0.0, 0.0, 0.34, 1.0])
    ax_bar = ax.inset_axes([0.40, 0.0, 0.60, 1.0], sharey=ax_mat)
    if brackets and connector == "tint":
        # No glyph at all: a pale band from the bar down to the UpSet, whose
        # right edge at pb_end is the whole signal.
        top = y_bar + 0.31
        ax_bar.add_patch(Rectangle((0, top), pb_end, -0.34 - top, facecolor="#E8EDF0",
                                   edgecolor="none", zorder=0))
    left = 0.0
    for _name, value, color in segs:
        ax_bar.barh(y_bar, value, left=left, height=0.62, color=color, linewidth=0)
        left += value
    if brackets and counts == "study":
        _square_bracket(ax_bar, 0, total, y_hi, TIP, "")
        # The stem sits at a fifth of the span, not the middle, so it passes
        # clear of the right-aligned overlap label on the line below it.
        x_stem = 0.22 * total
        ax_bar.plot([x_stem, x_stem], [y_hi, y_head + 0.06], color=_RULE_COLOR,
                    lw=0.7, clip_on=False, zorder=3)
        ax_bar.annotate(f"PacBio+ONT ({total:,})", xy=(total / 2, y_head),
                        ha="center", va="bottom", fontsize=MIN_FONT, color="#44525C")
        x_both = 0.5 * (ont_start + pb_end)
        ax_bar.plot([x_both, x_both], [y_bar - BAR_H / 2, y_both + 0.10],
                    color=f1._HIGH_PASS_COLOR, lw=0.7, clip_on=False, zorder=3)
        ax_bar.annotate(f"Both ({segs[1][1]:,})", xy=(total, y_both), ha="right",
                        va="bottom", fontsize=MIN_FONT, color=f1._HIGH_PASS_COLOR)
        _square_bracket(ax_bar, 0, pb_end, y_lo, -TIP, "")
        ax_bar.annotate(f"This study ({pb_end:,})", xy=(pb_end / 2, y_study),
                        ha="center", va="top", fontsize=MIN_FONT, color="#44525C")
        ax_bar.plot([pb_end / 2] * 2, [y_study + TXT + 0.06, y_study + TXT + 0.69],
                    color=_RULE_COLOR, lw=0.7, clip_on=False, zorder=3)
    elif brackets:
        # Two platform brackets above the bar. They overlap over the 876
        # participants sequenced on both, which is the point -- the platforms
        # are not a partition -- so they sit at two levels to stay legible,
        # with the labels on one baseline above both.
        pb_lab, ont_lab = "PacBio", "ONT"
        if counts == "in_brackets":
            # With a tip label the PacBio total is already stated at the point
            # of the brace, attached to the set it delivers to the UpSet, so
            # the bracket does not repeat it.
            pb_lab = "PacBio" if tip_label else f"PacBio\n{pb_end:,}"
            ont_lab = f"ONT\n{total - ont_start:,}"
        marks = [_square_bracket(ax_bar, ont_start, total, y_ont, TIP, ont_lab),
                 _square_bracket(ax_bar, 0, pb_end, y_pb, TIP, pb_lab)]
        # Labels are bottom-anchored on one baseline, so a one-line label next
        # to a two-line one would sit against the lower line. Lift it instead,
        # which keeps the platform names on a common first line.
        deep = max(lab.count("\n") for _x, lab in marks)
        for x_mid, label in marks:
            ax_bar.annotate(label, xy=(x_mid, y_label - TXT * (deep - label.count("\n"))),
                            xycoords=("data", "data"),
                            ha="center", va="bottom", fontsize=MIN_FONT, color="#44525C",
                            linespacing=1.1, multialignment="center")

    if y_counts is not None and counts not in {"in_brackets", "study"}:
        anchor = None
        for i, (_name, value, color) in enumerate(segs):
            text = f"{value:,}" if i == 0 else f" + {value:,}"
            shade = {f1._MID_PASS_COLOR: "#6F8697", f1._ONT_ONLY_COLOR: "#5C8FA8"}.get(color, color)
            if anchor is None:
                anchor = ax_bar.annotate(text, xy=(0, y_counts), xycoords=("data", "data"),
                                         ha="left", va="top", fontsize=MIN_FONT, color=shade)
            else:
                anchor = ax_bar.annotate(text, xy=(1, 0), xycoords=anchor, ha="left", va="bottom",
                                         fontsize=MIN_FONT, color=shade)
    ax_mat.annotate("All Phase 2", xy=(1.0, y_bar), xycoords=("axes fraction", "data"),
                    xytext=(-1, 0), textcoords="offset points", ha="right", va="center",
                    fontsize=MIN_FONT, color="#44525C")
    if brackets and counts == "study":
        pass
    elif brackets and connector == "brace":
        # The PacBio part of the bar is what the UpSet below partitions, so a
        # brace funnels that span into it. A curly, not a square bracket:
        # different glyph, different relationship.
        _curly_brace(ax_bar, 0, pb_end, y_flat, y_tip)
        if tip_label:
            # What the brace delivers to the UpSet, stated at the point where
            # it arrives: this is the denominator the four rows sum to, and
            # the one panel E also carries.
            # Left of the point, not right: to the right it runs past the end
            # of the bar and the row layout answers by narrowing the panel.
            ax_bar.annotate(tip_label, xy=(pb_end / 2, y_tip), xytext=(-6, 0),
                            textcoords="offset points", ha="right", va="center",
                            fontsize=MIN_FONT, color="#44525C")
    elif brackets and connector == "bracket":
        y_line = y_flat + 0.40
        _square_bracket(ax_bar, 0, pb_end, y_line, -0.40, "")
        ax_bar.plot([pb_end / 2] * 2, [y_line, y_tip], color=_RULE_COLOR, lw=0.7,
                    clip_on=False, zorder=3)

    ax_bar.barh(y, feat, height=0.68, color=f1._MULTIOMIC_COLOR, linewidth=0,
                label=_FREEZE_LABELS[0])
    ax_bar.barh(y, rest, left=feat, height=0.68, color=f1._MULTIOMIC_PALE,
                linewidth=0, label=_FREEZE_LABELS[1])
    for yi, total, f_, r_ in zip(y, totals, feat, rest):
        label = f"{f_:,}+{r_:,}" if f_ and r_ else f"{total:,}"
        ax_bar.annotate(label, xy=(total, yi), xytext=(3, 0), textcoords="offset points",
                        ha="left", va="center", fontsize=MIN_FONT)
    ax_bar.set_xlabel("Participants")
    ax_bar.set_ylim(len(combos) - 0.55, y_top)
    ax_bar.set_yticks([])
    ax_bar.set_xlim(0, sum(v for _n, v, _c in segs) * 1.02)
    ax_bar.xaxis.set_major_formatter(f1._THOUSANDS)
    ax_bar.set_xticks([0, 10000])
    for side in ("top", "right", "left"):
        ax_bar.spines[side].set_visible(False)
    handles, labels = ax_bar.get_legend_handles_labels()
    if legend_lift:
        # Pad with blank rows at the bottom so the two visible lines land on
        # panel B's top legend rows rather than its bottom ones. Every row
        # legend on this plate is bottom-anchored, so adding rows is what
        # moves text up, and it stays exact under any font change. The band
        # this frees under the legend is where the platform labels overflow
        # the axes.
        handles = handles + [Patch(facecolor="none", edgecolor="none")] * legend_lift
        labels = labels + [" "] * legend_lift
    ax_bar.legend(handles, labels, loc="lower right",
                  bbox_to_anchor=(1.0, f1._ROW_LEGEND_Y), ncol=1, **f1._ROW_LEGEND)

    names = ["5mC", "Olink", "RNA-seq"]
    for yi, combo in zip(y, combos):
        present = [i for i, v in enumerate(combo) if v]
        if len(present) > 1:
            ax_mat.plot([min(present), max(present)], [yi, yi], color=f1._MULTIOMIC_COLOR, lw=1.0)
        for xi, v in enumerate(combo):
            ax_mat.scatter(xi, yi, s=18, color=f1._MULTIOMIC_COLOR if v else "#D9DEE2",
                           edgecolors="none", zorder=2)
    ax_mat.set_xticks(range(len(names)))
    ax_mat.set_xticklabels(names, rotation=90, ha="center", va="top")
    ax_mat.set_xlim(-0.5, len(names) - 0.5)
    ax_mat.set_yticks([])
    ax_mat.tick_params(axis="x", length=0, pad=2)
    for spine in ax_mat.spines.values():
        spine.set_visible(False)


def _panel_platform_assays(ax, rows) -> None:
    """One UpSet over platform and assay, with 5mC dropped.

    Every releasable Phase 2 participant has 5mC, so that column is always
    filled and carries no information; the caption states it instead. Putting
    PacBio and ONT in its place makes every column informative and lets the
    panel partition all 13,530 participants - the 1,269 ONT-only ones
    included, who appear in no other panel.

    Platform columns sit slightly apart from the assay columns so the two
    kinds of thing read as two groups rather than four equivalents.
    """
    t = f1._truthy
    avail: dict[tuple[bool, ...], int] = {}
    featured: dict[tuple[bool, ...], int] = {}
    for row in rows:
        key = (True, t(row.get("has_ONT")), t(row.get("has_proteomics")), t(row.get("has_rna")))
        avail[key] = avail.get(key, 0) + 1
        rna_ok = (not key[3]) or t(row.get("in_eqtl")) or t(row.get("in_sqtl"))
        if ((not key[2]) or t(row.get("in_pqtl"))) and rna_ok:
            featured[key] = featured.get(key, 0) + 1
    n_ont_only = len(_load_ont_only_rows())
    ont_key = (False, True, False, False)
    avail[ont_key] = avail.get(ont_key, 0) + n_ont_only
    featured[ont_key] = featured.get(ont_key, 0) + n_ont_only

    combos = sorted(avail, key=lambda k: -avail[k])
    hidden = [k for k in combos if avail[k] < _MIN_DISCLOSED]
    combos = [k for k in combos if avail[k] >= _MIN_DISCLOSED]
    feat = [featured.get(k, 0) for k in combos]
    rest = [avail[k] - featured.get(k, 0) for k in combos]
    totals = [avail[k] for k in combos]
    y = np.arange(len(combos))

    ax.set_axis_off()
    ax_mat = ax.inset_axes([0.0, 0.0, 0.40, 1.0])
    ax_bar = ax.inset_axes([0.46, 0.0, 0.54, 1.0], sharey=ax_mat)
    ax_bar.barh(y, feat, height=0.68, color=f1._MULTIOMIC_COLOR, linewidth=0,
                label=_FREEZE_LABELS[0])
    ax_bar.barh(y, rest, left=feat, height=0.68, color=f1._MULTIOMIC_PALE,
                linewidth=0, label=_FREEZE_LABELS[1])
    ax_bar.set_xlabel("Participants")
    ax_bar.set_ylim(len(combos) - 0.55, -0.55)
    ax_bar.set_yticks([])
    ax_bar.set_xlim(0, max(totals) * 1.30)
    ax_bar.xaxis.set_major_formatter(f1._THOUSANDS)
    ax_bar.set_xticks([0, 5000])
    for yi, total, f_, r_ in zip(y, totals, feat, rest):
        label = f"{f_:,}+{r_:,}" if f_ and r_ else f"{total:,}"
        ax_bar.annotate(label, xy=(total, yi), xytext=(3, 0), textcoords="offset points",
                        ha="left", va="center", fontsize=MIN_FONT)
    for side in ("top", "right", "left"):
        ax_bar.spines[side].set_visible(False)
    ax_bar.legend(loc="lower right", bbox_to_anchor=(1.0, f1._ROW_LEGEND_Y), ncol=1,
                  **f1._ROW_LEGEND)
    if hidden:
        ax_bar.annotate(
            f"{len(hidden)} further combination{'s' if len(hidden) > 1 else ''} <20",
            xy=(1.0, -0.33), xycoords="axes fraction", ha="right", va="top",
            fontsize=MIN_FONT, color="#8A8A8A",
        )

    names = ["PacBio", "ONT", "Olink", "RNA-seq"]
    xs = [0.0, 1.0, 2.3, 3.3]
    on = [f1._HIGH_PASS_COLOR, f1._ONT_COLOR, f1._MULTIOMIC_COLOR, f1._MULTIOMIC_COLOR]
    for yi, combo in zip(y, combos):
        present = [i for i, v in enumerate(combo) if v]
        if len(present) > 1:
            ax_mat.plot([xs[min(present)], xs[max(present)]], [yi, yi],
                        color="#9AA6AC", lw=1.0, zorder=1)
        for xi, v in enumerate(combo):
            ax_mat.scatter(xs[xi], yi, s=16, color=on[xi] if v else "#D9DEE2",
                           edgecolors="none", zorder=2)
    ax_mat.set_xticks(xs)
    ax_mat.set_xticklabels(names, rotation=90, ha="center", va="top")
    for tick, color in zip(ax_mat.get_xticklabels(), on):
        tick.set_color(color if color != f1._MULTIOMIC_COLOR else "#111111")
    ax_mat.set_xlim(-0.5, xs[-1] + 0.5)
    ax_mat.set_yticks([])
    ax_mat.tick_params(axis="x", length=0, pad=2)
    for spine in ax_mat.spines.values():
        spine.set_visible(False)


def _panel_molecular_blocks(ax, rows, *, whole_cohort: bool = True,
                            with_pacbio_only: bool = False) -> None:
    """Assay UpSet over a platform block, both on one Participants axis.

    Every releasable Phase 2 participant has 5mC — the 1,269 ONT-only ones
    included, and they have no Olink or RNA-seq — so with ``whole_cohort``
    the assay block partitions all 13,530 and the 5mC-only row carries both
    platforms. Set it False to keep the assay block on the 12,261 PacBio
    cohort instead, which makes no claim that PacBio and ONT methylation are
    one resource.

    The platform block is text-labelled rather than a second dot matrix: at
    two or three rows a matrix costs a second column header the panel has no
    height for, and the rows are self-describing.
    """
    t = f1._truthy
    n_ont_only = _ont_counts(rows)[1]
    avail: dict[tuple[bool, bool, bool], int] = {}
    featured: dict[tuple[bool, bool, bool], int] = {}
    for row in rows:
        key = (t(row.get("has_methylation")), t(row.get("has_proteomics")), t(row.get("has_rna")))
        if not any(key):
            continue
        avail[key] = avail.get(key, 0) + 1
        rna_ok = (not key[2]) or t(row.get("in_eqtl")) or t(row.get("in_sqtl"))
        if ((not key[1]) or t(row.get("in_pqtl"))) and rna_ok:
            featured[key] = featured.get(key, 0) + 1
    meth_only = (True, False, False)
    if whole_cohort:  # ONT-only participants have 5mC and nothing else
        avail[meth_only] = avail.get(meth_only, 0) + n_ont_only
        featured[meth_only] = featured.get(meth_only, 0) + n_ont_only
    combos = sorted(avail, key=lambda key: -avail[key])
    feat = [featured.get(k, 0) for k in combos]
    rest = [avail[k] - featured.get(k, 0) for k in combos]
    totals = [avail[k] for k in combos]

    with_pb = _ont_counts(rows)[0]
    plats = [("PacBio + ONT", with_pb, f1._HIGH_PASS_COLOR),
             ("ONT only", n_ont_only, f1._ONT_ONLY_COLOR)]
    if with_pacbio_only:
        plats.insert(0, ("PacBio only", len(rows) - with_pb, f1._MID_PASS_COLOR))
    gap = 0.9
    y_plat = [len(combos) - 1 + gap + i for i in range(len(plats))]

    ax.set_axis_off()
    ax_mat = ax.inset_axes([0.0, 0.0, 0.34, 1.0])
    ax_bar = ax.inset_axes([0.40, 0.0, 0.60, 1.0], sharey=ax_mat)
    y = np.arange(len(combos))
    ax_bar.barh(y, feat, height=0.68, color=f1._MULTIOMIC_COLOR, linewidth=0,
                label=_FREEZE_LABELS[0])
    ax_bar.barh(y, rest, left=feat, height=0.68, color=f1._MULTIOMIC_PALE,
                linewidth=0, label=_FREEZE_LABELS[1])
    for yi, (name, value, color) in zip(y_plat, plats):
        ax_bar.barh(yi, value, height=0.68, color=color, linewidth=0)
        ax_bar.annotate(f"{value:,}", xy=(value, yi), xytext=(3, 0),
                        textcoords="offset points", ha="left", va="center",
                        fontsize=MIN_FONT)
        ax_mat.annotate(name, xy=(1.0, yi), xycoords=("axes fraction", "data"),
                        xytext=(-1, 0), textcoords="offset points", ha="right",
                        va="center", fontsize=MIN_FONT, color="#44525C")
    for yi, total, f_, r_ in zip(y, totals, feat, rest):
        label = f"{f_:,}+{r_:,}" if f_ and r_ else f"{total:,}"
        ax_bar.annotate(label, xy=(total, yi), xytext=(3, 0), textcoords="offset points",
                        ha="left", va="center", fontsize=MIN_FONT)
    ax_bar.axhline(len(combos) - 1 + gap / 2, color="#D0D5D8", lw=0.6, zorder=0)
    ax_bar.set_xlabel("Participants")
    ax_bar.set_ylim(y_plat[-1] + 0.6, -0.55)
    ax_bar.set_yticks([])
    ax_bar.set_xlim(0, max(totals + [v for _n, v, _c in plats]) * 1.26)
    ax_bar.xaxis.set_major_formatter(f1._THOUSANDS)
    for side in ("top", "right", "left"):
        ax_bar.spines[side].set_visible(False)
    ax_bar.legend(loc="lower right", bbox_to_anchor=(1.0, f1._ROW_LEGEND_Y), ncol=1,
                  **f1._ROW_LEGEND)

    names = ["5mC", "Olink", "RNA-seq"]
    for yi, combo in zip(y, combos):
        present = [i for i, v in enumerate(combo) if v]
        if len(present) > 1:
            ax_mat.plot([min(present), max(present)], [yi, yi], color=f1._MULTIOMIC_COLOR, lw=1.0)
        for xi, v in enumerate(combo):
            ax_mat.scatter(xi, yi, s=18, color=f1._MULTIOMIC_COLOR if v else "#D9DEE2",
                           edgecolors="none", zorder=2)
    # Column headers sit under both blocks, as in the single-block panel.
    # Above the dots they would run into this row's legends; the platform
    # rows in between are named in place and at a different x, so the headers
    # still read as belonging to the dot columns.
    for xi, name in enumerate(names):
        ax_mat.annotate(name, xy=(xi, y_plat[-1] + 0.5), ha="center", va="top",
                        rotation=90, fontsize=MIN_FONT, color="#111111",
                        annotation_clip=False)
    ax_mat.set_xlim(-0.5, len(names) - 0.5)
    ax_mat.set_xticks([])
    ax_mat.set_yticks([])
    for spine in ax_mat.spines.values():
        spine.set_visible(False)


def _panel_molecular_compact(ax, rows) -> None:
    """One bar per assay instead of the UpSet matrix.

    The four-panel cohort row leaves the UpSet too narrow for its matrix
    labels. Per-assay totals carry the same point in a third of the width;
    the intersections stay in the supplement.
    """
    # 5mC has no QTL set of its own, so its bar is assayed-only throughout.
    keys = (("has_methylation", "5mC", ()),
            ("has_proteomics", "Olink", ("in_pqtl",)),
            ("has_rna", "RNA-seq", ("in_eqtl", "in_sqtl")))
    names, feat_vals, other_vals = [], [], []
    for flag, name, qtl in keys:
        assayed = [r for r in rows if f1._truthy(r.get(flag))]
        feat = sum(any(f1._truthy(r.get(k)) for k in qtl) for r in assayed) if qtl else 0
        names.append(name)
        feat_vals.append(feat)
        other_vals.append(len(assayed) - feat)
    y = np.arange(len(names))
    ax.barh(y, feat_vals, height=0.62, color=f1._MULTIOMIC_COLOR, linewidth=0, label="In QTL analysis")
    ax.barh(y, other_vals, left=feat_vals, height=0.62, color=f1._MULTIOMIC_PALE,
            linewidth=0, label="Assayed only")
    ax.set_yticks(y)
    ax.set_yticklabels(names)
    ax.set_ylim(len(names) - 0.5, -0.95)  # blank row above for the all-three note
    ax.set_xlabel("Participants")
    ax.set_xlim(0, max(f + o for f, o in zip(feat_vals, other_vals)) * 1.06)
    ax.xaxis.set_major_formatter(f1._THOUSANDS)
    ax.set_xticks([0, 10000])
    despine(ax)
    triple = sum(
        f1._truthy(r.get("has_methylation")) and f1._truthy(r.get("has_proteomics"))
        and f1._truthy(r.get("has_rna")) for r in rows
    )
    ax.annotate(
        f"{triple:,} have all three", xy=(1.0, -0.72), xycoords=("axes fraction", "data"),
        ha="right", va="center", fontsize=MIN_FONT, color="#44525C",
    )
    ax.legend(loc="lower right", bbox_to_anchor=(1.0, f1._ROW_LEGEND_Y), ncol=1, **f1._ROW_LEGEND)


def _load_ont_coverage(rows) -> np.ndarray:
    """ONT depth for every releasable Phase 2 participant with ONT data.

    Participants with both platforms carry ``ont_coverage`` on their PacBio
    row; ONT-only participants are a separate technology in the covariate
    table and are not in ``rows`` at all.
    """
    import csv
    import gzip

    from .style import covariates_path

    vals = [
        float(r["ont_coverage"]) for r in rows
        if f1._truthy(r.get("has_ONT")) and r.get("ont_coverage") not in ("", None)
    ]
    with gzip.open(covariates_path(), "rt", newline="") as handle:
        for row in csv.DictReader(handle):
            if row.get("technology") != "ONT":
                continue
            if str(row.get("final_releasable_v9", "")).strip() not in {"True", "true"}:
                continue
            if row.get("ont_coverage") not in ("", None):
                vals.append(float(row["ont_coverage"]))
    return np.asarray(vals, dtype=float)


def _add_ont_coverage(ax, rows, bins) -> None:
    """Overlay the ONT depth distribution on the PacBio coverage histogram.

    Drawn as an outline rather than a fill: ONT sits on top of the high-pass
    PacBio mode, and two filled distributions there are unreadable. No mean
    line, because the ONT median is within ~1x of the high-pass mean and the
    two rules would print as one.
    """
    cov = _load_ont_coverage(rows)
    med = float(np.median(cov))
    # ONT reaches 65x; the axis stops at the PacBio range, so the 11 deepest
    # samples go into the last bin rather than being dropped from the count.
    drawn = np.clip(cov, bins[0], bins[-1] - 1e-9)
    ax.hist(drawn, bins=bins, color=f1._ONT_COLOR, alpha=0.10, linewidth=0, zorder=2)
    ax.hist(drawn, bins=bins, histtype="step", color=f1._ONT_COLOR, lw=1.0, zorder=3,
            label=f"ONT (n={cov.size:,})")
    counts, _edges = np.histogram(drawn, bins=bins)
    ax.annotate(
        f"{med:.0f}×", xy=(med, float(counts.max())), xytext=(3, 1),
        textcoords="offset points", ha="left", va="bottom", fontsize=MIN_FONT,
        color=f1._ONT_COLOR, zorder=6,
        path_effects=[patheffects.withStroke(linewidth=2.0, foreground="white")],
    )
    ax.set_xlabel("Coverage (×)")
    # _panel_coverage already reordered Phase 1 to the top; ONT goes last.
    handles, labels = ax.get_legend_handles_labels()
    keyed = dict(zip(labels, handles))
    want = [t for t in labels if t.startswith("Phase 1")]
    want += [t for t in labels if t.startswith("Mid-pass")]
    want += [t for t in labels if t.startswith("High-pass")]
    want += [t for t in labels if t.startswith("ONT")]
    ax.legend([keyed[t] for t in want], want, loc="lower right",
              bbox_to_anchor=(1.0, f1._ROW_LEGEND_Y), **f1._ROW_LEGEND)


def _coverage_by_center_rail(ax, rows, tag: str | None = None, phase1_cov=None) -> None:
    """Median and interquartile PacBio coverage per center, under the histogram.

    Drawn in axes coordinates below the histogram so it shares the coverage
    axis. Centers are alphabetical and carry no count, so the rail reads as
    batch structure rather than as volume.
    """
    cov = np.array(
        [float(r["coverage"]) if r.get("coverage") not in ("", None) else np.nan for r in rows],
        dtype=float,
    )
    gc = np.array([r.get("GC") or "" for r in rows])
    # One row per center. HudsonAlpha carries two marks, Phase 1 in the Phase 1
    # color and Phase 2 in the center gray; the legend above already keys that
    # color, so the row needs no extra label.
    phase1_center = "HA"
    rails: list[tuple[str, list]] = []
    for code in _CENTER_ALPHA:
        v = cov[(gc == code) & np.isfinite(cov)]
        if v.size < 20:
            continue
        marks = []
        if code == phase1_center and phase1_cov is not None and phase1_cov.size:
            marks.append((phase1_cov, f1._PHASE1, f1._PHASE1))
        marks.append((v, "#6B7A86", "#1F2D38"))
        rails.append((_CENTER_NAME[code], marks))
    n = len(rails)
    top = ax.get_ylim()[1]
    step = 500 if top <= 2000 else 1000
    ticks = list(np.arange(0, top, step))
    # Reserve a band below the histogram baseline, one row per rail.
    frac = 0.50
    ax.set_ylim(-top * frac / (1 - frac), top)
    ax.spines["left"].set_bounds(0, top)
    ax.set_yticks(ticks)
    span = -ax.get_ylim()[0]
    rows_n = n + 1.1 if tag else n + 0.6
    if tag:
        ax.annotate(
            tag, xy=(0.0, -span * (n + 0.85) / rows_n), xycoords=("axes fraction", "data"),
            xytext=(1, 0), textcoords="offset points", ha="left", va="center",
            fontsize=MIN_FONT, color="#44525C",
        )
    for i, (name, marks) in enumerate(rails):
        yi = -span * (i + 0.8) / rows_n
        right = 0.0
        for v, bar_color, tick_color in marks:
            q1, med, q3 = np.percentile(v, [25, 50, 75])
            ax.plot([q1, q3], [yi, yi], color=bar_color, lw=1.4, solid_capstyle="butt", zorder=4)
            ax.plot([med], [yi], marker="|", color=tick_color, ms=4.5, mew=1.1, zorder=5)
            right = max(right, q3)
        ax.annotate(
            name, xy=(right, yi), xytext=(3, 0),
            textcoords="offset points", ha="left", va="center",
            fontsize=MIN_FONT, color="#44525C", annotation_clip=False,
        )


# ---------------------------------------------------------------------------
# Plate
# ---------------------------------------------------------------------------


_LAYOUTS = {
    # cohort row height, catalog axes height, loci, spectrum share, row gap
    "band": dict(cohort_h=1.45, catalog_h=4.97, loci_n=4, spec=0.552, hspace=0.34, cat_gap=0.50),
    "tall": dict(cohort_h=1.85, catalog_h=4.42, loci_n=3, spec=0.552, hspace=0.34, cat_gap=0.50),
    "compact": dict(cohort_h=1.30, catalog_h=4.20, loci_n=4, spec=0.72, hspace=0.48, cat_gap=0.34),
}


def render(
    layout: str = "band",
    center: str = "named",
    *,
    assays: str = "upset",
    ont_coverage: bool = False,
    ont_row: bool = False,
    ont_outline: bool = False,
    connector: str = "brace",
    counts: str = "under_bar",
    tip_label: str | None = None,
    anc_labels: bool = True,
    strip_rotation: float = 45.0,
    strip_pad: float = 30,
    summary: bool = True,
    stem: str | None = None,
    rng_seed: int = 7,
) -> tuple[Path, Path]:
    """One merged plate.

    ``assays`` is ``upset`` (5mC / Olink / RNA-seq, Fig. 1E's shading rule),
    ``upset_any`` (the same panel, dark = in at least one QTL sample list),
    ``upset_freeze`` (dark = every assay the participant carries was in its
    companion analysis), ``upset_plain`` (availability only, no shading),
    ``ont`` (ONT as a fourth column), or ``compact`` (one bar per assay). ``summary`` keeps the
    repeat-allele carrier histogram above the allelic series.
    """
    if center not in {"named", "omit", "shares", "coverage"}:
        raise ValueError(f"unknown center mode {center!r}")
    if layout not in _LAYOUTS:
        raise ValueError(f"unknown layout {layout!r}")
    if assays not in {"upset", "upset_any", "upset_freeze", "upset_plain",
                      "blocks", "blocks_pb", "blocks_pacbio_only",
                      "platform", "bar_upset", "bracket", "ont", "compact"}:
        raise ValueError(f"unknown assay panel {assays!r}")

    spec = _LAYOUTS[layout]
    cohort_h, catalog_h, loci_n = spec["cohort_h"], spec["catalog_h"], spec["loci_n"]
    if not summary:
        # The carrier histogram and its gap come out of the lower catalog row.
        catalog_h -= 0.85

    rows = f1._load_model_rows()
    phase1_rows = f1._load_phase1_rows()
    illumina_f, illumina_m, illumina_n = f1._load_illumina_reference()
    ont_only = f1._load_ont_only_counts()

    block_guess = cohort_h + 1.05
    cat_gap = spec.get("cat_gap", _CAT_GAP)
    height = f1._MARGIN + block_guess + cat_gap + catalog_h + _CAT_FOOT
    fig = new_figure(height)

    # --- cohort row -------------------------------------------------------
    ax_pyr = fig.add_axes([0.08, 0.72, 0.34, 0.20])
    ax_ref, handles, labels = f1._panel_age_sex(
        ax_pyr, rows, illumina_f, illumina_m, illumina_n, f1._pyramid(),
        phase1_rows=phase1_rows,
    )
    # One narrow column. Side-by-side legends need more width than panel A
    # gets once the catalog tier is on the same plate, and the right-hand
    # Illumina axis is gray like its swatch, so the pairing still reads.
    if ont_outline:
        ont_handle = _add_ont_only_outline(ax_pyr, f1._pyramid(), _load_ont_only_rows())
        handles = [*handles[:3], ont_handle, handles[3]]
        labels = [*labels[:3], f"ONT only ({len(_load_ont_only_rows()):,})", labels[3]]
    compact = {
        labels[1]: labels[1].replace("Phase 2, with phenotypes", "Phase 2 + phenotypes"),
        labels[2]: labels[2].replace("Phase 2, without phenotypes", "Phase 2, no phenotypes"),
    }
    ax_pyr.legend(
        handles[:5] if ont_outline else handles[:4],
        [compact.get(t, t) for t in (labels[:5] if ont_outline else labels[:4])],
        loc="lower left", bbox_to_anchor=(0.0, 1.02), ncol=1,
        frameon=False, fontsize=MIN_FONT, handlelength=1.4, handleheight=0.7,
        handletextpad=0.4, borderaxespad=0.0, labelspacing=0.42,
        handler_map={tuple: HandlerTuple(ndivide=None, pad=0)},
    )
    # Four cohort panels leave panel A too narrow for a label every decade.
    ax_pyr.set_xticks([20, 40, 60, 90])
    ax_pyr.set_xticklabels(["20", "40", "60", "90+"])
    # A short cohort row thins the count axis to a single pair; the pyramid
    # reads better with the 50s kept. The formatter already prints magnitudes.
    ax_pyr.yaxis.set_major_locator(MultipleLocator(50))

    panels = [f1._Panel("A", ax_pyr, extra=(ax_ref,), weight=2.05)]

    if center in {"named", "shares"}:
        ax_gc = fig.add_axes([0.45, 0.72, 0.14, 0.20])
        if center == "named":
            f1._panel_center(ax_gc, rows, ont_only, phase1_n=len(phase1_rows))
        else:
            _panel_center_shares(ax_gc, rows, ont_only, phase1_n=len(phase1_rows))
        panels.append(f1._Panel(chr(ord("A") + len(panels)), ax_gc, weight=1.30))

    ax_cov = fig.add_axes([0.62, 0.72, 0.14, 0.20])
    cov_bins = np.arange(4, 53, 1)
    f1._panel_coverage(ax_cov, rows, phase1_cov=f1._phase1_coverage(phase1_rows))
    if ont_coverage:
        _add_ont_coverage(ax_cov, rows, cov_bins)
    if center == "coverage":
        _coverage_by_center_rail(ax_cov, rows, phase1_cov=f1._phase1_coverage(phase1_rows))
    panels.append(f1._Panel(chr(ord("A") + len(panels)), ax_cov, weight=1.30 if center == "coverage" else 1.15))

    ax_mol = fig.add_axes([0.80, 0.72, 0.14, 0.20])
    if len(panels) == 3 or assays == "compact":  # four-panel row: no space for a matrix
        _panel_molecular_compact(ax_mol, rows)
        panels.append(f1._Panel(chr(ord("A") + len(panels)), ax_mol, weight=1.05))
    elif assays == "ont":
        _panel_molecular_ont(ax_mol, rows, sum(ont_only.values()))
        panels.append(f1._Panel(chr(ord("A") + len(panels)), ax_mol, weight=1.35))
    elif assays == "bracket":
        _panel_cohort_bracket(ax_mol, rows)
        panels.append(f1._Panel(chr(ord("A") + len(panels)), ax_mol, weight=1.70))
    elif assays == "bar_upset":
        _panel_bar_over_upset(ax_mol, rows, connector=connector, counts=counts,
                              tip_label=tip_label)
        # The study header is wider than the bar, but it costs no width:
        # measured, the panel is identical with and without it. The UpSet's
        # count labels already set the row's right-hand ink, and the gap
        # between the matrix and bar insets absorbs the overhang on the left.
        panels.append(f1._Panel(chr(ord("A") + len(panels)), ax_mol, weight=1.3))
    elif assays == "platform":
        _panel_platform_assays(ax_mol, rows)
        panels.append(f1._Panel(chr(ord("A") + len(panels)), ax_mol, weight=1.35))
    elif assays.startswith("blocks"):
        _panel_molecular_blocks(
            ax_mol, rows,
            whole_cohort=assays != "blocks_pb",
            with_pacbio_only=assays == "blocks_pacbio_only",
        )
        panels.append(f1._Panel(chr(ord("A") + len(panels)), ax_mol, weight=1.3))
    elif assays in {"upset_any", "upset_freeze", "upset_plain"}:
        _panel_molecular_any_qtl(
            ax_mol, rows, split=assays != "upset_plain",
            rule="any" if assays == "upset_any" else "freeze",
            ont_row=ont_row,
        )
        panels.append(f1._Panel(chr(ord("A") + len(panels)), ax_mol, weight=1.25))
    else:
        f1._panel_molecular_data(ax_mol, rows)
        panels.append(f1._Panel(chr(ord("A") + len(panels)), ax_mol, weight=1.25))

    cohort_letters, floor = _settle_row(fig, panels, cohort_h, f1._MARGIN, inline_letters=True)

    # --- catalog block ----------------------------------------------------
    rng = np.random.default_rng(rng_seed)
    from .fig2_panel_c_options import load_chosen_loci, panel_rarefaction_examples

    data = f2.spectrum_data(rng, region_path=f2._region_bins_path())
    cat_top = (floor - cat_gap) / height
    cat_bottom = _CAT_FOOT / height
    gs = f2.GridSpec(
        2, 1, figure=fig, left=0.09, right=0.975, top=cat_top, bottom=cat_bottom,
        height_ratios=[spec["spec"], 1.2], hspace=spec["hspace"],
    )
    ax_spec = f2.draw_spectrum(
        fig, gs[0], data, tracks=(), p1=True, stack_ctx=True, bnd_right=True,
        scope=True, scope_brackets=False, bin_notes=True,
    )
    row2 = gs[1].subgridspec(1, 2, wspace=0.24, width_ratios=[1.0, 1.12])
    ax_disc = f2.panel_discovery_structure(
        fig, row2[0, 0], rng,
        strip_labels=ANC_LIKE if anc_labels else None,
        strip_label_rotation=strip_rotation if anc_labels else 0.0,
        strip_label_pad=strip_pad if anc_labels else 12,
    )
    ax_sum, ax_rows = panel_rarefaction_examples(
        fig, row2[0, 1], rng, examples=load_chosen_loci()[:loci_n],
        ratios=(0.48, 1.62), summary="carrier_unit" if summary else "none",
        row_scales="symlog",
    )
    if not summary:
        # The grey / blue key lived on the carrier histogram.
        from .fig2_panel_c_options import legend_handles

        ax_rows.legend(
            handles=legend_handles(), loc="lower right", bbox_to_anchor=(1.0, 1.005),
            ncol=2, frameon=False, fontsize=MIN_FONT, handlelength=0.9,
            handleheight=0.7, handletextpad=0.3, columnspacing=0.8, borderaxespad=0.0,
        )

    n_cohort = len(panels)
    cat_letters = [chr(ord("A") + n_cohort + i) for i in range(4)]
    fig_w, fig_h = fig.get_size_inches()
    # The cohort letters sit on the figure margin; the full-width spectrum and
    # the discovery panel below it share that column.
    letter_x = f1._MARGIN / fig_w
    # Without the carrier histogram the lower row sits closer to the spectrum's
    # x label, so its letters ride lower.
    row2_dy = 0.02 if summary else 0.010
    spec_letter = fig.text(
        letter_x, ax_spec.get_position().y1 + 0.024, cat_letters[0],
        fontsize=f1._LETTER_SIZE, fontweight="bold", ha="left", va="bottom",
    )
    fig.text(
        letter_x, ax_disc.get_position().y1 + row2_dy, cat_letters[1],
        fontsize=f1._LETTER_SIZE, fontweight="bold", ha="left", va="bottom",
    )
    f2.letter(fig, ax_sum, cat_letters[2], dx=0.075, dy=row2_dy)
    if summary:
        fig.text(
            ax_sum.get_position().x0 - 0.075, ax_rows.get_position().y1 + 0.012,
            cat_letters[3], fontsize=9, fontweight="bold", ha="left", va="bottom",
        )
    _align_letter_to_legend(fig, spec_letter, "US")
    f2._match_xlabel_baseline(fig, "Repeat units from reference", "Participants in discovery order")
    # Panel heights in inches, for tuning the cohort / spectrum balance.
    LAST_GEOMETRY.update(
        figure_h=fig_h, cohort=cohort_h, floor=floor,
        catalog_region=floor - cat_gap - _CAT_FOOT,
        spectrum=ax_spec.get_position().height * fig_h,
        discovery=ax_disc.get_position().height * fig_h,
    )

    fig_w, fig_h = fig.get_size_inches()
    for letter_, x, y in cohort_letters:
        fig.text(
            x / fig_w, y / fig_h, letter_, fontsize=f1._LETTER_SIZE,
            fontweight="bold", ha="left", va="bottom", color="#111111",
        )
    return save(fig, stem or f"fig12_{layout}_{center}")


VERSIONS = [
    ("band", "named"),
    ("band", "omit"),
    ("band", "shares"),
    ("band", "coverage"),
    ("tall", "named"),
    ("tall", "shares"),
]


def render_all() -> list[tuple[Path, Path]]:
    return [render(layout, center) for layout, center in VERSIONS]


if __name__ == "__main__":
    for pdf, png in render_all():
        print(png)
