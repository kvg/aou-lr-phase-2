"""Figure S9: Ebert discovery faceted by INS/DEL and unique vs repetitive sequence."""

from __future__ import annotations

from pathlib import Path

import numpy as np
from matplotlib.gridspec import GridSpec

from .figS8_ebert_discovery import (
    _ancestry_bar,
    _incremental,
    _p2_blocks,
    _stack,
    _strata,
)
from .style import (
    P2_N,
    SV_N,
    ax_panel,
    format_millions,
    new_figure,
    save,
    stamp_mockup,
)

# Placeholder splits of Table 2 ≥50 bp INS/DEL. Repetitive = rmsk ∪ simpleRepeat ∪ segdup.
FACETS = (
    ("INS, unique sequence", 0.22 * SV_N["INS_50"], 4_200, 2.05, 0.56),
    ("INS, repetitive", 0.78 * SV_N["INS_50"], 11_800, 1.40, 0.46),
    ("DEL, unique sequence", 0.38 * SV_N["DEL_50"], 5_400, 2.05, 0.54),
    ("DEL, repetitive", 0.62 * SV_N["DEL_50"], 6_200, 1.40, 0.44),
)


def render() -> tuple[Path, Path]:
    rng = np.random.default_rng(22)
    afr_start = next(s for a, s, _ in _p2_blocks() if a == "AFR")
    strata = []
    for _, target, per_genome, afr_mult, sing_max in FACETS:
        new = _incremental(P2_N, per_genome, float(target), afr_start, afr_mult, rng)
        strata.append(
            _strata(
                new,
                shared_inf=max(80.0, 0.004 * target),
                major_plat=max(400.0, 0.012 * target),
                major_tau=80,
                sing_tau=1500,
                sing_max=sing_max,
                afr_start=afr_start,
            )
        )

    fig = new_figure(7.55)
    gs = GridSpec(
        2,
        1,
        figure=fig,
        left=0.09,
        right=0.98,
        top=0.91,
        bottom=0.07,
        hspace=0.10,
        height_ratios=[1.0, 0.09],
    )
    top = gs[0].subgridspec(2, 2, wspace=0.22, hspace=0.38)
    axes = [fig.add_subplot(top[i, j]) for i in range(2) for j in range(2)]
    ax_bar = fig.add_subplot(gs[1])

    slugs = ["ins_unique", "ins_repetitive", "del_unique", "del_repetitive"]
    for i, (ax, st, (title, *_rest)) in enumerate(zip(axes, strata, FACETS)):
        _stack(ax, st, afr_start, ylabel=(i % 2 == 0), arrows=(i == 0))
        ax.set_title(title, loc="left", pad=8 if i < 2 else 2)
        ax.text(
            0.98,
            0.88,
            format_millions(st["total"][-1]),
            transform=ax.transAxes,
            ha="right",
            fontsize=6.0,
            color="#333333",
        )
        if i < 2:
            ax.tick_params(labelbottom=False)

    ymax_ins = max(strata[0]["total"][-1], strata[1]["total"][-1]) * 1.08
    ymax_del = max(strata[2]["total"][-1], strata[3]["total"][-1]) * 1.08
    for ax in axes[:2]:
        ax.set_ylim(0, ymax_ins)
    for ax in axes[2:]:
        ax.set_ylim(0, ymax_del)

    axes[0].legend(loc="upper left", fontsize=5.2, bbox_to_anchor=(0.0, 0.84))
    _ancestry_bar(ax_bar, P2_N)

    for i, (st, (title, *_rest), slug) in enumerate(zip(strata, FACETS, slugs)):
        def _facet(ax, st=st, title=title, i=i):
            _stack(ax, st, afr_start, ylabel=True, arrows=(i == 0))
            ax.set_title(title, loc="left", pad=2)
            ax.text(
                0.98,
                0.88,
                format_millions(st["total"][-1]),
                transform=ax.transAxes,
                ha="right",
                fontsize=6.0,
                color="#333333",
            )
            ymax = ymax_ins if i < 2 else ymax_del
            ax.set_ylim(0, ymax)
            if i == 0:
                ax.legend(loc="upper left", fontsize=5.2, bbox_to_anchor=(0.0, 0.84))
            ax.set_xlabel("Discovery sample")

        ax_panel(f"figS9_{slug}", 2.65, _facet, left=0.12, right=0.98, top=0.86, bottom=0.16)

    stamp_mockup(
        fig,
        extra="Phase 2 only. Repetitive = rmsk + simpleRepeat + segdup. Counts are placeholders",
    )
    return save(fig, "figS9_ebert_facets")
