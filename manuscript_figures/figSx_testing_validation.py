"""Supplementary figure mock-up: calibration and added value of dosage- and ancestry-aware testing.

Moved out of Fig. 3 when that figure became results-first. Numbers in B and C mark the loci shown in Fig. 3.
All panels are simulated placeholders until the scan is run.
"""

from __future__ import annotations

from . import fig3_scan_merged as M
from .fig3_scan import add, panel_ancestry, panel_dosage
from .style import new_figure, save, stamp_mockup

TOP, AX_H, H = 0.50, 1.35, 2.40


def render():
    fig = new_figure(H)
    for ch, ttl, fn, x_title, w in (  # each y-axis sits at the left edge of its panel title
        ("A", "Calibration on null phenotypes", M.panel_null_qq, 0.30, 2.08),
        ("B", "Dosage vs biallelic coding", lambda ax: panel_dosage(ax, M.SHOWN), 2.74, 1.93),
        ("C", "Ancestry-aware vs pooled", lambda ax: panel_ancestry(ax, M.SHOWN), 5.04, 1.93),
    ):
        ax = add(fig, x_title, TOP, w, AX_H)
        fn(ax)
        M._head(fig, ch, ttl + "*", x_title - M.HEAD_DX, TOP - 0.27)
    stamp_mockup(fig, "*Not real data yet: all panels are simulated. Numbers in B and C mark the loci shown in Fig. 3.")
    return save(fig, "figSx_testing_validation")
