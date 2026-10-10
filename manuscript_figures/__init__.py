"""Render Science-style mock-up figures for the AoU-LR Phase 2 manuscript.

Mock-up panels are written unlabeled to mockup_figures/panels/. Combined
plates with A/B/C letters are kept only for figures that already use real
data (Fig. 1 mosaic, Fig. S3, Fig. S16, Fig. S18).
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt

from . import fig1_null as fig1_cohort
from . import fig2_associations
from . import fig2_full_model
from . import fig3_mhc
from . import fig3_scan
from . import fig3_scan_chapter_cards
from . import fig3_scan_merged
from . import fig3_scan_results
from . import figSx_testing_validation
from . import fig4_mei
from . import fig5_mtdna
from . import fig6_imputation
from . import figS1_assembly
from . import figS2_callset
from . import figS3_methylation
from . import figS4_repeats
from . import figS5_qtl
from . import figS6_p1p2_sv
from . import figS7_size_spectrum
from . import figS8_ebert_discovery
from . import figS9_ebert_facets
from . import figS10_panel_novelty
from . import figS11_how_many_more
from . import figS12_new_alleles
from . import figS13_assoc_atlas
from . import figS14_phewas_coloc
from . import figS15_phewas_heatmap
from . import figS16_zip_map
from . import figS17_signed_length
from . import figS18_p1p2_counts
from .style import OUTDIR, apply_style

MAIN_FIGURES = [
    ("fig1_cohort", fig1_cohort.render),
    ("fig2_associations", fig2_associations.render),
    ("fig3_mhc", fig3_mhc.render),
    ("fig4_mei", fig4_mei.render),
    ("fig5_mtdna", fig5_mtdna.render),
    ("fig6_imputation", fig6_imputation.render),
]

# Layout exploration; not part of the six-figure contact sheet.
EXPLORATION_FIGURES = [
    ("fig2_full_model", fig2_full_model.render),
    ("fig3_scan", fig3_scan.render),
    ("fig3_scan_chapter_cards", fig3_scan_chapter_cards.render),
    ("fig3_scan_merged", fig3_scan_merged.render),
    ("fig3_scan_results", fig3_scan_results.render),
    ("figSx_testing_validation", figSx_testing_validation.render),
]

SUPP_FIGURES = [
    ("figS1_assembly", figS1_assembly.render),
    ("figS2_callset", figS2_callset.render),
    ("figS3_methylation", figS3_methylation.render),
    ("figS4_repeats", figS4_repeats.render),
    ("figS5_qtl", figS5_qtl.render),
    ("figS6_p1p2_sv", figS6_p1p2_sv.render),
    ("figS7_size_spectrum", figS7_size_spectrum.render),
    ("figS8_ebert_discovery", figS8_ebert_discovery.render),
    ("figS9_ebert_facets", figS9_ebert_facets.render),
    ("figS10_panel_novelty", figS10_panel_novelty.render),
    ("figS11_how_many_more", figS11_how_many_more.render),
    ("figS12_new_alleles", figS12_new_alleles.render),
    ("figS13_assoc_atlas", figS13_assoc_atlas.render),
    ("figS14_phewas_coloc", figS14_phewas_coloc.render),
    ("figS15_phewas_heatmap", figS15_phewas_heatmap.render),
    ("figS16_zip_map", figS16_zip_map.render),
    ("figS17_signed_length", figS17_signed_length.render),
    ("figS18_p1p2_counts", figS18_p1p2_counts.render),
]

FIGURES = MAIN_FIGURES + EXPLORATION_FIGURES + SUPP_FIGURES
_MAIN_NAMES = {name for name, _ in MAIN_FIGURES}
_SUPP_NAMES = {name for name, _ in SUPP_FIGURES}


def _contact_sheet(pngs: list[Path], titles: list[str], stem: str, heading: str, nrow: int, ncol: int) -> Path:
    apply_style()
    fig, axes = plt.subplots(nrow, ncol, figsize=(7.2, min(12.4, 2.55 * nrow + 1.0)))
    fig.subplots_adjust(left=0.03, right=0.99, top=0.96, bottom=0.03, wspace=0.04, hspace=0.10)
    axes_flat = np_axes(axes)
    for i, ax in enumerate(axes_flat):
        if i < len(pngs):
            img = plt.imread(pngs[i])
            ax.imshow(img)
            ax.set_title(titles[i], loc="left", fontsize=8, pad=2)
        ax.axis("off")
    fig.suptitle(heading, fontsize=10, fontweight="bold", y=0.99)
    OUTDIR.mkdir(parents=True, exist_ok=True)
    out = OUTDIR / f"{stem}.png"
    fig.savefig(out, dpi=220, bbox_inches="tight", pad_inches=0.04)
    plt.close(fig)
    return out


def np_axes(axes):
    import numpy as np

    return np.asarray(axes).ravel()


def render_all(which: list[str] | None = None) -> list[tuple[str, str, str]]:
    out = []
    main_pngs: list[Path] = []
    supp_pngs: list[Path] = []
    lookup = {name: fn for name, fn in FIGURES}
    selected = []
    if which:
        for name, fn in FIGURES:
            if name in which or any(name.startswith(w) or w in name for w in which):
                selected.append((name, fn))
    else:
        selected = list(FIGURES)
    for name, fn in selected:
        pdf, png = fn()
        out.append((name, str(pdf), str(png)))
        print(f"wrote {png}")
        if name in _SUPP_NAMES:
            supp_pngs.append(png)
        elif name in _MAIN_NAMES:
            main_pngs.append(png)
    if len(main_pngs) == 6:
        sheet = _contact_sheet(
            main_pngs,
            [
                "Fig. 1  Cohort and local ancestry",
                "Fig. 2  SV–EHR associations",
                "Fig. 3  MHC and KIR",
                "Fig. 4  MEIs and SVAs",
                "Fig. 5  mtDNA",
                "Fig. 6  Imputation",
            ],
            "fig0_contact_sheet",
            "AoU-LR Phase 2  ·  main-text figure mock-ups",
            3,
            2,
        )
        print(f"wrote {sheet}")
    if len(supp_pngs) == len(SUPP_FIGURES):
        import math

        n = len(supp_pngs)
        ncol = 2
        nrow = math.ceil(n / ncol)
        titles = [
            "Fig. S1  Assembly and pedigrees",
            "Fig. S2  SV callset",
            "Fig. S3  Methylation",
            "Fig. S4  Tandem repeats",
            "Fig. S5  Multi-omic QTLs",
            "Fig. S6  Phase 1 vs Phase 2 SVs",
            "Fig. S7  Callset length spectrum",
            "Fig. S8  Ebert-style SV discovery",
            "Fig. S9  Discovery by type and context",
            "Fig. S10  Novelty vs HPRC/HGSVC3",
            "Fig. S11  How many SVs remain",
            "Fig. S12  New alleles at old sites",
            "Fig. S13  SV–EHR association atlas",
            "Fig. S14  PheWAS and multi-omic coloc",
            "Fig. S15  SV–phenome Z heatmap",
            "Fig. S16  ZIP3 map by ancestry",
            "Fig. S17  Signed variant length",
            "Fig. S18  Phase 1 vs Phase 2 counts",
        ]
        sheet = _contact_sheet(
            supp_pngs,
            titles[:n],
            "figS0_contact_sheet",
            "AoU-LR Phase 2  ·  supplementary figure mock-ups",
            nrow,
            ncol,
        )
        print(f"wrote {sheet}")
    return out
