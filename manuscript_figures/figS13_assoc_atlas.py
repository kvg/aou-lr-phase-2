"""Figure S13: SV–EHR association atlas (scale + examples).

Distinct from Fig. 2 (design schematic, tagging, one Manhattan, coloc) and
Fig. 6 (imputation panel, imputed PheWAS). This plate is the many-G × many-P
highlight reel: how large the scan is, where hits sit in the phenome, one
variant's pleiotropy, and a few loci where the SV is the story.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
from matplotlib.gridspec import GridSpec
from matplotlib.lines import Line2D
from matplotlib.patches import FancyBboxPatch, Rectangle

from .style import (
    COL_IN,
    OKABE,
    P1_N,
    P2_EHR_N,
    P2_N,
    ax_panel,
    despine,
    new_figure,
    save,
    save_panel,
    stamp_mockup,
)

# Phenotype categories for the genome × phenome map (placeholder counts).
CATS = [
    "Hematology",
    "Lipids",
    "Renal",
    "Liver",
    "Metabolic",
    "Cardiovascular",
    "Autoimmune",
    "Pulmonary",
    "Anthropometry",
    "Inflammation",
    "Medications",
    "Other",
]

KIND_COLOR = {
    "sv": OKABE["vermillion"],
    "snv": "#9A9A9A",
    "anc": OKABE["blue"],
}

# One column per locus, genomic order. Hits: (category, log10p, kind).
# kind: sv = SV remains lead; snv = SNV remains lead; anc = ancestry-stratified.
LOCI: list[tuple[int, str, list[tuple[str, float, str]]]] = [
    (1, "AMY1", [("Metabolic", 8.4, "sv")]),
    (1, "GSTM1", [("Other", 8.0, "anc")]),
    (1, "RHD", [("Hematology", 8.7, "sv")]),
    (1, "NPHS2", [("Renal", 10.1, "sv")]),
    (1, "CFH", [("Autoimmune", 9.7, "sv"), ("Renal", 7.8, "sv")]),
    (2, "LCT", [("Metabolic", 7.9, "snv")]),
    (2, "UGT1A1", [("Liver", 12.4, "sv"), ("Other", 8.1, "sv")]),
    (3, "CCR5", [("Other", 8.3, "sv")]),
    (3, "GPX1", [("Inflammation", 7.7, "sv")]),
    (4, "HGFAC", [("Liver", 7.6, "sv")]),
    (4, "GYPB", [("Hematology", 9.4, "anc")]),
    (5, "SMN1", [("Other", 7.8, "sv")]),
    (5, "IRX1", [("Anthropometry", 7.6, "snv")]),
    (6, "MHC", [("Autoimmune", 15.6, "snv"), ("Inflammation", 9.2, "snv"), ("Medications", 7.8, "snv")]),
    (6, "C4", [("Autoimmune", 10.4, "sv")]),
    (6, "LPA", [("Lipids", 16.1, "sv"), ("Cardiovascular", 11.2, "sv")]),
    (7, "PMS2", [("Other", 7.5, "sv")]),
    (7, "CYP3A4", [("Medications", 8.6, "sv")]),
    (8, "DEFA1", [("Hematology", 8.8, "sv"), ("Inflammation", 7.9, "sv")]),
    (8, "FCGR1A", [("Hematology", 9.9, "anc")]),
    (9, "ABO", [("Hematology", 8.2, "snv"), ("Cardiovascular", 7.6, "snv")]),
    (10, "CYP2C19", [("Medications", 8.4, "sv")]),
    (11, "HBB", [("Hematology", 12.6, "sv")]),
    (11, "INS", [("Metabolic", 8.1, "snv")]),
    (12, "ALDH2", [("Liver", 8.2, "sv")]),
    (12, "HNF1A", [("Metabolic", 7.9, "snv")]),
    (14, "SERPINA1", [("Liver", 8.9, "sv"), ("Pulmonary", 8.0, "sv")]),
    (15, "CHRNA5", [("Pulmonary", 7.7, "snv")]),
    (16, "HBA", [("Hematology", 18.4, "sv"), ("Other", 8.1, "sv")]),
    (16, "PDXDC1", [("Anthropometry", 7.5, "sv")]),
    (17, "CCL3L1", [("Inflammation", 8.5, "sv"), ("Hematology", 7.8, "sv")]),
    (17, "NF1", [("Anthropometry", 7.6, "sv")]),
    (19, "LDLR", [("Lipids", 9.8, "snv"), ("Cardiovascular", 7.9, "snv")]),
    (19, "APOE", [("Lipids", 8.6, "snv")]),
    (19, "KIR", [("Autoimmune", 8.8, "sv")]),
    (21, "TMPRSS2", [("Pulmonary", 7.6, "sv")]),
    (22, "APOL1", [("Renal", 14.2, "anc"), ("Other", 7.6, "anc")]),
    (22, "IGL", [("Autoimmune", 7.9, "sv")]),
]


def _panel_hitmap(ax) -> None:
    """One column per locus, grouped by chromosome; gene names on the x-axis."""
    cat_y = {c: i for i, c in enumerate(CATS)}
    chrom_gap = 0.85
    xs: list[float] = []
    genes: list[str] = []
    x = 0.0
    prev = None
    spans: list[tuple[int, float, float]] = []  # chrom, x0, x1
    span_start = 0.0
    for chrom, gene, _hits in LOCI:
        if prev is not None and chrom != prev:
            spans.append((prev, span_start, x - 1.0))
            x += chrom_gap
            span_start = x
        elif prev is None:
            span_start = x
        xs.append(x)
        genes.append(gene)
        prev = chrom
        x += 1.0
    spans.append((prev, span_start, xs[-1]))

    # Alternating chromosome bands.
    ymax = len(CATS) - 0.35
    for i, (_c, x0, x1) in enumerate(spans):
        if i % 2 == 0:
            ax.add_patch(
                Rectangle(
                    (x0 - 0.45, -0.45),
                    x1 - x0 + 0.9,
                    ymax + 0.55,
                    facecolor="#F3F6F8",
                    edgecolor="none",
                    zorder=0,
                )
            )

    counts = {c: 0 for c in CATS}
    for xi, (_chrom, _gene, hits) in zip(xs, LOCI):
        for cat, logp, kind in hits:
            if cat not in cat_y:
                continue
            counts[cat] += 1
            ax.scatter(
                [xi],
                [cat_y[cat]],
                s=14 + 1.8 * max(0.0, logp - 7.3),
                c=KIND_COLOR[kind],
                linewidths=0.25,
                edgecolors="white",
                zorder=3,
            )

    ax.set_xticks(xs)
    ax.set_xticklabels(genes, fontsize=5.3, fontstyle="italic")
    ax.tick_params(axis="x", top=False, labeltop=False, bottom=True, labelbottom=True, length=2.0, pad=1.5)
    for tick in ax.get_xticklabels():
        tick.set_rotation(90)
        tick.set_ha("right")
        tick.set_va("center")
        tick.set_rotation_mode("anchor")

    ax.set_yticks(range(len(CATS)))
    ax.set_yticklabels(CATS, fontsize=6.2)
    xmax = xs[-1]
    ax.set_xlim(-0.7, xmax + 2.6)
    # Hematology at the top; room above for chromosome numbers.
    ax.set_ylim(len(CATS) - 0.40, -1.05)
    ax.set_title("Genome-wide SV–trait pairs, one column per locus", loc="left", pad=2)
    despine(ax, bottom=True)

    for chrom, x0, x1 in spans:
        ax.plot([x0 - 0.35, x1 + 0.35], [-0.62, -0.62], color="#888888", lw=0.55, clip_on=False, zorder=1)
        ax.text((x0 + x1) / 2, -0.92, str(chrom), ha="center", va="top", fontsize=6.0, color="#333333")

    for i, cat in enumerate(CATS):
        ax.text(xmax + 0.7, i, str(counts[cat]), ha="left", va="center", fontsize=5.8, color="#555555")
    ax.text(xmax + 0.7, -0.55, "n", ha="left", va="center", fontsize=5.2, color="#888888")

    handles = [
        Line2D([0], [0], marker="o", color="none", markerfacecolor=KIND_COLOR["sv"], markersize=5, label="SV is lead"),
        Line2D([0], [0], marker="o", color="none", markerfacecolor=KIND_COLOR["snv"], markersize=5, label="SNV remains lead"),
        Line2D([0], [0], marker="o", color="none", markerfacecolor=KIND_COLOR["anc"], markersize=5, label="Ancestry-stratified"),
    ]
    ax.legend(handles=handles, loc="lower left", bbox_to_anchor=(0.58, 1.02), ncol=3, fontsize=5.5, borderaxespad=0.0)


def _panel_scale(ax) -> None:
    """How large the scan is — not another pipeline cartoon."""
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    ax.set_title(
        f"Native LRS scan  ·  n={P2_N:,} PacBio  ·  EHR n={P2_EHR_N:,}  ·  7 ancestries, SAIGE + METAL",
        loc="left",
        pad=2,
    )
    items = [
        ("180k", "SVs tested\n(MAF and QC)"),
        ("1,186", "EHR traits\n(phecodes, labs)"),
        ("214M", "association\ntests"),
        ("847", "genome-wide\nhits"),
        ("94", "traits with\n≥1 hit"),
        ("61%", "loci where the\nSV is lead"),
    ]
    n = len(items)
    gap = 0.012
    w = (1.0 - gap * (n - 1)) / n
    for i, (num, lab) in enumerate(items):
        x0 = i * (w + gap)
        ax.add_patch(
            FancyBboxPatch(
                (x0, 0.04),
                w,
                0.90,
                boxstyle="round,pad=0.012,rounding_size=0.04",
                facecolor="#F4F4F4",
                edgecolor="#DDDDDD",
                linewidth=0.6,
                transform=ax.transAxes,
                clip_on=False,
            )
        )
        ax.text(x0 + w / 2, 0.62, num, ha="center", va="center", fontsize=11.5, fontweight="bold", transform=ax.transAxes)
        ax.text(x0 + w / 2, 0.26, lab, ha="center", va="center", fontsize=5.6, color="#444444", transform=ax.transAxes, linespacing=1.15)


def _panel_phewas(ax) -> None:
    """One SV × the phenome: HBA α-thal deletion (classic PheWAS)."""
    rows = [
        ("Mean corpuscular volume", 18.1, "heme"),
        ("Mean corpuscular hemoglobin", 16.4, "heme"),
        ("Red blood cell count", 12.8, "heme"),
        ("Hemoglobin", 9.6, "heme"),
        ("Thalassemia phecode", 8.9, "heme"),
        ("Anemia phecode", 7.4, "heme"),
        ("Hematocrit", 6.1, "heme"),
        ("Total bilirubin", 4.2, "liver"),
        ("Hemoglobinopathy phecode", 3.8, "inf"),
        ("Ferritin", 2.9, "heme"),
        ("Platelet count", 1.8, "heme"),
        ("LDL cholesterol", 1.4, "lipid"),
        ("eGFR", 1.2, "renal"),
        ("HbA1c", 1.1, "met"),
        ("BMI", 0.9, "anth"),
        ("C-reactive protein", 0.8, "infl"),
        ("Systolic BP", 0.7, "cv"),
        ("Type 2 diabetes", 0.6, "met"),
        ("ALT", 0.5, "liver"),
        ("White blood cell count", 0.4, "heme"),
    ]
    cat_col = {
        "heme": OKABE["vermillion"],
        "liver": OKABE["orange"],
        "inf": OKABE["blue"],
        "lipid": OKABE["green"],
        "renal": OKABE["sky"],
        "met": "#888888",
        "anth": "#888888",
        "infl": OKABE["pink"],
        "cv": "#888888",
    }
    y = np.arange(len(rows))[::-1]
    colors = [cat_col[c] for _, _, c in rows]
    logp = np.array([v for _, v, _ in rows])
    ax.barh(y, logp, color=colors, height=0.72, zorder=2)
    ax.axvline(7.5, color="#444444", ls="--", lw=0.6, zorder=3)
    ax.set_yticks(y)
    ax.set_yticklabels([n for n, *_ in rows], fontsize=5.6)
    ax.set_xlabel(r"$-\log_{10}(P)$")
    ax.set_xlim(0, 20.5)
    ax.set_title(r"PheWAS of one SV: HBA1/2 $\alpha$-thal deletion", loc="left", pad=2)
    despine(ax)
    ax.text(7.7, len(rows) - 0.3, "phenome-wide", fontsize=5.2, color="#444444", va="bottom")


def _panel_scale_vs_prior(ax) -> None:
    """Larger than Phase 1; imputed srWGS is the companion paper's full resource."""
    labels = [
        f"Phase 1 LRS\nn={P1_N:,} AFR",
        f"Phase 2 LRS\nn={P2_N:,}",
        "Imputed srWGS\nn≈245k",
    ]
    hits = np.array([38, 847, 4180], dtype=float)
    traits = np.array([12, 94, 210], dtype=float)
    x = np.arange(len(labels))
    w = 0.36
    ax.bar(x - w / 2, hits, width=w, color=OKABE["blue"], zorder=2, label="Genome-wide SV–trait pairs")
    ax.bar(x + w / 2, traits, width=w, color=OKABE["orange"], zorder=2, label="Traits with ≥1 hit")
    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=6.0)
    ax.set_yscale("log")
    ax.set_ylim(5, 1.2e4)
    ax.set_ylabel("Count")
    ax.set_title("Yield versus Phase 1 and versus imputed srWGS", loc="left", pad=2)
    despine(ax)
    ax.legend(loc="upper left", fontsize=5.5)
    for i, (h, t) in enumerate(zip(hits, traits)):
        ax.text(i - w / 2, h * 1.12, f"{int(h)}", ha="center", va="bottom", fontsize=5.8)
        ax.text(i + w / 2, t * 1.12, f"{int(t)}", ha="center", va="bottom", fontsize=5.8)
    ax.text(0.98, 0.06, "Imputed scan: companion manuscript", ha="right", va="bottom", fontsize=5.2, color="#666666", transform=ax.transAxes)


def _mini_locus(
    ax,
    rng: np.random.Generator,
    title: str,
    chrom_lab: str,
    sv_x: float,
    sv_h: float,
    snv_h: float,
    sv_lead: bool = True,
    ylabel: bool = False,
) -> None:
    n = 160
    x = rng.uniform(0.02, 0.98, n)
    u = rng.uniform(1e-3, 1.0, n)
    logp = np.clip(-np.log10(u), 0, 3.2)
    bump = snv_h * np.exp(-0.5 * ((x - sv_x) / 0.09) ** 2)
    logp = np.maximum(logp, bump + rng.normal(0, 0.28, n))
    ax.scatter(x, logp, s=3.5, c="#B8B8B8", linewidths=0, rasterized=True, zorder=2)
    # A few SNVs near the SV so the diamond is not floating in empty space.
    near = np.abs(x - sv_x) < 0.04
    ax.scatter(x[near], logp[near], s=4.5, c="#7A7A7A", linewidths=0, rasterized=True, zorder=3)
    mcolor = OKABE["vermillion"] if sv_lead else OKABE["blue"]
    ax.scatter([sv_x], [sv_h], marker="D", s=26, c=mcolor, edgecolors="white", linewidths=0.35, zorder=5)
    ax.axhline(7.5, color=OKABE["vermillion"], ls="--", lw=0.45, alpha=0.75, zorder=1)
    ymax = max(sv_h, float(logp.max())) + 2.4
    ax.set_xlim(0, 1)
    ax.set_ylim(0, ymax)
    ax.set_xticks([])
    ax.set_title(title, loc="left", pad=1, fontsize=6.3)
    ax.text(0.03, 0.90, chrom_lab, transform=ax.transAxes, fontsize=5.4, color="#666666", va="top")
    if ylabel:
        ax.set_ylabel(r"$-\log_{10}(P)$")
    else:
        ax.set_ylabel("")
    despine(ax)
    tag = "SV lead" if sv_lead else "SNV lead"
    ax.text(0.97, 0.90, tag, transform=ax.transAxes, fontsize=5.3, color=mcolor, ha="right", va="top")


def render() -> tuple[Path, Path]:
    rng = np.random.default_rng(23)
    locis = [
        ("HBA1/2  ·  MCV", "chr16p", 0.28, 14.8, 9.2, True),
        ("LPA KIV-2  ·  Lp(a)", "chr6q", 0.55, 15.4, 10.1, True),
        ("APOL1 tag SV  ·  eGFR", "chr22q", 0.62, 12.6, 8.4, True),
        ("UGT1A1 promoter  ·  bilirubin", "chr2q", 0.48, 11.9, 8.8, True),
        ("NPHS2 del  ·  urine ACR", "chr1q", 0.41, 10.4, 7.1, True),
        ("MHC  ·  rheumatoid arthritis", "chr6p", 0.50, 9.2, 14.6, False),
    ]
    ax_panel("figS13_scale", 1.55, _panel_scale, left=0.08, right=0.98, top=0.86, bottom=0.18)
    ax_panel("figS13_hitmap", 3.35, _panel_hitmap, left=0.10, right=0.98, top=0.88, bottom=0.14)
    ax_panel("figS13_phewas", 3.15, _panel_phewas, width=COL_IN * 1.25, left=0.28, right=0.96, top=0.88, bottom=0.16)
    ax_panel("figS13_scale_vs_prior", 2.85, _panel_scale_vs_prior, width=COL_IN * 1.12, left=0.16, right=0.96, top=0.88, bottom=0.16)

    fig_e = new_figure(4.15)
    gs_e = GridSpec(2, 3, figure=fig_e, left=0.08, right=0.98, top=0.90, bottom=0.12, wspace=0.28, hspace=0.45)
    for i, (title, chrom, sv_x, sv_h, snv_h, sv_lead) in enumerate(locis):
        ax = fig_e.add_subplot(gs_e[i // 3, i % 3])
        _mini_locus(ax, rng, title, chrom, sv_x, sv_h, snv_h, sv_lead=sv_lead, ylabel=(i % 3 == 0))
        if i >= 3:
            ax.set_xlabel("1.5 Mb window")
    save_panel(fig_e, "figS13_example_loci")

    fig = new_figure(9.15)
    gs = GridSpec(
        4,
        1,
        figure=fig,
        left=0.105,
        right=0.975,
        top=0.965,
        bottom=0.045,
        hspace=0.38,
        height_ratios=[0.32, 1.48, 1.10, 1.48],
    )
    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[1, 0])
    gs_cd = gs[2, 0].subgridspec(1, 2, wspace=0.42, width_ratios=[1.18, 0.82])
    ax_c = fig.add_subplot(gs_cd[0, 0])
    ax_d = fig.add_subplot(gs_cd[0, 1])
    gs_e = gs[3, 0].subgridspec(2, 3, wspace=0.22, hspace=0.42)

    _panel_scale(ax_a)
    _panel_hitmap(ax_b)
    _panel_phewas(ax_c)
    _panel_scale_vs_prior(ax_d)

    locis = [
        ("HBA1/2  ·  MCV", "chr16p", 0.28, 14.8, 9.2, True),
        ("LPA KIV-2  ·  Lp(a)", "chr6q", 0.55, 15.4, 10.1, True),
        ("APOL1 tag SV  ·  eGFR", "chr22q", 0.62, 12.6, 8.4, True),
        ("UGT1A1 promoter  ·  bilirubin", "chr2q", 0.48, 11.9, 8.8, True),
        ("NPHS2 del  ·  urine ACR", "chr1q", 0.41, 10.4, 7.1, True),
        ("MHC  ·  rheumatoid arthritis", "chr6p", 0.50, 9.2, 14.6, False),
    ]
    for i, (title, chrom, sv_x, sv_h, snv_h, sv_lead) in enumerate(locis):
        ax = fig.add_subplot(gs_e[i // 3, i % 3])
        _mini_locus(ax, rng, title, chrom, sv_x, sv_h, snv_h, sv_lead=sv_lead, ylabel=(i % 3 == 0))
        if i >= 3:
            ax.set_xlabel("1.5 Mb window")

    stamp_mockup(fig)
    return save(fig, "figS13_assoc_atlas")
