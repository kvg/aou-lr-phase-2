"""Figure S15: SV × phenome Z-score heatmap (pleiotropy catalog).

Layout after common PheWAS heatmaps: variants as columns, grouped traits as
rows, diverging Z, black outlines on phenome-wide hits. Columns here are
selected lead SVs from the native LRS scan, not SNPs.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
from matplotlib.gridspec import GridSpec
from matplotlib.patches import Rectangle

from .style import OKABE, new_figure, save, stamp_mockup

ZMAP = LinearSegmentedColormap.from_list(
    "z",
    [OKABE["blue"], "#D6EAF5", "#FFFFFF", "#F5D6C8", OKABE["vermillion"]],
)
ZMAP.set_bad("#C8C8C8")
ZNORM = TwoSlopeNorm(vmin=-5.0, vcenter=0.0, vmax=5.0)
Z_SIG = 4.4  # |Z| at or above this gets a black outline (placeholder)

# Columns: lead SVs, genomic order, label as in the source-style "id  GENE".
SVS = [
    "DEL 8.4 kb  NPHS2",
    "DUP 155 kb  CFH",
    "CNV  AMY1A",
    "STR  UGT1A1",
    "DEL 3.1 kb  GYPB",
    "Alu  HLA-DRB1",
    "SVA  C4",
    "KIV-2  LPA",
    "INS 312 bp  chr8",
    "DEL  HBB",
    "DEL 1.8 kb  HBA1/2",
    "DUP 42 kb  SMN1",
    "INS 6 kb L1  LDLR",
    "tag SV  APOL1",
    "DEL  NCF1",
    "VNTR  TERT",
]

# Rows: (phenotype, category). Category order is plot order, top to bottom.
PHENOS: list[tuple[str, str]] = [
    ("MCV", "Hematology"),
    ("MCH", "Hematology"),
    ("RBC count", "Hematology"),
    ("Hemoglobin", "Hematology"),
    ("Hematocrit", "Hematology"),
    ("Anemia", "Hematology"),
    ("Thalassemia", "Hematology"),
    ("Neutrophil count", "Hematology"),
    ("WBC count", "Hematology"),
    ("Platelet count", "Hematology"),
    ("Lp(a)", "Lipids"),
    ("LDL-C", "Lipids"),
    ("HDL-C", "Lipids"),
    ("Triglycerides", "Lipids"),
    ("Total cholesterol", "Lipids"),
    ("ApoB", "Lipids"),
    ("eGFR", "Renal"),
    ("Urine ACR", "Renal"),
    ("CKD", "Renal"),
    ("Proteinuria", "Renal"),
    ("ESRD", "Renal"),
    ("Total bilirubin", "Liver"),
    ("Direct bilirubin", "Liver"),
    ("ALT", "Liver"),
    ("AST", "Liver"),
    ("ALP", "Liver"),
    ("Type 2 diabetes", "Metabolic"),
    ("BMI", "Metabolic"),
    ("HbA1c", "Metabolic"),
    ("Fasting glucose", "Metabolic"),
    ("Gout", "Metabolic"),
    ("CAD", "Cardiovascular"),
    ("Myocardial infarction", "Cardiovascular"),
    ("Hypertension", "Cardiovascular"),
    ("Heart failure", "Cardiovascular"),
    ("Stroke", "Cardiovascular"),
    ("Rheumatoid arthritis", "Autoimmune"),
    ("SLE", "Autoimmune"),
    ("IBD", "Autoimmune"),
    ("Psoriasis", "Autoimmune"),
    ("Hypothyroidism", "Autoimmune"),
    ("Asthma", "Other"),
    ("COPD", "Other"),
    ("Serum amylase", "Other"),
    ("Statin use", "Other"),
    ("Gallstones", "Other"),
]

CAT_COLORS = {
    "Hematology": OKABE["vermillion"],
    "Lipids": OKABE["orange"],
    "Renal": OKABE["blue"],
    "Liver": OKABE["green"],
    "Metabolic": "#B8860B",
    "Cardiovascular": OKABE["pink"],
    "Autoimmune": OKABE["sky"],
    "Other": "#4A4A4A",
}

# Planted associations (phenotype, SV label) -> Z. Unsigned noise fills the rest.
PLANTED: dict[tuple[str, str], float] = {
    ("MCV", "DEL 1.8 kb  HBA1/2"): -5.9,
    ("MCH", "DEL 1.8 kb  HBA1/2"): -5.5,
    ("RBC count", "DEL 1.8 kb  HBA1/2"): 5.1,
    ("Hemoglobin", "DEL 1.8 kb  HBA1/2"): -3.6,
    ("Hematocrit", "DEL 1.8 kb  HBA1/2"): -3.1,
    ("Anemia", "DEL 1.8 kb  HBA1/2"): 4.6,
    ("Thalassemia", "DEL 1.8 kb  HBA1/2"): 5.2,
    ("MCV", "DEL  HBB"): -4.8,
    ("MCH", "DEL  HBB"): -4.5,
    ("RBC count", "DEL  HBB"): 4.2,
    ("Anemia", "DEL  HBB"): 3.9,
    ("Neutrophil count", "INS 312 bp  chr8"): 5.0,
    ("WBC count", "INS 312 bp  chr8"): 4.1,
    ("RBC count", "DEL 3.1 kb  GYPB"): 3.7,
    ("Lp(a)", "KIV-2  LPA"): 6.4,
    ("CAD", "KIV-2  LPA"): 4.9,
    ("Myocardial infarction", "KIV-2  LPA"): 4.5,
    ("LDL-C", "KIV-2  LPA"): 2.6,
    ("Statin use", "KIV-2  LPA"): 3.4,
    ("LDL-C", "INS 6 kb L1  LDLR"): 4.7,
    ("Total cholesterol", "INS 6 kb L1  LDLR"): 4.1,
    ("CAD", "INS 6 kb L1  LDLR"): 3.6,
    ("eGFR", "tag SV  APOL1"): -5.4,
    ("CKD", "tag SV  APOL1"): 5.0,
    ("Urine ACR", "tag SV  APOL1"): 4.6,
    ("ESRD", "tag SV  APOL1"): 4.3,
    ("Hypertension", "tag SV  APOL1"): 3.2,
    ("Urine ACR", "DEL 8.4 kb  NPHS2"): 5.1,
    ("Proteinuria", "DEL 8.4 kb  NPHS2"): 4.7,
    ("CKD", "DEL 8.4 kb  NPHS2"): 3.8,
    ("Total bilirubin", "STR  UGT1A1"): 5.6,
    ("Direct bilirubin", "STR  UGT1A1"): 4.8,
    ("Gallstones", "STR  UGT1A1"): 3.3,
    ("Rheumatoid arthritis", "Alu  HLA-DRB1"): 5.3,
    ("SLE", "Alu  HLA-DRB1"): 4.7,
    ("IBD", "Alu  HLA-DRB1"): 3.5,
    ("Psoriasis", "Alu  HLA-DRB1"): 4.0,
    ("Hypothyroidism", "Alu  HLA-DRB1"): 3.1,
    ("SLE", "SVA  C4"): 4.9,
    ("Rheumatoid arthritis", "SVA  C4"): 3.4,
    ("CKD", "DUP 155 kb  CFH"): 3.6,
    ("eGFR", "DUP 155 kb  CFH"): -3.2,
    ("Serum amylase", "CNV  AMY1A"): 5.2,
    ("BMI", "CNV  AMY1A"): -2.4,
    ("COPD", "DUP 42 kb  SMN1"): 2.8,
    ("Asthma", "DEL  NCF1"): 3.9,
    ("Rheumatoid arthritis", "DEL  NCF1"): 3.3,
    ("WBC count", "DEL  NCF1"): 2.9,
}


def _category_spans(names: list[str]) -> list[tuple[str, int, int]]:
    spans = []
    start = 0
    current = names[0]
    for i, name in enumerate(names + [None]):
        if name != current:
            spans.append((current, start, i - 1))
            start = i
            current = name
    return spans


def _matrix(rng: np.random.Generator) -> tuple[np.ndarray, np.ndarray]:
    n_p, n_s = len(PHENOS), len(SVS)
    z = rng.normal(0.0, 0.55, size=(n_p, n_s))
    z = np.clip(z, -2.2, 2.2)
    pheno_i = {p: i for i, (p, _) in enumerate(PHENOS)}
    sv_i = {s: i for i, s in enumerate(SVS)}
    for (ph, sv), val in PLANTED.items():
        z[pheno_i[ph], sv_i[sv]] = val
    missing = rng.random((n_p, n_s)) < 0.045
    for (ph, sv) in PLANTED:
        missing[pheno_i[ph], sv_i[sv]] = False
    z_m = np.ma.array(z, mask=missing)
    return z, z_m


def _style_xticklabels(ax) -> None:
    for tick in ax.get_xticklabels():
        tick.set_rotation(90)
        tick.set_ha("left")
        tick.set_va("bottom")
        tick.set_rotation_mode("anchor")
        tick.set_fontsize(5.6)


def render() -> tuple[Path, Path]:
    rng = np.random.default_rng(41)
    z, z_m = _matrix(rng)
    n_p, n_s = z.shape
    cats = [c for _, c in PHENOS]
    spans = _category_spans(cats)

    fig = new_figure(8.85)
    gs = GridSpec(
        1,
        2,
        figure=fig,
        left=0.07,
        right=0.80,
        top=0.78,
        bottom=0.035,
        wspace=0.012,
        width_ratios=[0.055, 1.0],
    )
    ax_cat = fig.add_subplot(gs[0, 0])
    ax = fig.add_subplot(gs[0, 1])

    ax.imshow(z_m, cmap=ZMAP, norm=ZNORM, aspect="auto", interpolation="nearest")
    for i in range(n_p):
        for j in range(n_s):
            if z_m.mask[i, j]:
                continue
            if abs(z[i, j]) >= Z_SIG:
                ax.add_patch(
                    Rectangle(
                        (j - 0.5, i - 0.5),
                        1,
                        1,
                        fill=False,
                        edgecolor="#111111",
                        lw=0.85,
                        zorder=3,
                    )
                )

    ax.set_xticks(range(n_s))
    ax.set_xticklabels(SVS, fontsize=5.6)
    ax.tick_params(axis="x", top=True, labeltop=True, bottom=False, labelbottom=False, length=2.2, pad=2)
    _style_xticklabels(ax)
    ax.set_yticks(range(n_p))
    ax.set_yticklabels([p for p, _ in PHENOS], fontsize=5.7)
    ax.yaxis.tick_right()
    ax.tick_params(axis="y", length=0, pad=3, right=True, labelright=True, left=False, labelleft=False)
    ax.set_xlim(-0.5, n_s - 0.5)
    ax.set_ylim(n_p - 0.5, -0.5)
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_linewidth(0.5)
        spine.set_color("#888888")
    ax.set_title("Lead SVs × EHR traits (placeholder Z from the native LRS scan)", loc="left", pad=8)

    for _, _s, end in spans[:-1]:
        ax.axhline(end + 0.5, color="white", lw=1.15, zorder=4)

    ax_cat.set_xlim(0, 1)
    ax_cat.set_ylim(n_p - 0.5, -0.5)
    ax_cat.axis("off")
    for name, start, end in spans:
        color = CAT_COLORS[name]
        ax_cat.add_patch(
            Rectangle((0.18, start - 0.5), 0.70, end - start + 1, facecolor=color, edgecolor="white", lw=0.6, clip_on=False)
        )
        text_color = "#111111" if name in {"Metabolic", "Lipids"} else "white"
        ax_cat.text(
            0.53,
            (start + end) / 2,
            name,
            rotation=90,
            ha="center",
            va="center",
            fontsize=6.0,
            color=text_color,
        )

    cax = fig.add_axes([0.125, 0.935, 0.22, 0.016])
    cb = fig.colorbar(ax.images[0], cax=cax, orientation="horizontal")
    cb.set_label("Z score", fontsize=6.5, labelpad=2)
    cb.set_ticks([-5, 0, 5])
    cb.ax.tick_params(labelsize=5.8, length=2.0, pad=1)
    cb.ax.xaxis.set_label_position("top")
    cb.ax.xaxis.set_ticks_position("top")
    fig.text(
        0.36,
        0.943,
        "Black outline: phenome-wide significant   ·   Grey: not tested",
        ha="left",
        va="center",
        fontsize=5.8,
        color="#444444",
    )

    stamp_mockup(fig)
    return save(fig, "figS15_phewas_heatmap", panel_stem="figS15_z_heatmap")
