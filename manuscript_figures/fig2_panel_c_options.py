"""Figure 2 panel C options: allelic series at repeat loci, read the way FELIX tests them.

Every option draws from one simulated set of repeat loci. Each haplotype carries a
locus dosage x_h = repeat units relative to REF, summed over the integrated-callset
records in the locus (the quantity `write_admixed_dosage_vcf.py` would test once
records are aggregated per locus). Phase 1 and Phase 2 are drawn from the same
locus truth with their own haplotype ancestry mixes.

    python -m manuscript_figures.fig2_panel_c_options   (from the repo root)
"""

from __future__ import annotations

import csv
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.ticker import MaxNLocator

from .fig2_options import LAI, MANUSCRIPT_ROOT, P1_HAP_ANC_P, P2_HAP_ANC_P
from .style import MIN_FONT, OUTDIR, P1_N, P2_N, despine, new_figure, save, stamp_mockup

FS = MIN_FONT
PANEL_W, PANEL_H = 3.55, 3.2
N_LOCI = 2000
MAC_MIN = 20
H1, H2 = 2 * P1_N, 2 * P2_N

C_P1 = "#7A7A7A"
C_P2_DS = "#8FA5B2"
C_P2 = "#22323F"
P1_FILL, P2_FILL = "#CFCFCF", "#BCCAD3"
C_RARE = "#D9D9D9"
C_TESTABLE = "#4F6878"
STAMP = "simulated repeat loci"
# Written by notebooks/terra/fig2_02_panel_c_repeat_loci.ipynb (pipeline repo).
RAREFACTION_PATH = MANUSCRIPT_ROOT / "data" / "fig2_panel_c_rarefaction.tsv"
EXAMPLES_PATH = MANUSCRIPT_ROOT / "data" / "fig2_panel_c_examples.tsv"
ALLELE_GAIN_PATH = MANUSCRIPT_ROOT / "data" / "fig2_panel_c_allele_gain.json"
ALTS_PATH = MANUSCRIPT_ROOT / "data" / "fig2_panel_c_alts.json"
# chr20:10–20 Mb, 200 people. Exact match and match within one motif, by the
# larger of the two TRGT haplotype length changes. From the TrgtLocusConcordance summary.
CONCORDANCE = (
    ("0", 99.84, 99.98),
    ("1–19", 58.08, 86.30),
    ("20–49", 20.49, 25.99),
    ("≥50", 8.10, 10.44),
)


def simulate(rng: np.random.Generator, n_loci: int = N_LOCI) -> dict:
    """Haplotype dosages per locus for Phase 1, Phase 2 and a Phase 2 order for rarefaction."""
    p1 = np.array([P1_HAP_ANC_P[a] for a in LAI])
    p2 = np.array([P2_HAP_ANC_P[a] for a in LAI])
    anc1 = rng.choice(len(LAI), size=H1, p=p1)
    anc2 = rng.choice(len(LAI), size=H2, p=p2)
    loci = []
    for _ in range(n_loci):
        sd = float(np.clip(rng.lognormal(np.log(0.55), 0.9), 0.12, 10.0))
        offset = rng.normal(0, 0.5 * sd)
        shift = rng.normal(0, 0.45 * sd, len(LAI))
        spread = sd * rng.uniform(0.85, 1.3, len(LAI))
        tail = rng.uniform(0.0005, 0.006) if rng.random() < 0.35 else 0.0

        def draw(anc):
            x = rng.normal(offset + shift[anc], spread[anc])
            if tail:
                hit = rng.random(anc.size) < tail
                x[hit] += rng.geometric(1 / (3 * sd + 2), hit.sum())
            return np.round(x).astype(np.int32)

        loci.append((draw(anc1), draw(anc2), sd))
    return {"loci": loci, "order": rng.permutation(H2)}


def alleles_per_locus(sim: dict) -> dict[str, np.ndarray]:
    sub = sim["order"][:H1]
    n1 = np.array([np.unique(u1).size for u1, _u2, _sd in sim["loci"]])
    n2 = np.array([np.unique(u2).size for _u1, u2, _sd in sim["loci"]])
    n2_ds = np.array([np.unique(u2[sub]).size for _u1, u2, _sd in sim["loci"]])
    return {"p1": n1, "p2_ds": n2_ds, "p2": n2}


def new_allele_share(sim: dict) -> float:
    new = total = 0
    for u1, u2, _sd in sim["loci"]:
        a2 = np.unique(u2)
        new += int((~np.isin(a2, u1)).sum())
        total += a2.size
    return new / total


def rarefaction(sim: dict, grid1: np.ndarray, grid2: np.ndarray) -> dict:
    def curve(get, grid):
        m = np.empty((len(sim["loci"]), grid.size))
        for i, loc in enumerate(sim["loci"]):
            u = get(loc)
            first = np.zeros(u.size, bool)
            _, idx = np.unique(u, return_index=True)
            first[idx] = True
            m[i] = np.cumsum(first)[grid - 1]
        return m.mean(0), np.percentile(m, 25, 0), np.percentile(m, 75, 0)

    order = sim["order"]
    return {
        "p1": curve(lambda loc: loc[0], grid1),
        "p2": curve(lambda loc: loc[1][order], grid2),
    }


def encoding_stats(sim: dict) -> dict:
    """Carrier haplotypes on alleles testable alone (MAC >= MAC_MIN) vs dosage."""
    out = {}
    for key, idx in (("p1", 0), ("p2", 1)):
        carriers = rare = n_alleles = rare_alleles = 0
        tests, dosage_only, dosage_ok = [], 0, 0
        for loc in sim["loci"]:
            u = loc[idx]
            _vals, cnt = np.unique(u[u != 0], return_counts=True)
            carriers += int(cnt.sum())
            rare += int(cnt[cnt < MAC_MIN].sum())
            n_alleles += cnt.size
            rare_alleles += int((cnt < MAC_MIN).sum())
            k = int((cnt >= MAC_MIN).sum())
            tests.append(k)
            if cnt.sum() >= MAC_MIN:
                dosage_ok += 1
                dosage_only += k == 0
        out[key] = {
            "rare_share": rare / carriers,
            "rare_allele_share": rare_alleles / n_alleles,
            "tests": np.array(tests),
            "dosage_ok": dosage_ok,
            "dosage_only": dosage_only,
        }
    return out


def _legend(ax, handles, **kw) -> None:
    opts = dict(frameon=False, fontsize=FS, handlelength=1.4, handletextpad=0.4, labelspacing=0.3, borderaxespad=0)
    ax.legend(handles=handles, **{**opts, **kw})


def draw_alleles_hist(ax, apl: dict, new_share: float) -> None:
    cap = 25
    edges = np.arange(0.5, cap + 1.5)
    series = (
        ("p1", C_P1, "-", f"Phase 1 (n = {P1_N:,})"),
        ("p2_ds", C_P2_DS, (0, (2.5, 1.2)), f"Phase 2, subsampled to {P1_N:,}"),
        ("p2", C_P2, "-", f"Phase 2 (n = {P2_N:,})"),
    )
    top = 0.0
    for key, color, ls, label in series:
        h, _ = np.histogram(np.minimum(apl[key], cap), bins=edges)
        frac = h / h.sum()
        top = max(top, frac.max())
        ax.stairs(frac, edges, color=color, lw=1.1, ls=ls, label=label)
    for key, color, *_ in series:
        med = float(np.median(apl[key]))
        ax.plot(med, top * 1.08, marker="v", ms=4, color=color, mec="none", clip_on=False)
    ax.set_xlim(0.5, cap + 0.5)
    ax.set_ylim(0, top * 1.15)
    ax.set_xticks([1, 5, 10, 15, 20, 25])
    ax.set_xticklabels(["1", "5", "10", "15", "20", f"≥{cap}"])
    ax.set_xlabel("Distinct alleles per repeat locus")
    ax.set_ylabel("Fraction of repeat loci")
    ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda v, _p: f"{v:.2f}".rstrip("0").rstrip(".")))
    _legend(ax, ax.get_legend_handles_labels()[0], loc="upper right")
    ax.text(0.98, 0.70, f"▼ median\n{new_share:.0%} of Phase 2 alleles\nare absent from Phase 1", transform=ax.transAxes, ha="right", va="top", fontsize=FS, color="#444444")
    despine(ax)


def draw_series_tail(ax) -> None:
    """Phase 2 histogram of distinct alleles per locus. The last bin is the tail."""
    d = json.loads(ALLELE_GAIN_PATH.read_text())
    last = 40
    h = d["h_n2"]
    total = sum(h)
    y = np.array([100.0 * h[k] / total for k in range(1, last)] + [100.0 * sum(h[last:]) / total])
    x = np.arange(1, last + 1)
    ax.bar(x, y, width=1.0, color=P2_FILL, lw=0, align="center", zorder=2)
    ax.stairs(y, np.arange(0.5, last + 1.5), color=C_P2, lw=0.7, zorder=3)
    ax.set_xlim(0.4, last + 0.6)
    ax.set_xticks([1, 10, 20, 30, last])
    ax.set_xticklabels(["1", "10", "20", "30", f"≥{last}"])
    ax.set_ylabel("Loci (%)")
    ax.set_xlabel("Distinct alleles per locus")
    _legend(
        ax,
        [Patch(facecolor=P2_FILL, edgecolor=C_P2, linewidth=0.7, label="Phase 2")],
        loc="upper right",
        handlelength=0.8,
        borderaxespad=0.1,
    )
    despine(ax)


def _load_alts() -> dict:
    return json.loads(ALTS_PATH.read_text())


# Fine carrier bins from the scan, collapsed so the labels fit this slot.
_CARRIER_SLICES = (slice(0, 1), slice(1, 4), slice(4, 5), slice(5, 7), slice(7, 8), slice(8, 10))
_CARRIER_LABELS = ["1", "2–9", "10–19", "20–99", "100–999", "≥1,000"]


def _carrier_counts(d: dict) -> tuple[np.ndarray, np.ndarray]:
    seen = d["alleles_seen_in_phase1"]
    new = d["alleles_new_in_phase2"]
    return (
        np.array([sum(seen[s]) for s in _CARRIER_SLICES], float),
        np.array([sum(new[s]) for s in _CARRIER_SLICES], float),
    )


def draw_carrier_spectrum(ax) -> None:
    """How many haplotypes carry each Phase 2 length. The cut is 20 haplotypes."""
    seen, new = _carrier_counts(_load_alts())
    y = 100.0 * (seen + new) / (seen + new).sum()
    x = np.arange(y.size)
    colors = [C_RARE] * 3 + [C_TESTABLE] * 3
    edges = ["#8A8A8A"] * 3 + [C_P2] * 3
    ax.bar(x, y, width=0.86, color=colors, edgecolor=edges, lw=0.6, zorder=2)
    ax.set_xticks(x)
    ax.set_xticklabels(_CARRIER_LABELS)
    ax.set_xlim(-0.6, y.size - 0.4)
    ax.set_ylabel("Alleles (%)")
    ax.set_xlabel("Haplotypes carrying the allele")
    _legend(
        ax,
        [
            Patch(facecolor=C_RARE, edgecolor="#8A8A8A", linewidth=0.6, label="< 20 haplotypes"),
            Patch(facecolor=C_TESTABLE, edgecolor=C_P2, linewidth=0.6, label="≥ 20 haplotypes"),
        ],
        loc="lower right",
        bbox_to_anchor=(1.0, 1.02),
        ncol=2,
        handlelength=0.8,
        borderaxespad=0,
        columnspacing=0.8,
    )
    despine(ax)


def draw_new_by_frequency(ax, *, shared_legend: bool = False) -> None:
    """Phase 2 lengths Phase 1 already had, against lengths it did not, by carrier count.

    With ``shared_legend`` the colors match panel B and the example rows below
    (gray, a Phase 2 length Phase 1 also had; slate, new) and one legend serves both.
    """
    seen, new = _carrier_counts(_load_alts())
    total = (seen + new).sum()
    ys, yn = 100.0 * seen / total, 100.0 * new / total
    x = np.arange(ys.size)
    if shared_legend:
        ax.bar(x, ys, width=0.86, color=C_IN_P1_FILL, edgecolor=C_IN_P1, lw=0.6, zorder=2)
        ax.bar(x, yn, width=0.86, bottom=ys, color=C_NEW, edgecolor=C_NEW, lw=0.4, zorder=3)
    else:
        ax.bar(x, ys, width=0.86, color=P1_FILL, edgecolor=C_P1, lw=0.6, zorder=2, label="Seen in Phase 1")
        ax.bar(x, yn, width=0.86, bottom=ys, color=C_P2, edgecolor=C_P2, lw=0.4, zorder=3, label="New in Phase 2")
    ax.set_xticks(x)
    ax.set_xticklabels(_CARRIER_LABELS)
    ax.set_xlim(-0.6, ys.size - 0.4)
    ax.set_ylabel("Phase 2 alleles (%)" if shared_legend else "Alleles (%)")
    ax.set_xlabel("Haplotypes carrying the allele")
    _legend(
        ax,
        legend_handles() if shared_legend else ax.get_legend_handles_labels()[0],
        loc="lower right",
        bbox_to_anchor=(1.0, 1.02),
        ncol=3 if shared_legend else 2,
        handlelength=0.8,
        borderaxespad=0,
        columnspacing=0.8,
    )
    despine(ax)


def draw_plvi_contrast(ax) -> None:
    """Share of loci reaching each allele count. Top 1% of PLVI within a motif group, against the other ranked loci."""
    d = _load_alts()
    cuts = (5, 10, 20, 40)
    labels = [f"≥{k}" for k in cuts]

    def tails(h: list[int]) -> np.ndarray:
        total = sum(h)
        return np.array([100.0 * sum(h[k:]) / total for k in cuts])

    rest = tails(d["alleles_per_locus_rest"])
    high = tails(d["alleles_per_locus_high_plvi"])
    x = np.arange(len(labels))
    w = 0.38
    ax.bar(x - w / 2, rest, width=w, color=P2_FILL, edgecolor=C_P2, lw=0.5, label="Other loci", zorder=2)
    ax.bar(x + w / 2, high, width=w, color=C_P2, edgecolor=C_P2, lw=0.4, label="Top 1% PLVI", zorder=3)
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.set_xlim(-0.55, len(labels) - 0.45)
    ax.set_ylim(0, 118)
    ax.set_ylabel("Loci (%)")
    ax.set_xlabel("Distinct alleles")
    _legend(ax, ax.get_legend_handles_labels()[0], loc="upper right", handlelength=0.8, borderaxespad=0.1)
    despine(ax)


def draw_concordance(ax) -> None:
    """Integrated genotypes against TRGT, split by how far TRGT moved from the reference length."""
    labels = [row[0] for row in CONCORDANCE]
    exact = np.array([row[1] for row in CONCORDANCE])
    within = np.array([row[2] for row in CONCORDANCE])
    x = np.arange(len(labels))
    w = 0.38
    ax.bar(x - w / 2, exact, width=w, color=C_P2, edgecolor=C_P2, lw=0.4, label="Exact", zorder=3)
    ax.bar(x + w / 2, within, width=w, color=P2_FILL, edgecolor=C_P2, lw=0.5, label="Within one motif", zorder=2)
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.set_xlim(-0.6, len(labels) - 0.4)
    ax.set_ylim(0, 112)
    ax.set_ylabel("Genotypes (%)")
    ax.set_xlabel("TRGT length change (bp), chr20, 200 people")
    _legend(ax, ax.get_legend_handles_labels()[0], loc="upper right", handlelength=0.8, borderaxespad=0.1)
    despine(ax)


def draw_rarefaction(ax, rf: dict, grid1: np.ndarray, grid2: np.ndarray, *, compact: bool = False, phase1_size: float | None = None) -> None:
    for key, grid, color, band in (("p2", grid2, C_P2, P2_FILL), ("p1", grid1, C_P1, P1_FILL)):
        mean, lo, hi = rf[key]
        ax.fill_between(grid, lo, hi, color=band, alpha=0.55, lw=0)
        ax.plot(grid, mean, color=color, lw=1.2)
    # The mean sits above the interquartile band when a few loci have very many alleles.
    hi = max(float(np.max(rf[k][i])) for k in ("p1", "p2") for i in (0, 2))
    y_top = hi * 1.18
    x1 = float(grid1[-1] if phase1_size is None else phase1_size)
    ax.axvline(x1, color="#999999", lw=0.6, ls=(0, (2, 1.5)), zorder=0)
    # The interquartile band reaches the bottom of the axes, so the marker label sits at the top.
    ax.annotate("Phase 1 size", xy=(x1, y_top), xytext=(3, -1), textcoords="offset points", fontsize=FS, color="#777777", ha="left", va="top")
    m1, m2 = float(rf["p1"][0][-1]), float(rf["p2"][0][-1])
    m2_at = float(np.interp(np.log(x1), np.log(grid2), rf["p2"][0]))
    ax.text(grid1[0] * 1.1, float(rf["p1"][0][0]) + 0.06 * y_top, "Phase 1", fontsize=FS, color=C_P1, ha="left", va="bottom")
    ax.text(grid2[-1], m2 - 0.04 * y_top, "Phase 2", fontsize=FS, color=C_P2, ha="right", va="top")
    ax.plot([x1], [m2_at], marker="o", ms=2.8, color=C_P2, zorder=4)
    ax.plot([grid1[-1]], [m1], marker="o", ms=2.8, color=C_P1, zorder=4)
    ax.set_xscale("log")
    ax.set_xlim(grid1[0], grid2[-1] * 1.05)
    ax.set_ylim(0, y_top)
    ticks = [100, 1000, 10000]
    ax.set_xticks(ticks)
    ax.set_xticklabels(["100", "1,000", "10,000"])
    ax.xaxis.set_minor_locator(plt.NullLocator())
    ax.set_xlabel("Haplotypes sampled")
    ax.set_ylabel("Alleles per locus" if compact else "Mean distinct alleles\nper repeat locus")
    if not compact:
        ax.text(0.02, 0.98, f"Bands: interquartile range across loci\nAt Phase 1 size: Phase 2 {m2_at:.1f}, Phase 1 {m1:.1f}\nFull Phase 2: {m2:.1f}", transform=ax.transAxes, fontsize=FS, color="#444444", va="top")
    despine(ax)


def draw_encoding(fig, cell, enc: dict) -> plt.Axes:
    sub = cell.subgridspec(2, 1, height_ratios=[0.55, 1.0], hspace=1.05)
    ax = fig.add_subplot(sub[0])
    rows = (("p2", "Phase 2"), ("p1", "Phase 1"))
    for y, (key, label) in enumerate(rows):
        r = enc[key]["rare_allele_share"]
        ax.barh(y, 1 - r, color=C_TESTABLE, height=0.62, lw=0)
        ax.barh(y, r, left=1 - r, color=C_RARE, height=0.62, lw=0)
        ax.text((1 - r) / 2, y, f"{1 - r:.0%}", ha="center", va="center", fontsize=FS, color="white")
        if r >= 0.12:
            ax.text(1 - r / 2, y, f"{r:.0%}", ha="center", va="center", fontsize=FS, color="#333333")
        else:
            ax.text(1.02, y, f"{r:.0%}", ha="left", va="center", fontsize=FS, color="#333333")
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels([lab for _k, lab in rows])
    ax.tick_params(axis="y", length=0)
    ax.set_xlim(0, 1)
    ax.set_xticks([0, 0.5, 1])
    ax.set_xticklabels(["0", "50%", "100%"])
    c1, c2 = enc["p1"]["rare_share"], enc["p2"]["rare_share"]
    ax.set_xlabel(f"Distinct non-reference alleles at repeat loci\n(they carry {c1:.0%} / {c2:.0%} of carriers; dosage tests all)")
    despine(ax, left=False)
    _legend(
        ax,
        [
            Patch(facecolor=C_TESTABLE, label=f"≥{MAC_MIN} copies (testable alone)"),
            Patch(facecolor=C_RARE, label=f"<{MAC_MIN} copies"),
        ],
        loc="lower left",
        bbox_to_anchor=(0, 1.04),
        ncol=2,
        columnspacing=0.8,
        handlelength=0.9,
    )
    ax2 = fig.add_subplot(sub[1])
    cap = 10
    edges = np.arange(-0.5, cap + 1.5)
    top = 0.0
    for key, color, label in (("p1", C_P1, "Phase 1"), ("p2", C_P2, "Phase 2")):
        h, _ = np.histogram(np.minimum(enc[key]["tests"], cap), bins=edges)
        frac = h / h.sum()
        top = max(top, frac.max())
        ax2.stairs(frac, edges, color=color, lw=1.1, label=label)
    ax2.set_ylim(0, top * 1.25)
    ax2.axvline(1, color="#B03A1A", lw=0.8, ls=(0, (2.5, 1.5)))
    ax2.text(1.2, top * 1.22, "Dosage: 1 test", color="#B03A1A", fontsize=FS, va="top")
    ax2.set_xlim(-0.5, cap + 0.5)
    ax2.set_xticks([0, 1, 2, 4, 6, 8, 10])
    ax2.set_xticklabels(["0", "1", "2", "4", "6", "8", f"≥{cap}"])
    ax2.set_xlabel(f"Biallelic tests per locus (alleles with ≥{MAC_MIN} copies)")
    ax2.set_ylabel("Fraction of loci")
    d1, d2 = enc["p1"]["dosage_only"], enc["p2"]["dosage_only"]
    ax2.text(0.98, 0.60, f"Loci testable only as dosage\n(0 alleles with ≥{MAC_MIN} copies):\nPhase 1 {d1:,}, Phase 2 {d2:,}", transform=ax2.transAxes, ha="right", va="top", fontsize=FS, color="#444444")
    _legend(ax2, ax2.get_legend_handles_labels()[0], loc="upper right")
    despine(ax2)
    return ax


def pick_examples(sim: dict, k: int = 4) -> list[dict]:
    """Two known pathogenic loci, then two high-PLVI stand-ins. Wider shapes sit lower."""
    # 135 is the rounded PLVI of the AAAAG intron in IQCB1 (measured 135.2).
    # 125 matches the width of another high value (the group-6 winner, 125.2) so the
    # second annotation can be judged; it is not a measurement at EP400.
    labels = (("FMR1", None), ("HTT", None), ("IQCB1", 135), ("EP400", 125))
    sds = np.array([sd for *_u, sd in sim["loci"]])
    order = np.argsort(sds)
    out = []
    for i, q in enumerate((0.45, 0.65, 0.82, 0.95)[:k]):
        u1, u2, _sd = sim["loci"][int(order[int(q * (order.size - 1))])]
        a1, c1 = np.unique(u1, return_counts=True)
        a2, c2 = np.unique(u2, return_counts=True)
        label, plvi = labels[i]
        out.append({
            "label": label,
            "italic": True,
            "plvi": plvi,
            "a1": a1, "c1": c1, "a2": a2, "c2": c2,
        })
    return out


def load_examples(path: Path = EXAMPLES_PATH) -> list[dict] | None:
    """PLVI-ranked loci from panel_c_repeat_alleles.py: allele (repeat units) and haplotype count per phase."""
    if not path.is_file():
        return None
    rows: dict[int, dict] = {}
    with path.open() as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            ex = rows.setdefault(int(r["rank"]), {
                "label": f"{r['chrom']}:{int(r['start']) + 1:,}\n({r['motif']})n",
                "phase1": ([], []), "phase2": ([], []),
            })
            vals, cnts = ex[r["phase"]]
            vals.append(float(r["dosage_units"]))
            cnts.append(int(r["n_hap_allele"]))
    out = []
    for _rank, ex in sorted(rows.items()):
        (a1, c1), (a2, c2) = ex["phase1"], ex["phase2"]
        out.append({"label": ex["label"], "a1": np.array(a1), "c1": np.array(c1), "a2": np.array(a2), "c2": np.array(c2)})
    return out


def load_phase1_size(path: Path = RAREFACTION_PATH) -> float | None:
    """Cohort haplotype count recorded next to the rarefaction table, when that run wrote one."""
    summary = path.with_name("fig2_panel_c.summary.json")
    if not summary.is_file():
        return None
    return float(json.loads(summary.read_text())["haplotypes"]["phase1"])


def load_rarefaction(path: Path = RAREFACTION_PATH) -> tuple[dict, np.ndarray, np.ndarray] | None:
    if not path.is_file():
        return None
    cols: dict[str, list[list[float]]] = {"phase1": [], "phase2": []}
    with path.open() as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            cols[r["phase"]].append([float(r[k]) for k in ("n_hap", "mean", "q25", "q75")])
    a1, a2 = np.array(cols["phase1"]), np.array(cols["phase2"])
    rf = {"p1": (a1[:, 1], a1[:, 2], a1[:, 3]), "p2": (a2[:, 1], a2[:, 2], a2[:, 3])}
    return rf, a1[:, 0], a2[:, 0]


def _weighted_quantile(vals: np.ndarray, weights: np.ndarray, q: float) -> float:
    order = np.argsort(vals)
    cum = np.cumsum(weights[order]) / weights.sum()
    return float(vals[order][np.searchsorted(cum, q)])


def mock_panel_examples(rng: np.random.Generator, n: int) -> list[dict]:
    """Layout stand-ins. Known loci carry a rare long tail; high-PLVI loci are broad.

    Four loci drop *C9orf72* and the third high-PLVI row. The high-PLVI names are
    stand-ins for the length of a gene label, not the loci that will be chosen.
    """
    specs = (
        ("FMR1", "known", 2.0, 1.4, 14),
        ("HTT", "known", 0.0, 1.8, 8),
        ("C9orf72", "known", 1.0, 1.5, 10),
        ("IQCB1", "plvi", 0.0, 3.6, 4),
        ("EP400", "plvi", 1.0, 2.8, 3),
        ("FAM193B", "plvi", -0.5, 4.2, 5),
    )
    if n == 4:
        specs = (specs[0], specs[1], specs[3], specs[4])
    else:
        specs = specs[:n]

    def draw(loc, scale, n_hap, n_tail):
        x = np.rint(rng.normal(loc, scale, n_hap)).astype(np.int32)
        if n_tail:
            far = np.rint(loc + 4 * scale + rng.exponential(max(n_tail, 1) / 3, size=n_tail)).astype(np.int32)
            x = np.concatenate([x, far])
        return np.unique(x, return_counts=True)

    out = []
    for label, group, loc, scale, tail in specs:
        a1, c1 = draw(loc, scale, H1, max(2, tail // 10))
        a2, c2 = draw(loc, scale * 1.05, H2, tail)
        out.append({"label": label, "italic": True, "group": group, "a1": a1, "c1": c1, "a2": a2, "c2": c2})
    return out


def draw_example_rows(fig, cell, examples: list[dict]) -> plt.Axes:
    sub = cell.subgridspec(1, 2, width_ratios=[1.0, 0.22], wspace=0.1)
    ax = fig.add_subplot(sub[0])
    axb = fig.add_subplot(sub[1], sharey=ax)
    n = len(examples)
    amp, tick = 0.36, 0.13
    half = 0.0
    for ex in examples:
        vals = np.abs(np.r_[ex["a1"], ex["a2"]])
        half = max(half, _weighted_quantile(vals, np.r_[ex["c1"], ex["c2"]].astype(float), 0.999))
    half = int(np.ceil(half / 5) * 5) + 2
    grid = np.arange(-half, half + 1)

    def density(a, c):
        d = np.zeros(grid.size)
        np.add.at(d, np.clip(np.round(a).astype(int), -half, half) + half, c)
        return d / c.sum()

    counts = []
    for i, ex in enumerate(examples):
        base = n - 1 - i
        d1, d2 = density(ex["a1"], ex["c1"]), density(ex["a2"], ex["c2"])
        s = amp / max(d1.max(), d2.max())
        ax.fill_between(grid, base, base + d1 * s, step="mid", color=P1_FILL, lw=0)
        ax.fill_between(grid, base, base - d2 * s, step="mid", color=P2_FILL, lw=0)
        a1, a2 = ex["a1"], ex["a2"]
        seen = np.isin(a2, a1)
        clip = lambda a: a[np.abs(a) <= half]
        ax.vlines(clip(a1), base, base + tick, color=C_P1, lw=0.5)
        ax.vlines(clip(a2[seen]), base - tick, base, color=C_P2_DS, lw=0.5)
        ax.vlines(clip(a2[~seen]), base - tick, base, color=C_P2, lw=0.7)
        ax.plot([-half, half], [base, base], color="#333333", lw=0.5)
        counts.append((a1.size, a2.size, int(seen.sum())))
        if i == 0:
            ax.text(half, base + 0.06, "Phase 1", ha="right", va="bottom", fontsize=FS, color="#555555")
            ax.text(half, base - 0.06, "Phase 2", ha="right", va="top", fontsize=FS)
    ax.axvline(0, color="#999999", lw=0.5, ls=(0, (2, 1.5)), zorder=0)
    ax.set_yticks([n - 1 - i for i in range(n)])
    ax.set_yticklabels([ex["label"] for ex in examples])
    ax.tick_params(axis="y", length=0)
    pad = ax.yaxis.get_major_ticks()[0].get_pad()
    for i, (t, ex) in enumerate(zip(ax.get_yticklabels(), examples)):
        if ex.get("italic"):
            t.set_style("italic")
        plvi = ex.get("plvi")
        if plvi is None:
            continue
        # Gene name sits on the row center; the PLVI hangs just beneath it,
        # right-aligned to the same pad as the tick label.
        t.set_verticalalignment("bottom")
        base = n - 1 - i
        ax.annotate(
            f"PLVI {plvi:g}",
            xy=(0, base),
            xycoords=ax.get_yaxis_transform(),
            xytext=(-pad, -1),
            textcoords="offset points",
            ha="right",
            va="top",
            fontsize=FS,
            color="#555555",
            annotation_clip=False,
        )
    ax.set_ylim(-0.55, n - 0.45)
    ax.set_xlim(-half - 0.5, half + 0.5)
    ax.set_xlabel("Repeat units relative to reference (tested dosage)")
    despine(ax, left=False)
    ax.spines["left"].set_visible(False)

    xmax = max(c for row in counts for c in row)
    for i, (n1, n2, n2s) in enumerate(counts):
        base = n - 1 - i
        axb.barh(base + 0.15, n1, height=0.2, color=C_P1, lw=0)
        axb.barh(base - 0.15, n2s, height=0.2, color=C_P2_DS, lw=0)
        axb.barh(base - 0.15, n2 - n2s, left=n2s, height=0.2, color=C_P2, lw=0)
        axb.text(n1 + 0.05 * xmax, base + 0.15, f"{n1}", va="center", fontsize=FS, color="#555555")
        axb.text(n2 + 0.05 * xmax, base - 0.15, f"{n2}", va="center", fontsize=FS)
    axb.set_xlim(0, xmax * 1.6)
    axb.axis("off")
    axb.text(0, n - 0.45, "Alleles", fontsize=FS, va="bottom")
    tick_kw = dict(marker="|", ls="none", ms=6, mew=1.0)
    _legend(
        ax,
        [
            Line2D([], [], color=C_P2_DS, label="Allele seen in Phase 1", **tick_kw),
            Line2D([], [], color=C_P2, label="New in Phase 2", **tick_kw),
        ],
        loc="lower left",
        bbox_to_anchor=(0, 1.0),
        ncol=2,
        columnspacing=0.9,
        handlelength=0.8,
    )
    return ax


def load_chosen_loci(path: Path | None = None) -> list[dict]:
    """Real panel C rows: own-scale repeat units, with a pathogenic line where one exists."""
    src = path or (MANUSCRIPT_ROOT / "data" / "fig2_panel_c_chosen_loci.json")
    out = []
    for ex in json.loads(src.read_text())["loci"]:
        out.append({
            "label": ex["label"],
            "italic": ex["italic"],
            "plvi": ex["plvi"],
            "threshold": ex["threshold"],
            "threshold_label": {"FMR1": ">200 CGG", "TCF4": "≥51 CTG", "DMPK": "≥50 CTG"}.get(ex["label"]),
            "outlier": ex.get("outlier"),
            "a1": np.asarray(ex["u1"], float),
            "c1": np.asarray(ex["n1"], float),
            "a2": np.asarray(ex["u2"], float),
            "c2": np.asarray(ex["n2"], float),
            "new2": np.asarray(ex["new2"], bool),
        })
    return out


# Repeat units are linear out to LINTHRESH, then logarithmic; the two
# regions get equal width per decade so bins stay the same drawn width.
LINTHRESH = 10.0


def _symlog_edges(lo: float, hi: float) -> np.ndarray:
    ratio = 10 ** 0.1
    inner = np.arange(-LINTHRESH - 0.5, LINTHRESH + 1.0, 1.0)
    up = [inner[-1]]
    while up[-1] < hi:
        up.append(up[-1] * ratio)
    down = [inner[0]]
    while down[-1] > lo:
        down.append(down[-1] * ratio)
    return np.r_[down[:0:-1], inner, up[1:]]


# Same classes as panel B's discovery bands.
C_IN_P1_FILL = "#D4D4D4"
C_IN_P1 = "#A6A6A6"
C_NEW = "#4F6878"
# Mark height runs from a singleton to 1,000 haplotypes; commoner lengths stop at the top.
MARK_MIN, MARK_MAX, MARK_TOP = 0.05, 0.35, 1000
SCALE_PAD_PT = 21


def legend_handles() -> list:
    return [
        Patch(facecolor=C_IN_P1_FILL, edgecolor=C_IN_P1, linewidth=0.6, label="In Phase 1"),
        Patch(facecolor=C_NEW, edgecolor=C_NEW, linewidth=0.6, label="New in Phase 2"),
    ]


def _mark_height(count: np.ndarray) -> np.ndarray:
    frac = np.log10(np.maximum(count, 1)) / np.log10(MARK_TOP)
    return MARK_MIN + (MARK_MAX - MARK_MIN) * np.clip(frac, 0, 1)


def draw_symlog_rows(fig, cell, examples: list[dict]) -> plt.Axes:
    """All rows on one axis: repeat units from the reference, linear near it and log in the tails."""
    sub = cell.subgridspec(1, 2, width_ratios=[1.0, 0.24], wspace=0.06)
    ax = fig.add_subplot(sub[0])
    axb = fig.add_subplot(sub[1], sharey=ax)
    n = len(examples)
    amp = MARK_MAX
    lo = min(float(min(ex["a1"].min(), ex["a2"].min())) for ex in examples)
    hi = max(float(max(ex["a1"].max(), ex["a2"].max())) for ex in examples)
    edges = _symlog_edges(lo, hi)
    # Room past the longest allele for the Phase 1 / Phase 2 labels.
    x_right = edges[-1] * 2.0
    bar_max = max(max(ex["a1"].size, int(ex["new2"].sum())) for ex in examples)
    for i, ex in enumerate(examples):
        base = n - 1 - i
        a1, c1, a2, c2 = ex["a1"], ex["c1"], ex["a2"], ex["c2"]
        new = ex["new2"]
        ax.vlines(a1, base, base + _mark_height(c1), color=C_IN_P1, lw=0.5, zorder=3)
        ax.vlines(a2[~new], base - _mark_height(c2[~new]), base, color=C_IN_P1, lw=0.5, zorder=3)
        ax.vlines(a2[new], base - _mark_height(c2[new]), base, color=C_NEW, lw=0.6, zorder=4)
        ax.plot([edges[0], x_right], [base, base], color="#333333", lw=0.5, zorder=2)
        threshold = ex.get("threshold")
        if threshold is not None:
            ax.plot([threshold, threshold], [base - amp, base + amp], color="#B03A1A", lw=0.8, zorder=5)
            ax.text(threshold / 1.08, base + amp, ex.get("threshold_label") or "", color="#B03A1A",
                    fontsize=FS, ha="right", va="top")
        outlier = ex.get("outlier")
        if outlier is not None:
            ax.plot([outlier, outlier], [base - amp, base + amp], color="#555555", lw=0.8,
                    ls=(0, (2.5, 1.5)), zorder=5)
            ax.text(outlier / 1.08, base + amp, "99th pct", color="#555555",
                    fontsize=FS, ha="right", va="top")
        rows = (
            (base + 0.15, a1.size, C_IN_P1_FILL, C_IN_P1, f"{a1.size:,}"),
            (base - 0.15, int(new.sum()), C_NEW, C_NEW, f"+{int(new.sum()):,}"),
        )
        for y, count, face, edge, text in rows:
            axb.barh(y, count, height=0.24, color=face, edgecolor=edge, lw=0.5)
            axb.text(count + 0.05 * bar_max, y, text, va="center", fontsize=FS, color="#333333")
        # Haplotype scale at the left edge of each row, mirrored like the marks.
        yt = ax.get_yaxis_transform()
        ax.plot([0, 0], [base - amp, base + amp], color="#333333", lw=0.5, transform=yt, clip_on=False)
        for k in (0, 1, 2, 3):
            hk = float(_mark_height(np.array([10 ** k]))[0])
            for y in (base + hk, base - hk):
                ax.annotate("", xy=(0, y), xycoords=yt, xytext=(-2, 0), textcoords="offset points",
                            arrowprops=dict(arrowstyle="-", lw=0.5, color="#333333", shrinkA=0, shrinkB=0))
            if k == 3:
                for y in (base + hk, base - hk):
                    ax.annotate("$10^3$", xy=(0, y), xycoords=yt, xytext=(-3, 0), textcoords="offset points",
                                fontsize=FS, ha="right", va="center", annotation_clip=False)
        ax.annotate("1", xy=(0, base), xycoords=yt, xytext=(-3, 0), textcoords="offset points",
                    fontsize=FS, ha="right", va="center", annotation_clip=False)
        if i == 0:
            right = x_right
            ax.text(right, base + 0.06, "Phase 1", ha="right", va="bottom", fontsize=FS, color="#555555")
            ax.text(right, base - 0.06, "Phase 2", ha="right", va="top", fontsize=FS)
    ax.axvline(0, color="#999999", lw=0.5, ls=(0, (2, 1.5)), zorder=0)
    ax.set_xscale("symlog", linthresh=LINTHRESH, linscale=1.0)
    ax.set_xlim(edges[0], x_right)
    ticks = [t for t in (-100, -10, 0, 10, 100, 1000) if edges[0] <= t <= x_right]
    ax.set_xticks(ticks)
    ax.set_xticklabels([f"{t:,}".replace("-", "−") for t in ticks])
    ax.xaxis.set_minor_locator(plt.NullLocator())
    ax.set_xlabel("Repeat units from reference")
    ax.set_yticks([n - 1 - i for i in range(n)])
    ax.set_yticklabels([ex["label"] for ex in examples])
    ax.tick_params(axis="y", length=0, pad=SCALE_PAD_PT)
    pad = SCALE_PAD_PT
    for i, (t, ex) in enumerate(zip(ax.get_yticklabels(), examples)):
        if ex.get("italic"):
            t.set_style("italic")
        plvi = ex.get("plvi")
        if plvi is None:
            continue
        t.set_verticalalignment("bottom")
        ax.annotate(
            f"PLVI {plvi:.0f}",
            xy=(0, n - 1 - i),
            xycoords=ax.get_yaxis_transform(),
            xytext=(-pad, -1),
            textcoords="offset points",
            ha="right",
            va="top",
            fontsize=FS,
            color="#555555",
            annotation_clip=False,
        )
    ax.set_ylim(-0.55, n - 0.45)
    despine(ax, left=False)
    ax.spines["left"].set_visible(False)
    axb.set_xlim(0, bar_max * 1.55)
    axb.axis("off")
    axb.text(0, n - 0.45, "Alleles", fontsize=FS, va="bottom")
    # Narrow the rows from the left so the gene and PLVI labels clear the haplotype scale.
    box = ax.get_position()
    inset = 8 / (fig.get_figwidth() * 72)
    ax.set_position([box.x0 + inset, box.y0, box.width - inset, box.height])
    ax.annotate("Haplotypes", xy=(0, n - 0.45), xycoords=ax.get_yaxis_transform(), xytext=(0, 1),
                textcoords="offset points", fontsize=FS, ha="center", va="bottom", annotation_clip=False)
    return ax


def panel_rarefaction_examples(
    fig, cell, rng: np.random.Generator, examples: list[dict] | None = None, ratios: tuple[float, float] = (0.55, 1.50),
    summary: str = "new", row_scales: str = "shared",
) -> tuple[plt.Axes, plt.Axes]:
    """A summary above the example loci, or the example loci alone.

    Example rows stay simulated until a list is passed in. ``summary`` is
    ``new`` (Phase 2 lengths by carrier count, seen in Phase 1 against new),
    ``carriers``, ``alleles``, ``plvi``, ``concordance``, or ``none``.
    """
    needs_curve = summary == "alleles" and not ALLELE_GAIN_PATH.is_file()
    real_rf = load_rarefaction() if needs_curve else None
    sim = None if examples is not None else simulate(rng)
    chosen = examples if examples is not None else pick_examples(sim)
    draw_rows = draw_symlog_rows if row_scales == "symlog" else draw_example_rows
    if summary == "none":
        ax = draw_rows(fig, cell, chosen)
        return ax, ax
    # The shared legend sits on the summary, so the rows need no room above them.
    gs = cell.subgridspec(2, 1, height_ratios=list(ratios), hspace=0.40 if row_scales == "symlog" else 0.62)
    ax = fig.add_subplot(gs[0])
    if summary == "alleles":
        if ALLELE_GAIN_PATH.is_file():
            draw_series_tail(ax)
        else:
            if real_rf is not None:
                rf, grid1, grid2 = real_rf
            else:
                grid1 = np.unique(np.geomspace(50, H1, 40).astype(int))
                grid2 = np.unique(np.geomspace(50, H2, 60).astype(int))
                rf = rarefaction(sim if sim is not None else simulate(rng), grid1, grid2)
            draw_rarefaction(ax, rf, grid1, grid2, compact=True, phase1_size=load_phase1_size() if real_rf is not None else None)
    elif summary == "carriers":
        draw_carrier_spectrum(ax)
    elif summary == "new":
        draw_new_by_frequency(ax, shared_legend=row_scales == "symlog")
    elif summary == "plvi":
        draw_plvi_contrast(ax)
    elif summary == "concordance":
        draw_concordance(ax)
    else:
        raise ValueError(summary)
    ax_rows = draw_rows(fig, gs[1], chosen)
    return ax, ax_rows


def _panel(stem: str, draw) -> Path:
    fig = new_figure(PANEL_H, PANEL_W)
    draw(fig)
    stamp_mockup(fig, extra=STAMP)
    return save(fig, stem)[1]


def render(seed: int = 7) -> list[Path]:
    rng = np.random.default_rng(seed)
    sim = simulate(rng)
    apl = alleles_per_locus(sim)
    share = new_allele_share(sim)
    grid1 = np.unique(np.geomspace(50, H1, 40).astype(int))
    grid2 = np.unique(np.geomspace(50, H2, 60).astype(int))
    rf = rarefaction(sim, grid1, grid2)
    enc = encoding_stats(sim)

    def c1(fig):
        ax = fig.add_axes([0.16, 0.16, 0.8, 0.76])
        draw_alleles_hist(ax, apl, share)

    def c2(fig):
        ax = fig.add_axes([0.18, 0.16, 0.78, 0.76])
        draw_rarefaction(ax, rf, grid1, grid2)

    def c3(fig):
        gs = fig.add_gridspec(1, 1, left=0.16, right=0.96, top=0.86, bottom=0.14)
        draw_encoding(fig, gs[0], enc)

    def c4(fig):
        gs = fig.add_gridspec(2, 1, left=0.16, right=0.97, top=0.95, bottom=0.15, height_ratios=[0.8, 1.25], hspace=0.62)
        draw_rarefaction(fig.add_subplot(gs[0]), rf, grid1, grid2, compact=True)
        draw_example_rows(fig, gs[1], pick_examples(sim))

    pngs = [
        _panel("fig2c_opt1_alleles_per_locus", c1),
        _panel("fig2c_opt2_rarefaction", c2),
        _panel("fig2c_opt3_encoding", c3),
        _panel("fig2c_opt4_rarefaction_examples", c4),
    ]
    _contact_sheet(pngs)
    return pngs


def _contact_sheet(pngs: list[Path]) -> Path:
    titles = (
        "1 · Alleles per locus",
        "2 · Rarefaction",
        "3 · Biallelic vs dosage",
        "4 · Rarefaction + example loci",
    )
    fig, axes = plt.subplots(2, 2, figsize=(9, 8.4))
    fig.subplots_adjust(left=0.01, right=0.99, top=0.96, bottom=0.01, wspace=0.03, hspace=0.08)
    for ax, png, title in zip(axes.flat, pngs, titles):
        ax.imshow(plt.imread(png))
        ax.set_title(title, loc="left", fontsize=10)
        ax.axis("off")
    out = OUTDIR / "fig2c_options_contact_sheet.png"
    fig.savefig(out, dpi=200)
    plt.close(fig)
    return out


if __name__ == "__main__":
    for p in render():
        print(p)
