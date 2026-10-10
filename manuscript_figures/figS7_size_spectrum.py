"""Figure S7: Variant length spectrum across the four callset partitions."""

from __future__ import annotations

from pathlib import Path

import numpy as np
from matplotlib.gridspec import GridSpec
from matplotlib.patches import Rectangle

from .style import (
    COL_IN,
    OKABE,
    SV_N,
    ax_panel,
    despine,
    format_millions,
    new_figure,
    save,
    save_panel,
    stamp_mockup,
    subplot_panel,
)

# Partition colors (not ancestry, not INS/DEL).
C_DV = OKABE["sky"]
C_MAIN = OKABE["blue"]
C_LARGE = OKABE["vermillion"]
C_BND = OKABE["grey"]

# Placeholder catalog sizes. Main INS/DEL ≥20 bp match Table 2.
N_SNV = 95_000_000
N_DV_INS = 5_600_000  # <50 bp
N_DV_DEL = 7_300_000
N_MAIN_INS = SV_N["INS_20"]
N_MAIN_DEL = SV_N["DEL_20"]
N_MAIN_INV = SV_N["INV_20"]
N_LARGE = 225_000  # WDL comment: ~225k ultralong intervals
N_BND = 85_000
N_BND_INTRA = 62_000
N_BND_INTER = N_BND - N_BND_INTRA

BINS = np.logspace(0, 7.0, 92)  # 1 bp → 10 Mb
ZOOM_BINS = np.linspace(20, 50, 16)


def _clip(x: np.ndarray, lo: float, hi: float) -> np.ndarray:
    return x[(x >= lo) & (x < hi)]


def _lognorm(n: int, mean: float, sigma: float, rng: np.random.Generator) -> np.ndarray:
    return rng.lognormal(np.log(mean), sigma, n)


def _geom_lengths(n: int, p: float, rng: np.random.Generator, lo: int, hi: int) -> np.ndarray:
    x = rng.geometric(p, n) + (lo - 1)
    return np.clip(x, lo, hi).astype(float)


def _simulate(rng: np.random.Generator) -> dict[str, np.ndarray]:
    """Placeholder length draws. Peaks at 1 bp, Alu, SVA, L1 are layout cues."""

    def mix_ins_small(n: int) -> np.ndarray:
        w = rng.random(n)
        out = np.empty(n)
        m1 = w < 0.62
        m2 = (w >= 0.62) & (w < 0.88)
        m3 = w >= 0.88
        out[m1] = _geom_lengths(m1.sum(), 0.55, rng, 1, 8)
        out[m2] = rng.integers(8, 20, m2.sum())
        out[m3] = rng.integers(20, 50, m3.sum())
        return out

    def mix_del_small(n: int) -> np.ndarray:
        w = rng.random(n)
        out = np.empty(n)
        m1 = w < 0.58
        m2 = (w >= 0.58) & (w < 0.86)
        m3 = w >= 0.86
        out[m1] = _geom_lengths(m1.sum(), 0.50, rng, 1, 8)
        out[m2] = rng.integers(8, 20, m2.sum())
        out[m3] = rng.integers(20, 50, m3.sum())
        return out

    def mix_sv(n: int, ins: bool) -> np.ndarray:
        w = rng.random(n)
        out = np.empty(n)
        # Small / VNTR, Alu, SVA, L1, background.
        cuts = (0.48, 0.72, 0.82, 0.93) if ins else (0.52, 0.78, 0.86, 0.94)
        b0 = w < cuts[0]
        b1 = (w >= cuts[0]) & (w < cuts[1])
        b2 = (w >= cuts[1]) & (w < cuts[2])
        b3 = (w >= cuts[2]) & (w < cuts[3])
        b4 = w >= cuts[3]
        out[b0] = np.clip(_lognorm(b0.sum(), 42 if ins else 55, 0.70, rng), 20, 180)
        out[b1] = np.clip(_lognorm(b1.sum(), 310, 0.11, rng), 250, 380)
        out[b2] = np.clip(_lognorm(b2.sum(), 1400, 0.32, rng), 700, 2800)
        out[b3] = np.clip(_lognorm(b3.sum(), 6000, 0.12, rng), 4500, 7500)
        out[b4] = np.clip(_lognorm(b4.sum(), 400 if ins else 800, 1.05, rng), 20, 9999)
        return out

    def mix_large(n: int) -> np.ndarray:
        return 10 ** rng.uniform(4.0, 6.8, n)

    n_plot = 40_000  # subsample for speed; hist is scaled to catalog n
    dv_ins = mix_ins_small(n_plot)
    dv_del = mix_del_small(n_plot)
    main_ins = mix_sv(n_plot, ins=True)
    main_del = mix_sv(n_plot, ins=False)
    large = mix_large(n_plot)
    large_ins = large[rng.random(n_plot) < 0.18]
    large_del = large[rng.random(n_plot) < 0.62]
    bnd_span = 10 ** rng.uniform(2.8, 7.7, int(n_plot * 0.45))
    return {
        "dv_ins": dv_ins,
        "dv_del": dv_del,
        "main_ins": main_ins,
        "main_del": main_del,
        "large_ins": large_ins,
        "large_del": large_del,
        "bnd_span": bnd_span,
    }


def _stairs(ax, values: np.ndarray, n_true: int, bins: np.ndarray, color: str, **kw):
    """Histogram scaled so the area matches n_true catalog sites."""
    raw, edges = np.histogram(values, bins=bins)
    scale = n_true / max(values.size, 1)
    ax.stairs(raw * scale / 1e3, edges, fill=True, color=color, **kw)


def _panel_rails(ax) -> None:
    ax.set_xscale("log")
    ax.set_xlim(0.85, 1.15e7)
    ax.set_ylim(0.35, 4.15)
    ax.set_yticks([])
    ax.tick_params(labelbottom=False, bottom=False)
    for sp in ax.spines.values():
        sp.set_visible(False)
    ax.set_title("Callset partitions by resolved length", loc="left", pad=2)

    bands = [
        (3.10, 1.0, 50, C_DV, "DeepVariant", "SNVs + indels <50 bp"),
        (2.10, 20, 10_000, C_MAIN, "Integrated SV", "20 bp – 10 kb  ·  Kanpig"),
        (1.10, 10_000, 1.0e7, C_LARGE, "Large events", ">10 kb  ·  not Kanpig-genotyped"),
    ]
    for y, lo, hi, color, name, note in bands:
        ax.add_patch(Rectangle((lo, y), hi - lo, 0.82, facecolor=color, edgecolor="none", alpha=0.90, zorder=2))
        mid = np.sqrt(lo * hi)
        ax.text(mid, y + 0.52, name, ha="center", va="center", fontsize=6.4, color="white", fontweight="bold", zorder=3)
        ax.text(mid, y + 0.22, note, ha="center", va="center", fontsize=5.2, color="white", zorder=3)

    for y in (3.10, 2.10):
        ax.add_patch(
            Rectangle((20, y), 30, 0.82, facecolor="none", edgecolor="#111111", hatch="////", linewidth=0.45, alpha=0.40, zorder=4)
        )
    ax.annotate(
        "20–49 bp in both",
        xy=(35, 3.92),
        xytext=(80, 3.92),
        fontsize=5.4,
        color="#333333",
        va="center",
        arrowprops=dict(arrowstyle="-|>", color="#555555", lw=0.5, mutation_scale=6),
        zorder=5,
    )
    for x in (1, 20, 50, 10_000):
        ax.axvline(x, color="#888888", lw=0.5, ls=":", zorder=1)


def _style_logx(ax, labels: bool) -> None:
    ticks = [1, 20, 50, 100, 300, 1e3, 1e4, 1e5, 1e6, 1e7]
    labs = ["1", "20", "50", "100", "300", "1 kb", "10 kb", "100 kb", "1 Mb", "10 Mb"]
    ax.set_xscale("log")
    ax.set_xlim(0.85, 1.15e7)
    ax.set_xticks(ticks)
    if labels:
        ax.set_xticklabels(labs, fontsize=5.5)
    else:
        ax.tick_params(labelbottom=False)
    ax.tick_params(axis="x", which="minor", bottom=False)
    for x in (20, 50, 10_000):
        ax.axvline(x, color="#888888", lw=0.5, ls=":", zorder=1)


def _landmarks(ax) -> None:
    for x, lab in [(310, "Alu"), (1400, "SVA"), (6000, "L1")]:
        ax.axvline(x, color="#B0B0B0", lw=0.45, ls="--", zorder=1)
        ax.text(x, 0.90, lab, transform=ax.get_xaxis_transform(), ha="center", va="top", fontsize=5.3, color="#666666")


def _panel_ins(ax, data: dict[str, np.ndarray]) -> None:
    _stairs(ax, data["dv_ins"], N_DV_INS, BINS, C_DV, alpha=0.55, zorder=3)
    _stairs(ax, data["main_ins"], N_MAIN_INS, BINS, C_MAIN, alpha=0.50, zorder=2)
    _stairs(ax, data["large_ins"], int(N_LARGE * 0.18), BINS, C_LARGE, alpha=0.75, zorder=4)
    _style_logx(ax, labels=False)
    ax.set_yscale("log")
    ax.set_ylim(0.7, 900)
    ax.set_ylabel("Insertions\n(thousands / bin)")
    ax.set_title("Length spectrum, insertions up / deletions down", loc="left", pad=2)
    despine(ax)
    _landmarks(ax)
    ax.text(0.02, 0.88, f"SNVs  {format_millions(N_SNV)}  (off scale)", transform=ax.transAxes, fontsize=5.6, color="#444444")
    ax.plot([], [], color=C_DV, lw=6, alpha=0.7, label="DeepVariant")
    ax.plot([], [], color=C_MAIN, lw=6, alpha=0.7, label="Integrated SV")
    ax.plot([], [], color=C_LARGE, lw=6, alpha=0.8, label="Large events")
    ax.legend(loc="upper right", fontsize=5.5, ncol=3)


def _panel_del(ax, data: dict[str, np.ndarray]) -> None:
    _stairs(ax, data["dv_del"], N_DV_DEL, BINS, C_DV, alpha=0.55, zorder=3)
    _stairs(ax, data["main_del"], N_MAIN_DEL, BINS, C_MAIN, alpha=0.50, zorder=2)
    _stairs(ax, data["large_del"], int(N_LARGE * 0.62), BINS, C_LARGE, alpha=0.75, zorder=4)
    _style_logx(ax, labels=True)
    ax.set_yscale("log")
    ax.set_ylim(0.7, 900)
    ax.invert_yaxis()
    ax.set_ylabel("Deletions\n(thousands / bin)")
    ax.set_xlabel("Resolved variant length")
    despine(ax)
    ax.text(
        0.98,
        0.12,
        f"INV in main callset  n = {N_MAIN_INV}  (not shown)",
        transform=ax.transAxes,
        fontsize=5.4,
        color="#666666",
        ha="right",
    )


def _panel_bnd(ax, data: dict[str, np.ndarray]):
    fig = ax.figure
    spec = ax.get_subplotspec()
    ax.remove()
    inner = spec.subgridspec(2, 1, height_ratios=[1.05, 1.15], hspace=0.42)
    ax_top = fig.add_subplot(inner[0])
    ax_span = fig.add_subplot(inner[1])

    ax_top.set_xlim(0, 1.0)
    ax_top.set_ylim(-0.15, 2.55)
    ax_top.axis("off")
    ax_top.set_title("Breakends have no SVLEN", loc="left", pad=2)
    ax_top.add_patch(Rectangle((0.0, 1.55), 1.0, 0.78, facecolor=C_BND, edgecolor="none", zorder=2))
    ax_top.text(0.50, 1.94, f"{format_millions(N_BND)} sites", ha="center", va="center", fontsize=6.6, color="white", fontweight="bold")
    ax_top.text(0.50, 1.70, "unresolved / inter-break", ha="center", va="center", fontsize=5.3, color="white")

    labs = ["Intra-chromosomal", "Inter-chromosomal"]
    vals = np.array([N_BND_INTRA, N_BND_INTER], dtype=float)
    y = np.array([0.85, 0.18])
    ax_top.barh(y, vals / vals.max() * 0.72, height=0.48, color=[OKABE["blue"], "#B0B0B0"], zorder=2, left=0.0)
    for yi, v, lab in zip(y, vals, labs):
        ax_top.text(0.0, yi + 0.28, lab, fontsize=5.4, color="#444444", va="bottom")
        ax_top.text(0.98, yi, format_millions(v), va="center", ha="right", fontsize=5.6, color="#222222")

    raw, edges = np.histogram(data["bnd_span"], bins=np.logspace(2.6, 7.8, 32))
    scale = N_BND_INTRA / max(data["bnd_span"].size, 1)
    ax_span.stairs(raw * scale / 1e3, edges, fill=True, color=C_BND, alpha=0.85)
    ax_span.set_xscale("log")
    ax_span.set_xlim(500, 8e7)
    ax_span.set_xlabel("Inferred intra-chrom. span")
    ax_span.set_ylabel("Thousands / bin")
    ax_span.set_title("When both breaks are on one chromosome", loc="left", pad=2)
    despine(ax_span)
    ax_span.axvline(10_000, color="#888888", lw=0.5, ls=":")
    return ax_top


def _panel_zoom(ax, data: dict[str, np.ndarray]) -> None:
    dv = _clip(data["dv_ins"], 20, 50)
    main = _clip(data["main_ins"], 20, 50)
    # Scale by the 20–49 fraction of each catalog (placeholder).
    n_dv = int(N_DV_INS * 0.12)
    n_main = int(N_MAIN_INS - SV_N["INS_50"])
    _stairs(ax, dv, n_dv, ZOOM_BINS, C_DV, alpha=0.55)
    _stairs(ax, main, n_main, ZOOM_BINS, C_MAIN, alpha=0.50)
    ax.set_xlim(20, 50)
    ax.set_xlabel("Insertion length (bp)")
    ax.set_ylabel("Thousands / bin")
    ax.set_title("20–49 bp is in DeepVariant and the SV callset", loc="left", pad=2)
    despine(ax)
    ax.axvline(50, color="#888888", lw=0.6, ls=":")
    ax.text(0.04, 0.88, "Not the same sites — two callers, one size window", transform=ax.transAxes, fontsize=5.5, color="#444444")
    ax.plot([], [], color=C_DV, lw=6, alpha=0.7, label="DeepVariant")
    ax.plot([], [], color=C_MAIN, lw=6, alpha=0.7, label="Integrated SV")
    ax.legend(loc="upper right", fontsize=5.4)


def _panel_counts(ax) -> None:
    labels = [
        "DV SNVs",
        "DV indels <50 bp",
        "SV 20 bp–10 kb",
        "Large >10 kb",
        "Breakends",
    ]
    vals = np.array(
        [N_SNV, N_DV_INS + N_DV_DEL, N_MAIN_INS + N_MAIN_DEL + N_MAIN_INV, N_LARGE, N_BND],
        dtype=float,
    )
    colors = [C_DV, C_DV, C_MAIN, C_LARGE, C_BND]
    y = np.arange(len(labels))
    ax.barh(y, vals, color=colors, height=0.68, zorder=2)
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=6.0)
    ax.invert_yaxis()
    ax.set_xscale("log")
    ax.set_xlim(4e4, 2.2e8)
    ax.set_xlabel("Catalog sites")
    ax.set_title("Partition sizes (placeholder)", loc="left", pad=2)
    despine(ax)
    for yi, v in zip(y, vals):
        ax.text(v * 1.12, yi, format_millions(v), va="center", fontsize=5.6, color="#333333")


def render() -> tuple[Path, Path]:
    rng = np.random.default_rng(7)
    data = _simulate(rng)

    fig_spec = new_figure(4.35)
    gs_spec = GridSpec(
        3,
        1,
        figure=fig_spec,
        left=0.11,
        right=0.98,
        top=0.92,
        bottom=0.10,
        hspace=0.08,
        height_ratios=[0.55, 1.0, 1.0],
    )
    ax_rails = fig_spec.add_subplot(gs_spec[0])
    ax_ins = fig_spec.add_subplot(gs_spec[1], sharex=ax_rails)
    ax_del = fig_spec.add_subplot(gs_spec[2], sharex=ax_rails)
    _panel_rails(ax_rails)
    _panel_ins(ax_ins, data)
    _panel_del(ax_del, data)
    ax_ins.tick_params(labelbottom=False)
    save_panel(fig_spec, "figS7_length_spectrum")

    subplot_panel("figS7_bnd", 3.55, lambda ax: _panel_bnd(ax, data), width=COL_IN * 1.15, left=0.16, right=0.96, top=0.90, bottom=0.12)
    ax_panel("figS7_zoom", 2.65, lambda ax: _panel_zoom(ax, data), left=0.12, right=0.98, top=0.86, bottom=0.16)
    ax_panel("figS7_counts", 2.65, _panel_counts, width=COL_IN * 1.15, left=0.28, right=0.96, top=0.88, bottom=0.16)

    fig = new_figure(7.90)
    gs = GridSpec(
        2,
        1,
        figure=fig,
        left=0.11,
        right=0.98,
        top=0.95,
        bottom=0.07,
        hspace=0.22,
        height_ratios=[2.55, 1.18],
    )
    top = gs[0].subgridspec(3, 2, height_ratios=[0.70, 1.05, 1.05], width_ratios=[3.20, 1.05], hspace=0.06, wspace=0.28)
    bot = gs[1].subgridspec(1, 2, width_ratios=[3.20, 1.05], wspace=0.28)

    ax_rails = fig.add_subplot(top[0, 0])
    ax_ins = fig.add_subplot(top[1, 0], sharex=ax_rails)
    ax_del = fig.add_subplot(top[2, 0], sharex=ax_rails)
    ax_bnd_host = fig.add_subplot(top[:, 1])
    ax_zoom = fig.add_subplot(bot[0, 0])
    ax_counts = fig.add_subplot(bot[0, 1])

    _panel_rails(ax_rails)
    _panel_ins(ax_ins, data)
    _panel_del(ax_del, data)
    _panel_bnd(ax_bnd_host, data)
    _panel_zoom(ax_zoom, data)
    _panel_counts(ax_counts)

    stamp_mockup(fig, extra="Length draws and partition counts are layout placeholders")
    return save(fig, "figS7_size_spectrum")
