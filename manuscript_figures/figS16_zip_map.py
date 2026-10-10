"""Figure S16: Phase 2 participants on a US map, 3-digit ZIP, colored by ancestry.

Points are jittered ZIP3 centroids from covariates.v6, not addresses.
ZIP3 000 (withheld / unknown) is omitted.
"""

from __future__ import annotations

import csv
import gzip
import json
from pathlib import Path

import numpy as np
from matplotlib.lines import Line2D

from .style import (
    ANC_COLORS,
    ANC_LABELS,
    ANC_ORDER,
    FIG_W,
    ROOT,
    covariates_path,
    new_figure,
    save,
)

GEO_DIR = ROOT / "data" / "figS16"
CONUS = (-125.0, -66.5, 24.5, 49.4)
# Include the Aleutians (GeoJSON uses lon < -180). lat_ref 58 is a visual
# compromise between the mainland and the chain so the inset is not too tall.
AK = (-188.5, -129.5, 51.25, 71.55)
HI = (-160.55, -154.70, 18.85, 22.30)
LAT_CONUS = 39.0
LAT_AK = 58.0
LAT_HI = 20.5


def _truthy(v) -> bool:
    return str(v).strip().lower() in {"true", "1", "yes"}


def _load_centroids() -> dict[str, tuple[float, float]]:
    path = GEO_DIR / "zip3_centroids.tsv"
    out: dict[str, tuple[float, float]] = {}
    with path.open() as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            out[row["zip3"]] = (float(row["lon"]), float(row["lat"]))
    return out


def _load_states() -> list[dict]:
    return json.loads((GEO_DIR / "us-states.json").read_text())["features"]


def _rings(geom: dict):
    t = geom["type"]
    coords = geom["coordinates"]
    if t == "Polygon":
        yield coords[0]
    elif t == "MultiPolygon":
        for poly in coords:
            yield poly[0]


def _lon_scale(lon, lat_ref: float):
    return np.asarray(lon, dtype=float) * np.cos(np.radians(lat_ref))


def _ring_in_view(ring, xlim, ylim) -> bool:
    xs = [p[0] for p in ring]
    ys = [p[1] for p in ring]
    return not (max(xs) < xlim[0] or min(xs) > xlim[1] or max(ys) < ylim[0] or min(ys) > ylim[1])


def _draw_named(
    ax,
    features: list[dict],
    names: set[str],
    xlim,
    ylim,
    lat_ref: float,
    fc="#F2F2F2",
    ec="#B5B5B5",
    lw=0.35,
) -> None:
    for feat in features:
        if feat["properties"]["name"] not in names:
            continue
        for ring in _rings(feat["geometry"]):
            if not _ring_in_view(ring, xlim, ylim):
                continue
            xs = _lon_scale([p[0] for p in ring], lat_ref)
            ys = [p[1] for p in ring]
            ax.fill(xs, ys, facecolor=fc, edgecolor=ec, lw=lw, zorder=0, closed=True, clip_on=True)


def _load_people() -> tuple[np.ndarray, np.ndarray, np.ndarray, int, int]:
    """Return lon, lat, ancestry for mappable people, plus n_total and n_withheld."""
    path = covariates_path()
    cents = _load_centroids()
    lon, lat, anc = [], [], []
    n_total = 0
    n_withheld = 0
    with gzip.open(path, "rt") as fh:
        for row in csv.DictReader(fh):
            if not _truthy(row.get("final_releasable_v9", "")):
                continue
            if row.get("technology") != "PacBio":
                continue
            n_total += 1
            a = (row.get("ancestry_pred_other") or "oth").upper()
            if a not in ANC_COLORS:
                a = "OTH"
            z = str(row.get("zip3") or "").strip().split(".")[0]
            if not z or z.lower() in {"nan", "none", "na"}:
                n_withheld += 1
                continue
            z = z.zfill(3)[:3]
            if z == "000" or z not in cents:
                n_withheld += 1
                continue
            x, y = cents[z]
            lon.append(x)
            lat.append(y)
            anc.append(a)
    return np.asarray(lon), np.asarray(lat), np.asarray(anc), n_total, n_withheld


def _jitter(lon: np.ndarray, lat: np.ndarray, rng: np.random.Generator) -> tuple[np.ndarray, np.ndarray]:
    # ~10–15 km; ZIP3 is already coarse. Reproducible; not an address.
    return lon + rng.normal(0, 0.12, lon.size), lat + rng.normal(0, 0.10, lat.size)


def _regions(lon: np.ndarray, lat: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    conus = (lon > -126) & (lon < -66) & (lat > 24) & (lat < 50)
    hi = (lon > -162) & (lon < -154) & (lat > 18) & (lat < 23)
    ak = lat > 50
    return conus, hi, ak


def _scatter(ax, lon, lat, anc, rng, lat_ref: float, s=7.5, alpha=0.62) -> None:
    order = rng.permutation(lon.size)
    ax.scatter(
        _lon_scale(lon[order], lat_ref),
        lat[order],
        c=[ANC_COLORS[a] for a in anc[order]],
        s=s,
        alpha=alpha,
        linewidths=0,
        rasterized=True,
        zorder=2,
        clip_on=True,
    )


def _geo_aspect(xlim: tuple[float, float], ylim: tuple[float, float], lat_ref: float) -> float:
    """Axes box width / height after cosine lon scaling at lat_ref."""
    return ((xlim[1] - xlim[0]) * np.cos(np.radians(lat_ref))) / (ylim[1] - ylim[0])


def _set_geo(ax, xlim, ylim, lat_ref: float, *, frame: bool) -> None:
    ax.set_xlim(_lon_scale(xlim[0], lat_ref), _lon_scale(xlim[1], lat_ref))
    ax.set_ylim(*ylim)
    ax.set_aspect("equal", adjustable="box", anchor="SW")
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_facecolor("white")
    if frame:
        for sp in ax.spines.values():
            sp.set_visible(True)
            sp.set_linewidth(0.45)
            sp.set_color("#888888")
    else:
        ax.axis("off")


def _add_map_axes(fig, left_in, bottom_in, width_in, xlim, ylim, lat_ref: float):
    fig_w, fig_h = fig.get_size_inches()
    height_in = width_in / _geo_aspect(xlim, ylim, lat_ref)
    ax = fig.add_axes(
        [left_in / fig_w, bottom_in / fig_h, width_in / fig_w, height_in / fig_h]
    )
    return ax, height_in


def _ancestry_legend(ax, anc: np.ndarray) -> None:
    handles = []
    for a in ANC_ORDER:
        n = int((anc == a).sum())
        handles.append(
            Line2D(
                [0],
                [0],
                marker="o",
                color="none",
                markerfacecolor=ANC_COLORS[a],
                markersize=5.2,
                label=f"{ANC_LABELS[a]}  {n:,}",
            )
        )
    ax.legend(
        handles=handles,
        loc="lower right",
        bbox_to_anchor=(0.995, 0.02),
        fontsize=5.5,
        frameon=False,
        borderpad=0.12,
        handletextpad=0.32,
        labelspacing=0.28,
        title="Ancestry",
        title_fontsize=6.0,
    )


def render() -> tuple[Path, Path]:
    rng = np.random.default_rng(12)
    features = _load_states()
    lon, lat, anc, _n_total, _n_withheld = _load_people()
    jlon, jlat = _jitter(lon, lat, rng)
    conus, hi, ak = _regions(jlon, jlat)
    conus_names = {
        f["properties"]["name"]
        for f in features
        if f["properties"]["name"] not in {"Alaska", "Hawaii", "Puerto Rico"}
    }

    pad = 0.04
    conus_w = FIG_W - 2 * pad
    conus_h = conus_w / _geo_aspect(CONUS[:2], CONUS[2:], LAT_CONUS)
    fig_h = conus_h + 2 * pad
    fig = new_figure(fig_h)
    ax, _ = _add_map_axes(fig, pad, pad, conus_w, CONUS[:2], CONUS[2:], LAT_CONUS)
    _draw_named(ax, features, conus_names, CONUS[:2], CONUS[2:], LAT_CONUS)
    _scatter(ax, jlon[conus], jlat[conus], anc[conus], rng, LAT_CONUS, s=6.5, alpha=0.58)
    _set_geo(ax, CONUS[:2], CONUS[2:], LAT_CONUS, frame=False)
    ax.set_zorder(1)

    ax_ak, _ = _add_map_axes(fig, pad + 0.04, pad + 0.12, 1.42, AK[:2], AK[2:], LAT_AK)
    _draw_named(ax_ak, features, {"Alaska"}, AK[:2], AK[2:], LAT_AK, lw=0.4)
    if ak.any():
        _scatter(ax_ak, jlon[ak], jlat[ak], anc[ak], rng, LAT_AK, s=12, alpha=0.95)
    _set_geo(ax_ak, AK[:2], AK[2:], LAT_AK, frame=True)
    ax_ak.text(0.035, 0.96, "Alaska", transform=ax_ak.transAxes, fontsize=5.6, ha="left", va="top", color="#444444")
    ax_ak.set_zorder(3)

    hi_w = 0.82
    ax_hi, _ = _add_map_axes(fig, pad + 0.04 + 1.42 + 0.07, pad + 0.12, hi_w, HI[:2], HI[2:], LAT_HI)
    _draw_named(ax_hi, features, {"Hawaii"}, HI[:2], HI[2:], LAT_HI, lw=0.4)
    if hi.any():
        _scatter(ax_hi, jlon[hi], jlat[hi], anc[hi], rng, LAT_HI, s=12, alpha=0.95)
    _set_geo(ax_hi, HI[:2], HI[2:], LAT_HI, frame=True)
    ax_hi.text(0.04, 0.96, "Hawaii", transform=ax_hi.transAxes, fontsize=5.6, ha="left", va="top", color="#444444")
    ax_hi.set_zorder(3)

    _ancestry_legend(ax, anc)
    return save(fig, "figS16_zip_map")
