#!/usr/bin/env python3
"""Per-locus dossier for every hit the FELIX screen retains: the diagnostics that decide whether a hit is real.

This is the data side of the one-page-per-locus review plate (manuscript_figures/figS_hit_dossier.py). It takes the
whole-scan hit table written by felix_08_hit_sweep (felix_sweep_all_hits.tsv), reduces it to loci, and assembles for
each locus everything a reviewer needs short of the read pileups:

======================================  ======================================================  ===================
output                                  holds                                                   cost
======================================  ======================================================  ===================
felix_dossier_loci.tsv                  one row per locus: identity, lead, arms, window, the    free (hit table plus
                                        phenotype's case count and per-arm lambda               calibration_summary)
felix_dossier_leads.tsv                 one row per locus-arm: the lead marker's FULL           free (harvest cache
                                        per-ancestry block, genotype counts, missingness,       top_hits.tsv)
                                        SPA status, P_het / P_hom / P_cct
felix_dossier_context.tsv               one row per locus: cytoband, segmental duplications,    UCSC hg38 API,
                                        tandem repeat, gene and nearest gene                    ~4 requests a locus
felix_dossier_genes.tsv                 one row per gene in a locus window (the gene track)      with the above
felix_dossier_region.<id>.<arm>.tsv.gz  thinned regional support around the locus, per arm       ONE STREAM PER RUN
======================================  ======================================================  ===================

The split matters: everything except the regional files comes out of the harvest cache felix_08 already builds, so the
expensive part (streaming each run's multi-GB merged.felix.tsv) is isolated in one step and reads every window of a
run in a single pass.

There is no QQ export. The harvest cache keeps only calibration_summary.tsv and top_hits.tsv, which hold no full
p-value list, so the page reports the phenotype's lambda per arm rather than drawing its QQ.

``dossier_loci`` must agree with manuscript_figures.fig3_real_data.Real.loci() -- the plate reviews exactly the hits
Figure 3 counts -- which test_felix_hit_dossier.py checks against the real table in the manuscript repo's data/.

Definitions follow felix_phenotype_report: a hit is retained when ``verdict == 'keep'`` (none of spa_fallback,
no_case_carriers, impossible_p, excess_het fired) and it has at least MIN_CASE_CARRIERS case carriers. Hits of one
phenotype within LOCUS_WINDOW on a chromosome are one locus.
"""

from __future__ import annotations

import math
import re
from pathlib import Path

import numpy as np
import pandas as pd

LOCUS_WINDOW = 500_000  # retained hits of one phenotype closer than this are one locus (fig3_real_data.LOCUS_WINDOW)
MIN_CASE_CARRIERS = 10  # a hit with fewer case carriers is not reviewed (fig3_real_data.MIN_CASE_CARRIERS)
MIN_STRATUM_CARRIERS = 10  # a PER-ANCESTRY beta resting on fewer case carriers than this is not shown as an estimate.
# Sparse-data bias pushes |beta| away from zero as the carrier count falls, measured over the 144 estimable strata of
# the 49 structural loci: median |beta| 0.70 and max 1.49 at 20+ case carriers, against median 1.79 and max 6.15 at
# 2-4. The floor matches MIN_CASE_CARRIERS, which is the project's locus-level convention. It is a floor, not a
# guarantee -- the 10-19 band still holds one |beta| of 5.09 -- so the carrier count is printed beside every stratum.
SITE_WINDOW = 1_000  # a repeat-arm and a main-arm hit this close are the same site
REGION_HALF = 250_000  # the regional plot spans the lead +/- this
ARMS = ("main", "ru", "ru_biallelic")
FLARE_ORDER = ("EAS", "AMR", "EUR", "AFR", "SAS")  # FELIX *_anc1..5, the FLARE reference-panel order
STRATA = tuple(f"anc{k}" for k in range(1, len(FLARE_ORDER) + 1)) + ("ancALL",)
GW_P = 5e-8
P_COL = "P_cct_admixed_c"  # the headline joint p (felix_phenotype_report.hit_report)

# Per-ancestry columns copied verbatim into felix_dossier_leads.tsv, for each stratum in STRATA.
LEAD_STRATUM_COLS = (
    "BETA_c", "SE_c", "p.value_c", "Tstat_c", "var_c", "Is.SPA",
    "AF_Allele2", "AF_case", "AF_ctrl", "MissingRate",
    "N_case", "N_ctrl", "N_case_het", "N_case_hom", "N_ctrl_het", "N_ctrl_hom",
)
LEAD_SCALAR_COLS = ("CHR", "POS", "MarkerID", "Allele1", "Allele2",
                    "P_het_admixed_c", "P_hom_admixed_c", "P_cct_admixed_c",
                    "P_het_admixed", "P_hom_admixed", "P_cct_admixed")

LOCI_FILE = "felix_dossier_loci.tsv"
LEADS_FILE = "felix_dossier_leads.tsv"
REGION_FMT = "felix_dossier_region.{locus_id}.{arm}.tsv.gz"
CALIB_FILE = "felix_dossier_calibration.tsv"  # run x arm lambda_gc for the WHOLE sweep, not just the dossier runs:
# the chip reports this run against the others, so the comparison set has to be every run that was scanned.
LAMBDA_MIN = 0.90  # genomic-control lambda below this is deflated, above LAMBDA_MAX inflated: either way the
LAMBDA_MAX = 1.20  # phenotype's whole scan is suspect and no single locus from it can be read at face value.


def region_path(directory, locus_id: str, arm: str) -> Path:
    """Where the compact region file of one locus and arm lives."""
    return Path(directory) / REGION_FMT.format(locus_id=locus_id, arm=arm)


# --------------------------------------------------------------------------------------------------- locus list


def _chrom_int(v) -> int | None:
    s = str(v).strip().lower().removeprefix("chr")
    return int(s) if s.isdigit() and 1 <= int(s) <= 22 else None


def phecode(run: str) -> str:
    """felix.prod.BI_174_2 -> BI_174.2 (the PhecodeX form), matching fig3_real_data.phecode."""
    name = re.sub(r"^felix\.prod\.", "", str(run))
    m = re.match(r"^([A-Za-z]{2})_(\d+)(?:[_.](\d+))?$", name)
    if not m:
        return name
    base = f"{m.group(1).upper()}_{m.group(2)}"
    return base if m.group(3) is None else f"{base}.{m.group(3)}"


def locus_id(run: str, chrom: int, pos: int) -> str:
    """Stable, filename-safe name for a locus: BI_174_2__chr9_8991058."""
    return f"{str(run).removeprefix('felix.prod.')}__chr{chrom}_{int(pos)}"


def retained(hits: pd.DataFrame, min_case_carriers: int = MIN_CASE_CARRIERS) -> pd.DataFrame:
    """The hits the plate reviews: verdict 'keep' and at least ``min_case_carriers`` case carriers."""
    need = {"run", "arm", "CHR", "POS", "p", "case_carriers", "verdict"}
    missing = need - set(hits.columns)
    if missing:
        raise KeyError(f"hit table lacks column(s) {sorted(missing)}")
    out = hits[hits["verdict"] == "keep"].copy()
    out["chrom"] = out["CHR"].map(_chrom_int)
    out = out[out["chrom"].notna() & (out["case_carriers"] >= min_case_carriers)]
    out["chrom"] = out["chrom"].astype(int)
    return out.reset_index(drop=True)


def dossier_loci(hits: pd.DataFrame, window: int = LOCUS_WINDOW,
                 min_case_carriers: int = MIN_CASE_CARRIERS) -> pd.DataFrame:
    """Retained hits merged into loci, one row each, in genomic order.

    Hits of one phenotype within ``window`` on a chromosome are one locus. The lead is the hit with the smallest p over
    every arm; ``arms`` lists the arms that retained the locus. Mirrors fig3_real_data.Real.loci().
    """
    keep = retained(hits, min_case_carriers)
    rows = []
    for (run, chrom), g in keep.groupby(["run", "chrom"], sort=True):
        g = g.sort_values("POS")
        group: list[pd.Series] = []
        for _, h in g.iterrows():
            if group and h["POS"] - group[-1]["POS"] > window:
                rows.append(_locus_row(run, chrom, group))
                group = []
            group.append(h)
        if group:
            rows.append(_locus_row(run, chrom, group))
    out = pd.DataFrame(rows)
    if out.empty:
        return out
    return out.sort_values(["chrom", "pos"]).reset_index(drop=True)


def _locus_row(run: str, chrom: int, group: list[pd.Series]) -> dict:
    best = min(group, key=lambda h: h["p"])
    arms = sorted({h["arm"] for h in group})
    pos = int(best["POS"])
    return dict(
        locus_id=locus_id(run, chrom, pos),
        run=run,
        phenotype=phecode(run),
        chrom=chrom,
        pos=pos,
        lead_arm=best["arm"],
        lead_marker=best.get("MarkerID"),
        p=float(best["p"]),
        case_carriers=int(best["case_carriers"]),
        min_case_carriers=int(min(h["case_carriers"] for h in group)),
        arms=",".join(arms),
        n_arms=len(arms),
        n_markers=len({(h["arm"], int(h["POS"])) for h in group}),
        ctrl_het_frac=_first_float(group, "ctrl_het_frac"),
        ctrl_F=_first_float(group, "ctrl_F"),
        window_start=max(1, pos - REGION_HALF),
        window_end=pos + REGION_HALF,
    )


def _first_float(group: list[pd.Series], col: str) -> float:
    """The lead hit's value of ``col``, or the first non-null among the locus' hits; nan when the column is absent."""
    best = min(group, key=lambda h: h["p"])
    for h in [best] + list(group):
        v = h.get(col)
        if v is not None and not (isinstance(v, float) and math.isnan(v)) and str(v) not in ("", "nan", "None"):
            try:
                return float(v)
            except (TypeError, ValueError):
                return math.nan
    return math.nan


# -------------------------------------------------------------------------------------------- lead sites

_ALLELE = re.compile(r"^(?:[ACGTNacgtn]+|<[A-Za-z:_]+>|\*)$")


def split_alleles(marker: str | None) -> tuple[str, str] | None:
    """(REF, ALT) of a main-arm marker, or None when the string does not end in two alleles.

    felixla writes ``<ID>_<REF>_<ALT>`` and the ID may itself be a ';'-joined list of ``chr_pos_ref_alt``, so the split
    allele is always the LAST two fields. Records with no ID are written ``chr:pos:ref:alt``. A duplicate of
    manuscript_figures.fig3_real_data.split_alleles, which the renderer cannot share because this module has to import
    on the Terra VM with no manuscript repo; test_split_alleles_agrees_with_the_figure_loader keeps the two in step.
    """
    if not marker or str(marker) in ("nan", "None"):
        return None
    m = str(marker).strip()
    if m.count(":") >= 3 and "_" not in m.rsplit(":", 2)[0].split(";")[-1][-1:]:
        f = m.rsplit(":", 2)
        return (f[1], f[2]) if _ALLELE.match(f[1]) and _ALLELE.match(f[2]) else None
    f = m.rsplit("_", 2)
    if len(f) == 3 and _ALLELE.match(f[1]) and _ALLELE.match(f[2]):
        return f[1], f[2]
    return None


SV_MIN = 50  # a length change of at least this many bp makes the record an SV (fig3_real_data.SV_MIN, aggregate_repeat_loci --split-bp)
STRUCTURAL_CLASSES = ("SV", "repeat")  # the allele-based half of the structural rule; the other half is the repeat context


def variant_class(marker: str | None, arm: str = "main") -> str | None:
    """'SV', 'indel' or 'SNV' for a main-arm marker, 'repeat' for the repeat arms, None when the alleles cannot be read.

    Mirrors fig3_real_data.variant_class, which the dossier cannot import: this module runs on Terra with only
    felix/scripts bootstrapped. test_the_variant_classes_match_fig3 pins the two together.
    """
    if arm in ("ru", "ru_biallelic"):
        return "repeat"
    al = split_alleles(marker)
    if not al:
        return None
    ref, alt = al
    if ref.startswith("<") or alt.startswith("<"):
        return "SV"
    dlen = abs(len(alt) - len(ref))
    if dlen >= SV_MIN or (dlen == 0 and len(ref) >= SV_MIN):
        return "SV"
    return "indel" if dlen >= 1 else "SNV"


def length_change(marker: str | None, arm: str = "main") -> int | None:
    """Signed len(ALT) - len(REF) of a main-arm marker; None for the repeat arms and for symbolic or unreadable alleles."""
    if arm in ("ru", "ru_biallelic"):
        return None
    al = split_alleles(marker)
    if not al or al[0].startswith("<") or al[1].startswith("<"):
        return None
    return len(al[1]) - len(al[0])


def is_structural(loci: pd.DataFrame, context: pd.DataFrame | None = None) -> pd.Series:
    """Which loci are the long-read-discoverable class: an SV or repeat-arm lead, or any lead inside a tandem repeat.

    The allele string alone undercounts badly. On the 227-locus dossier only 8 leads are an SV or a repeat-arm hit,
    yet 43 sit inside a UCSC simpleRepeat interval -- a 4 bp TG contraction at a dinucleotide repeat is a repeat-length
    change by mechanism while its alleles read as an ordinary indel (GI_542.11 chr1:26,664,446 is exactly that). The
    repeat context is therefore part of the rule, not a refinement of it.

    ``context`` is the felix_dossier_context.tsv frame; without it only the allele-based half can be applied, so pass
    it whenever it exists or the filter will silently keep too little.
    """
    cls = pd.Series([variant_class(m, a) for m, a in zip(loci["lead_marker"], loci["lead_arm"])], index=loci.index)
    keep = cls.isin(STRUCTURAL_CLASSES)
    if context is not None and len(context) and "in_tandem_repeat" in context.columns:
        tr = context.set_index("locus_id")["in_tandem_repeat"] if "locus_id" in context.columns else context["in_tandem_repeat"]
        keep = keep | loci["locus_id"].map(lambda i: bool(tr.get(i, False))).astype(bool)
    return keep


def _ordinal(n: float) -> str:
    """'1st', '2nd', '3rd', '4th' ... '11th', '12th', '13th' -- the teens take 'th' whatever their last digit."""
    i = int(round(n))
    suffix = "th" if 11 <= i % 100 <= 13 else {1: "st", 2: "nd", 3: "rd"}.get(i % 10, "th")
    return f"{i}{suffix}"


def lambda_context(calib: pd.DataFrame | None, run: str, arm: str = "main") -> dict:
    """This run's genomic-control lambda for one arm, and where it sits among every run of the sweep.

    A lambda is not interpretable on its own: it tracks case count, so a 99-case phenotype and a 5,000-case phenotype
    have different expected lambdas under the same null. ``pct`` is the percentile of this run's lambda among the
    other runs of the same arm, which is the comparison the page actually needs, and ``n_runs`` says how many that is.
    """
    out = dict(available=False, lam=math.nan, pct=math.nan, n_runs=0, median=math.nan,
               n_variants=math.nan, n_gw_sig=math.nan, n_cases=math.nan)
    if calib is None or len(calib) == 0 or "lambda_gc" not in calib.columns:
        return out
    same = calib[(calib["arm"] == arm)].copy()
    same["lambda_gc"] = pd.to_numeric(same["lambda_gc"], errors="coerce")
    same = same[same["lambda_gc"].notna()]
    mine = same[same["run"] == run]
    if same.empty or mine.empty:
        return out
    lam = float(mine["lambda_gc"].iloc[0])
    out.update(available=True, lam=lam, n_runs=len(same), median=float(same["lambda_gc"].median()),
               pct=float((same["lambda_gc"] <= lam).mean() * 100.0),
               n_variants=_f(mine["n_variants"].iloc[0]) if "n_variants" in mine else math.nan,
               n_gw_sig=_f(mine["n_gw_sig"].iloc[0]) if "n_gw_sig" in mine else math.nan,
               n_cases=_f(mine["n_cases"].iloc[0]) if "n_cases" in mine else math.nan)
    return out


def lead_sites(loci: pd.DataFrame) -> pd.DataFrame:
    """The site list for pull_site_fields: locus_id, run, chrom, pos, ref, alt (ref/alt empty where the arm has no alleles).

    ``run`` travels with the site because pheno_cov is a per-run input: most runs share one file, but a sex-restricted
    phenotype has its own, and mapping such a locus with the shared file would put excluded participants back into the
    control set.

    Only the lead position of each locus, because that is all the call-rate question needs -- 90 index seeks rather
    than a window read. A repeat-arm lead carries a TRID rather than alleles, so its ref/alt are left blank and every
    record at the position comes back; the caller decides what to do with that.
    """
    rows = []
    for _, L in loci.iterrows():
        al = split_alleles(L.get("lead_marker")) if L.get("lead_arm") == "main" else None
        rows.append(dict(locus_id=L["locus_id"], run=L["run"], chrom=int(L["chrom"]), pos=int(L["pos"]),
                         phenotype=L["phenotype"], lead_arm=L["lead_arm"],
                         ref=(al[0] if al else ""), alt=(al[1] if al else "")))
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------------------- lead per-ancestry block


def lead_rows(loci: pd.DataFrame, top_hits: dict[tuple[str, str], pd.DataFrame],
              p_col: str = P_COL) -> pd.DataFrame:
    """One row per locus-arm with the lead marker's full per-ancestry block, from the harvest cache's top_hits tables.

    ``top_hits`` maps (run, arm) to that arm's top_hits.tsv. For every arm that retained a locus, the lead is the row
    of that arm with the smallest ``p_col`` inside the locus window -- matched on (POS, MarkerID) where the hit table
    gives a marker, because several markers can share one POS. An arm with no matching row is skipped and reported in
    the ``note`` column of the locus it belongs to.
    """
    out = []
    for _, L in loci.iterrows():
        for arm in str(L["arms"]).split(","):
            df = top_hits.get((L["run"], arm))
            if df is None or df.empty:
                out.append(dict(locus_id=L["locus_id"], run=L["run"], arm=arm, found=False,
                                note=f"no top_hits for {arm}"))
                continue
            sel = df[(df["POS"] >= L["window_start"]) & (df["POS"] <= L["window_end"])]
            if p_col in sel.columns:
                sel = sel[pd.to_numeric(sel[p_col], errors="coerce").notna()]
            if sel.empty:
                out.append(dict(locus_id=L["locus_id"], run=L["run"], arm=arm, found=False,
                                note=f"{arm}: no top_hits row in the window"))
                continue
            row = sel.loc[pd.to_numeric(sel[p_col], errors="coerce").idxmin()]
            rec = dict(locus_id=L["locus_id"], run=L["run"], phenotype=L["phenotype"], arm=arm, found=True, note="")
            for c in LEAD_SCALAR_COLS:
                rec[c] = row.get(c)
            for s in STRATA:
                for c in LEAD_STRATUM_COLS:
                    rec[f"{c}_{s}"] = row.get(f"{c}_{s}")
            out.append(rec)
    return pd.DataFrame(out)


def ancestry_block(lead: pd.Series | dict, order: tuple[str, ...] = FLARE_ORDER) -> pd.DataFrame:
    """The lead's per-ancestry statistics as one tidy frame: ancestry, beta, se, p, z, af_case, af_ctrl, carriers, missing.

    ``ancestry`` is the FLARE label for anc1..K plus 'ALL' for the pooled stratum. Rows whose beta and se are both
    missing are kept with nan so the forest plot can show the stratum as untested rather than drop it.
    """
    get = (lambda c: lead.get(c)) if isinstance(lead, dict) else (lambda c: lead.get(c))
    rows = []
    for i, s in enumerate(STRATA):
        label = order[i] if i < len(order) else "ALL"
        beta, se = _f(get(f"BETA_c_{s}")), _f(get(f"SE_c_{s}"))
        tstat, var = _f(get(f"Tstat_c_{s}")), _f(get(f"var_c_{s}"))
        n_case_het, n_case_hom = _f(get(f"N_case_het_{s}")), _f(get(f"N_case_hom_{s}"))
        n_ctrl_het, n_ctrl_hom = _f(get(f"N_ctrl_het_{s}")), _f(get(f"N_ctrl_hom_{s}"))
        rows.append(dict(
            ancestry=label, stratum=s,
            beta=beta, se=se, p=_f(get(f"p.value_c_{s}")),
            z=(beta / se) if (se and se > 0 and not math.isnan(beta)) else math.nan,
            score_z=(tstat / math.sqrt(var)) if (var and var > 0 and not math.isnan(tstat)) else math.nan,
            spa=_bool(get(f"Is.SPA_{s}")),
            af=_f(get(f"AF_Allele2_{s}")), af_case=_f(get(f"AF_case_{s}")), af_ctrl=_f(get(f"AF_ctrl_{s}")),
            missing=_f(get(f"MissingRate_{s}")),
            n_case=_f(get(f"N_case_{s}")), n_ctrl=_f(get(f"N_ctrl_{s}")),
            case_carriers=_nansum(n_case_het, n_case_hom), ctrl_carriers=_nansum(n_ctrl_het, n_ctrl_hom),
            n_case_het=n_case_het, n_case_hom=n_case_hom, n_ctrl_het=n_ctrl_het, n_ctrl_hom=n_ctrl_hom,
        ))
    out = pd.DataFrame(rows)
    # 'thin' marks a stratum whose effect estimate rests on too few case carriers to believe. The case cell is the
    # limiting one in 99.3% of the 144 estimable strata across the 49 structural loci (controls outnumber cases about
    # 9:1), so the rule keys on case carriers alone.
    out["thin"] = (out.case_carriers < MIN_STRATUM_CARRIERS) & (out.stratum != "ancALL")
    return out


def _f(v) -> float:
    if v is None:
        return math.nan
    try:
        x = float(v)
    except (TypeError, ValueError):
        return math.nan
    return x


def _bool(v):
    s = str(v).strip().lower()
    if s in ("true", "1", "1.0", "t", "yes"):
        return True
    if s in ("false", "0", "0.0", "f", "no"):
        return False
    return None


def _nansum(*vals) -> float:
    vs = [v for v in vals if not math.isnan(v)]
    return float(sum(vs)) if vs else math.nan


# ------------------------------------------------------------------------------------------- genomic context


UCSC_API = "https://api.genome.ucsc.edu/getData/track"
CONTEXT_FILE = "felix_dossier_context.tsv"
GENES_FILE = "felix_dossier_genes.tsv"
SEGDUP_FLANK = 10_000  # a lead this close to a segmental duplication is called segdup-adjacent
GENE_TRACK = "ncbiRefSeqCurated"


def ucsc_track(track: str, chrom: str, start: int, end: int, timeout: int = 60) -> list[dict]:
    """Records of one UCSC hg38 track in a window. Returns [] when the track has nothing there.

    The API nests the records under the track name, except that a track with no records in the window may omit the key
    altogether, and some tracks return a dict keyed by chromosome.
    """
    import json
    import urllib.request

    url = f"{UCSC_API}?genome=hg38;track={track};chrom={chrom};start={int(start)};end={int(end)}"
    with urllib.request.urlopen(url, timeout=timeout) as fh:
        payload = json.load(fh)
    got = payload.get(track, [])
    if isinstance(got, dict):  # keyed by chromosome
        got = got.get(chrom, [])
    return got if isinstance(got, list) else []


def _overlap(a0: int, a1: int, b0: int, b1: int) -> int:
    return max(0, min(a1, b1) - max(a0, b0))


def locus_context(chrom: int, pos: int, half: int = REGION_HALF, fetch=ucsc_track) -> tuple[dict, list[dict]]:
    """(context, genes) for one locus from the UCSC hg38 tracks.

    ``context`` carries the cytoband and its stain, whether the lead sits in or beside a segmental duplication and how
    much of its +/-SEGDUP_FLANK neighbourhood is duplicated, whether the lead sits inside a tandem repeat (with the
    repeat's period and copy number), the gene the lead is in, and the nearest gene with its distance. ``genes`` is one
    record per gene in the +/-``half`` window, for the regional plot's gene track. ``fetch`` is injectable for tests.
    """
    c = f"chr{chrom}"
    w0, w1 = max(0, pos - half), pos + half
    band = fetch("cytoBand", c, max(0, pos - 1), pos + 1)
    dups = fetch("genomicSuperDups", c, max(0, pos - SEGDUP_FLANK), pos + SEGDUP_FLANK)
    reps = fetch("simpleRepeat", c, max(0, pos - 1_000), pos + 1_000)
    tx = fetch(GENE_TRACK, c, w0, w1)

    in_dup = any(d["chromStart"] <= pos <= d["chromEnd"] for d in dups)
    dup_bases = 0
    for d in dups:
        dup_bases += _overlap(d["chromStart"], d["chromEnd"], max(0, pos - SEGDUP_FLANK), pos + SEGDUP_FLANK)
    here = [r for r in reps if r["chromStart"] <= pos <= r["chromEnd"]]
    here.sort(key=lambda r: -(r.get("score") or 0))

    genes: dict[str, dict] = {}
    for t in tx:
        name = t.get("name2") or t.get("name")
        g = genes.get(name)
        if g is None:
            genes[name] = dict(gene=name, start=int(t["txStart"]), end=int(t["txEnd"]), strand=t.get("strand", ""))
        else:  # several transcripts of one gene: take the widest span
            g["start"], g["end"] = min(g["start"], int(t["txStart"])), max(g["end"], int(t["txEnd"]))
    glist = sorted(genes.values(), key=lambda g: g["start"])
    inside = [g["gene"] for g in glist if g["start"] <= pos <= g["end"]]
    nearest, nearest_bp = "", math.nan
    if glist:
        dist = {g["gene"]: 0 if g["start"] <= pos <= g["end"] else min(abs(pos - g["start"]), abs(pos - g["end"]))
                for g in glist}
        nearest = min(dist, key=dist.get)
        nearest_bp = dist[nearest]

    ctx = dict(
        chrom=chrom, pos=pos,
        cytoband=(f"{chrom}{band[0]['name']}" if band else ""),
        stain=(band[0]["gieStain"] if band else ""),
        pericentromeric=bool(band and band[0]["gieStain"] in ("acen", "gvar", "stalk")),
        in_segdup=in_dup,
        segdup_frac=round(min(1.0, dup_bases / (2 * SEGDUP_FLANK)), 4),
        n_segdups=len(dups),
        in_tandem_repeat=bool(here),
        repeat_period=(here[0].get("period") if here else None),
        repeat_copies=(here[0].get("copyNum") if here else None),
        gene=";".join(inside),
        nearest_gene=nearest,
        nearest_gene_bp=nearest_bp,
        n_genes_in_window=len(glist),
    )
    return ctx, glist


TRACKS_FILE = "felix_dossier_tracks.tsv"
CALLRATE_FILE = "felix_dossier_callrate.tsv"
TRACK_MIN_PERIOD = 10  # tandem repeats shorter than this are everywhere and would fill the strip with noise
TRACK_MAX_PER_KIND = 60  # keep the widest features of each kind, so a repeat-dense window stays drawable


def locus_tracks(chrom: int, pos: int, half: int = REGION_HALF, fetch=ucsc_track) -> list[dict]:
    """Segmental duplications and tandem repeats across a locus window, as a drawable strip under the regional plot.

    Only repeats with a unit of at least TRACK_MIN_PERIOD bp are kept -- microsatellites are dense enough to fill the
    strip with ink that says nothing -- and only the TRACK_MAX_PER_KIND widest features of each kind.
    """
    c = f"chr{chrom}"
    w0, w1 = max(0, pos - half), pos + half
    out = []
    for kind, track in (("segdup", "genomicSuperDups"), ("repeat", "simpleRepeat")):
        feats = []
        for r in fetch(track, c, w0, w1):
            if kind == "repeat" and (r.get("period") or 0) < TRACK_MIN_PERIOD:
                continue
            feats.append(dict(kind=kind, start=int(r["chromStart"]), end=int(r["chromEnd"]),
                              period=r.get("period"), copies=r.get("copyNum")))
        feats.sort(key=lambda f: f["start"] - f["end"])  # widest first
        out += feats[:TRACK_MAX_PER_KIND]
    return sorted(out, key=lambda f: (f["kind"], f["start"]))


def annotate_context(loci: pd.DataFrame, half: int = REGION_HALF, fetch=ucsc_track,
                     pause: float = 0.1) -> tuple[pd.DataFrame, pd.DataFrame]:
    """``locus_context`` over every locus: (context frame, gene frame, track frame), all keyed by locus_id.

    One failed lookup does not stop the pass -- the locus gets an ``error`` note and empty context -- because the plate
    must still render for the other loci.
    """
    import time

    ctxs, genes, tracks = [], [], []
    for _, L in loci.iterrows():
        try:
            ctx, glist = locus_context(int(L["chrom"]), int(L["pos"]), half=half, fetch=fetch)
            tlist = locus_tracks(int(L["chrom"]), int(L["pos"]), half=half, fetch=fetch)
            ctx["error"] = ""
        except Exception as exc:  # noqa: BLE001 - one locus must not take down the pass
            ctx = dict(chrom=int(L["chrom"]), pos=int(L["pos"]), error=f"{type(exc).__name__}: {exc}")
            glist, tlist = [], []
        ctx["locus_id"] = L["locus_id"]
        ctxs.append(ctx)
        for g in glist:
            genes.append(dict(locus_id=L["locus_id"], **g))
        for t in tlist:
            tracks.append(dict(locus_id=L["locus_id"], **t))
        if pause:
            time.sleep(pause)
    cols = ["locus_id"] + [c for c in ctxs[0] if c != "locus_id"] if ctxs else ["locus_id"]
    return pd.DataFrame(ctxs)[cols], pd.DataFrame(genes), pd.DataFrame(tracks)


def _thin_region(df: pd.DataFrame, p_col: str, keep_below: float, budget: int, seed: int = 0) -> pd.DataFrame:
    """Every row with p below ``keep_below``, plus a random sample of the rest up to ``budget`` rows."""
    p = pd.to_numeric(df[p_col], errors="coerce")
    strong = df[p < keep_below]
    rest = df[~(p < keep_below)]
    if len(rest) > budget:
        rest = rest.sample(budget, random_state=seed)
    return pd.concat([strong, rest]).sort_values("POS")


# -------------------------------------------------------------------------------------- compact region files

# What the regional panel of the plate reads. A full FELIX row is 161 columns and a 1 Mb window is 8-14k rows, so
# shipping the raw extract for 90 loci would be about 1 GB; these columns plus the thinning below bring it to a few MB.
REGION_COLS = (
    "CHR", "POS", "MarkerID", "Allele1", "Allele2",
    "P_cct_admixed_c", "P_het_admixed_c", "P_hom_admixed_c",
    "BETA_c_ancALL", "SE_c_ancALL", "Is.SPA_ancALL", "MissingRate_ancALL",
    "AF_Allele2_ancALL", "AF_case_ancALL", "AF_ctrl_ancALL",
    "N_case_ancALL", "N_ctrl_ancALL",
    "N_case_het_ancALL", "N_case_hom_ancALL", "N_ctrl_het_ancALL", "N_ctrl_hom_ancALL",
) + tuple(f"p.value_c_anc{k}" for k in range(1, len(FLARE_ORDER) + 1))

REGION_KEEP_BELOW = 1e-3  # every marker at least this suggestive is kept, so the peak and its shoulders are exact
REGION_BUDGET = 3_000  # plus this many of the rest, sampled, which is enough for the background cloud of a 1 Mb window


def compact_region(src, out, locus: pd.Series | dict | None = None, keep_below: float = REGION_KEEP_BELOW,
                   budget: int = REGION_BUDGET, seed: int = 0) -> int:
    """Reduce one raw region extract to the columns and rows the plate needs; returns the row count written.

    Keeps every marker with ``P_cct_admixed_c`` below ``keep_below`` and a reproducible sample of ``budget`` others, so
    the peak is exact and the background is only indicative. Writes gzipped TSV. When ``locus`` is given, rows outside
    its window are dropped first, which matters when one extract covers several of a run's loci.
    """
    df = pd.read_csv(src, sep="\t", low_memory=False)
    if locus is not None:
        df = df[(df.POS >= int(locus["window_start"])) & (df.POS <= int(locus["window_end"]))]
    have = [c for c in REGION_COLS if c in df.columns]
    if P_COL not in df.columns:
        raise KeyError(f"{src}: no {P_COL} column; got {list(df.columns)[:8]}...")
    out_df = _thin_region(df[have], P_COL, keep_below, budget, seed=seed)
    out_df.to_csv(out, sep="\t", index=False, compression="gzip")
    return len(out_df)


def region_plan(loci: pd.DataFrame) -> dict[str, list[dict]]:
    """The windows to cut per run: {run: [{locus_id, chrom, start, end}, ...]}, merged where they overlap.

    One streaming pass per run and arm then serves every locus of that run, which is what makes the regional step
    affordable -- a run's merged.felix.tsv is 0.5-27 GB and must be read once, not once per locus.

    Spans are merged only where they actually overlap, which with REGION_HALF below half LOCUS_WINDOW cannot happen
    for two loci of one phenotype: widen REGION_HALF and it starts to.
    """
    plan: dict[str, list[dict]] = {}
    for run, g in loci.groupby("run", sort=True):
        spans: list[dict] = []
        for _, L in g.sort_values(["chrom", "window_start"]).iterrows():
            last = spans[-1] if spans else None
            if last and last["chrom"] == L["chrom"] and L["window_start"] <= last["end"]:
                last["end"] = max(last["end"], int(L["window_end"]))
                last["locus_id"] = f"{last['locus_id']},{L['locus_id']}"
            else:
                spans.append(dict(locus_id=L["locus_id"], chrom=int(L["chrom"]),
                                  start=int(L["window_start"]), end=int(L["window_end"])))
        plan[run] = spans
    return plan


# ------------------------------------------------------------------------------------- the checks that decide

# Thresholds. Each is a measured property of this scan, not a convention, and the comment says where it comes from.
CTRL_F_LOW = -0.3  # excess heterozygosity in controls: felix_phenotype_report.EXCESS_HET_F, the valley of the F distribution
CTRL_F_HIGH = 0.3  # heterozygote deficit. Measured over felix_sweep_all_hits.tsv (329-run sweep): of the 193 kept hits
# carrying an F, 44 exceed 0.05 and 4 exceed 0.3 -- DE_660.12 at chr7:129,570,069 (0.580), BI_181 chr18:70,042,366
# (0.549), CV_440 chr22:21,651,719 (0.447) and GE_972 chr19:54,111,497 (0.307); EM_232.2's lead sits just under at 0.299.
# Admixture alone (Wahlund) produces positive F, so this is a prompt to test F within local-ancestry strata, not a verdict.
COHERENCE_MAX = 5.0  # the K-df chi2 implied by P_het divided by the sum of squared per-ancestry z. Measured over the 67
# rows with P_het < 1e-3 in the five 1-Mb region extracts in data/ (BI_174, CA_103_22, CV_438, CV_446_2, EM_232_2): median
# 1.23, 95th percentile 3.29, 99th 4.66. Exactly one row exceeds 5 -- CV_446.2's lead at chr13:63,057,660, ratio 5.17.
# Read off that same extract, its P_cct of 1.33e-25 (P_het 6.66e-26) comes with P_hom 0.908 and per-ancestry p of 0.0043
# (SAS), 0.0066 (AMR), 0.0095 (EUR), 0.174 (AFR), 0.435 (EAS) -- which do not point the same way, let alone to 1e-25.
# Recheck the cut once the dossier export gives per-ancestry statistics for all 90 loci rather than five windows.
MISSING_MAX = 0.05  # per-ancestry genotype missingness above this is worth seeing
CALLRATE_DIFF_MAX = 0.02  # case minus control call rate at the lead, from the PRE-phasing callset
IMPUTED_DIFF_MAX = 0.02  # case minus control share of genotypes phasing imputed rather than observed
SCORE_P_MAX_LOG10 = 1.0  # orders of magnitude the reported p may differ from the normal tail of its own Tstat/sqrt(var)
HET_Q_P = 0.01  # Cochran's Q across local ancestries below this means the effect is not one effect
LD_PARTNER_AF = 0.05  # above this allele frequency a real signal should bring neighbours, so a lonely peak is a flag
LD_PARTNER_WINDOW = 100_000  # distance within which a supporting neighbour is looked for
LD_PARTNER_ORDERS = 3.0  # a neighbour counts as support when its p is within this many orders of magnitude of the lead
FIXED_HET_FRAC = 0.90  # control heterozygous fraction at or above which the record is a collapsed paralogue, not a variant
ALLELIC_MIN_TESTS = 3  # fewer split tests than this at one record is not an allelic series and the check does not apply
ALLELIC_ORDERS = 3.0  # a second allele of the series within this many orders of the lead counts as a graded effect


def _norm_sf(z: float) -> float:
    """Upper tail of the standard normal."""
    from math import erfc, sqrt

    return 0.5 * erfc(abs(z) / sqrt(2.0))


def chi2_sf(x: float, df: int) -> float:
    """Upper tail of a chi-square with ``df`` degrees of freedom; scipy when it is installed, else a series.

    The Terra notebook VM has no guaranteed scipy, and this is used for Hardy-Weinberg and Cochran's Q, so the
    fallback is the regularised upper incomplete gamma Q(df/2, x/2) by its standard series / continued fraction.
    """
    if x is None or (isinstance(x, float) and math.isnan(x)) or x < 0 or df < 1:
        return math.nan
    try:
        from scipy.stats import chi2 as _chi2

        return float(_chi2.sf(x, df))
    except ImportError:
        pass
    a, z = df / 2.0, x / 2.0
    if z == 0:
        return 1.0
    lg = math.lgamma(a)
    if z < a + 1.0:  # series for the lower tail P(a, z)
        term = total = 1.0 / a
        for n in range(1, 1000):
            term *= z / (a + n)
            total += term
            if abs(term) < abs(total) * 1e-15:
                break
        return max(0.0, 1.0 - total * math.exp(-z + a * math.log(z) - lg))
    b, c = z + 1.0 - a, 1e300  # continued fraction for the upper tail Q(a, z)
    d = 1.0 / b
    hh = d
    for i in range(1, 1000):
        an = -i * (i - a)
        b += 2.0
        d = an * d + b
        if abs(d) < 1e-300:
            d = 1e-300
        c = b + an / c
        if abs(c) < 1e-300:
            c = 1e-300
        d = 1.0 / d
        delta = d * c
        hh *= delta
        if abs(delta - 1.0) < 1e-15:
            break
    return min(1.0, max(0.0, math.exp(-z + a * math.log(z) - lg) * hh))


def hardy_weinberg(n_het: float, n_hom: float, n_hap: float) -> dict:
    """Hardy-Weinberg in one stratum from FELIX's counts: samples, allele frequency, observed and expected heterozygotes, F, p.

    ``n_het``/``n_hom`` count SAMPLES and ``n_hap`` counts HAPLOTYPES (2 x samples), which is how FELIX writes them.
    ``F = 1 - het_obs / (2pq)``: negative means more heterozygotes than Hardy-Weinberg allows (a collapsed paralogue
    at the extreme), positive means fewer (population structure, or a dropped allele). ``p`` is the 1-df chi-square
    goodness-of-fit over the three genotype classes, not the exact test, so read it as a scale rather than a verdict.
    """
    out = dict(n_samples=math.nan, af=math.nan, het_obs=math.nan, het_exp=math.nan, hom_obs=math.nan,
               hom_exp=math.nan, F=math.nan, chi2=math.nan, p=math.nan, het_frac=math.nan)
    n_het, n_hom, n_hap = _f(n_het), _f(n_hom), _f(n_hap)
    if math.isnan(n_het) or math.isnan(n_hom) or math.isnan(n_hap) or n_hap <= 0:
        return out
    n = n_hap / 2.0
    if n <= 0:
        return out
    p = (n_het + 2 * n_hom) / n_hap
    q = 1.0 - p
    out.update(n_samples=n, af=p, het_obs=n_het, hom_obs=n_hom, het_frac=n_het / n,
               het_exp=2 * p * q * n, hom_exp=p * p * n)
    if p <= 0 or q <= 0:
        return out
    out["F"] = 1.0 - (n_het / n) / (2 * p * q)
    exp = (q * q * n, 2 * p * q * n, p * p * n)
    obs = (n - n_het - n_hom, n_het, n_hom)
    if min(exp) > 0:
        out["chi2"] = sum((o - e) ** 2 / e for o, e in zip(obs, exp))
        out["p"] = chi2_sf(out["chi2"], 1)
    return out


def joint_coherence(lead: pd.Series | dict, num_ancs: int = len(FLARE_ORDER)) -> dict:
    """Whether the joint p is supported by the per-ancestry statistics it is built from.

    FELIX's ``P_het_admixed_c`` is a K-df joint score test over the local-ancestry strata, so the chi-square it implies
    should be close to the sum of squared per-ancestry z. ``ratio`` is implied / observed; above COHERENCE_MAX the
    headline p rests on something the strata do not show, which is how CV_446.2's 1.3e-25 lead was caught.
    ``n_strata`` counts the strata that actually contributed a z.
    """
    blk = ancestry_block(lead)
    strata = blk[blk.stratum != "ancALL"]
    zs = [z for z in strata["z"] if not math.isnan(z)]
    z2 = float(sum(z * z for z in zs)) if zs else math.nan
    p_het = _f(lead.get("P_het_admixed_c") if hasattr(lead, "get") else None)
    implied = math.nan
    if not math.isnan(p_het) and 0 < p_het < 1:
        implied = _chi2_isf(p_het, num_ancs)
    return dict(
        z2_sum=z2, n_strata=len(zs), p_het=p_het,
        p_hom=_f(lead.get("P_hom_admixed_c")), p_cct=_f(lead.get("P_cct_admixed_c")),
        implied_chi2=implied,
        ratio=(implied / z2) if (z2 and z2 > 0 and not math.isnan(implied)) else math.nan,
        min_anc_p=(min([p for p in strata["p"] if not math.isnan(p)], default=math.nan)),
    )


def _chi2_isf(p: float, df: int) -> float:
    """The chi-square quantile with upper tail ``p`` -- bisection on chi2_sf, which is enough for a ratio."""
    if not (0 < p < 1):
        return math.nan
    lo, hi = 0.0, 10.0
    while chi2_sf(hi, df) > p and hi < 1e7:
        hi *= 2
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if chi2_sf(mid, df) > p:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def score_p_coherence(lead: pd.Series | dict) -> dict:
    """How far each stratum's reported p is from the normal tail of its own score statistic Tstat/sqrt(var).

    A file-level check in the EasyQC tradition: p and the statistic it came from must agree. Saddlepoint correction
    legitimately moves them apart -- that is its purpose -- so the gap is only diagnostic in strata where SPA was NOT
    applied, and those are what ``max_log10_gap`` is taken over. ``n_spa`` and ``n_plain`` say how the strata split, so
    a locus whose every stratum was saddlepoint-corrected reports no gap rather than a spurious one.
    """
    blk = ancestry_block(lead)
    worst, worst_gap = "", math.nan
    n_spa = n_plain = 0
    for _, r in blk.iterrows():
        if math.isnan(r["score_z"]) or math.isnan(r["p"]) or r["p"] <= 0:
            continue
        if r["spa"]:
            n_spa += 1
            continue
        n_plain += 1
        implied = max(_norm_sf(r["score_z"]) * 2, 1e-300)
        gap = abs(math.log10(r["p"] / implied))
        if math.isnan(worst_gap) or gap > worst_gap:
            worst, worst_gap = r["ancestry"], gap
    return dict(worst_stratum=worst, max_log10_gap=worst_gap, n_spa=n_spa, n_plain=n_plain)


def effect_consistency(lead: pd.Series | dict) -> dict:
    """Cochran's Q over the local-ancestry effect estimates: one effect, or several?

    Inverse-variance weighted mean of BETA_c across the strata with a usable SE, Q on (n-1) df, and the share of
    strata whose effect points the same way as that mean. This is independent of FELIX's own P_het, which is a score
    test rather than a heterogeneity test.
    """
    blk = ancestry_block(lead)
    estimable = blk[(blk.stratum != "ancALL") & blk.beta.notna() & blk.se.notna() & (blk.se > 0)]
    # Thin strata are excluded rather than down-weighted: inverse-variance weighting already down-weights them, but
    # sparse-data bias inflates the POINT estimate too, so a thin stratum contributes a spuriously large (b - b_fe)
    # and manufactures heterogeneity. MS_723 chr3:134,785,621 is the case in point -- its EAS beta of -0.81 rests on
    # 2 case haplotypes and pulls the forest the opposite way from every other stratum.
    use = estimable[~estimable.thin]
    n_thin = int(estimable.thin.sum())
    if use.empty:
        return dict(beta_fe=math.nan, se_fe=math.nan, Q=math.nan, df=0, p=math.nan, same_sign=math.nan, n_strata=0,
                    n_thin=n_thin, n_estimable=len(estimable))
    w = 1.0 / use.se.to_numpy() ** 2
    b = use.beta.to_numpy()
    beta_fe = float((w * b).sum() / w.sum())
    q = float((w * (b - beta_fe) ** 2).sum())
    df = len(use) - 1
    return dict(beta_fe=beta_fe, se_fe=float(math.sqrt(1.0 / w.sum())), Q=q, df=df,
                p=(chi2_sf(q, df) if df >= 1 else math.nan),
                same_sign=float((np.sign(b) == np.sign(beta_fe)).mean()), n_strata=len(use),
                n_thin=n_thin, n_estimable=len(estimable))


def regional_support(region: pd.DataFrame | None, pos: int, lead_p: float, p_col: str = P_COL,
                     window: int = LD_PARTNER_WINDOW, orders: float = LD_PARTNER_ORDERS) -> dict:
    """Does the peak have neighbours? Counted from the regional extract, which is the only input that can answer it.

    ``n_support`` counts markers within ``window`` bp (excluding the lead's own position) whose p is within ``orders``
    of magnitude of the lead. ``best_neighbour_p`` and ``best_neighbour_bp`` describe the strongest of them.
    ``available`` is False when no extract has been produced for the locus, and also when the extract on hand does not
    reach the lead -- an extract cut for some other window says nothing about this one, and must not be read as
    'no support'.
    """
    out = dict(available=False, n_in_window=0, n_support=0, n_gw_sig=0, best_neighbour_p=math.nan,
               best_neighbour_bp=math.nan, n_markers=0)
    if region is None or len(region) == 0 or p_col not in region.columns:
        return out
    r = region.copy()
    r["_p"] = pd.to_numeric(r[p_col], errors="coerce")
    r = r[r["_p"].notna() & (r["_p"] > 0)]
    near = r[(r.POS - pos).abs() <= window]
    if near.empty:  # the extract does not cover this locus
        out["n_markers"] = len(r)
        return out
    others = near[near.POS != pos]
    # Within `orders` of the lead OR genome-wide significant in its own right, whichever is more permissive. The
    # relative rule alone fails for a very strong lead: at p 7e-17 it demands neighbours below 7e-14, so a locus with
    # hundreds of markers at 1e-12 would be called a lone peak (GE_983.2 in the HBB region did read that way).
    cut = max(lead_p * (10.0 ** orders), GW_P)
    sup = others[others["_p"] <= cut]
    best = sup.loc[sup["_p"].idxmin()] if len(sup) else None
    out.update(available=True, n_markers=len(r), n_in_window=len(near), n_support=len(sup),
               n_gw_sig=int((others["_p"] < GW_P).sum()),
               best_neighbour_p=(float(best["_p"]) if best is not None else math.nan),
               best_neighbour_bp=(int(abs(best["POS"] - pos)) if best is not None else math.nan))
    return out


def repeat_tract_window(context: pd.Series | dict | None) -> float:
    """Half-width in bp over which split tests are treated as one repeat tract: the tract's own length.

    locus_context exports the repeat's PERIOD and COPIES but not its interval, so the extent has to be approximated.
    Period x copies is the tract length, and the lead sits somewhere inside it, so +/- that length is a safe bound on
    where the same tract can reach. Exporting the simpleRepeat start/end in felix_09 section 6 would make this exact
    and is a one-line change on the next VM run; until then this is deliberately generous and paired with the
    length-changing-allele filter in ``allelic_series``, which is what keeps the generosity from pulling in SNVs.
    """
    if context is None:
        return 0.0
    if not _bool(context.get("in_tandem_repeat")):
        return 0.0
    period, copies = _f(context.get("repeat_period")), _f(context.get("repeat_copies"))
    if math.isnan(period) or math.isnan(copies):
        return 0.0
    return float(period * copies)


def allelic_series(region: pd.DataFrame | None, pos: int, lead_marker: str | None = None, p_col: str = P_COL,
                   min_tests: int = ALLELIC_MIN_TESTS, orders: float = ALLELIC_ORDERS,
                   window: float = 0.0) -> dict:
    """Does the effect grade across the allele series of a merged multiallelic record, or sit on one allele?

    felixla splits a merged record into one REF/ALT test per ALT, and every split test keeps the record's POS, so the
    series is the set of rows at the lead's own position in the regional extract. This is the within-site counterpart
    of ``regional_support``, and unlike it the reasoning survives a hypermutable repeat: recurrent mutation erodes LD
    with the flanking SNPs, so an empty 100 kb window around a short tandem repeat is uninformative, but it says
    nothing about why the neighbouring *lengths* of the same repeat should be flat. If repeat length drives risk, a
    25 bp allele being strongly associated while 23 bp and 27 bp are null is not what the mechanism predicts.

    Reports the size of the series, the lead's own length change, the second-best p within it, and how many alleles
    fall within ``orders`` of the lead. ``n_merged`` is the number of alleles in the lead's MarkerID, which is the
    number of contrasts the site offered -- a multiplicity the single headline p does not carry.
    """
    out = dict(available=False, n_tests=0, n_merged=0, lead_delta_bp=math.nan, second_p=math.nan,
               second_delta_bp=math.nan, n_within_orders=0)
    if lead_marker:
        out["n_merged"] = str(lead_marker).count(";") + 1
    if region is None or len(region) == 0 or p_col not in region.columns or "POS" not in region.columns:
        return out
    r = region.copy()
    r["_p"] = pd.to_numeric(r[p_col], errors="coerce")
    r = r[r["_p"].notna() & (r["_p"] > 0)]
    if window and window > 0:
        # A caller splits one tandem-repeat tract into records at adjacent coordinates, so grouping on exact POS
        # under-counts: at MS_723 chr3:134,785,621 the +24 bp allele is recorded one base along at 134,785,622 and
        # was missed entirely. Within a tract only LENGTH-CHANGING alleles belong to the series -- the SNVs that
        # share the window are not alternative lengths of anything -- and that filter is what keeps a 1.1 kb tract
        # from sweeping in dozens of unrelated markers.
        site = r[(r.POS - pos).abs() <= window]
        keep = [(lambda d: d is not None and d != 0)(length_change(m)) for m in site["MarkerID"]]
        site = site[pd.Series(keep, index=site.index)]
    else:
        site = r[r.POS == pos]
    if len(site) < min_tests:
        out["available"] = len(site) > 0
        out["n_tests"] = len(site)
        return out
    site = site.sort_values("_p")
    lead_row, rest = site.iloc[0], site.iloc[1:]

    def delta(marker) -> float:
        al = split_alleles(marker)
        return float(len(al[1]) - len(al[0])) if al else math.nan

    second = rest.iloc[0]
    out.update(available=True, n_tests=len(site), lead_delta_bp=delta(lead_row.get("MarkerID")),
               second_p=float(second["_p"]), second_delta_bp=delta(second.get("MarkerID")),
               n_within_orders=int((rest["_p"] <= float(lead_row["_p"]) * (10.0 ** orders)).sum()))
    if lead_marker:  # the lead row of the extract is the same record, so prefer the caller's full MarkerID
        out["n_merged"] = str(lead_marker).count(";") + 1
    elif isinstance(lead_row.get("MarkerID"), str):
        out["n_merged"] = lead_row["MarkerID"].count(";") + 1
    return out


def checks(locus: pd.Series | dict, lead: pd.Series | dict | None, context: pd.Series | dict | None,
           region: pd.DataFrame | None = None, callrate: pd.Series | dict | None = None,
           repeat_is_context: bool = False, calib: pd.DataFrame | None = None) -> list[dict]:
    """The flag chips for one locus: name, status ('pass', 'warn', 'fail' or 'n/a') and a one-line reading.

    'fail' is reserved for a record that cannot be a variant as reported (fixed heterozygosity, a joint p its own
    strata do not support). 'warn' marks something a reviewer must look at before believing the locus. 'n/a' means the
    input is not in data/ yet, and is deliberately distinct from 'pass' -- an absent regional extract is not support.
    """
    out: list[dict] = []

    def add(name, status, detail):
        out.append(dict(check=name, status=status, detail=detail))

    # First, because it is about the scan rather than the variant: if the phenotype's own null is miscalibrated then
    # no p-value it produced can be read at face value, however clean the locus looks. The lambda is the lead arm's,
    # not the main arm's -- the repeat arms test ~1.6M markers against ~25M and their lambdas are not interchangeable.
    arm = str(locus.get("lead_arm", "main"))
    lam = _f(locus.get(f"lambda_gc_{arm}"))
    lc = lambda_context(calib, str(locus.get("run", "")), arm)
    where = f" ({_ordinal(lc['pct'])} percentile of {lc['n_runs']} {arm}-arm runs)" if lc["available"] else ""
    if math.isnan(lam):
        add("Phenotype calibration", "n/a", f"no lambda_gc_{arm} for this run")
    elif lam < LAMBDA_MIN:
        add("Phenotype calibration", "warn",
            f"lambda {lam:.3f} is deflated below {LAMBDA_MIN:.2f}{where}: the test is conservative, so the p is not "
            f"calibrated either")
    elif lam > LAMBDA_MAX:
        add("Phenotype calibration", "warn",
            f"lambda {lam:.3f} is inflated above {LAMBDA_MAX:.2f}{where}: the whole scan is suspect, not this locus")
    else:
        add("Phenotype calibration", "pass", f"lambda {lam:.3f}{where}")

    ctrl_f = _f(locus.get("ctrl_F"))
    het_frac = _f(locus.get("ctrl_het_frac"))
    if math.isnan(ctrl_f):
        add("Hardy-Weinberg", "n/a", "no control genotype counts (dosage arm, or MAF too low to measure)")
    elif het_frac >= FIXED_HET_FRAC:
        add("Hardy-Weinberg", "fail", f"{het_frac:.1%} of controls heterozygous: a collapsed paralogue, not a variant")
    elif ctrl_f < CTRL_F_LOW:
        add("Hardy-Weinberg", "fail", f"control F {ctrl_f:+.3f}: more heterozygotes than Hardy-Weinberg allows")
    elif ctrl_f > CTRL_F_HIGH:
        add("Hardy-Weinberg", "warn", f"control F {ctrl_f:+.3f}: heterozygote deficit; test F within ancestry strata")
    else:
        add("Hardy-Weinberg", "pass", f"control F {ctrl_f:+.3f}")

    if lead is None:
        add("Joint-test coherence", "n/a", "lead per-ancestry statistics not exported yet")
        add("Effect consistency", "n/a", "lead per-ancestry statistics not exported yet")
        add("p vs score statistic", "n/a", "lead per-ancestry statistics not exported yet")
    else:
        co = joint_coherence(lead)
        if math.isnan(co["ratio"]):
            add("Joint-test coherence", "n/a", "no per-ancestry z or no P_het")
        elif co["ratio"] > COHERENCE_MAX:
            add("Joint-test coherence", "fail",
                f"P_het implies chi2 {co['implied_chi2']:.0f} but the strata sum to {co['z2_sum']:.0f} "
                f"(ratio {co['ratio']:.1f})")
        else:
            add("Joint-test coherence", "pass", f"implied / observed chi2 = {co['ratio']:.2f}")

        ec = effect_consistency(lead)
        thin = (f"; {ec.get('n_thin', 0)} stratum/strata set aside under {MIN_STRATUM_CARRIERS} case carriers"
                .replace("stratum/strata", "stratum" if ec.get("n_thin") == 1 else "strata")) if ec.get("n_thin") else ""
        if ec["n_strata"] < 2:
            add("Effect consistency", "n/a",
                f"fewer than two strata with a believable effect estimate{thin}")
        elif not math.isnan(ec["p"]) and ec["p"] < HET_Q_P:
            add("Effect consistency", "warn",
                f"Cochran Q {ec['Q']:.1f} on {ec['df']} df (p {ec['p']:.1e}): the strata do not share one effect{thin}")
        else:
            add("Effect consistency", "pass",
                f"Q p {ec['p']:.2f}, {ec['same_sign']:.0%} of strata agree in sign{thin}")

        sc = score_p_coherence(lead)
        if sc["n_plain"] == 0:
            add("p vs score statistic", "pass" if sc["n_spa"] else "n/a",
                f"saddlepoint applied in all {sc['n_spa']} strata, so p is expected to differ from the normal tail"
                if sc["n_spa"] else "no Tstat / var")
        elif sc["max_log10_gap"] > SCORE_P_MAX_LOG10:
            add("p vs score statistic", "warn",
                f"{sc['worst_stratum']}: p differs from its own score statistic by 10^{sc['max_log10_gap']:.1f} "
                f"with no saddlepoint correction")
        else:
            add("p vs score statistic", "pass",
                f"agrees to 10^{sc['max_log10_gap']:.1f} in the {sc['n_plain']} stratum/strata without correction".replace("stratum/strata", "stratum" if sc["n_plain"] == 1 else "strata"))

    add("Call rate", *_callrate_chip(callrate, lead))

    rs = regional_support(region, int(locus["pos"]), float(locus["p"]))
    af = math.nan
    if lead is not None:
        af = _f(lead.get("AF_Allele2_ancALL"))
    if not rs["available"]:
        add("Regional support", "n/a", "no region extract for this locus yet")
    elif rs["n_support"] > 0:
        gw = f", {rs['n_gw_sig']} of them genome-wide significant" if rs.get("n_gw_sig") else ""
        add("Regional support", "pass",
            f"{rs['n_support']} supporting neighbour(s) within {LD_PARTNER_WINDOW // 1000} kb"
            f"{gw}; best p {rs['best_neighbour_p']:.1e}")
    elif not math.isnan(af) and min(af, 1 - af) >= LD_PARTNER_AF:
        add("Regional support", "warn",
            f"no neighbour within {LD_PARTNER_ORDERS:.0f} orders of the lead, at allele frequency {af:.2f}")
    else:
        add("Regional support", "pass", "lone peak, as expected for a rare allele")

    al = allelic_series(region, int(locus["pos"]),
                        lead_marker=(lead.get("MarkerID") if lead is not None else locus.get("lead_marker")),
                        window=repeat_tract_window(context))
    if not al["available"]:
        add("Allelic series", "n/a", "no region extract for this locus yet")
    elif al["n_tests"] < ALLELIC_MIN_TESTS:
        add("Allelic series", "n/a",
            f"{al['n_tests']} allele(s) tested at this record: no series to grade the effect across")
    elif al["n_within_orders"] > 0:
        add("Allelic series", "pass",
            f"{al['n_tests']} alleles tested here, {al['n_within_orders']} within {ALLELIC_ORDERS:.0f} orders of the "
            f"lead: the effect grades across the series")
    else:
        d = al["lead_delta_bp"]
        lead_bit = f"the {d:+.0f} bp allele" if not math.isnan(d) else "one allele"
        add("Allelic series", "warn",
            f"{al['n_tests']} alleles tested here but only {lead_bit} is associated "
            f"(next best p {al['second_p']:.1e}): a length effect should not skip its own neighbours")

    if context is None or str(context.get("cytoband", "")) == "":
        add("Genomic context", "n/a", "context not annotated")
    else:
        flagging, noted = [], []
        if context.get("pericentromeric"):
            flagging.append(f"{context.get('stain')} band")
        if context.get("in_segdup"):
            flagging.append("inside a segmental duplication")
        elif _f(context.get("segdup_frac")) > 0:
            flagging.append(f"{_f(context.get('segdup_frac')):.0%} of the flank duplicated")
        if context.get("in_tandem_repeat"):
            # In a document restricted to SVs and repeats the tandem repeat is the INCLUSION CRITERION, so warning on
            # it says nothing the reviewer did not choose -- it fired on 44 of the 49 structural loci. The period is
            # still worth printing, so it moves to the detail text and stops driving the status.
            period = _f(context.get("repeat_period"))
            (noted if repeat_is_context else flagging).append(
                "in a tandem repeat" if math.isnan(period) else f"in a {period:.0f} bp tandem repeat")
        bits = flagging + noted
        add("Genomic context", "warn" if flagging else "pass",
            "; ".join(bits) if bits else f"unique sequence in {context.get('cytoband')}")

    return out


def _callrate_chip(callrate: pd.Series | dict | None, lead: pd.Series | dict | None = None) -> tuple[str, str]:
    """(status, detail) for the call-rate chip from a felix_dossier_callrate.tsv row.

    This check does not depend on the per-ancestry export: FELIX reads phase_vcfs = flare_lai_prod.gt_vcf, the ligated
    phased callset, where phasing has filled every genotype, so its MissingRate is 0 at every record by construction --
    measured as exactly 0 in all 331,518 MissingRate cells of the five region extracts in data/. Reporting that as a
    pass would claim a check that did not run, so the real numbers come from pull_site_fields reading the pre-phasing
    callset at the lead position (felix_09, section 8).

    Reads the PRE-phasing call rate, because the phased callset's is 1.0 by construction, and the imputed fraction --
    the share of the genotypes the scan used that phasing filled in. A case/control gap in either is what matters:
    CALLRATE_DIFF_MAX on the call rate, IMPUTED_DIFF_MAX on the imputed share. ``lead`` is read only so that a
    genuinely non-zero MissingRate, which would mean the scan's input changed, is reported rather than ignored.
    """
    if callrate is None:
        if lead is not None:
            miss = ancestry_block(lead).missing.max()
            if not math.isnan(miss) and miss > 0:
                return ("warn" if miss > MISSING_MAX else "pass",
                        f"FELIX reports up to {miss:.2%} missing, so its input is no longer fully phased")
        return "n/a", "phased callset has no no-calls; needs the pre-phasing VCF (felix_09 section 8)"
    cr_case, cr_ctrl = _f(callrate.get("pre_call_rate_case")), _f(callrate.get("pre_call_rate_control"))
    imp_case, imp_ctrl = _f(callrate.get("imputed_frac_case")), _f(callrate.get("imputed_frac_control"))
    if math.isnan(cr_case) and math.isnan(imp_case):
        note = str(callrate.get("note") or "").strip()
        return "n/a", note or "no pre-phasing genotypes pulled for this locus"
    bits, worst = [], "pass"
    if not (math.isnan(cr_case) or math.isnan(cr_ctrl)):
        diff = cr_case - cr_ctrl
        bits.append(f"call rate {cr_case:.3f} in cases vs {cr_ctrl:.3f} in controls")
        if abs(diff) > CALLRATE_DIFF_MAX:
            worst = "warn"
            bits[-1] += f" (differs by {diff:+.3f})"
    if not (math.isnan(imp_case) or math.isnan(imp_ctrl)):
        diff = imp_case - imp_ctrl
        bits.append(f"{imp_case:.1%} of case and {imp_ctrl:.1%} of control genotypes imputed by phasing")
        if abs(diff) > IMPUTED_DIFF_MAX:
            worst = "warn"
    return worst, "; ".join(bits)


def verdict(flags: list[dict]) -> str:
    """'fail' if any check failed, else 'warn' if any warned, else 'incomplete' if any is n/a, else 'pass'."""
    st = {f["status"] for f in flags}
    for level in ("fail", "warn"):
        if level in st:
            return level
    return "incomplete" if "n/a" in st else "pass"
