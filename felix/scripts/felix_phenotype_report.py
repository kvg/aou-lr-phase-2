"""Per-phenotype FELIX result review — analysis functions.

Reads one FelixGenome `run_summary_tar` bundle and produces:
  * a calibration table across the three arms (main / ru / ru_biallelic)
  * lambda recomputed independently from the arm-comparison p-values
  * the SPA-fallback screen applied to every genome-wide hit
  * case-carrier counts at each hit
  * a side-by-side QQ (as reported vs screened) for a chosen arm
  * repeat-genotyping QC rolled up across chromosomes

Everything works off the bundle alone. Where the bundle cannot answer a
question (it carries only the top 50 loci per arm) the report says so.
"""
from __future__ import annotations

import io
import json
import os
import re
import subprocess
import tarfile
from pathlib import Path

import numpy as np
import pandas as pd

ARMS = ("main", "ru", "ru_biallelic")
SPA_CUTOFF = 2.0
P_GW = 5e-8
CHI2_MEDIAN = 0.4549364231195733


# ---------------------------------------------------------------- fetching

def fetch(uri: str, cache_dir="bundle_cache") -> Path:
    """gs:// URI -> local path, cached. A local path is returned unchanged."""
    if not str(uri).startswith("gs://"):
        return Path(uri)
    Path(cache_dir).mkdir(parents=True, exist_ok=True)
    dest = Path(cache_dir) / uri.rsplit("/", 1)[-1].replace("/", "_")
    if dest.exists() and dest.stat().st_size > 0:
        return dest
    cmd = ["gsutil"]
    proj = os.environ.get("GOOGLE_PROJECT")
    if proj:                      # AoU buckets are requester-pays
        cmd += ["-u", proj]
    cmd += ["cp", uri, str(dest)]
    subprocess.run(cmd, check=True, capture_output=True)
    return dest


def load_bundle(uri_or_path, cache_dir="bundle_cache") -> dict:
    """Parse the members we need out of felix_run.summary.tar.gz."""
    path = fetch(uri_or_path, cache_dir)
    out = {"arms": {}, "repeat_qc": {}, "arm_compare": {}, "source": str(uri_or_path)}
    with tarfile.open(path) as t:
        names = t.getnames()
        for arm in ARMS:
            a = {}
            for member, key in ((f"./{arm}/calibration_summary.tsv", "calibration"),
                                (f"./{arm}/results_manifest.tsv", "manifest")):
                if member in names:
                    a[key] = pd.read_csv(t.extractfile(member), sep="\t")
            hits = [n for n in names if n.startswith(f"./{arm}/qc/") and n.endswith("top_hits.tsv")]
            if hits:
                a["top_hits"] = pd.read_csv(t.extractfile(hits[0]), sep="\t", low_memory=False)
                a["phenotype"] = hits[0].split("/")[3]
            a["figures"] = {Path(n).stem: n for n in names
                            if n.startswith(f"./{arm}/qc/") and n.endswith(".png")}
            out["arms"][arm] = a
        for member, key in (("./arm_compare/arms_summary.json", "summary"),
                            ("./arm_compare/arms_table.tsv.gz", "table_member")):
            if member in names:
                out["arm_compare"][key] = (json.load(t.extractfile(member))
                                           if key == "summary" else member)
        for n in names:
            m = re.match(r"\./repeat_qc/chr(\w+)\.(loci\.summary|ru\.admixed\.stats|"
                         r"ru_biallelic\.admixed\.stats)\.json$", n)
            if m:
                out["repeat_qc"].setdefault(m.group(2), {})[m.group(1)] = json.load(t.extractfile(n))
    out["_tar"] = str(path)
    return out


def arms_table(bundle: dict) -> pd.DataFrame:
    """The full four-arm comparison table (1.6M rows) — loaded on demand."""
    member = bundle["arm_compare"].get("table_member")
    if member is None:
        return pd.DataFrame()
    with tarfile.open(bundle["_tar"]) as t:
        comp = "gzip" if member.endswith(".gz") else None
        return pd.read_csv(t.extractfile(member), sep="\t",
                           compression=comp, low_memory=False)


# ---------------------------------------------------------------- screens

def as_bool(s: pd.Series) -> pd.Series:
    if s.dtype == bool:
        return s
    return s.astype(str).str.strip().str.lower().map(
        {"true": True, "false": False, "1": True, "0": False})


def _strata(num_ancs: int):
    return [f"anc{k}" for k in range(1, num_ancs + 1)] + ["ancALL"]


def spa_fallback(df: pd.DataFrame, num_ancs=5, cutoff=SPA_CUTOFF):
    """Tests where SPA was attempted but the uncorrected p was reported.

    Two detectors, OR-ed:
      unconditioned  Is.SPA == false      AND |Tstat|/sqrt(var)   > cutoff
      conditioned    p.value_c == p.value.NA_c AND |Tstat_c|/sqrt(var_c) > cutoff

    The conditioned one matters because the documented headline column is
    P_cct_admixed_c and FELIX emits no Is.SPA_c flag for it.
    """
    mask = pd.Series(False, index=df.index)
    which = pd.Series("", index=df.index, dtype=object)
    for a in _strata(num_ancs):
        tc, vc, sc = f"Tstat_{a}", f"var_{a}", f"Is.SPA_{a}"
        if {tc, vc, sc} <= set(df.columns):
            t = pd.to_numeric(df[tc], errors="coerce")
            v = pd.to_numeric(df[vc], errors="coerce")
            std = t.abs() / np.sqrt(v.where(v > 0))
            hit = ((as_bool(df[sc]) == False) & (std > cutoff)).fillna(False)  # noqa: E712
            mask |= hit
            which = which.where(~hit | which.str.contains("uncond"), which + "uncond;")
        tcc, vcc = f"Tstat_c_{a}", f"var_c_{a}"
        pc, pcn = f"p.value_c_{a}", f"p.value.NA_c_{a}"
        if {tcc, vcc, pc, pcn} <= set(df.columns):
            t = pd.to_numeric(df[tcc], errors="coerce")
            v = pd.to_numeric(df[vcc], errors="coerce")
            a1 = pd.to_numeric(df[pc], errors="coerce")
            a2 = pd.to_numeric(df[pcn], errors="coerce")
            std = t.abs() / np.sqrt(v.where(v > 0))
            hit = ((a1 == a2) & a1.notna() & (std > cutoff)).fillna(False)
            mask |= hit
            which = which.where(~hit | which.str.contains("cond"), which + "cond;")
    return mask, which.str.rstrip(";").replace("", ".")


def case_carriers(df: pd.DataFrame) -> pd.Series:
    het = pd.to_numeric(df.get("N_case_het_ancALL"), errors="coerce").fillna(0)
    hom = pd.to_numeric(df.get("N_case_hom_ancALL"), errors="coerce").fillna(0)
    return (het + hom).astype(int)


def lambda_gc(p) -> float:
    from scipy.stats import chi2
    p = np.asarray(pd.to_numeric(pd.Series(p), errors="coerce"), dtype=float)
    p = p[np.isfinite(p) & (p > 0) & (p <= 1)]
    return float(np.median(chi2.isf(p, 1)) / CHI2_MEDIAN) if p.size else float("nan")


# ---------------------------------------------------------------- reports

def calibration_report(bundle: dict) -> pd.DataFrame:
    rows = []
    for arm in ARMS:
        a = bundle["arms"].get(arm, {})
        cal = a.get("calibration")
        if cal is None or not len(cal):
            rows.append(dict(arm=arm, status="missing"))
            continue
        c = cal.iloc[0]
        rows.append(dict(arm=arm, status="ok", phenotype=c.get("phenotype"),
                         n_samples=int(c["n_samples"]), n_cases=int(c["n_cases"]),
                         n_variants=int(c["n_variants"]),
                         lambda_gc=round(float(c["lambda_gc_cct"]), 4),
                         n_gw_sig=int(c["n_gw_sig"]),
                         N_eff=int(round(4 * c["n_samples"] *
                                         (c["n_cases"] / c["n_samples"]) *
                                         (1 - c["n_cases"] / c["n_samples"])))))
    return pd.DataFrame(rows)


def hit_report(bundle: dict, num_ancs=5, p_col="P_cct_admixed_c") -> pd.DataFrame:
    rows = []
    for arm in ARMS:
        a = bundle["arms"].get(arm, {})
        th = a.get("top_hits")
        if th is None or p_col not in th.columns:
            continue
        p = pd.to_numeric(th[p_col], errors="coerce")
        mask, which = spa_fallback(th, num_ancs)
        cc = case_carriers(th)
        for i in th.index[p < P_GW]:
            rows.append(dict(arm=arm, MarkerID=th.at[i, "MarkerID"], CHR=th.at[i, "CHR"],
                             POS=th.at[i, "POS"], p=p[i], case_carriers=int(cc[i]),
                             spa_fallback=bool(mask[i]), detector=which[i],
                             verdict="ARTIFACT" if mask[i] else "keep"))
    return pd.DataFrame(rows)


def coverage_note(bundle: dict, p_col="P_cct_admixed_c") -> list:
    """Hits the bundle cannot screen: top_hits holds only the strongest N loci.

    Compares like with like — `n_gw_sig` in calibration_summary is counted on
    that arm's own `joint_p_col` (usually P_cct_admixed), which is NOT always
    the column used for screening (P_cct_admixed_c), so both are reported.
    """
    notes = []
    for arm in ARMS:
        a = bundle["arms"].get(arm, {})
        cal, th = a.get("calibration"), a.get("top_hits")
        if cal is None or th is None or not len(cal):
            continue
        c = cal.iloc[0]
        jcol = str(c.get("joint_p_col", "P_cct_admixed"))
        reported = int(c["n_gw_sig"])
        in_file_j = int((pd.to_numeric(th.get(jcol), errors="coerce") < P_GW).sum()) \
            if jcol in th.columns else -1
        if reported > in_file_j >= 0:
            notes.append(
                f"{arm}: {reported} hits on {jcol} but only {in_file_j} are in the "
                f"bundle's top-{len(th)} file — {reported - in_file_j} NOT screened "
                f"(needs the merged TSV)")
        if jcol != p_col and p_col in th.columns:
            n_alt = int((pd.to_numeric(th[p_col], errors="coerce") < P_GW).sum())
            if n_alt != in_file_j:
                notes.append(f"{arm}: {in_file_j} hits on {jcol} vs {n_alt} on {p_col} "
                             f"(screening used {p_col})")
    return notes


def qq_screen_figure(bundle: dict, arm="ru", p_col="P_la_dosage", num_ancs=5, out_png=None):
    """Side-by-side QQ: as reported, and with SPA-fallback tests removed.

    Uses the genome-wide arm-comparison table for the p-values and attaches the
    SPA columns from top_hits, which exist only for the strongest loci — so the
    screen is applied where the bundle can support it and declines elsewhere.
    Returns (fig, summary_dict).
    """
    import matplotlib.pyplot as plt
    at = arms_table(bundle)
    th = bundle["arms"].get(arm, {}).get("top_hits")
    if at.empty or th is None or p_col not in at.columns:
        raise ValueError(f"cannot build QQ: need {p_col} in arms_table and {arm} top_hits")
    keep = ["MarkerID"] + [f"{s}_{a}" for a in _strata(num_ancs)
                           for s in ("Tstat", "var", "Is.SPA", "Tstat_c", "var_c",
                                     "p.value_c", "p.value.NA_c") if f"{s}_{a}" in th.columns]
    merged = at[["MarkerID", p_col]].merge(th[keep], on="MarkerID", how="left")
    mask, _ = spa_fallback(merged, num_ancs)

    p = pd.to_numeric(merged[p_col], errors="coerce")
    ok = np.isfinite(p) & (p > 0) & (p <= 1)
    pv, mk = p[ok].to_numpy(float), mask[ok].to_numpy(bool)
    order = np.argsort(pv)
    obs, flag = -np.log10(pv[order]), mk[order]
    n = obs.size
    exp = -np.log10((np.arange(1, n + 1) - 0.5) / n)
    kept = np.sort(pv[~mk])
    obs_k = -np.log10(kept)
    exp_k = -np.log10((np.arange(1, kept.size + 1) - 0.5) / kept.size)

    def thin(e, o, keep_above=2.0, budget=20000):
        it, ib = np.flatnonzero(o > keep_above), np.flatnonzero(o <= keep_above)
        b = max(budget - it.size, 1)
        if ib.size > b:
            step = -(-ib.size // b)
            sel = ib[::step]
            if sel[-1] != ib[-1]:
                sel = np.append(sel, ib[-1])
            ib = sel
        i = np.union1d(it, ib)
        return e[i], o[i]

    lam_a, lam_k = lambda_gc(pv), lambda_gc(kept)
    full = max(obs.max(), exp.max()) * 1.05
    zoom = max(np.percentile(obs, 99.999), obs_k.max() if obs_k.size else 0, exp.max()) * 1.10
    fig, ax = plt.subplots(2, 2, figsize=(7.2, 6.4))
    for row, lim in ((0, full), (1, zoom)):
        for col, (e, o, f) in enumerate(((exp, obs, flag), (exp_k, obs_k, None))):
            A = ax[row][col]
            A.plot([0, lim], [0, lim], ls="--", lw=0.8, color="#888888")
            if f is None:
                ee, oo = thin(e, o)
                A.scatter(ee, oo, s=5, color="#3b6ea5", linewidths=0)
            else:
                ee, oo = thin(e[~f], o[~f])
                A.scatter(ee, oo, s=5, color="#3b6ea5", linewidths=0)
                A.scatter(e[f], o[f], s=18, facecolors="none", edgecolors="#d1495b",
                          linewidths=0.9,
                          label="SPA correction failed" if row == 0 else None)
                if row == 0:
                    A.legend(frameon=False, loc="upper left", handletextpad=0.3)
            A.set_xlim(0, max(e.max(), 1.0) * 1.08)
            A.set_ylim(0, lim)
            if row == 1:
                A.set_xlabel(r"Expected $-\log_{10}(p)$")
            if col == 0:
                A.set_ylabel(r"Observed $-\log_{10}(p)$")
    ax[0][0].set_title("As reported by FELIX")
    ax[0][1].set_title("SPA-fallback tests removed")
    pheno = bundle["arms"][arm].get("phenotype", "")
    fig.suptitle(f"{pheno} — {arm} arm ({p_col})", x=0.01, ha="left")
    fig.tight_layout()
    if out_png:
        fig.savefig(out_png, dpi=200)
    return fig, dict(arm=arm, p_col=p_col, n_tested=int(n), n_flagged=int(mk.sum()),
                     lambda_as_reported=lam_a, lambda_screened=lam_k,
                     n_sig_as_reported=int((pv < P_GW).sum()),
                     n_sig_screened=int((kept < P_GW).sum()),
                     n_rows_with_spa_cols=int(merged["Is.SPA_ancALL"].notna().sum())
                     if "Is.SPA_ancALL" in merged.columns else 0)


def repeat_qc_report(bundle: dict) -> pd.DataFrame:
    loci = bundle["repeat_qc"].get("loci.summary", {})
    if not loci:
        return pd.DataFrame()
    rows = []
    for chrom, d in loci.items():
        c = d["counts"]
        rows.append(dict(chrom=f"chr{chrom}", n_samples=d["n_samples"],
                         loci_written=c["loci_written"],
                         unassigned=c["records_unassigned"],
                         hap_missing_frac=d["hap_missing_frac"],
                         dup_frac=d["hap_possible_duplicate_frac"],
                         anc_inconsistent_frac=d["hap_ancestry_inconsistent_frac"]))
    df = pd.DataFrame(rows)
    df["_k"] = df.chrom.str.replace("chr", "", regex=False).astype(int)
    return df.sort_values("_k").drop(columns="_k").reset_index(drop=True)
