"""Per-phenotype FELIX result review — analysis functions.

Reads one FelixGenome `run_summary_tar` bundle and produces:
  * a calibration table across the three arms (main / ru / ru_biallelic)
  * lambda recomputed independently from the arm-comparison p-values
  * the SPA-fallback screen applied to every genome-wide hit
  * case-carrier counts at each hit
  * the excess-heterozygosity screen (controls vs Hardy-Weinberg) at each hit
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


# Arms whose N_*_het / N_*_hom are genotype counts. The repeat-dosage arm ("ru") reports dosage sums there (its
# AF columns are mean signed repeat units, negative for contractions), so Hardy-Weinberg does not apply to it.
HET_ARMS = ("main", "ru_biallelic")
# Controls' inbreeding coefficient below which a hit is excluded. F = 1 - (observed heterozygotes) / (2pq), so
# -0.3 means 30% more heterozygotes than Hardy-Weinberg expects. Chosen from the distribution of F over the
# 5,430 common (MAF >= 1/3) main-arm markers in five 1-Mb windows around sweep hits: the bulk lies in (0, 0.5], a
# separate cluster of fixed-heterozygous sites lies in (-1, -0.5], and the valley between them (-0.4 to -0.1)
# holds 1.3% of the markers. In the two windows with no such cluster the lowest F is -0.20 and -0.06. Real
# heterozygote excess of this size does not arise from population structure (admixture raises homozygosity,
# F > 0), so a cut anywhere in the valley gives the same result for the hits seen so far. Revisit it against the
# ctrl_F column of a full sweep (felix_08, "Genotype check").
EXCESS_HET_F = -0.3


def excess_het(df: pd.DataFrame, arm: str = "main", f_cut: float = EXCESS_HET_F):
    """Controls with more heterozygotes than Hardy-Weinberg allows: a record every sample carries once.

    Returns (flag, F, het_frac), each aligned to ``df.index``. ``het_frac`` is the share of controls that are
    heterozygous; ``F = 1 - het_frac / (2pq)`` with p the control allele frequency rebuilt from the same het and hom
    counts. ``flag`` is ``F < f_cut``.

    Why it is a screen. Under Hardy-Weinberg at most half the controls are heterozygous, whatever the allele
    frequency. A record in which 98-100% of cases and 99% of controls are heterozygous, allele frequency 0.5 in
    every ancestry, is a collapsed paralogue or a reference artefact, not a polymorphism. It also breaks the test:
    every sample then carries exactly one alternate allele, the per-ancestry dosages sum to a constant, and FELIX's
    K-df joint score (P_het, hence P_cct) is singular. Its p-value can be 1e-75 with per-ancestry statistics that
    support nothing like it.

    Scope. Only arms in ``HET_ARMS`` (the dosage arm has no genotype counts), and only where the control counts are
    in the table; otherwise nothing is flagged and F is NaN, so the caller can see the screen did not run. The
    bound het <= 2p means F < -0.3 needs MAF >= 0.23: a low-frequency record cannot show excess heterozygosity by
    this measure, and this screen does not claim to catch one.
    """
    off = pd.Series(False, index=df.index)
    nan = pd.Series(np.nan, index=df.index)
    need = {"N_ctrl_ancALL", "N_ctrl_het_ancALL", "N_ctrl_hom_ancALL"}
    if arm not in HET_ARMS or not need <= set(df.columns):
        return off, nan, nan
    n = (pd.to_numeric(df["N_ctrl_ancALL"], errors="coerce") / 2.0).where(lambda s: s > 0)  # haplotypes -> samples
    het = pd.to_numeric(df["N_ctrl_het_ancALL"], errors="coerce")
    hom = pd.to_numeric(df["N_ctrl_hom_ancALL"], errors="coerce")
    p = (het + 2 * hom) / (2 * n)
    expected = (2 * p * (1 - p)).where(lambda s: s > 0)
    het_frac = het / n
    f = 1 - het_frac / expected
    return (f < f_cut).fillna(False), f, het_frac


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
    """Every genome-wide hit, with the reasons it cannot be believed.

    Four independent exclusions, OR-ed into the verdict:

    spa_fallback      SPA was attempted and bailed, so the reported p is the
                      uncorrected normal-approximation score p.
    no_case_carriers  N_case_het + N_case_hom == 0. No case carries the allele,
                      so there is nothing for the phenotype to associate with,
                      whatever the p-value says. This needs no statistics.
    impossible_p      p is exactly 0. Underflowing a double takes chi2 ~ 1425
                      at 1 df; a scan this size cannot produce that, so it
                      marks a broken statistic rather than a strong hit.
    excess_het        the controls are more heterozygous than Hardy-Weinberg
                      allows (F < EXCESS_HET_F, see excess_het). A record every
                      sample carries once is a paralogue or reference artefact and
                      its joint test is singular. main and ru_biallelic only.

    Earlier versions reported case_carriers as a column but let only
    spa_fallback decide. Over 70 phenotypes that left 10 zero-carrier loci and
    2 exact-zero p-values among 34 "retained" hits, so all three now gate it.
    """
    rows = []
    for arm in ARMS:
        a = bundle["arms"].get(arm, {})
        th = a.get("top_hits")
        if th is None or p_col not in th.columns:
            continue
        p = pd.to_numeric(th[p_col], errors="coerce")
        mask, which = spa_fallback(th, num_ancs)
        cc = case_carriers(th)
        het_flag, ctrl_f, het_frac = excess_het(th, arm)
        for i in th.index[p < P_GW]:
            reasons = []
            if mask[i]:
                reasons.append("spa_fallback")
            if int(cc[i]) == 0:
                reasons.append("no_case_carriers")
            if p[i] == 0:
                reasons.append("impossible_p")
            if het_flag[i]:
                reasons.append("excess_het")
            rows.append(dict(arm=arm, MarkerID=th.at[i, "MarkerID"], CHR=th.at[i, "CHR"],
                             POS=th.at[i, "POS"], p=p[i], case_carriers=int(cc[i]),
                             spa_fallback=bool(mask[i]), detector=which[i],
                             ctrl_het_frac=het_frac[i], ctrl_F=ctrl_f[i],
                             reasons=",".join(reasons) or ".",
                             verdict="EXCLUDE" if reasons else "keep"))
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


# ---------------------------------------------------------------- sweep

SMALL_MEMBERS = ("./*/calibration_summary.tsv", "./*/qc/*/top_hits.tsv")


def harvest(uri: str, run_id: str, store="harvest_cache") -> Path:
    """Extract only the small members of a run bundle into a compact cache.

    Streams the object through tar so the 125 MB tarball never lands on disk.
    Returns the per-run cache directory; re-harvesting is a no-op.
    """
    out = Path(store) / run_id
    if (out / ".done").exists():
        return out
    out.mkdir(parents=True, exist_ok=True)

    streamed = False
    if str(uri).startswith("gs://"):
        # Stream through tar so the 125 MB tarball never lands on disk.
        try:
            proj = os.environ.get("GOOGLE_PROJECT")
            cat = ["gsutil"] + (["-u", proj] if proj else []) + ["cat", uri]
            tar = ["tar", "xz", "--wildcards", "-C", str(out)] + list(SMALL_MEMBERS)
            p1 = subprocess.Popen(cat, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
            p2 = subprocess.Popen(tar, stdin=p1.stdout, stderr=subprocess.DEVNULL)
            p1.stdout.close()
            p2.communicate()
            p1.wait()
            streamed = (p2.returncode == 0 and any(out.rglob("calibration_summary.tsv")))
        except (OSError, subprocess.SubprocessError):
            streamed = False      # no gsutil, no GNU tar, or a broken pipe

    if not streamed:
        # Local path, or the stream did not work: extract members directly.
        local = fetch(uri, cache_dir=str(Path(store) / "_tmp"))
        with tarfile.open(local) as t:
            want = [m for m in t.getmembers()
                    if m.name.endswith(("calibration_summary.tsv", "top_hits.tsv"))]
            t.extractall(out, members=want, filter="data")
        if str(uri).startswith("gs://"):
            local.unlink(missing_ok=True)

    (out / ".done").write_text("ok")
    return out


def read_harvest(run_dir: Path) -> dict:
    """Rebuild the subset of a bundle dict that the screens need."""
    b = {"arms": {}, "repeat_qc": {}, "arm_compare": {}, "source": str(run_dir)}
    for arm in ARMS:
        a = {}
        cal = run_dir / arm / "calibration_summary.tsv"
        if cal.exists():
            a["calibration"] = pd.read_csv(cal, sep="\t")
        th = list((run_dir / arm).rglob("top_hits.tsv")) if (run_dir / arm).is_dir() else []
        if th:
            a["top_hits"] = pd.read_csv(th[0], sep="\t", low_memory=False)
            a["phenotype"] = th[0].parent.name
        b["arms"][arm] = a
    return b


def sweep(table: pd.DataFrame, id_col: str, tar_col="run_summary_tar",
          store="harvest_cache", num_ancs=5, workers=8, progress=True):
    """Harvest and screen every run with a bundle.

    Returns (per_run_calibration, all_hits, failures).
    """
    from concurrent.futures import ThreadPoolExecutor, as_completed
    rows = [(r[id_col], r[tar_col]) for _, r in table.iterrows()
            if isinstance(r.get(tar_col), str) and r[tar_col]]
    failures, dirs = [], {}
    with ThreadPoolExecutor(max_workers=workers) as ex:
        futs = {ex.submit(harvest, uri, rid, store): rid for rid, uri in rows}
        for i, f in enumerate(as_completed(futs), 1):
            rid = futs[f]
            try:
                dirs[rid] = f.result()
            except Exception as exc:
                failures.append(dict(run=rid, error=str(exc)[:160]))
            if progress and i % 10 == 0:
                print(f"  harvested {i}/{len(rows)}", flush=True)

    cal_rows, hit_rows = [], []
    for rid, d in sorted(dirs.items()):
        b = read_harvest(d)
        c = calibration_report(b)
        c.insert(0, "run", rid)
        cal_rows.append(c)
        h = hit_report(b, num_ancs=num_ancs)
        if len(h):
            h.insert(0, "run", rid)
            pheno = next((b["arms"][a].get("phenotype") for a in ARMS
                          if b["arms"].get(a, {}).get("phenotype")), rid)
            h.insert(1, "phenotype", pheno)
            hit_rows.append(h)
    cal = pd.concat(cal_rows, ignore_index=True) if cal_rows else pd.DataFrame()
    hits = pd.concat(hit_rows, ignore_index=True) if hit_rows else pd.DataFrame()
    return cal, hits, pd.DataFrame(failures)


def cluster_positions(hits: pd.DataFrame, window=1000) -> pd.Series:
    """Label hits that sit within `window` bp of each other on a chromosome.

    Needed because MarkerID is not comparable across arms — the SNV arm writes
    `chr9_8991058_T_TA;...` where the repeat arms write `9-8991058-8991069-A`
    for the same site — so cross-arm agreement has to be matched on position.
    """
    if not len(hits):
        return pd.Series(dtype=object)
    d = hits[["CHR", "POS"]].copy()
    d["POS"] = pd.to_numeric(d["POS"], errors="coerce")
    d["_i"] = np.arange(len(d))
    lab = pd.Series(index=hits.index, dtype=object)
    for chrom, g in d.groupby("CHR", sort=False):
        g = g.sort_values("POS")
        new = (g.POS.diff() > window).fillna(True).cumsum()
        for idx, c in zip(g.index, new):
            lab.loc[idx] = f"{chrom}:{int(g.loc[g.index[new.eq(c).values.argmax()], 'POS'])}"
    return lab


def recurrence(hits: pd.DataFrame, min_phenotypes=2, window=1000) -> pd.DataFrame:
    """Loci that are hits in several phenotypes, matched by position.

    A locus surfacing across unrelated phenotypes is more often a property of
    the locus (mapping, genotyping, a degenerate dosage) than real pleiotropy,
    so this is a screen in its own right.
    """
    if not len(hits):
        return pd.DataFrame()
    h = hits.copy()
    h["locus"] = cluster_positions(h, window)
    g = (h.groupby("locus")
          .agg(n_phenotypes=("phenotype", "nunique"), n_hits=("p", "size"),
               arms=("arm", lambda s: ",".join(sorted(set(s)))),
               kept=("verdict", lambda s: int((s == "keep").sum())),
               excluded=("verdict", lambda s: int((s == "EXCLUDE").sum())),
               min_p=("p", "min"),
               phenotypes=("phenotype", lambda s: ",".join(sorted(set(s)))))
          .reset_index())
    return (g[g.n_phenotypes >= min_phenotypes]
            .sort_values(["n_phenotypes", "n_hits"], ascending=False)
            .reset_index(drop=True))


def corroborated(hits: pd.DataFrame, window=1000) -> pd.DataFrame:
    """Retained hits seen in more than one arm of the same phenotype.

    Matched by position, so a repeat locus and the SNV/indel records over the
    same tandem array are recognised as the same finding.
    """
    if not len(hits):
        return pd.DataFrame()
    k = hits[hits.verdict == "keep"].copy()
    if not len(k):
        return pd.DataFrame()
    k["locus"] = cluster_positions(k, window)
    g = (k.groupby(["phenotype", "locus"])
          .agg(n_arms=("arm", "nunique"), arms=("arm", lambda s: ",".join(sorted(set(s)))),
               n_markers=("MarkerID", "nunique"), min_p=("p", "min"),
               case_carriers=("case_carriers", "max"))
          .reset_index())
    return g[g.n_arms >= 2].sort_values("min_p").reset_index(drop=True)
