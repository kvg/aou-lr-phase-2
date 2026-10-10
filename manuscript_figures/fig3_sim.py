"""Shared inputs and simulation for the Figure 3 mock-ups (all-by-all SV scan).

Decisions (2026-10-08): PheCodeX phenotypes only; case threshold 50 (pooled);
sex-specific phecodes tested only in the matching sex, and only participants
with a defined sex_at_birth are analysed; significance line 5e-8 per test with
no correction for the number of phenotypes (the All by All convention).

Ancestry axis (ANC_MODE)
------------------------
"group" (default): the seven participant genetic-ancestry groups of Figs. 1-2. A matrix cell is the
number of phecodes with at least MIN_CASES cases within the group; a card effect is drawn only for a
group with at least MIN_CASES cases of that phecode (group-stratified estimates, which need within-group
analyses). "local": the K = 5 FLARE local ancestries that FELIX's per-ancestry tests report (MID and
OTH have no column); cells and cards use expected case haplotypes with a MIN_HAP floor.

Observed inputs
---------------
* PheCodeX chapters and phecode counts: data/phecodeX_info.csv
* data/fig3_cohort_summary.json: analysis cohort (long-read, releasable, in the
  phenotype table, defined sex_at_birth, complete age/PC1-10/coverage; n = 9,375)
  and its expected haplotypes per local ancestry.
* data/fig3_phecode_case_counts.tsv, one row per phecode: pooled cases and
  eligible participants in that cohort, plus expected case haplotypes per local
  ancestry (hap_AFR..hap_SAS). Built from
  tractor_mix/resources/AoU_Phase2_Phenotype.csv.gz and the chr1 FLARE global
  ancestry proportions (aou_lr_phase2_v1.chr1.global.anc.gz). A haplotype count
  is 2 x the sum of a participant's ancestry proportion over cases, so it is an
  expectation from chr1 global ancestry, not a count of genome-wide local calls.
  K = 5 local ancestries (FLARE panels); K = 6 awaits the nanc 5 vs 6 decision.
  If either file is absent, case counts are simulated (Coverage.simulated).

Simulated (placeholder) until the Terra scan is run
---------------------------------------------------
* All association statistics (hit map, QQ, scatters, card effects and P values).
  Card standard errors are scaled from the observed case haplotypes, and an
  ancestry is drawn as not estimable below MIN_HAP; effect sizes are simulated.
  Locus choices on the cards are illustrative literature priors, not results.
"""

from __future__ import annotations

import csv
import json
import math
from dataclasses import dataclass

import numpy as np

from .style import HG38_MB, ROOT

_erfc = np.vectorize(math.erfc, otypes=[float])  # numpy-only: the repo .venv has no scipy
CHI2_MEDIAN_1DF = 0.454936423119572

MIN_CASES = 50  # decided: minimum pooled cases for a phecode to be tested
LOCAL_ANC = ["AFR", "AMR", "EAS", "EUR", "SAS"]  # K = 5 FLARE reference panels (display order)
MIN_HAP = 2 * MIN_CASES  # case haplotypes of one local ancestry = the same 50-case equivalent
GROUP_ANC = ["AFR", "AMR", "EAS", "EUR", "MID", "SAS", "OTH"]  # participant groups, as in Figs. 1-2
ANC_MODE = "group"  # "group" (7 participant groups) or "local" (K=5 FLARE local ancestries)
GW_LOGP = -np.log10(5e-8)  # decided: 5e-8 per test, no correction for phenotype count (7.30)
HET_LOGP = 3.0  # placeholder: P_het <= 1e-3 counts as ancestry-specific
SEED = 20261008

CHAPTER_FIX = {"Muscloskeletal": "Musculoskeletal"}  # typo in phecodeX_info.csv
# Log10 prevalence centre per chapter for the simulated case counts.
CH_MEAN = {
    "Symptoms": -2.00, "Cardiovascular": -2.05, "Endocrine/Metab": -2.10,
    "Musculoskeletal": -2.15, "Gastrointestinal": -2.25, "Genitourinary": -2.25,
    "Respiratory": -2.25, "Sense organs": -2.25, "Mental": -2.25,
    "Neurological": -2.55, "Dermatological": -2.45, "Infections": -2.45,
    "Neoplasms": -2.75, "Blood/Immune": -2.55, "Pregnancy": -2.45,
    "Congenital": -3.40, "Genetic": -3.70, "Neonatal": -4.00,
}
CH_SIGMA = 0.62
CH_SHIFT = -0.22  # global shift so ~1 in 4 phecodes reaches MIN_CASES pooled

_FALLBACK_SUMMARY = {  # used only when data/fig3_cohort_summary.json is absent
    "n": 9375, "n_female": 5261, "n_male": 4114,
    "n_group": {"AFR": 2106, "AMR": 1272, "EAS": 1029, "EUR": 1576, "MID": 321, "SAS": 1094, "OTH": 1977},
    "hap_total": {"AFR": 4786.0, "AMR": 988.0, "EAS": 2451.0, "EUR": 7708.0, "SAS": 2817.0},
}


def load_phecodes() -> list[dict]:
    path = ROOT / "data" / "phecodeX_info.csv"
    with open(path, encoding="latin-1", newline="") as fh:
        rows = list(csv.DictReader(fh))
    for r in rows:
        r["category"] = CHAPTER_FIX.get(r["category"], r["category"])
    return rows


def load_cohort_summary() -> dict:
    path = ROOT / "data" / "fig3_cohort_summary.json"
    if path.is_file():
        with open(path) as fh:
            return json.load(fh)
    return dict(_FALLBACK_SUMMARY)


def load_case_counts():
    """phecode -> {n_cases, n_eligible, cases: {group: n}, hap: {local ancestry: expected case haplotypes}}."""
    path = ROOT / "data" / "fig3_phecode_case_counts.tsv"
    if not path.is_file():
        return None
    out = {}
    with open(path, newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            out[r["phecode"]] = {
                "n_cases": int(r["n_cases"]),
                "n_eligible": int(r["n_eligible"]),
                "cases": {g: int(r[f"cases_{g}"]) for g in GROUP_ANC},
                "hap": {k: float(r[f"hap_{k}"]) for k in LOCAL_ANC},
            }
    return out


def ancestry_keys(mode: str | None = None) -> list[str]:
    return GROUP_ANC if (mode or ANC_MODE) == "group" else LOCAL_ANC


@dataclass
class Coverage:
    chapters: list[str]  # display order (most testable phecodes first)
    n_in_chapter: dict[str, int]
    testable: dict[str, dict[str, int]]  # chapter -> {ancestry key | ALL: n testable phecodes}
    mode: str  # "group" or "local"
    keys: list[str]  # ancestry columns in display order
    totals: dict[str, float]  # participants per group, or expected haplotypes per local ancestry
    n_total: int  # participants in the analysis cohort
    simulated: bool = True


def build_coverage(order: str = "testable", mode: str | None = None) -> Coverage:
    mode = mode or ANC_MODE
    rng = np.random.default_rng(SEED)
    phe = load_phecodes()
    summ = load_cohort_summary()
    keys = ancestry_keys(mode)
    chapters = sorted({r["category"] for r in phe})
    n_in = {c: sum(r["category"] == c for r in phe) for c in chapters}
    testable = {c: {k: 0 for k in keys + ["ALL"]} for c in chapters}
    real = load_case_counts()
    if real is not None:
        for r in phe:
            row = real.get(r["phecode"])
            if row is None:  # not in the phenotype table
                continue
            c = r["category"]
            testable[c]["ALL"] += int(row["n_cases"] >= MIN_CASES)
            for k in keys:
                if mode == "group":
                    testable[c][k] += int(row["cases"][k] >= MIN_CASES)
                else:
                    testable[c][k] += int(row["hap"][k] >= MIN_HAP)
    else:
        n_el = {"Both": summ["n"], "Female": summ["n_female"], "Male": summ["n_male"]}
        for r in phe:
            c = r["category"]
            sex = r["sex"] if r["sex"] in n_el else "Both"
            lp = rng.normal(CH_MEAN[c] + CH_SHIFT, CH_SIGMA)
            k_cases = rng.binomial(n_el[sex], min(10 ** lp, 0.5))
            testable[c]["ALL"] += int(k_cases >= MIN_CASES)
            for k in keys:
                if mode == "group":
                    testable[c][k] += int(k_cases * summ["n_group"][k] / summ["n"] >= MIN_CASES)
                else:
                    testable[c][k] += int(2 * k_cases * summ["hap_total"][k] / (2 * summ["n"]) >= MIN_HAP)
    if order == "testable":
        chapters = sorted(chapters, key=lambda c: (-testable[c]["ALL"], c))
    else:  # canonical PheCodeX order
        num = {r["category"]: r["category_num"] for r in phe}
        chapters = sorted(chapters, key=lambda c: num[c])
    totals = dict(summ["n_group"]) if mode == "group" else dict(summ["hap_total"])
    return Coverage(chapters, n_in, testable, mode, keys, totals, summ["n"], simulated=real is None)


@dataclass
class Locus:
    num: int
    gene: str
    chrom: int
    pos_mb: float
    chapter: str
    trait: str
    phecode: str  # PheCodeX code whose observed case haplotypes set each ancestry's power
    kind: str  # "repeat" | "del" | "dup"
    lead: str  # "SV" | "SNV"
    logp: float  # simulated headline -log10 P of the shared-effect signal
    unit: str  # effect unit for the lollipop strip
    beta: dict[str, float]  # simulated effect per local ancestry


NAN = float("nan")
LOCI = [
    Locus(1, "LPA KIV-2", 6, 160.6, "Cardiovascular", "Coronary atherosclerosis", "CV_404.2", "repeat", "SV", 14.8,
          "log-OR per repeat", {"AFR": -0.020, "AMR": -0.032, "EAS": -0.044, "EUR": -0.041, "SAS": -0.036, "MID": -0.030, "OTH": -0.028}),
    Locus(2, "APOL1 G1/G2 SV", 22, 36.25, "Genitourinary", "Chronic kidney disease", "GU_582.2", "del", "SV", 11.6,
          "log-OR per allele", {"AFR": 0.62, "AMR": 0.30, "EAS": 0.05, "EUR": 0.02, "SAS": 0.05, "MID": 0.10, "OTH": 0.41}),
    Locus(3, "HBA1/2 deletion", 16, 0.18, "Blood/Immune", "Anemia", "BI_164", "del", "SV", 13.1,
          "log-OR per allele", {"AFR": 0.46, "AMR": 0.28, "EAS": 0.71, "EUR": 0.10, "SAS": 0.64, "MID": 0.52, "OTH": 0.40}),
    Locus(4, "TCF4 CTG18.1", 18, 55.59, "Sense organs", "Disorders of the cornea", "SO_369", "repeat", "SV", 9.9,
          "log-OR per repeat", {"AFR": 0.004, "AMR": 0.012, "EAS": 0.006, "EUR": 0.031, "SAS": 0.010, "MID": 0.010, "OTH": 0.018}),
    Locus(5, "C4A/C4B", 6, 31.95, "Musculoskeletal", "Systemic lupus erythematosus", "MS_700.11", "dup", "SV", 9.4,
          "log-OR per copy", {"AFR": 0.52, "AMR": 0.35, "EAS": 0.30, "EUR": 0.45, "SAS": 0.35, "MID": 0.40, "OTH": 0.44}),
    Locus(6, "CFHR3/1 deletion", 1, 196.80, "Sense organs", "Macular degeneration", "SO_374.5", "del", "SNV", 8.1,
          "log-OR per allele", {"AFR": -0.10, "AMR": -0.22, "EAS": -0.20, "EUR": -0.41, "SAS": -0.25, "MID": -0.20, "OTH": -0.30}),
    Locus(7, "LCE3B/C deletion", 1, 152.57, "Dermatological", "Psoriasis", "DE_664.4", "del", "SV", 8.9,
          "log-OR per allele", {"AFR": 0.08, "AMR": 0.21, "EAS": 0.33, "EUR": 0.29, "SAS": 0.25, "MID": 0.20, "OTH": 0.22}),
]


def _z_from_logp(logp: float) -> float:
    """Two-sided z whose P equals 10**-logp (bisection on erfc)."""
    target = 10.0 ** (-logp)
    lo, hi = 0.0, 40.0
    for _ in range(80):
        mid = (lo + hi) / 2
        lo, hi = (mid, hi) if math.erfc(mid / math.sqrt(2)) > target else (lo, mid)
    return (lo + hi) / 2


def chi2_sf_small_df(x: float, df: int) -> float:
    """Chi-square survival function for integer df >= 1 (closed forms; math only)."""
    h = x / 2
    if df % 2 == 0:
        return math.exp(-h) * sum(h ** j / math.factorial(j) for j in range(df // 2))
    return math.erfc(math.sqrt(h)) + math.exp(-h) * sum(h ** (j - 0.5) / math.gamma(j + 0.5) for j in range(1, (df - 1) // 2 + 1))


def card_effects(locus: Locus, mode: str | None = None):
    """Per-ancestry (beta, se) for a locus card, plus the heterogeneity P among estimable ancestries.

    mode "group": size = cases of the phecode in each participant group, floor MIN_CASES.
    mode "local": size = expected case haplotypes of each local ancestry, floor MIN_HAP.
    The shared ('ALL') effect is the size-weighted mean of the simulated per-ancestry effects, with its SE
    set so its z matches the card's headline P; each ancestry's SE then scales with the square root of its
    observed size (se_k = se_all * sqrt(size_all / size_k)). Below the floor an ancestry is not estimable.
    Returns (effects, n_estimable, p_het or None).
    """
    mode = mode or ANC_MODE
    keys = ancestry_keys(mode)
    floor = MIN_CASES if mode == "group" else MIN_HAP
    real = load_case_counts()
    row = real.get(locus.phecode) if real else None
    size = ({k: float(row["cases"][k]) for k in keys} if mode == "group" else dict(row["hap"])) if row else {k: 1000.0 for k in keys}
    size_all = sum(size.values())
    b_all = sum(size[k] * locus.beta[k] for k in keys) / size_all
    se_all = abs(b_all) / _z_from_logp(locus.logp)
    out = {"ALL": (b_all, se_all)}
    est = []
    for k in keys:
        if size[k] >= floor:
            se = se_all * math.sqrt(size_all / size[k])
            out[k] = (locus.beta[k], se)
            est.append((locus.beta[k], se))
        else:
            out[k] = (NAN, NAN)
    p_het = None
    if len(est) >= 2:
        w = [1 / se ** 2 for _, se in est]
        bbar = sum(wi * b for wi, (b, _) in zip(w, est)) / sum(w)
        q = sum(wi * (b - bbar) ** 2 for wi, (b, _) in zip(w, est))
        p_het = chi2_sf_small_df(q, len(est) - 1)
    return out, len(est), p_het


def regional(locus: Locus, rng: np.random.Generator, n: int = 170):
    """Simulated 1 Mb regional plot: SNV positions (Mb offset), -log10 P, SV logp."""
    x = rng.uniform(-0.5, 0.5, n)
    h = locus.logp if locus.lead == "SV" else locus.logp + 1.6
    ld = np.exp(-np.abs(x) / rng.uniform(0.07, 0.13))
    r2 = np.clip(ld * rng.uniform(0.35, 1.0, n), 0, 1)
    y = np.clip(h * (r2 ** 1.25) * rng.uniform(0.55, 1.0, n), 0, None) + rng.exponential(0.5, n)
    if locus.lead == "SV":  # the SV is the top signal; the best SNVs trail it
        y = np.minimum(y, locus.logp - rng.uniform(0.3, 1.8, n))
    sv_logp = locus.logp if locus.lead == "SV" else locus.logp - 1.4
    return x, np.clip(y, 0.05, None), sv_logp


# --------------------------------------------------------------------------
# Hit map: one mark per locus x chapter, genome-wide significant.
# --------------------------------------------------------------------------
@dataclass
class Hit:
    chrom: int
    pos_mb: float
    chapter: str
    logp: float
    repeat: bool
    num: int = 0


def build_hits(cov: Coverage) -> list[Hit]:
    rng = np.random.default_rng(SEED + 1)
    hits = [Hit(l.chrom, l.pos_mb, l.chapter, l.logp, l.kind == "repeat", l.num) for l in LOCI]
    w = np.array([cov.testable[c]["ALL"] for c in cov.chapters], float)
    w = np.where(w < 10, 0.0, w) ** 1.2  # no hits where almost nothing is testable
    w /= w.sum()
    chrom_w = np.array([HG38_MB[c] for c in HG38_MB])
    chrom_w = chrom_w / chrom_w.sum()
    for _ in range(78):
        ch = int(rng.choice(list(HG38_MB), p=chrom_w))
        pos = float(rng.uniform(1, HG38_MB[ch] - 1))
        k = int(rng.choice([1, 2, 3], p=[0.62, 0.28, 0.10]))
        chs = rng.choice(len(cov.chapters), size=k, replace=False, p=w)
        rep = bool(rng.random() < 0.28)
        for ci in chs:
            hits.append(Hit(ch, pos, cov.chapters[ci], GW_LOGP + rng.exponential(1.7), rep))
    # an MHC-region cluster across immune-mediated chapters
    for c in ("Blood/Immune", "Musculoskeletal", "Endocrine/Metab", "Dermatological"):
        if c in cov.chapters:
            hits.append(Hit(6, 31.8 + rng.uniform(-1.5, 1.5), c, GW_LOGP + rng.exponential(3.0), False))
    return hits


# --------------------------------------------------------------------------
# Calibration and method-comparison scatters
# --------------------------------------------------------------------------
def chi2_sf_1df(x):
    return _erfc(np.sqrt(np.asarray(x, float) / 2))


def chi2_sf_5df(x):
    x = np.asarray(x, float)
    return _erfc(np.sqrt(x / 2)) + np.exp(-x / 2) * np.sqrt(2 * x / np.pi) * (1 + x / 3)


def neglog10_two_sided(z):
    return -np.log10(np.clip(_erfc(np.abs(np.asarray(z, float)) / np.sqrt(2)), 1e-300, 1))


def sim_qq(n: int = 60000):
    rng = np.random.default_rng(SEED + 2)
    arms = [("Pooled SAIGE", 1.005), ("FELIX, biallelic", 1.012), ("FELIX, repeat dosage", 1.018)]
    out = []
    for name, lam in arms:
        z2 = rng.chisquare(1, n) * lam
        p = chi2_sf_1df(z2)
        k = 320  # true signals, from weak to genome-wide significant
        p[:k] = 10 ** (-(2.5 + rng.exponential(3.4, k)))
        lam_gc = float(np.median(z2[k:]) / CHI2_MEDIAN_1DF)  # on the null part
        out.append((name, np.sort(p), lam_gc))
    return out


def cct(p1: np.ndarray, p2: np.ndarray) -> np.ndarray:
    """Cauchy combination of two p-values (equal weights); stable for tiny p."""
    t = 0.5 * (np.tan((0.5 - p1) * np.pi) + np.tan((0.5 - p2) * np.pi))
    big = t > 1e15
    return np.where(big, 1.0 / (np.pi * np.where(big, t, 1.0)), 0.5 - np.arctan(t) / np.pi)


def sim_dosage_vs_biallelic(n: int = 720):
    """-log10 P at repeat loci: biallelic-split coding vs repeat-length dosage."""
    rng = np.random.default_rng(SEED + 3)
    z0 = rng.normal(size=n)
    zs = 0.55 * z0 + np.sqrt(1 - 0.55 ** 2) * rng.normal(size=n)
    sig = np.zeros(n)
    idx = rng.choice(n, 46, replace=False)
    sig[idx] = rng.uniform(3.8, 8.8, idx.size) * rng.choice([1, -1], idx.size)
    z_dos = z0 + sig
    shrink = np.where(rng.random(n) < 0.12, rng.uniform(1.0, 1.25, n), rng.uniform(0.35, 0.8, n))
    z_split = zs + sig * shrink
    return neglog10_two_sided(z_split), neglog10_two_sided(z_dos)


def sim_ancestry_vs_pooled(n: int = 4200):
    """-log10 P: pooled test vs FELIX Cauchy combination of shared and heterogeneity tests."""
    rng = np.random.default_rng(SEED + 4)
    z = rng.normal(size=n)
    sig = np.zeros(n)
    idx = rng.choice(n, 150, replace=False)
    sig[idx] = rng.uniform(2.5, 7.5, idx.size)
    z_pool = z + sig
    het = np.zeros(n, bool)
    het_idx = rng.choice(idx, 52, replace=False)  # ancestry-specific: pooled test loses power
    het[het_idx] = True
    z_pool[het_idx] = z[het_idx] + sig[het_idx] * rng.uniform(0.25, 0.6, het_idx.size)
    nc = np.where(het, (sig * rng.uniform(0.8, 1.3, n)) ** 2, 0.0)
    chi_het = rng.noncentral_chisquare(5, np.maximum(nc, 1e-12))
    p_het = np.clip(chi2_sf_5df(chi_het), 1e-300, 1)
    p_pool = np.clip(_erfc(np.abs(z_pool) / np.sqrt(2)), 1e-300, 1)
    p_cct = np.clip(cct(p_pool, p_het), 1e-300, 1)
    return -np.log10(p_pool), -np.log10(p_cct), -np.log10(p_het)
