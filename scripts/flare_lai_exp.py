#!/usr/bin/env python3
"""Fetch and summarize the Terra ``flare_lai_exp`` method-grid table."""

from __future__ import annotations

import argparse
import json
from io import StringIO
from pathlib import Path
from typing import Any, Optional

DEFAULT_NAMESPACE = "allofus-drc-wgs-LR-prodData"
DEFAULT_WORKSPACE = "AoU_DRC_LongReads_PhaseTwo_Storage"
DEFAULT_ENTITY_TYPE = "flare_lai_exp"
DEFAULT_ID_COLUMN = "flare_lai_exp_id"


def _entities_to_rows(entities: list[dict[str, Any]], *, id_column: str) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for ent in entities:
        rid = str(ent.get("name") or "")
        attrs = ent.get("attributes") or {}
        flat = {str(k): "" if v is None else str(v) for k, v in attrs.items()}
        # Terra Array[File] attributes often arrive as JSON-ish lists; keep as text.
        flat[id_column] = rid
        rows.append(flat)
    return rows


def _read_entity_tsv(text: str, *, id_column: str) -> "Any":
    import pandas as pd

    df = pd.read_csv(StringIO(text), sep="\t", dtype=str)
    if df.empty:
        raise RuntimeError("Terra entity TSV was empty")
    if id_column not in df.columns:
        candidates = [
            c for c in df.columns if str(c).startswith("entity:") or str(c).endswith("_id")
        ]
        if not candidates:
            raise ValueError(f"No entity id column in Terra TSV; columns={list(df.columns)}")
        df = df.rename(columns={candidates[0]: id_column})
    return df.fillna("")


def _fetch_entities_tsv(fapi: Any, namespace: str, workspace: str, entity_type: str) -> Any:
    try:
        return fapi.get_entities_tsv(namespace, workspace, entity_type, model="flexible")
    except TypeError:
        return fapi.get_entities_tsv(namespace, workspace, entity_type)


def fetch_lai_exp_table(
    *,
    tsv: Path | str | None = None,
    from_firecloud: bool = False,
    namespace: str = DEFAULT_NAMESPACE,
    workspace: str = DEFAULT_WORKSPACE,
    entity_type: str = DEFAULT_ENTITY_TYPE,
    id_column: str = DEFAULT_ID_COLUMN,
) -> "Any":
    """Load ``flare_lai_exp`` from a TSV export or the Firecloud/Terra API (FISS)."""
    import pandas as pd

    if tsv:
        df = pd.read_csv(tsv, sep="\t", dtype=str).fillna("")
        if id_column not in df.columns:
            candidates = [c for c in df.columns if c.startswith("entity:") or c.endswith("_id")]
            if not candidates:
                raise ValueError(f"No entity id column in {tsv}; columns={list(df.columns)}")
            df = df.rename(columns={candidates[0]: id_column})
        return df
    if not from_firecloud:
        raise ValueError("fetch_lai_exp_table requires tsv= or from_firecloud=True")

    from firecloud import api as fapi

    where = f"{namespace}/{workspace}/{entity_type}"
    tsv_error = ""
    if hasattr(fapi, "get_entities_tsv"):
        resp = _fetch_entities_tsv(fapi, namespace, workspace, entity_type)
        if resp.status_code == 200 and str(getattr(resp, "text", "")).strip():
            return _read_entity_tsv(resp.text, id_column=id_column)
        tsv_error = (
            f"get_entities_tsv {resp.status_code}: {str(getattr(resp, 'text', ''))[:300]}"
        )

    resp = fapi.get_entities(namespace, workspace, entity_type)
    if resp.status_code != 200:
        detail = tsv_error + ("; " if tsv_error else "") + resp.text[:500]
        raise RuntimeError(
            f"firecloud get_entities failed ({resp.status_code}) for {where}: {detail}"
        )
    entities = resp.json()
    if not isinstance(entities, list):
        raise RuntimeError(f"Unexpected firecloud payload type: {type(entities)}")
    rows = _entities_to_rows(entities, id_column=id_column)
    if not rows:
        raise RuntimeError(f"No entities returned for {where}")
    return pd.DataFrame(rows, dtype=str).fillna("")


def _nonempty_gs(val: Any) -> bool:
    s = "" if val is None else str(val).strip()
    return bool(s) and s.lower() != "nan" and s.startswith("gs://")


def row_is_complete(row: "Any") -> bool:
    """True when MergePopulationFlare outputs are written back to the table."""
    return _nonempty_gs(row.get("anc_vcf")) and _nonempty_gs(row.get("models_tsv"))


def status_frame(df: "Any", *, id_column: str = DEFAULT_ID_COLUMN) -> "Any":
    import pandas as pd

    rows = []
    for _, row in df.iterrows():
        done = row_is_complete(row)
        rows.append(
            {
                id_column: row.get(id_column, ""),
                "notes": row.get("notes", ""),
                "em": row.get("em", ""),
                "gen": row.get("gen", ""),
                "gen_by_pop": row.get("gen_by_pop", ""),
                "include_pops": row.get("include_pops", ""),
                "min_mac": row.get("min_mac", ""),
                "min_maf": row.get("min_maf", ""),
                "region": row.get("region", ""),
                "complete": done,
                "anc_vcf": row.get("anc_vcf", "") if done else "",
                "models_tsv": row.get("models_tsv", "") if done else "",
            }
        )
    out = pd.DataFrame(rows)
    return out.sort_values(["complete", id_column], ascending=[False, True]).reset_index(drop=True)


def parse_gen_by_pop(raw: str) -> dict[str, float]:
    out: dict[str, float] = {}
    text = (raw or "").strip()
    if not text:
        return out
    for part in text.split(","):
        part = part.strip()
        if not part:
            continue
        if ":" not in part:
            raise ValueError(f"gen_by_pop entry {part!r} must look like POP:T")
        key, val = part.split(":", 1)
        out[key.strip()] = float(val.strip())
    return out


def pinned_t_for_pop(row: "Any", population: str) -> Optional[float]:
    overrides = parse_gen_by_pop(str(row.get("gen_by_pop") or ""))
    if population in overrides:
        return overrides[population]
    gen = str(row.get("gen") or "").strip()
    if not gen:
        return None
    try:
        return float(gen)
    except ValueError:
        return None


def load_models_tsv(path: Path) -> "Any":
    import pandas as pd

    return pd.read_csv(path, sep="\t")


def summarize_models(models: "Any", row: "Any") -> list[dict[str, Any]]:
    from flare_switch_qc import expected_switches_per_hap_per_mb, het_admixture_factor

    out = []
    if models is None or models.empty:
        return out
    prop_cols = [c for c in models.columns if str(c).startswith("prop_")]
    for _, m in models.iterrows():
        pop = str(m.get("population", ""))
        t_gen = float(m["t_gen"]) if "t_gen" in m and pd_notna(m["t_gen"]) else None
        pinned = pinned_t_for_pop(row, pop)
        em = str(row.get("em", "")).strip().lower() in {"true", "1", "yes"}
        props = None
        if prop_cols:
            props = [
                float(m[c]) for c in prop_cols if pd_notna(m[c])
            ]
            if len(props) != len(prop_cols):
                props = None
        het = het_admixture_factor(props) if props else None
        out.append(
            {
                "population": pop,
                "model_t_gen": t_gen,
                "pinned_gen": pinned,
                "em": em,
                "props": props,
                "het_factor": het,
                "expected_rate_model_t": expected_switches_per_hap_per_mb(t_gen, props)
                if t_gen
                else None,
                "expected_rate_pinned": expected_switches_per_hap_per_mb(pinned, props)
                if pinned
                else None,
                "prop_afr": float(m["prop_afr"]) if "prop_afr" in m and pd_notna(m["prop_afr"]) else None,
                "prop_eur": float(m["prop_eur"]) if "prop_eur" in m and pd_notna(m["prop_eur"]) else None,
                "prop_amr": float(m["prop_amr"]) if "prop_amr" in m and pd_notna(m["prop_amr"]) else None,
                "mu_afr": float(m["mu_afr"]) if "mu_afr" in m and pd_notna(m["mu_afr"]) else None,
                "mu_eur": float(m["mu_eur"]) if "mu_eur" in m and pd_notna(m["mu_eur"]) else None,
                "mu_amr": float(m["mu_amr"]) if "mu_amr" in m and pd_notna(m["mu_amr"]) else None,
            }
        )
    return out


def pd_notna(val: Any) -> bool:
    try:
        import pandas as pd

        return bool(pd.notna(val))
    except Exception:
        return val is not None and str(val).strip() not in {"", "nan", "None"}


def tract_metrics_from_summary(
    summary: dict[str, Any],
    *,
    n_markers: Optional[float] = None,
) -> dict[str, Any]:
    tr = summary.get("tracts") or {}
    span_mb = tr.get("span_mb")
    n_switches = summary.get("n_switches")
    n_haps = tr.get("n_haps")
    switches_per_hap_per_marker = None
    markers_per_mb = None
    if n_markers and span_mb and span_mb > 0:
        markers_per_mb = float(n_markers) / float(span_mb)
    if (
        n_markers
        and n_switches is not None
        and n_haps
        and float(n_markers) > 0
        and float(n_haps) > 0
    ):
        switches_per_hap_per_marker = float(n_switches) / float(n_haps) / float(n_markers)
    out = {
        "n_samples": summary.get("n_samples"),
        "n_switches": n_switches,
        "n_flicker": summary.get("n_flicker_switches"),
        "span_mb": span_mb,
        "n_markers": n_markers,
        "markers_per_mb": markers_per_mb,
        "switches_per_hap_per_mb": tr.get("switches_per_hap_per_mb"),
        "switches_per_hap_per_marker": switches_per_hap_per_marker,
        "flicker_frac": tr.get("flicker_frac"),
        "flicker_per_hap_per_mb": tr.get("flicker_per_hap_per_mb"),
        "frac_haps_with_switch": tr.get("frac_haps_with_switch"),
        "implied_t_gen": tr.get("implied_t_gen"),
        "median_length_bp": tr.get("median_length_bp"),
        "median_complete_length_bp": tr.get("median_complete_length_bp"),
        "gq_dp_available": summary.get("gq_dp_available"),
    }
    cqf = summary.get("call_qc_filtered") or {}
    if cqf.get("n_switches_all") is not None:
        out["call_qc_frac_fail"] = cqf.get("frac_switches_fail_call_qc")
        out["call_qc_switches_per_hap_per_mb"] = cqf.get("switches_per_hap_per_mb")
        out["call_qc_implied_t_gen"] = cqf.get("implied_t_gen")
        n_pass = cqf.get("n_switches_pass_call_qc")
        if (
            n_markers
            and n_pass is not None
            and n_haps
            and float(n_markers) > 0
            and float(n_haps) > 0
        ):
            out["call_qc_switches_per_hap_per_marker"] = (
                float(n_pass) / float(n_haps) / float(n_markers)
            )
    return out


def filter_provenance(row: "Any") -> tuple[str, str, str]:
    """Return (biallelic_snvs_only, include_sites basename, exclude_regions basename)."""

    def _basename(val: Any) -> str:
        s = "" if val is None else str(val).strip()
        if not s or s.lower() in {"nan", "none", "null", "[]"}:
            return ""
        # Terra may store gs://.../file.bed
        return Path(s.split(",")[0].strip().strip('"')).name

    bial = str(row.get("biallelic_snvs_only", "")).strip().lower()
    if bial in {"true", "1", "yes"}:
        bial_s = "true"
    elif bial in {"false", "0", "no", ""}:
        bial_s = "false"
    else:
        bial_s = bial
    return (bial_s, _basename(row.get("include_sites")), _basename(row.get("exclude_regions")))


def warn_mixed_filter_provenance(rows: list[dict[str, Any]], *, id_key: str = "experiment") -> list[str]:
    """Warn when a comparison set spans more than one filter-provenance group."""
    groups: dict[tuple[str, str, str], list[str]] = {}
    for r in rows:
        key = (
            str(r.get("biallelic_snvs_only", "")),
            str(r.get("include_sites_basename", "")),
            str(r.get("exclude_regions_basename", "")),
        )
        if "filter_provenance" in r and isinstance(r["filter_provenance"], (list, tuple)):
            key = tuple(str(x) for x in r["filter_provenance"])  # type: ignore[assignment]
        groups.setdefault(key, []).append(str(r.get(id_key, "")))
    warnings: list[str] = []
    if len(groups) > 1:
        detail = "; ".join(f"{k}->{v}" for k, v in groups.items())
        warnings.append(
            "REFUSING apples-to-apples claim across different site filters: " + detail
        )
    return warnings


def enrich_compare_row(
    row: "Any",
    summary: dict[str, Any],
    *,
    n_markers: Optional[float] = None,
) -> dict[str, Any]:
    """Merge tract metrics + filter provenance for Part 5 scoring."""
    prov = filter_provenance(row)
    metrics = tract_metrics_from_summary(summary, n_markers=n_markers)
    metrics.update(
        {
            "biallelic_snvs_only": prov[0],
            "include_sites_basename": prov[1],
            "exclude_regions_basename": prov[2],
            "filter_provenance": prov,
        }
    )
    return metrics


def enrich_association_scores(
    row: "Any",
    *,
    allele_json: dict[str, Any] | Path | str | None = None,
    mendel_json: dict[str, Any] | Path | str | None = None,
    lambda_json: dict[str, Any] | Path | str | None = None,
) -> dict[str, Any]:
    """Attach association-facing scores to a compare row.

    Selection v2 (``select_recipe``): allele concordance is a gate, Mendelian
    violation rate ranks. Null-λ fields are carried as diagnostics only.
    Switch/tract metrics remain diagnostics.
    """

    def _load(obj: dict[str, Any] | Path | str | None) -> dict[str, Any]:
        if obj is None:
            return {}
        if isinstance(obj, dict):
            return obj
        return json.loads(Path(obj).read_text())

    allele = _load(allele_json)
    mendel = _load(mendel_json)
    lam = _load(lambda_json)
    out: dict[str, Any] = {
        "experiment": str(
            row.get("flare_lai_exp_id") or row.get("experiment") or allele.get("experiment") or ""
        ),
    }
    if allele:
        out.update(
            {
                "mean_ll": allele.get("mean_ll"),
                "mean_brier": allele.get("mean_brier"),
                "n_scored_haps": allele.get("n_scored_haps"),
                "n_markers_used": allele.get("n_markers_used"),
                "frac_carried_forward": allele.get("frac_carried_forward"),
            }
        )
    if mendel:
        out.update(
            {
                "violations_per_informative_locus": mendel.get(
                    "violations_per_informative_locus"
                ),
                "n_trios_scored": mendel.get("n_trios_scored"),
                "excess_recomb_over_expected": mendel.get("excess_recomb_over_expected"),
                "expected_crossovers_per_trio": mendel.get("expected_crossovers_per_trio"),
                # select_recipe refuses to rank recipes scored on different
                # ancestry alphabets, so these must reach it. Dropping them here
                # made that check inert: every recipe looked alphabet-less.
                "ancestry_alphabet": mendel.get("ancestry_alphabet"),
                "projected_to": mendel.get("projected_to"),
                "n_ancestries": mendel.get("n_ancestries"),
            }
        )
        ci = trio_bootstrap_violation_ci(mendel.get("trios") or [])
        if ci:
            out.update(
                {
                    "violations_ci_low": ci["ci_low"],
                    "violations_ci_high": ci["ci_high"],
                }
            )
    if lam:
        out.update(
            {
                "mean_lambda": lam.get("mean_lambda"),
                "abs_lambda_dev": lam.get("abs_lambda_dev"),
                "abs_lambda_ci_low": (lam.get("abs_lambda_ci") or {}).get("ci_low"),
                "abs_lambda_ci_high": (lam.get("abs_lambda_ci") or {}).get("ci_high"),
                "n_tested_sites": lam.get("n_tested_sites"),
                "unstable_lambda": lam.get("unstable_lambda"),
            }
        )
    return out


def trio_bootstrap_violation_ci(
    trios: list[dict[str, Any]],
    *,
    n_boot: int = 2000,
    alpha: float = 0.05,
    seed: int = 1,
) -> dict[str, float]:
    """Percentile CI for pooled violations / informative loci, resampling trios.

    Trios (not loci) are the independent unit: loci within a trio share tracts.
    """
    import random

    pairs = [
        (float(t.get("n_hard_violations") or 0), float(t.get("n_informative") or 0))
        for t in trios
    ]
    pairs = [(v, n) for v, n in pairs if n > 0]
    if len(pairs) < 2:
        return {}
    rng = random.Random(seed)
    rates = []
    for _ in range(n_boot):
        draw = [pairs[rng.randrange(len(pairs))] for _ in pairs]
        den = sum(n for _, n in draw)
        rates.append(sum(v for v, _ in draw) / den)
    rates.sort()
    lo = rates[int((alpha / 2) * (n_boot - 1))]
    hi = rates[int((1 - alpha / 2) * (n_boot - 1))]
    return {"ci_low": lo, "ci_high": hi, "n_trios": float(len(pairs))}


def _bcftools_any_record(vcf: str, region: str, expr: str) -> bool:
    """True if any record in ``region`` satisfies the bcftools filter ``expr``.

    Stops at the first hit. An empty result only counts as "none" when bcftools
    exited cleanly, so a failed read cannot pass for a clean scan."""
    import subprocess

    proc = subprocess.Popen(
        ["bcftools", "view", "-H", "-i", expr, "-r", region, vcf],
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True,
    )
    hit = bool(proc.stdout.readline())
    if hit:
        proc.terminate()
        proc.wait()
        return True
    err = proc.stderr.read()
    if proc.wait() != 0:
        raise RuntimeError(f"bcftools failed on {vcf} {region}: {err.strip()[-300:]}")
    return False


def _panel_windows(panel: "Path | str", n_windows: int, window_bp: int) -> list[str]:
    """Evenly spaced ``chrom:start-end`` windows anchored on marker-panel positions."""
    import flare_score_mendelian_lai as scorer

    pts = sorted(scorer.load_panel_positions(Path(panel)))
    if not pts:
        raise ValueError(f"{panel} has no marker positions to anchor sampling windows on")
    n = min(n_windows, len(pts))
    picks = [pts[round(k * (len(pts) - 1) / max(n - 1, 1))] for k in range(n)]
    return [f"{chrom}:{pos}-{pos + window_bp}" for chrom, pos in picks]


def backfill_mendelian_alphabet(
    mendel: dict[str, Any],
    anc_vcf: str,
    *,
    panel: "Path | str | None" = None,
    n_windows: int = 10,
    window_bp: int = 100_000,
) -> dict[str, Any]:
    """Give a Mendelian score written by an older scorer its ancestry alphabet.

    Scorers before 767872d hard-coded ancestry codes 0-4 and silently dropped
    every locus carrying any other code, so such a JSON is valid only if the
    VCF cannot contain codes outside 0-4. Returns ``mendel`` unchanged when it
    already carries an alphabet; raises when the cached rate cannot be trusted.

    * ``##ANCESTRY`` header present: read from it (header read only).
    * Header absent (the original-FLARE pin VCF has none): needs ``panel`` and
      checks the data. In ``n_windows`` windows spread over the marker panel it
      requires that no AN1/AN2 call exceeds 4 and that each of the five codes
      0-4 occurs. This is a sample, not an exhaustive scan; it is recorded as
      such in ``ancestry_alphabet_source``. Re-scoring with the current scorer
      (which refuses out-of-range codes on every locus) is the exhaustive check.
    """
    if mendel.get("ancestry_alphabet"):
        return mendel
    import flare_score_mendelian_lai as scorer

    legacy = set(scorer.ANCESTRY)
    names = scorer.parse_ancestry_header(scorer.vcf_header_text(anc_vcf))
    if names is not None:
        extra = sorted(set(names) - legacy)
        if extra:
            raise ValueError(
                f"{anc_vcf} declares ancestry codes {extra} outside 0-4; the older scorer "
                "dropped those loci, so its cached rate is not comparable. Delete its "
                "mendelian_lai_score.json and re-score."
            )
        alphabet = sorted(names)
        source = "##ANCESTRY header (read from the VCF; score from an older scorer)"
    else:
        if panel is None:
            raise ValueError(
                f"{anc_vcf} has no ##ANCESTRY header and no marker panel was given to sample "
                "its AN1/AN2 codes, so the alphabet of a score written by an older scorer "
                "cannot be confirmed. Pass panel=, or delete its mendelian_lai_score.json "
                "and re-score."
            )
        windows = _panel_windows(panel, n_windows, window_bp)
        top = max(legacy)
        for w in windows:
            if _bcftools_any_record(anc_vcf, w, f"FMT/AN1>{top} || FMT/AN2>{top}"):
                raise ValueError(
                    f"{anc_vcf} has AN1/AN2 calls above {top} in {w}; the older scorer dropped "
                    "those loci, so its cached rate is not comparable. Delete its "
                    "mendelian_lai_score.json and re-score."
                )
        absent = [
            c for c in sorted(legacy)
            if not any(_bcftools_any_record(anc_vcf, w, f"FMT/AN1=={c} || FMT/AN2=={c}")
                       for w in windows)
        ]
        if absent:
            raise ValueError(
                f"{anc_vcf}: ancestry code(s) {absent} never occur in the {len(windows)} sampled "
                "windows, so the five-panel alphabet is not confirmed. Re-score with the "
                "current scorer (it refuses out-of-range codes on every locus)."
            )
        alphabet = sorted(legacy)
        source = (f"sampled AN1/AN2 in {len(windows)} windows of {window_bp} bp over the marker "
                  "panel: codes 0-4 all present, none above 4 (no ##ANCESTRY header; "
                  "score from an older scorer)")
    return {
        **mendel,
        "ancestry_alphabet": alphabet,
        "ancestry_alphabet_source": source,
        "n_ancestries": len(alphabet),
        "projected_to": mendel.get("projected_to"),
        "alphabet_backfilled": True,
    }


def _num(val: Any) -> Optional[float]:
    try:
        x = float(val)
    except (TypeError, ValueError):
        return None
    return None if x != x else x


# Relative excess of the Mendelian violation rate on a recipe scored with ONE MORE
# ancestry label, at identical painting accuracy. Two synthetic measurements:
#   +9.1% to +11.2%  uniform label frequencies, K=4 -> 5, per-call error 0.5-10%
#                    (flare/eval/mendelian_k_bias/k_bias_summary.json)
#   +23.5% to +26.2% one pin-like VCF scored with and without folding a rare label
#                    (SAS, 5% of ancestry) into EUR (scripts/test_flare_alphabet_pipeline.py)
# The effect depends on label frequencies and on which labels get confused, not only
# on the label count, so the larger value is used. It only flags results that cannot
# be told apart from the label-count effect; it never adjusts a rate. The only exact
# resolution is to score both recipes on one alphabet (--project-labels).
LABEL_COUNT_BIAS_PER_LABEL = 0.26


def select_recipe(
    rows: list[dict[str, Any]],
    *,
    negative_controls: dict[str, str],
    gates: Optional[dict[str, Any]] = None,
    allow_mixed_alphabets: bool = False,
) -> dict[str, Any]:
    """Selection v2: pick a production LAI recipe (time-boxed rule).

    ``rows`` are enriched score rows for candidates **and** their negative
    controls, all scored on one sample list and one region.
    ``negative_controls`` maps candidate experiment → its within-population
    permuted control (``flare_make_negative_control.py``).

    1. Metric check: a metric is usable only if every candidate beats its own
       control on it (higher ``mean_ll``, lower Mendelian violation rate).
       If the Mendelian metric fails this, no winner is declared.
    2. Gates (fail closed on missing values): beat own control on both metrics;
       optional ``min_mean_ll`` / ``max_violations_per_informative_locus``.
    3. Rank survivors by violation rate. Candidates whose trio-bootstrap CI
       overlaps the best one's are reported as ``tied``; ties go to the human
       (tie-break suggestion: smaller |excess_recomb_over_expected|).
    4. Ancestry alphabets: the rate is only comparable between recipes scored
       on one alphabet, so by default mixed alphabets are refused. With
       ``allow_mixed_alphabets=True`` the recipes are ranked anyway, the
       override and each recipe's label count are recorded, and any recipe
       whose margin over a runner-up with MORE labels is within the simulated
       per-label bias is reported as tied rather than as a winner.
    """
    gates = {k: v for k, v in (gates or {}).items() if not str(k).startswith("_") and v is not None}
    by_id = {str(r.get("experiment")): r for r in rows}
    checks: list[dict[str, Any]] = []
    survivors: list[dict[str, Any]] = []
    eliminated: list[dict[str, Any]] = []
    mendel_separates = True
    ll_separates = True
    for cand, ctl in negative_controls.items():
        r = by_id.get(cand)
        c = by_id.get(ctl)
        if r is None:
            continue
        reasons: list[str] = []
        v_r = _num(r.get("violations_per_informative_locus"))
        v_c = _num(c.get("violations_per_informative_locus")) if c else None
        ll_r = _num(r.get("mean_ll"))
        ll_c = _num(c.get("mean_ll")) if c else None
        if c is None:
            reasons.append("negative_control_missing")
        if v_r is None:
            reasons.append("mendelian_missing")
        if ll_r is None:
            reasons.append("mean_ll_missing")
        sep_v = v_r is not None and v_c is not None and v_r < v_c
        sep_ll = ll_r is not None and ll_c is not None and ll_r > ll_c
        if c is not None and v_r is not None and v_c is not None and not sep_v:
            reasons.append("mendelian_not_better_than_control")
            mendel_separates = False
        if c is not None and ll_r is not None and ll_c is not None and not sep_ll:
            reasons.append("mean_ll_not_better_than_control")
            ll_separates = False
        min_ll = gates.get("min_mean_ll")
        if min_ll is not None and ll_r is not None and ll_r < float(min_ll):
            reasons.append("mean_ll_below_gate")
        max_v = gates.get("max_violations_per_informative_locus")
        if max_v is not None and v_r is not None and v_r > float(max_v):
            reasons.append("mendelian_above_gate")
        checks.append(
            {
                "experiment": cand,
                "negative_control": ctl,
                "violations": v_r,
                "violations_control": v_c,
                "mean_ll": ll_r,
                "mean_ll_control": ll_c,
            }
        )
        (eliminated if reasons else survivors).append({**r, "gate_fail": reasons})

    decision: dict[str, Any] = {
        "rule": "selection_v2",
        "metric_checks": checks,
        "mendelian_separates_controls": mendel_separates,
        "mean_ll_separates_controls": ll_separates,
        "gates": gates,
        "survivors": [r["experiment"] for r in survivors],
        "eliminated": {r["experiment"]: r["gate_fail"] for r in eliminated},
        "winner": None,
        "tied": [],
        "status": "",
    }
    # The violation rate is only comparable between recipes scored on the same
    # ancestry alphabet: at fixed painting accuracy it rises with the number of
    # labels, because label compatibility survives coarsening but not
    # refinement (flare/eval/mendelian_k_bias/). Ranking nanc 5 against nanc 6
    # on the raw rate therefore favours nanc 5 whatever the painting quality.
    # Score with flare_score_mendelian_lai.py --project-labels to put every
    # recipe in one alphabet before ranking.
    def _alphabet(r):
        for key in ("projected_to", "ancestry_alphabet"):
            v = r.get(key)
            if isinstance(v, (list, tuple)) and len(v):
                return tuple(int(x) for x in v)
        return None

    alphabets = {str(r.get("experiment")): _alphabet(r) for r in survivors}
    known = {a for a in alphabets.values() if a is not None}
    unknown = sorted(k for k, a in alphabets.items() if a is None)
    decision["ancestry_alphabets"] = {k: (list(v) if v else None) for k, v in alphabets.items()}
    decision["alphabets_unknown"] = unknown
    # Checkable only when at least one recipe reports its alphabet. A recipe
    # that does not (scored by an older scorer) cannot be assumed to match.
    decision["alphabets_checked"] = bool(known)
    distinct = known if not (known and unknown) else known | {None}
    decision["alphabets_comparable"] = len(distinct) <= 1
    if not mendel_separates:
        decision["status"] = "metric_invalid_mendelian_does_not_beat_controls"
        return decision
    if not survivors:
        decision["status"] = "no_survivors"
        return decision
    decision["allow_mixed_alphabets"] = bool(allow_mixed_alphabets)
    mixed = len(distinct) > 1
    if mixed and not allow_mixed_alphabets:
        decision["status"] = "metric_incomparable_mixed_ancestry_alphabets"
        decision["winner"] = None
        return decision
    decision["mixed_alphabets_override"] = mixed
    ranked = sorted(survivors, key=lambda r: _num(r.get("violations_per_informative_locus")))
    best = ranked[0]
    best_hi = _num(best.get("violations_ci_high"))
    tied = []
    for r in ranked[1:]:
        lo = _num(r.get("violations_ci_low"))
        if best_hi is None or lo is None or lo <= best_hi:
            tied.append(r["experiment"])
    decision["ranking"] = [
        {
            "experiment": r["experiment"],
            "violations_per_informative_locus": _num(r.get("violations_per_informative_locus")),
            "ci": [_num(r.get("violations_ci_low")), _num(r.get("violations_ci_high"))],
            "excess_recomb_over_expected": _num(r.get("excess_recomb_over_expected")),
            "mean_ll": _num(r.get("mean_ll")),
        }
        for r in ranked
    ]
    if mixed:
        sizes = {k: (len(v) if v else None) for k, v in alphabets.items()}
        decision["alphabet_sizes"] = sizes
        decision["label_count_bias_per_label"] = LABEL_COUNT_BIAS_PER_LABEL
        k_best = sizes.get(str(best["experiment"]))
        best_v = _num(best.get("violations_per_informative_locus"))
        pairwise: dict[str, Any] = {}
        within_bias: list[str] = []
        for r in ranked[1:]:
            name = str(r["experiment"])
            k_c = sizes.get(name)
            v_c = _num(r.get("violations_per_informative_locus"))
            gap = (v_c / best_v - 1.0) if (best_v and v_c is not None) else None
            if k_best is None or k_c is None:
                bias, allowance = "unknown", None
            elif k_best < k_c:
                # best is scored on fewer labels: its rate is biased LOW.
                bias = "favours_best"
                allowance = (1.0 + LABEL_COUNT_BIAS_PER_LABEL) ** (k_c - k_best) - 1.0
            elif k_best > k_c:
                bias, allowance = "conservative_for_best", 0.0
            else:
                bias, allowance = "none", 0.0
            ambiguous = bias == "unknown" or (
                bias == "favours_best" and (gap is None or gap <= allowance)
            )
            pairwise[name] = {
                "labels_best": k_best, "labels_competitor": k_c,
                "relative_gap": gap, "bias": bias,
                "bias_allowance": allowance, "within_bias": bool(ambiguous),
            }
            if ambiguous:
                within_bias.append(name)
        decision["label_count_comparison_vs_best"] = pairwise
        decision["caveat"] = (
            "Recipes were ranked across different ancestry alphabets "
            "(allow_mixed_alphabets=True). At identical accuracy the violation "
            "rate is 9-26% higher per extra label in synthetic checks (flare/eval/"
            "mendelian_k_bias; test_flare_alphabet_pipeline), so a recipe on FEWER "
            "labels is favoured and one on MORE labels is disadvantaged. A winner is "
            "declared only when it is not within that allowance of a runner-up with "
            "more labels. To settle a tie exactly, score the pin on the FLARE2 "
            "alphabet (--project-labels) and rank without the override."
        )
        tied = sorted(set(tied) | set(within_bias), key=lambda n: [str(x["experiment"]) for x in ranked].index(n))
        if within_bias:
            decision["tie_reason"] = "within_label_count_bias"
        # A lead taken by a recipe on FEWER labels (or on an unknown count) cannot be
        # certified by the allowance alone: the label-count effect depends on label
        # frequencies and was larger than the allowance in one synthetic check. Only a
        # lead held against the bias (best on MORE labels) is final.
        provisional = [n for n, c in pairwise.items() if c["bias"] in ("favours_best", "unknown")]
        decision["provisional"] = bool(provisional)
        decision["provisional_against"] = provisional
        if provisional:
            decision["provisional_reason"] = (
                "best recipe is scored on fewer ancestry labels than "
                + ", ".join(provisional)
                + "; confirm by scoring the pin on the common alphabet (PIN_SAS_TO) "
                "and ranking without the override"
            )
    decision["tied"] = tied
    if tied:
        decision["status"] = "tie_human_decision"
        decision["winner"] = None
        decision["best_point_estimate"] = best["experiment"]
    else:
        decision["status"] = "winner"
        decision["winner"] = best["experiment"]
    return decision


def apply_eval_gates(
    rows: list[dict[str, Any]],
    gates: dict[str, Any] | Path | str,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    """Split rows into (survivors, eliminated) using empirical ``eval_gates.json``.

    Expected keys (set after first notebook pass; absent key ⇒ no filter)::

        {
          "min_mean_ll": -0.8,
          "max_violations_per_informative_locus": 0.05,
          "max_excess_recomb_over_expected": 2.0
        }
    """
    if not isinstance(gates, dict):
        gates = json.loads(Path(gates).read_text())
    survivors: list[dict[str, Any]] = []
    eliminated: list[dict[str, Any]] = []
    for r in rows:
        reasons: list[str] = []
        min_ll = gates.get("min_mean_ll")
        if min_ll is not None and r.get("mean_ll") is not None:
            try:
                if float(r["mean_ll"]) < float(min_ll):
                    reasons.append("mean_ll")
            except (TypeError, ValueError):
                reasons.append("mean_ll_missing")
        max_viol = gates.get("max_violations_per_informative_locus")
        if max_viol is not None and r.get("violations_per_informative_locus") is not None:
            try:
                if float(r["violations_per_informative_locus"]) > float(max_viol):
                    reasons.append("mendelian_violations")
            except (TypeError, ValueError):
                reasons.append("mendelian_missing")
        max_xo = gates.get("max_excess_recomb_over_expected")
        if max_xo is not None and r.get("excess_recomb_over_expected") is not None:
            try:
                if float(r["excess_recomb_over_expected"]) > float(max_xo):
                    reasons.append("excess_recomb")
            except (TypeError, ValueError):
                reasons.append("excess_recomb_missing")
        if reasons:
            eliminated.append({**r, "gate_fail": reasons})
        else:
            survivors.append(r)
    return survivors, eliminated


def main(argv: Optional[list[str]] = None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--tsv", type=Path, default=None)
    p.add_argument("--from-firecloud", action="store_true")
    p.add_argument("--namespace", default=DEFAULT_NAMESPACE)
    p.add_argument("--workspace", default=DEFAULT_WORKSPACE)
    p.add_argument("--entity-type", default=DEFAULT_ENTITY_TYPE)
    p.add_argument("--out-tsv", type=Path, default=None)
    p.add_argument("--out-status", type=Path, default=None)
    args = p.parse_args(argv)

    df = fetch_lai_exp_table(
        tsv=args.tsv,
        from_firecloud=args.from_firecloud or args.tsv is None,
        namespace=args.namespace,
        workspace=args.workspace,
        entity_type=args.entity_type,
    )
    status = status_frame(df)
    print(status.to_string(index=False))
    if args.out_tsv:
        df.to_csv(args.out_tsv, sep="\t", index=False)
        print("wrote", args.out_tsv)
    if args.out_status:
        status.to_csv(args.out_status, sep="\t", index=False)
        print("wrote", args.out_status)
    n_done = int(status["complete"].sum())
    print(f"complete={n_done}/{len(status)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
