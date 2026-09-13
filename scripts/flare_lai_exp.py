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
    """Attach Part 2–4 association-facing scores to a compare row.

    Screens (allele / Mendelian) are gates; null-λ is decisive. Switch/tract
    metrics remain diagnostics only.
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
