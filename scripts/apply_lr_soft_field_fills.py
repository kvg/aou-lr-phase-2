#!/usr/bin/env python3
"""Apply soft joint-callset fills for ``population`` / ``sex_at_birth``.

Replays ``tractor_mix/reference_controls/lr_soft_field_fills.tsv`` (and optionally
discovers any new joint-callset gaps with the same rules):

  - missing ``population`` + ``lr_pop_population`` → uppercase continental code
  - missing ``sex_at_birth`` + diploid ``inferred_sex`` ``XX``/``XY`` → Female/Male

Never overwrites non-missing ``sex_at_birth`` (including sentinels like
``PMI: Skip``). Does not invent labels without a source column.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

SEX_FROM_INFERRED = {"XX": "Female", "XY": "Male"}
ALLOWED_FIELDS = frozenset({"population", "sex_at_birth"})


def is_missing(value) -> bool:
    if value is None or (isinstance(value, float) and pd.isna(value)):
        return True
    s = str(value).strip()
    return s == "" or s.upper() in {"NA", "NAN", "NONE"}


def apply_audit(
    cov: pd.DataFrame,
    audit: pd.DataFrame,
    *,
    force: bool = False,
) -> tuple[pd.DataFrame, list[dict]]:
    """Apply audit rows. Returns (updated_frame, per-row status dicts)."""
    out = cov.copy()
    out["research_id"] = out["research_id"].astype(str)
    id_to_idx = {rid: i for i, rid in enumerate(out["research_id"])}
    statuses: list[dict] = []

    required = {"research_id", "field", "after"}
    missing_cols = required - set(audit.columns)
    if missing_cols:
        raise ValueError(f"audit TSV missing columns: {sorted(missing_cols)}")

    for _, row in audit.iterrows():
        rid = str(row["research_id"]).strip()
        field = str(row["field"]).strip()
        after = row["after"]
        method = str(row.get("method", "")).strip() or "audit"
        if field not in ALLOWED_FIELDS:
            raise ValueError(f"unsupported audit field {field!r} for {rid}")
        if field not in out.columns:
            raise ValueError(f"covariates missing column {field!r}")
        if rid not in id_to_idx:
            statuses.append(
                {
                    "research_id": rid,
                    "field": field,
                    "status": "missing_sample",
                    "method": method,
                    "after": after,
                }
            )
            continue

        idx = out.index[id_to_idx[rid]]
        before = out.at[idx, field]
        if not force and not is_missing(before):
            same = str(before).strip() == str(after).strip()
            statuses.append(
                {
                    "research_id": rid,
                    "field": field,
                    "status": "already_set" if same else "skipped_non_missing",
                    "method": method,
                    "before": before,
                    "after": after,
                }
            )
            continue

        out.at[idx, field] = after
        statuses.append(
            {
                "research_id": rid,
                "field": field,
                "status": "applied",
                "method": method,
                "before": before if not is_missing(before) else pd.NA,
                "after": after,
            }
        )
    return out, statuses


def discover_soft_fills(cov: pd.DataFrame) -> pd.DataFrame:
    """Find joint-callset gaps fillable by lr_pop_population / inferred_sex."""
    rows: list[dict] = []
    jc = cov["has_lr_pcs"] == True if "has_lr_pcs" in cov.columns else pd.Series(False, index=cov.index)

    for idx in cov.index[jc]:
        rid = str(cov.at[idx, "research_id"])

        if "population" in cov.columns and is_missing(cov.at[idx, "population"]):
            lr_pop = cov.at[idx, "lr_pop_population"] if "lr_pop_population" in cov.columns else pd.NA
            if not is_missing(lr_pop):
                after = str(lr_pop).strip().upper()
                rows.append(
                    {
                        "research_id": rid,
                        "field": "population",
                        "before": pd.NA,
                        "after": after,
                        "method": "from_lr_pop_population",
                        "source_value": str(lr_pop).strip().lower(),
                        "notes": "Joint-callset sample missing population; copy uppercase lr_pop_population",
                    }
                )

        if "sex_at_birth" in cov.columns and is_missing(cov.at[idx, "sex_at_birth"]):
            inferred = cov.at[idx, "inferred_sex"] if "inferred_sex" in cov.columns else pd.NA
            if is_missing(inferred):
                continue
            key = str(inferred).strip().upper()
            if key not in SEX_FROM_INFERRED:
                continue
            after = SEX_FROM_INFERRED[key]
            rows.append(
                {
                    "research_id": rid,
                    "field": "sex_at_birth",
                    "before": pd.NA,
                    "after": after,
                    "method": "from_inferred_sex",
                    "source_value": key,
                    "notes": (
                        f"Map inferred_sex {key} → {after}; "
                        "did not overwrite PMI:Skip/sentinels"
                    ),
                }
            )

    return pd.DataFrame(rows)


def apply_lr_soft_field_fills(
    cov: pd.DataFrame,
    audit: pd.DataFrame | None = None,
    *,
    discover: bool = False,
    force: bool = False,
) -> tuple[pd.DataFrame, pd.DataFrame, list[dict]]:
    """Return (updated_covariates, audit_applied_or_merged, status_rows)."""
    out = cov.copy()
    statuses: list[dict] = []
    audit_parts: list[pd.DataFrame] = []

    if audit is not None and len(audit):
        audit = audit.copy()
        audit["research_id"] = audit["research_id"].astype(str)
        out, st = apply_audit(out, audit, force=force)
        statuses.extend(st)
        audit_parts.append(audit)

    if discover:
        found = discover_soft_fills(out)
        if len(found):
            out, st = apply_audit(out, found, force=force)
            statuses.extend(st)
            audit_parts.append(found)

    if audit_parts:
        merged = pd.concat(audit_parts, ignore_index=True)
        merged = merged.drop_duplicates(subset=["research_id", "field"], keep="last")
    else:
        merged = pd.DataFrame(
            columns=[
                "research_id",
                "field",
                "before",
                "after",
                "method",
                "source_value",
                "notes",
            ]
        )
    return out, merged, statuses


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--covariates",
        type=Path,
        default=Path("tractor_mix/covariates.source_rebuilt.csv.gz"),
    )
    parser.add_argument(
        "--audit",
        type=Path,
        default=Path("tractor_mix/reference_controls/lr_soft_field_fills.tsv"),
        help="Audit TSV to replay (created if --discover finds new rows).",
    )
    parser.add_argument(
        "--out-covariates",
        type=Path,
        default=None,
        help="Defaults to --covariates (in-place).",
    )
    parser.add_argument(
        "--out-audit",
        type=Path,
        default=None,
        help="Defaults to --audit.",
    )
    parser.add_argument(
        "--discover",
        action="store_true",
        help="Also fill any remaining has_lr_pcs gaps with the same rules.",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Overwrite non-missing field values from the audit.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Report actions without writing covariates.",
    )
    args = parser.parse_args()
    out_cov = args.out_covariates or args.covariates
    out_audit = args.out_audit or args.audit

    cov = pd.read_csv(args.covariates, dtype={"research_id": str}, low_memory=False)
    audit = None
    if args.audit.is_file():
        audit = pd.read_csv(args.audit, sep="\t", dtype={"research_id": str})
    elif not args.discover:
        raise FileNotFoundError(
            f"audit not found: {args.audit} (pass --discover to find fills without an audit)"
        )

    out, merged_audit, statuses = apply_lr_soft_field_fills(
        cov,
        audit,
        discover=args.discover,
        force=args.force,
    )

    status_df = pd.DataFrame(statuses)
    applied = int((status_df["status"] == "applied").sum()) if len(status_df) else 0
    already = int((status_df["status"] == "already_set").sum()) if len(status_df) else 0
    skipped = int((status_df["status"] == "skipped_non_missing").sum()) if len(status_df) else 0
    missing = int((status_df["status"] == "missing_sample").sum()) if len(status_df) else 0

    print(f"Soft fills: applied={applied} already_set={already} skipped_non_missing={skipped} missing_sample={missing}")
    if len(status_df):
        print(status_df.to_string(index=False))

    if args.dry_run:
        print("Dry run: no files written")
        return

    out_audit.parent.mkdir(parents=True, exist_ok=True)
    merged_audit.to_csv(out_audit, sep="\t", index=False)
    out.to_csv(out_cov, index=False, compression="gzip")
    print(f"Wrote {out_cov}")
    print(f"Wrote {out_audit}")


if __name__ == "__main__":
    main()
