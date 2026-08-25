#!/usr/bin/env python3
"""Fill missing ancestry / population for long-read samples.

Priority:
  1. If continental ``population`` is present but ancestry is missing, copy it into
     ``ancestry_pred`` / ``ancestry_pred_other`` (lowercase).
  2. Else if ``has_lr_pcs`` and ancestry+population are missing, infer from
     ``lr_PC1``–``lr_PC10`` via nearest centroid + kNN (k=15).

Writes an audit TSV of applied fills. Does not overwrite existing non-missing
ancestry or population values.

Notebook entry point: ``fill_lr_ancestry``.
"""

from __future__ import annotations

import argparse
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd

LR_PC_INFER = [f"lr_PC{i}" for i in range(1, 11)]
K = 15
HIGH_CONF_KNN_FRAC = 0.8
HIGH_CONF_MARGIN = 0.5
APPLIED_METHODS = frozenset({"from_population", "lr_pc_knn"})


def is_missing(s: pd.Series) -> pd.Series:
    return s.isna() | (s.astype(str).str.strip() == "") | (s.astype(str).str.upper() == "NA")


def is_long_read(df: pd.DataFrame) -> pd.Series:
    return (
        (df["has_lr_pcs"] == True)
        | (df["has_lr_tech"] == True)
        | (df["lr_meet_qc"] == True)
        | (df["final_releasable_v9"] == True)
        | (~is_missing(df["lr_phase"]))
        | (df.get("is_reference_control", False) == True)
    )


def ancestry_gap_mask(df: pd.DataFrame) -> pd.Series:
    return is_missing(df["ancestry_pred"]) | is_missing(df["ancestry_pred_other"])


def fill_from_population(df: pd.DataFrame) -> list[dict]:
    """Copy continental population → ancestry_* when ancestry is missing."""
    miss_anc = ancestry_gap_mask(df)
    has_pop = ~is_missing(df["population"])
    mask = is_long_read(df) & miss_anc & has_pop
    audits: list[dict] = []
    for idx in df.index[mask]:
        pop = str(df.at[idx, "population"]).strip().upper()
        label = pop.lower()
        df.at[idx, "ancestry_pred"] = label
        df.at[idx, "ancestry_pred_other"] = label
        df.at[idx, "has_ancestry_annotation"] = True
        audits.append(
            {
                "research_id": df.at[idx, "research_id"],
                "method": "from_population",
                "suggested": label,
                "centroid": pd.NA,
                "population_before": pop,
                "population_after": pop,
                "high_confidence": True,
                "knn15_frac": pd.NA,
                "centroid_margin": pd.NA,
                "knn15_counts": pd.NA,
                "lr_PC1": df.at[idx, "lr_PC1"] if "lr_PC1" in df.columns else pd.NA,
                "lr_PC2": df.at[idx, "lr_PC2"] if "lr_PC2" in df.columns else pd.NA,
            }
        )
    return audits


def infer_knn(df: pd.DataFrame, *, fill_low_confidence: bool) -> list[dict]:
    """Infer ancestry+population from lr_PCs when both are missing."""
    miss_anc = ancestry_gap_mask(df)
    miss_pop = is_missing(df["population"])
    mask = is_long_read(df) & (df["has_lr_pcs"] == True) & miss_anc & miss_pop
    gap_idx = df.index[mask]
    if len(gap_idx) == 0:
        return []

    missing_pcs = [c for c in LR_PC_INFER if c not in df.columns]
    if missing_pcs:
        raise ValueError(
            "covariates missing lr_PC columns required for kNN "
            f"(run merge_lr_global_pcs_into_covariates.py first): {missing_pcs[:3]}"
        )

    ref_mask = (df["has_lr_pcs"] == True) & ~is_missing(df["ancestry_pred_other"])
    ref = df.loc[ref_mask]
    if len(ref) < K:
        raise ValueError(f"need at least {K} labeled lr_PC references; found {len(ref)}")

    X_ref = ref[LR_PC_INFER].to_numpy(dtype=float)
    y_ref = ref["ancestry_pred_other"].astype(str).str.lower().to_numpy()

    mu = X_ref.mean(axis=0)
    sd = X_ref.std(axis=0)
    sd[sd == 0] = 1.0
    X_ref_z = (X_ref - mu) / sd

    centroids = {
        lab: X_ref_z[y_ref == lab].mean(axis=0) for lab in sorted(set(y_ref))
    }

    audits: list[dict] = []
    for idx in gap_idx:
        x = (df.loc[idx, LR_PC_INFER].to_numpy(dtype=float) - mu) / sd

        dists = {lab: float(np.linalg.norm(x - c)) for lab, c in centroids.items()}
        ordered = sorted(dists.items(), key=lambda kv: kv[1])
        lab_c, d1 = ordered[0]
        margin = ordered[1][1] - d1

        d = np.linalg.norm(X_ref_z - x, axis=1)
        nn = np.argpartition(d, K)[:K]
        counts = Counter(y_ref[nn])
        lab_k, n_top = counts.most_common(1)[0]
        frac = n_top / K

        agree = lab_c == lab_k
        high_conf = bool(agree and (frac >= HIGH_CONF_KNN_FRAC or margin >= HIGH_CONF_MARGIN))
        if not high_conf and not fill_low_confidence:
            audits.append(
                {
                    "research_id": df.at[idx, "research_id"],
                    "method": "lr_pc_knn_skipped_low_conf",
                    "suggested": lab_k,
                    "centroid": lab_c,
                    "population_before": pd.NA,
                    "population_after": pd.NA,
                    "high_confidence": False,
                    "knn15_frac": frac,
                    "centroid_margin": margin,
                    "knn15_counts": dict(counts),
                    "lr_PC1": float(df.at[idx, "lr_PC1"]),
                    "lr_PC2": float(df.at[idx, "lr_PC2"]),
                }
            )
            continue

        label = lab_k  # prefer kNN vote
        df.at[idx, "ancestry_pred"] = label
        df.at[idx, "ancestry_pred_other"] = label
        df.at[idx, "has_ancestry_annotation"] = True
        pop_before = df.at[idx, "population"]
        pop_missing = pd.isna(pop_before) or str(pop_before).strip() in {"", "nan", "NA", "None"}
        if pop_missing:
            df.at[idx, "population"] = label.upper()
            pop_after = label.upper()
            pop_before_out = pd.NA
        else:
            pop_after = str(pop_before).strip().upper()
            pop_before_out = pop_after
        audits.append(
            {
                "research_id": df.at[idx, "research_id"],
                "method": "lr_pc_knn",
                "suggested": label,
                "centroid": lab_c,
                "population_before": pop_before_out,
                "population_after": pop_after,
                "high_confidence": high_conf,
                "knn15_frac": frac,
                "centroid_margin": margin,
                "knn15_counts": dict(counts),
                "lr_PC1": float(df.at[idx, "lr_PC1"]),
                "lr_PC2": float(df.at[idx, "lr_PC2"]),
            }
        )
    return audits


def fill_lr_ancestry(
    df: pd.DataFrame,
    *,
    fill_low_confidence: bool = True,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Return (updated_covariates, audit_table).

    ``df`` is copied; the input frame is not modified.
    """
    out = df.copy()
    if "research_id" not in out.columns:
        raise ValueError("covariates missing research_id")
    out["research_id"] = out["research_id"].astype(str)

    audits: list[dict] = []
    audits.extend(fill_from_population(out))
    audits.extend(infer_knn(out, fill_low_confidence=fill_low_confidence))
    audit_df = pd.DataFrame(audits)
    if len(audit_df) and "research_id" in audit_df.columns:
        audit_df["research_id"] = audit_df["research_id"].astype(str)
    return out, audit_df


def summarize_gaps(df: pd.DataFrame) -> dict[str, int]:
    lr = is_long_read(df)
    miss = ancestry_gap_mask(df)
    return {
        "n_rows": int(len(df)),
        "n_long_read": int(lr.sum()),
        "lr_missing_ancestry": int((lr & miss).sum()),
        "has_lr_pcs_missing_ancestry": int(((df["has_lr_pcs"] == True) & miss).sum()),
        "final_releasable_missing_ancestry": int(((df["final_releasable_v9"] == True) & miss).sum()),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--covariates",
        type=Path,
        default=Path("tractor_mix/covariates.source_rebuilt.csv.gz"),
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
        default=Path("tractor_mix/reference_controls/lr_ancestry_knn_fills.tsv"),
    )
    parser.add_argument(
        "--high-confidence-only",
        action="store_true",
        help="Skip low-confidence kNN fills.",
    )
    args = parser.parse_args()
    out_cov = args.out_covariates or args.covariates

    df = pd.read_csv(args.covariates, dtype={"research_id": str}, low_memory=False)
    before = summarize_gaps(df)
    out, audit_df = fill_lr_ancestry(df, fill_low_confidence=not args.high_confidence_only)
    after = summarize_gaps(out)

    args.out_audit.parent.mkdir(parents=True, exist_ok=True)
    audit_df.to_csv(args.out_audit, sep="\t", index=False)
    out.to_csv(out_cov, index=False, compression="gzip")

    applied = audit_df.loc[audit_df["method"].isin(APPLIED_METHODS)] if len(audit_df) else audit_df
    skipped = (
        audit_df.loc[audit_df["method"] == "lr_pc_knn_skipped_low_conf"]
        if len(audit_df)
        else audit_df
    )

    print(f"Wrote {out_cov}")
    print(f"Wrote {args.out_audit}")
    print(f"LR ancestry gaps: {before['lr_missing_ancestry']} → {after['lr_missing_ancestry']}")
    print(f"Fills applied: {len(applied)}")
    if len(applied):
        print(applied["method"].value_counts().to_string())
        print("suggested labels:")
        print(applied["suggested"].value_counts().to_string())
    if len(skipped):
        print(f"Low-confidence kNN skipped: {len(skipped)}")


if __name__ == "__main__":
    main()
