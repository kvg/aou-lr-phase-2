#!/usr/bin/env python3
"""Fill missing ancestry / population for long-read samples.

Priority:
  1. If continental ``population`` is present but ancestry is missing, copy it into
     ``ancestry_pred_other`` (lowercase). Copy into ``ancestry_pred`` only when the
     population is a hard continental label (never ``OTH`` → ``oth``).
  2. Else if ``has_lr_pcs`` and ancestry+population are missing, infer from
     ``lr_PC1``–``lr_PC10`` via nearest centroid + kNN (k=15). ``ancestry_pred``
     is taken from hard-labeled neighbors only; ``ancestry_pred_other`` may be
     ``oth``.
  3. If ``ancestry_pred`` is still missing or ``oth`` but ``has_lr_pcs``, assign a
     hard continental label from hard-labeled neighbors (does not overwrite
     ``ancestry_pred_other`` or ``population``).

Writes an audit TSV of applied fills. Does not overwrite existing non-missing
hard ``ancestry_pred``, non-missing ``ancestry_pred_other``, or ``population``
except for rule 3's hard ``ancestry_pred`` remediation of ``oth``.

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
# Hard six-class labels for ancestry_pred (external pred without oth).
HARD_ANCESTRY = frozenset({"afr", "amr", "eas", "eur", "sas", "mid"})
APPLIED_METHODS = frozenset(
    {"from_population", "lr_pc_knn", "lr_pc_knn_hard_pred"}
)


def is_missing(s: pd.Series) -> pd.Series:
    return s.isna() | (s.astype(str).str.strip() == "") | (s.astype(str).str.upper() == "NA")


def is_missing_value(value) -> bool:
    if pd.isna(value):
        return True
    text = str(value).strip()
    return text == "" or text.upper() == "NA"


def _norm_label(value) -> str:
    if is_missing_value(value):
        return ""
    return str(value).strip().lower()


def is_hard_ancestry(label: str) -> bool:
    return _norm_label(label) in HARD_ANCESTRY


def needs_hard_ancestry_pred(s: pd.Series) -> pd.Series:
    """Missing or soft ``oth`` — both need a hard continental ancestry_pred."""
    lab = s.astype(str).str.strip().str.lower()
    return is_missing(s) | (lab == "oth")


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
    """Copy continental population → ancestry_* when ancestry is missing.

    ``OTH`` fills ``ancestry_pred_other`` only; never writes ``oth`` into
    ``ancestry_pred`` (hard labels come from rule B/C).
    """
    miss_other = is_missing(df["ancestry_pred_other"])
    miss_pred = is_missing(df["ancestry_pred"])
    has_pop = ~is_missing(df["population"])
    mask = is_long_read(df) & has_pop & (miss_other | miss_pred)
    audits: list[dict] = []
    for idx in df.index[mask]:
        pop = str(df.at[idx, "population"]).strip().upper()
        label = pop.lower()
        wrote_other = False
        wrote_pred = False
        if is_missing_value(df.at[idx, "ancestry_pred_other"]):
            df.at[idx, "ancestry_pred_other"] = label
            wrote_other = True
        if is_missing_value(df.at[idx, "ancestry_pred"]) and is_hard_ancestry(label):
            df.at[idx, "ancestry_pred"] = label
            wrote_pred = True
        if not (wrote_other or wrote_pred):
            continue
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
                "ancestry_pred_after": df.at[idx, "ancestry_pred"],
                "ancestry_pred_other_after": df.at[idx, "ancestry_pred_other"],
            }
        )
    return audits


def _knn_vote(
    x: np.ndarray,
    X_ref_z: np.ndarray,
    y_ref: np.ndarray,
    centroids: dict[str, np.ndarray],
) -> tuple[str, str, float, float, Counter, bool]:
    """Return (knn_label, centroid_label, knn_frac, margin, counts, high_conf)."""
    dists = {lab: float(np.linalg.norm(x - c)) for lab, c in centroids.items()}
    ordered = sorted(dists.items(), key=lambda kv: kv[1])
    lab_c, d1 = ordered[0]
    margin = ordered[1][1] - d1 if len(ordered) > 1 else float("inf")

    d = np.linalg.norm(X_ref_z - x, axis=1)
    nn = np.argpartition(d, K)[:K]
    counts = Counter(y_ref[nn])
    lab_k, n_top = counts.most_common(1)[0]
    frac = n_top / K

    agree = lab_c == lab_k
    high_conf = bool(agree and (frac >= HIGH_CONF_KNN_FRAC or margin >= HIGH_CONF_MARGIN))
    return lab_k, lab_c, frac, margin, counts, high_conf


def _build_ref(
    df: pd.DataFrame, *, labels: pd.Series, allowed: frozenset[str] | None
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, dict[str, np.ndarray]]:
    mask = (df["has_lr_pcs"] == True) & ~is_missing(labels)
    if allowed is not None:
        mask = mask & labels.astype(str).str.strip().str.lower().isin(allowed)
    ref = df.loc[mask]
    if len(ref) < K:
        raise ValueError(
            f"need at least {K} labeled lr_PC references; found {len(ref)}"
        )
    X_ref = ref[LR_PC_INFER].to_numpy(dtype=float)
    y_ref = labels.loc[ref.index].astype(str).str.lower().to_numpy()
    mu = X_ref.mean(axis=0)
    sd = X_ref.std(axis=0)
    sd[sd == 0] = 1.0
    X_ref_z = (X_ref - mu) / sd
    centroids = {lab: X_ref_z[y_ref == lab].mean(axis=0) for lab in sorted(set(y_ref))}
    return X_ref_z, y_ref, mu, sd, centroids


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

    # Soft labels (may include oth) for ancestry_pred_other / population.
    X_soft, y_soft, mu_s, sd_s, cent_soft = _build_ref(
        df, labels=df["ancestry_pred_other"], allowed=None
    )
    # Hard labels only for ancestry_pred.
    hard_labels = df["ancestry_pred"].where(
        df["ancestry_pred"].astype(str).str.strip().str.lower().isin(HARD_ANCESTRY),
        other=pd.NA,
    )
    # Prefer existing hard ancestry_pred; fall back to hard ancestry_pred_other.
    hard_src = hard_labels.fillna(
        df["ancestry_pred_other"].where(
            df["ancestry_pred_other"].astype(str).str.strip().str.lower().isin(HARD_ANCESTRY),
            other=pd.NA,
        )
    )
    X_hard, y_hard, mu_h, sd_h, cent_hard = _build_ref(
        df, labels=hard_src, allowed=HARD_ANCESTRY
    )

    audits: list[dict] = []
    for idx in gap_idx:
        x_soft = (df.loc[idx, LR_PC_INFER].to_numpy(dtype=float) - mu_s) / sd_s
        lab_soft, lab_c_soft, frac_s, margin_s, counts_s, high_s = _knn_vote(
            x_soft, X_soft, y_soft, cent_soft
        )
        x_hard = (df.loc[idx, LR_PC_INFER].to_numpy(dtype=float) - mu_h) / sd_h
        lab_hard, lab_c_hard, frac_h, margin_h, counts_h, high_h = _knn_vote(
            x_hard, X_hard, y_hard, cent_hard
        )

        high_conf = high_h  # gate on hard vote quality
        if not high_conf and not fill_low_confidence:
            audits.append(
                {
                    "research_id": df.at[idx, "research_id"],
                    "method": "lr_pc_knn_skipped_low_conf",
                    "suggested": lab_hard,
                    "centroid": lab_c_hard,
                    "population_before": pd.NA,
                    "population_after": pd.NA,
                    "high_confidence": False,
                    "knn15_frac": frac_h,
                    "centroid_margin": margin_h,
                    "knn15_counts": dict(counts_h),
                    "lr_PC1": float(df.at[idx, "lr_PC1"]),
                    "lr_PC2": float(df.at[idx, "lr_PC2"]),
                    "ancestry_pred_after": pd.NA,
                    "ancestry_pred_other_after": pd.NA,
                }
            )
            continue

        df.at[idx, "ancestry_pred"] = lab_hard
        df.at[idx, "ancestry_pred_other"] = lab_soft
        df.at[idx, "has_ancestry_annotation"] = True
        pop_before = df.at[idx, "population"]
        pop_missing = is_missing_value(pop_before)
        if pop_missing:
            df.at[idx, "population"] = lab_soft.upper()
            pop_after = lab_soft.upper()
            pop_before_out = pd.NA
        else:
            pop_after = str(pop_before).strip().upper()
            pop_before_out = pop_after
        audits.append(
            {
                "research_id": df.at[idx, "research_id"],
                "method": "lr_pc_knn",
                "suggested": lab_hard,
                "centroid": lab_c_hard,
                "population_before": pop_before_out,
                "population_after": pop_after,
                "high_confidence": high_conf,
                "knn15_frac": frac_h,
                "centroid_margin": margin_h,
                "knn15_counts": dict(counts_h),
                "lr_PC1": float(df.at[idx, "lr_PC1"]),
                "lr_PC2": float(df.at[idx, "lr_PC2"]),
                "ancestry_pred_after": lab_hard,
                "ancestry_pred_other_after": lab_soft,
                "soft_suggested": lab_soft,
                "soft_centroid": lab_c_soft,
                "soft_knn15_frac": frac_s,
                "soft_centroid_margin": margin_s,
                "soft_knn15_counts": dict(counts_s),
                "soft_high_confidence": high_s,
            }
        )
    return audits


def fill_hard_ancestry_pred(df: pd.DataFrame, *, fill_low_confidence: bool) -> list[dict]:
    """Assign hard ancestry_pred when missing or oth (AoU + controls with lr_PCs).

    Does not modify ancestry_pred_other or population. Skips reference controls
    only when they lack lr_PCs; HG002-style oth controls with PCs are left as-is
    if ``is_reference_control`` (curated labels).
    """
    miss_hard = needs_hard_ancestry_pred(df["ancestry_pred"])
    # Only remediate non-controls; curated control oth (e.g. HG002) stays.
    is_ctrl = df.get("is_reference_control", False) == True
    mask = (
        is_long_read(df)
        & (df["has_lr_pcs"] == True)
        & miss_hard
        & ~is_ctrl
    )
    gap_idx = df.index[mask]
    if len(gap_idx) == 0:
        return []

    missing_pcs = [c for c in LR_PC_INFER if c not in df.columns]
    if missing_pcs:
        raise ValueError(
            "covariates missing lr_PC columns required for kNN "
            f"(run merge_lr_global_pcs_into_covariates.py first): {missing_pcs[:3]}"
        )

    hard_labels = df["ancestry_pred"].where(
        df["ancestry_pred"].astype(str).str.strip().str.lower().isin(HARD_ANCESTRY),
        other=pd.NA,
    )
    hard_src = hard_labels.fillna(
        df["ancestry_pred_other"].where(
            df["ancestry_pred_other"].astype(str).str.strip().str.lower().isin(HARD_ANCESTRY),
            other=pd.NA,
        )
    )
    X_hard, y_hard, mu_h, sd_h, cent_hard = _build_ref(
        df, labels=hard_src, allowed=HARD_ANCESTRY
    )

    audits: list[dict] = []
    for idx in gap_idx:
        before = _norm_label(df.at[idx, "ancestry_pred"]) or pd.NA
        x_hard = (df.loc[idx, LR_PC_INFER].to_numpy(dtype=float) - mu_h) / sd_h
        lab_hard, lab_c, frac, margin, counts, high_conf = _knn_vote(
            x_hard, X_hard, y_hard, cent_hard
        )
        if not high_conf and not fill_low_confidence:
            audits.append(
                {
                    "research_id": df.at[idx, "research_id"],
                    "method": "lr_pc_knn_hard_pred_skipped_low_conf",
                    "suggested": lab_hard,
                    "centroid": lab_c,
                    "population_before": str(df.at[idx, "population"]).strip().upper()
                    if not is_missing_value(df.at[idx, "population"])
                    else pd.NA,
                    "population_after": pd.NA,
                    "high_confidence": False,
                    "knn15_frac": frac,
                    "centroid_margin": margin,
                    "knn15_counts": dict(counts),
                    "lr_PC1": float(df.at[idx, "lr_PC1"]),
                    "lr_PC2": float(df.at[idx, "lr_PC2"]),
                    "ancestry_pred_before": before,
                    "ancestry_pred_after": pd.NA,
                    "ancestry_pred_other_after": df.at[idx, "ancestry_pred_other"],
                }
            )
            continue

        df.at[idx, "ancestry_pred"] = lab_hard
        df.at[idx, "has_ancestry_annotation"] = True
        pop = df.at[idx, "population"]
        pop_out = str(pop).strip().upper() if not is_missing_value(pop) else pd.NA
        audits.append(
            {
                "research_id": df.at[idx, "research_id"],
                "method": "lr_pc_knn_hard_pred",
                "suggested": lab_hard,
                "centroid": lab_c,
                "population_before": pop_out,
                "population_after": pop_out,
                "high_confidence": high_conf,
                "knn15_frac": frac,
                "centroid_margin": margin,
                "knn15_counts": dict(counts),
                "lr_PC1": float(df.at[idx, "lr_PC1"]),
                "lr_PC2": float(df.at[idx, "lr_PC2"]),
                "ancestry_pred_before": before,
                "ancestry_pred_after": lab_hard,
                "ancestry_pred_other_after": df.at[idx, "ancestry_pred_other"],
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
    audits.extend(fill_hard_ancestry_pred(out, fill_low_confidence=fill_low_confidence))
    audit_df = pd.DataFrame(audits)
    if len(audit_df) and "research_id" in audit_df.columns:
        audit_df["research_id"] = audit_df["research_id"].astype(str)
    return out, audit_df


def merge_audit(previous: pd.DataFrame | None, new: pd.DataFrame) -> pd.DataFrame:
    """Keep historical fills; replace prior hard-pred remediation rows for same IDs."""
    if previous is None or len(previous) == 0:
        return new
    if len(new) == 0:
        return previous.copy()
    prev = previous.copy()
    prev["research_id"] = prev["research_id"].astype(str)
    hard_ids = set(
        new.loc[new["method"] == "lr_pc_knn_hard_pred", "research_id"].astype(str)
    )
    if hard_ids and "method" in prev.columns:
        prev = prev[
            ~(
                (prev["method"] == "lr_pc_knn_hard_pred")
                & (prev["research_id"].astype(str).isin(hard_ids))
            )
        ]
    return pd.concat([prev, new], ignore_index=True)


def summarize_gaps(df: pd.DataFrame) -> dict[str, int]:
    lr = is_long_read(df)
    miss = ancestry_gap_mask(df)
    need_hard = needs_hard_ancestry_pred(df["ancestry_pred"])
    aou = df.get("is_reference_control", False) != True
    return {
        "n_rows": int(len(df)),
        "n_long_read": int(lr.sum()),
        "lr_missing_ancestry": int((lr & miss).sum()),
        "has_lr_pcs_missing_ancestry": int(((df["has_lr_pcs"] == True) & miss).sum()),
        "final_releasable_missing_ancestry": int(
            ((df["final_releasable_v9"] == True) & miss).sum()
        ),
        "aou_ancestry_pred_oth_or_missing": int((aou & need_hard).sum()),
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
    parser.add_argument(
        "--replace-audit",
        action="store_true",
        help="Overwrite audit TSV instead of merging with any existing file.",
    )
    args = parser.parse_args()
    out_cov = args.out_covariates or args.covariates

    df = pd.read_csv(args.covariates, dtype={"research_id": str}, low_memory=False)
    before = summarize_gaps(df)
    out, audit_df = fill_lr_ancestry(df, fill_low_confidence=not args.high_confidence_only)
    after = summarize_gaps(out)

    args.out_audit.parent.mkdir(parents=True, exist_ok=True)
    if args.replace_audit or not args.out_audit.exists():
        merged_audit = audit_df
    else:
        prev = pd.read_csv(args.out_audit, sep="\t", dtype={"research_id": str})
        merged_audit = merge_audit(prev, audit_df)
    merged_audit.to_csv(args.out_audit, sep="\t", index=False)
    out.to_csv(out_cov, index=False, compression="gzip")

    applied = (
        audit_df.loc[audit_df["method"].isin(APPLIED_METHODS)] if len(audit_df) else audit_df
    )
    skipped = (
        audit_df.loc[
            audit_df["method"].isin(
                {"lr_pc_knn_skipped_low_conf", "lr_pc_knn_hard_pred_skipped_low_conf"}
            )
        ]
        if len(audit_df)
        else audit_df
    )

    print(f"Wrote {out_cov}")
    print(f"Wrote {args.out_audit} ({len(merged_audit)} rows)")
    print(f"LR ancestry gaps: {before['lr_missing_ancestry']} → {after['lr_missing_ancestry']}")
    print(
        "AoU ancestry_pred oth/missing: "
        f"{before['aou_ancestry_pred_oth_or_missing']} → "
        f"{after['aou_ancestry_pred_oth_or_missing']}"
    )
    print(f"Fills applied this run: {len(applied)}")
    if len(applied):
        print(applied["method"].value_counts().to_string())
        print("suggested labels:")
        print(applied["suggested"].value_counts().to_string())
    if len(skipped):
        print(f"Low-confidence kNN skipped: {len(skipped)}")


if __name__ == "__main__":
    main()
