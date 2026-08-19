#!/usr/bin/env python3
"""Prepare matched association inputs from source-rebuilt covariates.

Writes a shared recommended-model-complete phenotype/covariate table, explicit
limited/full covariate matrices, and their column lists for matched Tractor-Mix
and SAIGE runs. ``full`` means the recommended production model: limited
covariates plus standardized coverage and genome-center dummies. It deliberately
does not include extraction method, platform, methylation caller, or SV counts.
"""

from __future__ import annotations

import argparse
import gzip
import re
from pathlib import Path
from typing import Iterable

import pandas as pd


TRUE_VALUES = {"TRUE", "True", "true", "1", "1.0", "T", "t", "YES", "Yes", "yes"}
FALSE_VALUES = {"FALSE", "False", "false", "0", "0.0", "F", "f", "NO", "No", "no"}
SEX_FEMALE = {
    "F",
    "f",
    "female",
    "Female",
    "2",
    "XX",
    "xx",
}
SEX_MALE = {
    "M",
    "m",
    "male",
    "Male",
    "1",
    "XY",
    "xy",
}
MISSING_STR = {"", "nan", "None", "<NA>", "NA", "NaN", "null", "NULL"}


def open_text(path: Path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "rt")


def read_id_list(path: Path) -> list[str]:
    ids: list[str] = []
    with open_text(path) as fh:
        for line in fh:
            s = line.strip()
            if s and not s.startswith("#"):
                ids.append(s.split()[0])
    return ids


def encode_binary(series: pd.Series | Iterable) -> pd.Series:
    series = pd.Series(series)
    out = pd.Series(pd.NA, index=series.index, dtype="Float64")
    as_str = series.astype(str).str.strip()
    out[as_str.isin(TRUE_VALUES)] = 1.0
    out[as_str.isin(FALSE_VALUES)] = 0.0
    out[as_str.isin(MISSING_STR)] = pd.NA
    return out


def encode_sex(series: pd.Series | Iterable) -> pd.Series:
    series = pd.Series(series)
    out = pd.Series(pd.NA, index=series.index, dtype="Float64")
    as_str = series.astype(str).str.strip()
    out[as_str.isin(SEX_FEMALE)] = 0.0
    out[as_str.isin(SEX_MALE)] = 1.0
    return out


def pick_pc_columns(columns: Iterable[str], n_pcs: int) -> list[str]:
    preferred = [f"PC{i}" for i in range(1, n_pcs + 1)]
    present = [c for c in preferred if c in columns]
    if len(present) < n_pcs:
        raise SystemExit(
            f"Need PC1..PC{n_pcs} in covariates; found {present}. "
            f"Available PC-like columns: {[c for c in columns if c.upper().startswith('PC')][:40]}"
        )
    return present


def sanitize_token(value: str) -> str:
    tok = re.sub(r"[^A-Za-z0-9]+", "_", value.strip())
    tok = re.sub(r"_+", "_", tok).strip("_")
    return tok or "UNK"


def one_hot_categorical(
    series: pd.Series,
    *,
    prefix: str,
) -> tuple[pd.DataFrame, str, list[str], pd.DataFrame]:
    """Return dummy frame, reference level, non-ref columns, and documentation rows."""
    clean = series.astype(str).str.strip()
    clean = clean.mask(clean.isin(MISSING_STR), other=pd.NA)
    counts = clean.dropna().value_counts()
    if counts.empty:
        raise SystemExit(f"No non-missing values for categorical covariate {prefix}")
    # Deterministic reference = most common level; ties broken alphabetically
    max_count = int(counts.max())
    ref_candidates = sorted(counts[counts == max_count].index.tolist())
    ref = ref_candidates[0]
    levels = sorted(counts.index.tolist())
    dummy_cols: list[str] = []
    out = pd.DataFrame(index=series.index)
    docs = []
    for level in levels:
        col = f"{prefix}_{sanitize_token(str(level))}"
        is_ref = level == ref
        if is_ref:
            docs.append(
                {
                    "source_column": prefix,
                    "level": level,
                    "dummy_column": col,
                    "is_reference": True,
                    "n": int(counts[level]),
                    "included_in_model": False,
                }
            )
            continue
        out[col] = (clean == level).astype("Float64")
        out.loc[clean.isna(), col] = pd.NA
        dummy_cols.append(col)
        docs.append(
            {
                "source_column": prefix,
                "level": level,
                "dummy_column": col,
                "is_reference": False,
                "n": int(counts[level]),
                "included_in_model": True,
            }
        )
    return out, str(ref), dummy_cols, pd.DataFrame(docs)


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--phenotype-csv", required=True, type=Path)
    p.add_argument("--covariates-csv", required=True, type=Path)
    p.add_argument(
        "--vcf-samples",
        required=True,
        type=Path,
        help="Sample IDs present in the FLARE VCF (bcftools query -l), in VCF order.",
    )
    p.add_argument("--out-dir", required=True, type=Path)
    p.add_argument("--n-phenotypes", type=int, default=10)
    p.add_argument("--min-cases", type=int, default=100)
    p.add_argument("--n-pcs", type=int, default=10)
    p.add_argument(
        "--covariate-id-col",
        default="research_id",
        help="ID column in covariates (fallback: sample_id).",
    )
    p.add_argument("--phenotype-id-col", default="person_id")
    p.add_argument(
        "--keep-withdrawn",
        action="store_true",
        help="Keep covariates rows with withdrawn==True (default: drop them).",
    )
    args = p.parse_args()
    args.drop_withdrawn = not args.keep_withdrawn
    args.out_dir.mkdir(parents=True, exist_ok=True)

    vcf_samples = read_id_list(args.vcf_samples)
    vcf_set = set(vcf_samples)
    print(f"VCF samples: {len(vcf_samples):,}")

    pheno = pd.read_csv(args.phenotype_csv, dtype=str, low_memory=False)
    if args.phenotype_id_col not in pheno.columns:
        raise SystemExit(f"Missing phenotype id column {args.phenotype_id_col}")
    pheno[args.phenotype_id_col] = pheno[args.phenotype_id_col].astype(str)
    drop_cols = [c for c in pheno.columns if c.startswith("Unnamed") or c == ""]
    pheno = pheno.drop(columns=drop_cols, errors="ignore")

    cov = pd.read_csv(args.covariates_csv, dtype=str, low_memory=False)
    id_col = args.covariate_id_col
    if id_col not in cov.columns:
        if "sample_id" in cov.columns:
            id_col = "sample_id"
        else:
            raise SystemExit(
                f"Covariates missing {args.covariate_id_col} and sample_id; "
                f"have {list(cov.columns)[:20]}"
            )
    cov[id_col] = cov[id_col].astype(str)

    if args.drop_withdrawn and "withdrawn" in cov.columns:
        withdrawn_mask = cov["withdrawn"].astype(str).str.lower().isin(
            {"true", "1", "1.0", "t", "yes"}
        )
        n_w = int(withdrawn_mask.sum())
        cov = cov.loc[~withdrawn_mask].copy()
        print(f"Dropped withdrawn covariates rows: {n_w}")

    for required in ["GC", "coverage"]:
        if required not in cov.columns:
            raise SystemExit(f"Covariates missing required technical column: {required}")

    pc_cols = pick_pc_columns(cov.columns, args.n_pcs)

    sex_candidates = []
    # Prefer rebuilt-table demographics / inferred sex before phenotype sex
    # (phenotype "sex" is sometimes present but unusable).
    if "sex_at_birth" in cov.columns:
        sex_candidates.append(("cov", "sex_at_birth"))
    if "inferred_sex" in cov.columns:
        sex_candidates.append(("cov", "inferred_sex"))
    if "sex" in pheno.columns:
        sex_candidates.append(("pheno", "sex"))
    if not sex_candidates:
        raise SystemExit("No sex column found in phenotypes or covariates")

    # Prefer the source that encodes to mostly non-missing on the VCF∩pheno∩cov set.
    pheno_ids = set(pheno[args.phenotype_id_col])
    cov_ids = set(cov[id_col])
    keep_ids = [s for s in vcf_samples if s in pheno_ids and s in cov_ids]
    print(
        f"Intersected analysis samples (VCF order): {len(keep_ids):,} "
        f"(vcf={len(vcf_set):,}, pheno={len(pheno_ids):,}, cov={len(cov_ids):,})"
    )
    if not keep_ids:
        raise SystemExit("No overlapping sample IDs across VCF, phenotypes, and covariates")

    pheno_i = pheno.set_index(args.phenotype_id_col).loc[keep_ids]
    cov_i = cov.drop_duplicates(subset=[id_col]).set_index(id_col).loc[keep_ids]

    best = None
    for src, col in sex_candidates:
        series = pheno_i[col] if src == "pheno" else cov_i[col]
        enc = encode_sex(series)
        n_ok = int(enc.notna().sum())
        print(f"Sex source candidate {src}.{col}: {n_ok:,}/{len(keep_ids):,} encoded")
        if best is None or n_ok > best[0]:
            best = (n_ok, src, col, enc)
    if best is None or best[0] == 0:
        raise SystemExit(
            "Could not encode sex from any candidate column; "
            f"tried {sex_candidates}. Check value spellings."
        )
    sex_source = (best[1], best[2])
    print(f"Using sex from {sex_source[0]}.{sex_source[1]}")

    age_col = "age" if "age" in cov.columns else None
    skip = {args.phenotype_id_col, "sex"}
    candidate_cols = [c for c in pheno.columns if c not in skip]

    # Build covariate columns on the intersected set first
    table = pd.DataFrame({"ID": keep_ids})
    table["sex"] = best[3].to_numpy()

    if age_col:
        table["age"] = pd.to_numeric(cov_i[age_col], errors="coerce").to_numpy()

    for pc in pc_cols:
        table[pc] = pd.to_numeric(cov_i[pc], errors="coerce").to_numpy()

    coverage = pd.to_numeric(cov_i["coverage"], errors="coerce")
    table["coverage_raw"] = coverage.to_numpy()

    gc_dummies, gc_ref, gc_cols, gc_docs = one_hot_categorical(cov_i["GC"], prefix="GC")
    for col in gc_cols:
        table[col] = gc_dummies[col].to_numpy()

    limited_cols = ["sex"]
    if age_col:
        limited_cols.append("age")
    limited_cols.extend(pc_cols)

    # Standardize coverage on samples that will remain after full completeness
    # First drop incomplete limited + technical covariates (except standardized coverage)
    tech_required = ["coverage_raw"] + gc_cols
    full_complete_cols = ["sex"] + pc_cols + tech_required
    if age_col:
        n_age_ok = int(table["age"].notna().sum())
        if n_age_ok == 0:
            print(
                "WARNING: age is entirely missing on the intersected cohort; "
                "omitting age from covariate lists"
            )
            table = table.drop(columns=["age"])
            age_col = None
            limited_cols = [c for c in limited_cols if c != "age"]
        else:
            full_complete_cols.append("age")

    print("Missingness before full-model filter:")
    for col in full_complete_cols:
        n_miss = int(pd.isna(table[col]).sum())
        print(f"  {col}: {n_miss:,} / {len(table):,} missing ({100.0 * n_miss / len(table):.1f}%)")

    before = len(table)
    table = table.dropna(subset=full_complete_cols).copy()
    print(
        f"Dropped {before - len(table)} samples missing full-model covariates; "
        f"remaining {len(table):,}"
    )
    if table.empty:
        raise SystemExit(
            "No samples remain after requiring full-model covariate completeness. "
            "See missingness counts above (often sex encoding, age, or coverage)."
        )

    # Standardize coverage on the shared full-complete cohort
    cov_mean = float(table["coverage_raw"].mean())
    cov_sd = float(table["coverage_raw"].std(ddof=0))
    if cov_sd == 0:
        raise SystemExit("coverage has zero variance in the analysis cohort")
    table["coverage"] = (table["coverage_raw"] - cov_mean) / cov_sd
    table = table.drop(columns=["coverage_raw"])

    full_cols = limited_cols + ["coverage"] + gc_cols

    # Phenotype selection on the shared full-complete cohort (VCF order preserved)
    analysis_ids = table["ID"].tolist()
    pheno_complete = pheno_i.loc[analysis_ids]

    ranks = []
    for col in candidate_cols:
        enc = encode_binary(pheno_complete[col])
        n_case = int((enc == 1.0).sum())
        n_ctrl = int((enc == 0.0).sum())
        if n_case < args.min_cases:
            continue
        balance = min(n_case, n_ctrl)
        ranks.append((balance, n_case, n_ctrl, col))
    ranks.sort(reverse=True)
    selected = [r[3] for r in ranks[: args.n_phenotypes]]
    if len(selected) < args.n_phenotypes:
        print(
            f"WARNING: only {len(selected)} phenotypes met min_cases={args.min_cases}; "
            f"requested {args.n_phenotypes}"
        )
    if not selected:
        raise SystemExit("No phenotypes passed selection criteria")

    selected_stats = []
    for col in selected:
        enc = encode_binary(pheno_complete[col])
        safe = col.replace(".", "_").replace("-", "_").replace(" ", "_")
        table[safe] = enc.to_numpy()
        selected_stats.append(
            {
                "phenotype": col,
                "safe_name": safe,
                "n_cases": int((enc == 1.0).sum()),
                "n_controls": int((enc == 0.0).sum()),
                "n_missing": int(enc.isna().sum()),
            }
        )

    out_samples = args.out_dir / "analysis_samples.txt"
    out_pheno = args.out_dir / "pheno_cov.tsv"
    out_limited_matrix = args.out_dir / "covariates_limited.tsv"
    out_full_matrix = args.out_dir / "covariates_full.tsv"
    out_selected = args.out_dir / "selected_phenotypes.txt"
    out_stats = args.out_dir / "selected_phenotype_stats.tsv"
    out_limited = args.out_dir / "covariate_columns_limited.txt"
    out_full = args.out_dir / "covariate_columns_full.txt"
    out_legacy = args.out_dir / "covariate_columns.txt"
    out_levels = args.out_dir / "technical_covariate_levels.tsv"

    if not analysis_ids:
        raise SystemExit("Refusing to write empty analysis_samples.txt")
    Path(out_samples).write_text("\n".join(analysis_ids) + "\n")
    table.to_csv(out_pheno, sep="\t", index=False, na_rep="NA")
    table[["ID", *limited_cols]].to_csv(out_limited_matrix, sep="\t", index=False, na_rep="NA")
    table[["ID", *full_cols]].to_csv(out_full_matrix, sep="\t", index=False, na_rep="NA")
    Path(out_selected).write_text("\n".join(s["safe_name"] for s in selected_stats) + "\n")
    pd.DataFrame(selected_stats).to_csv(out_stats, sep="\t", index=False)
    Path(out_limited).write_text("\n".join(limited_cols) + "\n")
    Path(out_full).write_text("\n".join(full_cols) + "\n")
    # Back-compat default points at limited covariates
    Path(out_legacy).write_text("\n".join(limited_cols) + "\n")

    level_docs = gc_docs.copy()
    level_docs["coverage_mean"] = cov_mean
    level_docs["coverage_sd"] = cov_sd
    level_docs["n_analysis_samples"] = len(analysis_ids)
    level_docs["gc_reference"] = gc_ref
    level_docs.to_csv(out_levels, sep="\t", index=False)

    print("Wrote:")
    for path in [
        out_samples,
        out_pheno,
        out_limited_matrix,
        out_full_matrix,
        out_selected,
        out_stats,
        out_limited,
        out_full,
        out_legacy,
        out_levels,
    ]:
        print(f"  {path}")
    print(f"Limited covariates: {limited_cols}")
    print(f"Full covariates: {full_cols}")
    print(f"GC reference={gc_ref}")
    print(f"coverage standardized: mean={cov_mean:.4f}, sd={cov_sd:.4f}")
    print("Selected phenotypes:")
    for s in selected_stats:
        print(
            f"  {s['phenotype']} -> {s['safe_name']}: "
            f"cases={s['n_cases']}, controls={s['n_controls']}, missing={s['n_missing']}"
        )


if __name__ == "__main__":
    main()
