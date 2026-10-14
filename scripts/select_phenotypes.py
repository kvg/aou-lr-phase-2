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


def phecode_of(column: str) -> str:
    """Phecode from a phenotype column name, in the definitions table's dialect.

    PhecodeX keeps its chapter and uses a dot before the sub-code:
    ``CV_401_1`` and ``CV_401.1`` both give ``CV_401.1``. The chapter is part of
    the identifier, not a prefix to strip -- ``CM`` is one of the 18 PhecodeX
    chapters, so ``CM_761`` is a code in its own right.

    Classic phecodes are bare numbers, optionally behind ``PH_``/``PHECODE_``:
    ``PH_761_13`` and ``761.13`` both give ``761.13``.

    Returns "" for anything else, e.g. a covariate column.
    """
    s = column.strip()
    # PH_ / PHECODE_ first: neither is a PhecodeX chapter, so an explicit
    # classic prefix is unambiguous and wins over the chapter reading.
    body = re.sub(r"^(PH|PHECODE)[_.]", "", s, flags=re.IGNORECASE)
    if body == s:
        m = re.match(r"^([A-Za-z]{2})[_.](\d+)(?:[._](\d+))?$", s)
        if m:
            chapter, num, dec = m.groups()
            base = f"{chapter.upper()}_{num}"
            return base if dec is None else f"{base}.{dec}"
    m = re.match(r"^(\d+)(?:[._](\d+))?$", body)
    if not m:
        return ""
    return m.group(1) if m.group(2) is None else f"{m.group(1)}.{m.group(2)}"


def phecode_parent(code: str) -> str:
    """Parent in the phecode hierarchy: 761.13 -> 761.1 -> 761 -> ''."""
    if not code or "." not in code:
        return ""
    head, dec = code.split(".", 1)
    return f"{head}.{dec[:-1]}" if len(dec) > 1 else head


def read_phecode_sex(path: Path) -> dict[str, str]:
    """phecode -> Male|Female|Both from a phecode definitions table.

    The published PhecodeX table is Latin-1, not UTF-8 (it carries names such
    as Meniere's with an accented e), so decoding is attempted both ways before
    giving up.
    """
    for encoding in ("utf-8", "latin-1"):
        try:
            df = pd.read_csv(path, sep=None, engine="python", dtype=str,
                             encoding=encoding)
            break
        except UnicodeDecodeError:
            continue
    else:
        raise SystemExit(f"{path}: could not decode as UTF-8 or Latin-1")
    cols = {c.lower().strip(): c for c in df.columns}
    code_col = next((cols[c] for c in ("phecode", "phenotype", "code") if c in cols), None)
    sex_col = next((cols[c] for c in ("sex", "sex_specific") if c in cols), None)
    if code_col is None or sex_col is None:
        raise SystemExit(
            f"{path}: need a phecode/phenotype column and a sex column; have {list(df.columns)}"
        )
    out: dict[str, str] = {}
    for code, sex in zip(df[code_col], df[sex_col]):
        if pd.isna(code) or pd.isna(sex):
            continue
        norm = phecode_of(str(code)) or str(code).strip()
        s = str(sex).strip().lower()
        out[norm] = "Female" if s.startswith("f") else "Male" if s.startswith("m") else "Both"
    return out


def observed_case_sex(enc: pd.Series, sex: pd.Series) -> tuple[float, str]:
    """Fraction of cases in the majority sex, and which sex that is.

    sex is encoded 0 = female, 1 = male (``encode_sex``).
    """
    case_sex = sex[enc == 1.0].dropna()
    if case_sex.empty:
        return float("nan"), ""
    frac_male = float((case_sex == 1.0).mean())
    return (frac_male, "Male") if frac_male >= 0.5 else (1.0 - frac_male, "Female")


def prune_correlated(
    encoded: pd.DataFrame, order: list[str], max_r2: float
) -> dict[str, tuple[str, float]]:
    """Greedy prune in ``order`` (best first). Returns dropped -> (kept, r2).

    Pearson r on 0/1 vectors with pairwise-complete observations; r2 above the
    cutoff means the two phecodes carry nearly the same case set, which is what
    a parent/child pair or a synonym pair looks like.
    """
    dropped: dict[str, tuple[str, float]] = {}
    kept: list[str] = []
    corr = encoded[order].corr(min_periods=50)
    for col in order:
        hit = ""
        best = 0.0
        for k in kept:
            r = corr.at[col, k]
            if pd.isna(r):
                continue
            r2 = float(r) ** 2
            if r2 > best:
                best, hit = r2, k
        if hit and best >= max_r2:
            dropped[col] = (hit, best)
        else:
            kept.append(col)
    return dropped


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
    p.add_argument(
        "--n-phenotypes", type=int, default=10,
        help="Number to keep after filtering, ranked by min(n_cases, n_controls). "
             "0 keeps every phenotype that passes (use for the full phecode list).",
    )
    p.add_argument("--min-cases", type=int, default=100)
    p.add_argument(
        "--phecode-sex-tsv", type=Path, default=None,
        help="Phecode definitions table with a sex column (Male/Female/Both). Authoritative "
             "for --sex-specific; without it the call is made from the observed case sex ratio",
    )
    p.add_argument(
        "--sex-specific", choices=("flag", "restrict"), default="flag",
        help="flag (default): report sex-specific phecodes in the stats table and leave the "
             "phenotype alone. restrict: set the other sex to NA so they are not counted as "
             "controls. fit_felix_null.R drops `sex` from the model where it goes constant",
    )
    p.add_argument(
        "--sex-case-fraction", type=float, default=0.99,
        help="A phecode whose cases are at least this fraction one sex is called sex-specific "
             "when no definitions table is given [default 0.99]",
    )
    p.add_argument(
        "--nesting", choices=("keep", "leaves", "roots"), default="keep",
        help="keep (default): parent and child phecodes are separate phenotypes. leaves: drop a "
             "phecode that has a descendant in the set. roots: keep only integer-level codes. "
             "Counts are always of phecodes, not independent traits",
    )
    p.add_argument(
        "--max-r2", type=float, default=None,
        help="Greedily drop a phenotype whose r2 with an already-kept, better-powered one is at "
             "least this (All of Us All by All used 0.5). Off by default",
    )
    p.add_argument(
        "--tier-sizes", default="",
        help="Comma-separated prefix sizes, e.g. 10,100: also writes "
             "selected_phenotypes.top10.txt etc. from the head of the final ranking",
    )
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

    sex_by_id = pd.Series(table["sex"].to_numpy(), index=analysis_ids)
    phecode_sex = read_phecode_sex(args.phecode_sex_tsv) if args.phecode_sex_tsv else {}
    if args.phecode_sex_tsv:
        print(f"Read sex restrictions for {len(phecode_sex):,} phecodes from {args.phecode_sex_tsv}")

    encoded = pd.DataFrame(index=pheno_complete.index)
    audit: dict[str, dict] = {}
    for col in candidate_cols:
        enc = encode_binary(pheno_complete[col])
        code = phecode_of(col)
        frac, majority = observed_case_sex(enc, sex_by_id)
        from_table = phecode_sex.get(code, "")
        # A definitions label may only make the restriction STRONGER, never
        # weaker. "Both" used to short-circuit the observed check entirely, so a
        # phecode the table calls Both whose cases are in fact one sex kept sex
        # in the covariate set -- complete separation, and SAIGE step 1 diverges
        # while still writing a .rda (observed on CA_142.2 "Abnormal mammogram",
        # labelled Both, 360/360 cases female: theta 2.0e13, intercept -1.7e15,
        # every fitted value 2.22e-16, and every step-2 shard for it then died
        # with a gamma_q_inv overflow). Trust the data over the label.
        observed_specific = bool(majority) and frac >= args.sex_case_fraction
        if from_table in ("Male", "Female"):
            sex_specific, restrict_to, basis = True, from_table, "definitions"
        elif observed_specific:
            sex_specific, restrict_to = True, majority
            basis = "observed (definitions said Both)" if from_table == "Both" else "observed"
        else:
            sex_specific, restrict_to = False, ""
            basis = "definitions" if from_table else ""
        row = {
            "phenotype": col,
            "phecode": code,
            "parent_phecode": phecode_parent(code),
            "n_cases_raw": int((enc == 1.0).sum()),
            "n_controls_raw": int((enc == 0.0).sum()),
            "case_sex_fraction": None if pd.isna(frac) else round(frac, 4),
            "case_sex_majority": majority,
            "sex_specific": sex_specific,
            "sex_specific_basis": basis,
            "sex_restricted_to": "",
            "excluded_reason": "",
            "r2_with": "",
            "r2": None,
        }
        if sex_specific and args.sex_specific == "restrict":
            other = 0.0 if restrict_to == "Male" else 1.0
            enc = enc.mask(sex_by_id.reindex(enc.index) == other, pd.NA)
            row["sex_restricted_to"] = restrict_to
        encoded[col] = enc
        row["n_cases"] = int((enc == 1.0).sum())
        row["n_controls"] = int((enc == 0.0).sum())
        row["n_missing"] = int(enc.isna().sum())
        audit[col] = row

    n_sex_specific = sum(1 for r in audit.values() if r["sex_specific"])
    print(
        f"{n_sex_specific:,} of {len(audit):,} candidate phenotypes called sex-specific "
        f"({args.sex_specific})"
    )
    if phecode_sex:
        # A definitions table in a different dialect from the phenotype columns
        # joins on nothing, every call silently falls back to the observed case
        # sex ratio, and the output looks exactly like a successful run.
        n_matched = sum(1 for r in audit.values() if r["sex_specific_basis"] == "definitions")
        n_by_ratio = sum(
            1 for r in audit.values() if r["sex_specific"] and r["sex_specific_basis"] == "observed"
        )
        print(
            f"  {n_matched:,} of {len(audit):,} matched the definitions table "
            f"({len(phecode_sex):,} codes); {n_by_ratio:,} called from the case sex ratio"
        )
        if n_matched == 0:
            raise SystemExit(
                f"--phecode-sex-tsv matched none of the {len(audit):,} phenotype columns. "
                f"Example column {candidate_cols[0]!r} parses to "
                f"{phecode_of(candidate_cols[0])!r}; the table holds codes like "
                f"{sorted(phecode_sex)[:3]}. The two are in different dialects, so no "
                "phenotype would be restricted and the run would look fine."
            )

    passing = [c for c in candidate_cols if audit[c]["n_cases"] >= args.min_cases]
    for col in candidate_cols:
        if col not in passing:
            audit[col]["excluded_reason"] = f"n_cases<{args.min_cases}"

    if args.nesting != "keep":
        codes = {audit[c]["phecode"]: c for c in passing if audit[c]["phecode"]}
        drop = set()
        for code, col in codes.items():
            if args.nesting == "roots" and "." in code:
                drop.add(col)
            elif args.nesting == "leaves" and any(
                other != code and other.startswith(code + ".") for other in codes
            ):
                drop.add(col)
        for col in drop:
            audit[col]["excluded_reason"] = f"nesting={args.nesting}"
        passing = [c for c in passing if c not in drop]
        print(f"Nesting filter '{args.nesting}' dropped {len(drop):,} phecodes")

    # Rank by the smaller arm: that is what drives power for a binary trait.
    # Name is the last key so the order is total and reproducible across runs.
    passing.sort(
        key=lambda c: (
            min(audit[c]["n_cases"], audit[c]["n_controls"]),
            audit[c]["n_cases"],
            audit[c]["n_controls"],
            c,
        ),
        reverse=True,
    )

    if args.max_r2 is not None:
        dropped = prune_correlated(encoded, passing, args.max_r2)
        for col, (kept, r2) in dropped.items():
            audit[col]["excluded_reason"] = f"r2>={args.max_r2}"
            audit[col]["r2_with"] = kept
            audit[col]["r2"] = round(r2, 4)
        passing = [c for c in passing if c not in dropped]
        print(f"Correlation pruning at r2>={args.max_r2} dropped {len(dropped):,} phenotypes")

    selected = passing if args.n_phenotypes <= 0 else passing[: args.n_phenotypes]
    for col in passing[len(selected):]:
        audit[col]["excluded_reason"] = f"rank>{args.n_phenotypes}"
    if args.n_phenotypes > 0 and len(selected) < args.n_phenotypes:
        print(
            f"WARNING: only {len(selected)} phenotypes met the filters; "
            f"requested {args.n_phenotypes}"
        )
    if not selected:
        raise SystemExit("No phenotypes passed selection criteria")

    selected_stats = []
    for col in selected:
        enc = encoded[col]
        safe = col.replace(".", "_").replace("-", "_").replace(" ", "_")
        table[safe] = enc.to_numpy()
        row = audit[col]
        selected_stats.append(
            {
                "phenotype": col,
                "safe_name": safe,
                "n_cases": row["n_cases"],
                "n_controls": row["n_controls"],
                "n_missing": row["n_missing"],
                "phecode": row["phecode"],
                "parent_phecode": row["parent_phecode"],
                "sex_specific": row["sex_specific"],
                "sex_specific_basis": row["sex_specific_basis"],
                "sex_restricted_to": row["sex_restricted_to"],
                "case_sex_fraction": row["case_sex_fraction"],
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

    # Every candidate with the reason it is in or out, so a phenotype count in
    # the manuscript can be traced to a filter rather than re-derived.
    out_audit = args.out_dir / "phenotype_selection_audit.tsv"
    selected_set = {s["phenotype"] for s in selected_stats}
    audit_rows = []
    for col in candidate_cols:
        row = dict(audit[col])
        row["selected"] = col in selected_set
        audit_rows.append(row)
    audit_df = pd.DataFrame(audit_rows).sort_values(
        ["selected", "n_cases"], ascending=[False, False]
    )
    audit_df.to_csv(out_audit, sep="\t", index=False)

    tier_paths = []
    tier_sizes = [int(t) for t in args.tier_sizes.split(",") if t.strip()]
    for n in tier_sizes:
        if n <= 0:
            continue
        head = [s["safe_name"] for s in selected_stats[:n]]
        if len(head) < n:
            print(f"WARNING: tier {n} requested but only {len(head)} phenotypes selected")
        tier_path = args.out_dir / f"selected_phenotypes.top{n}.txt"
        Path(tier_path).write_text("\n".join(head) + "\n")
        tier_paths.append(tier_path)
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
        out_audit,
        *tier_paths,
    ]:
        print(f"  {path}")
    print(f"Limited covariates: {limited_cols}")
    print(f"Full covariates: {full_cols}")
    print(f"GC reference={gc_ref}")
    print(f"coverage standardized: mean={cov_mean:.4f}, sd={cov_sd:.4f}")
    print(f"Selected {len(selected_stats):,} phenotypes (phecodes, not independent traits):")
    for s in selected_stats[:20]:
        print(
            f"  {s['phenotype']} -> {s['safe_name']}: "
            f"cases={s['n_cases']}, controls={s['n_controls']}, missing={s['n_missing']}"
            + (f", sex-restricted to {s['sex_restricted_to']}" if s["sex_restricted_to"] else "")
        )
    if len(selected_stats) > 20:
        print(f"  ... {len(selected_stats) - 20:,} more; see {out_stats}")


if __name__ == "__main__":
    main()
