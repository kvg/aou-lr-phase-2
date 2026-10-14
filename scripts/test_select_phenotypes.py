#!/usr/bin/env python3
"""Tests for scripts/select_phenotypes.py — the parts that matter at full scale.

The pilot path (top-10 by balance, min 100 cases) is unchanged; these cover the
filters that only bite once the list is the whole phecode set: sex-specific
phecodes, parent/child nesting, and correlation pruning.
"""

from __future__ import annotations

import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from select_phenotypes import (  # noqa: E402
    encode_binary,
    observed_case_sex,
    phecode_of,
    phecode_parent,
    prune_correlated,
    read_phecode_sex,
)

SCRIPT = ROOT / "scripts" / "select_phenotypes.py"
N = 400


def test_phecode_parsing_and_hierarchy():
    # PhecodeX: the two-letter chapter is part of the code, and the sub-code is
    # written with a dot in the definitions table but an underscore in a column
    # name. CM is a chapter, not a prefix to strip.
    assert phecode_of("CM_761_13") == "CM_761.13"
    assert phecode_of("CM_761.13") == "CM_761.13"
    assert phecode_of("CM_761") == "CM_761"
    assert phecode_of("CV_401_1") == "CV_401.1"
    assert phecode_parent("CM_761.13") == "CM_761.1"
    assert phecode_parent("CM_761.1") == "CM_761"
    assert phecode_parent("CM_761") == ""

    # Classic phecodes are bare numbers, optionally behind PH_ / PHECODE_.
    assert phecode_of("PH_761_13") == "761.13"
    assert phecode_of("PHECODE_761") == "761"
    assert phecode_of("761.13") == "761.13"
    assert phecode_parent("761.13") == "761.1"
    assert phecode_parent("761") == ""

    # Covariate columns are not phecodes.
    for col in ("coverage", "sex", "age", "PC1", "GC_BCM"):
        assert phecode_of(col) == "", col


def test_observed_case_sex_calls_the_majority_sex():
    enc = pd.Series([1.0] * 10 + [0.0] * 10)
    sex = pd.Series([1.0] * 10 + [0.0] * 10)  # every case male
    frac, majority = observed_case_sex(enc, sex)
    assert (frac, majority) == (1.0, "Male")
    enc2 = pd.Series([1.0] * 10 + [0.0] * 10)
    sex2 = pd.Series([0.0] * 9 + [1.0] + [0.0] * 10)
    frac2, majority2 = observed_case_sex(enc2, sex2)
    assert majority2 == "Female" and round(frac2, 2) == 0.9


def test_prune_correlated_keeps_the_first_of_a_duplicated_pair():
    rng = np.random.default_rng(0)
    a = pd.Series(rng.integers(0, 2, 200).astype(float))
    df = pd.DataFrame({"a": a, "a_copy": a, "b": rng.integers(0, 2, 200).astype(float)})
    dropped = prune_correlated(df, ["a", "a_copy", "b"], 0.5)
    assert set(dropped) == {"a_copy"}
    assert dropped["a_copy"][0] == "a" and dropped["a_copy"][1] > 0.99


def _inputs(tmp: Path) -> tuple[Path, Path, Path]:
    """A cohort with a male-only phecode, a parent/child pair and a duplicate."""
    rng = np.random.default_rng(7)
    ids = [f"S{i:04d}" for i in range(N)]
    sex = np.array([1.0 if i % 2 else 0.0 for i in range(N)])  # 1 = male

    child = np.zeros(N)
    child[rng.choice(N, 120, replace=False)] = 1.0
    parent = child.copy()
    parent[rng.choice(np.flatnonzero(child == 0), 40, replace=False)] = 1.0

    male_only = np.zeros(N)
    male_idx = np.flatnonzero(sex == 1.0)
    male_only[rng.choice(male_idx, 110, replace=False)] = 1.0

    other = np.zeros(N)
    other[rng.choice(N, 150, replace=False)] = 1.0

    pheno = pd.DataFrame({
        "person_id": ids,
        "CM_761": parent,
        "CM_761_1": child,
        "CM_185": male_only,       # prostate-like: cases are male only
        "CM_250": other,
        "CM_250_2": other,         # exact duplicate of CM_250
    })
    cov = pd.DataFrame({
        "research_id": ids,
        "sex_at_birth": np.where(sex == 1.0, "Male", "Female"),
        "age": rng.normal(55, 10, N).round(1),
        "coverage": rng.normal(30, 3, N).round(2),
        "GC": rng.choice(["bi", "uw", "bcm"], N),
    })
    for i in range(1, 11):
        cov[f"PC{i}"] = rng.normal(size=N).round(4)
    pheno_path, cov_path = tmp / "pheno.csv", tmp / "cov.csv"
    samples_path = tmp / "vcf_samples.txt"
    pheno.to_csv(pheno_path, index=False)
    cov.to_csv(cov_path, index=False)
    samples_path.write_text("\n".join(ids) + "\n")
    return pheno_path, cov_path, samples_path


def _run(tmp: Path, *extra: str) -> tuple[pd.DataFrame, pd.DataFrame, list[str]]:
    pheno_path, cov_path, samples_path = _inputs(tmp)
    out = tmp / "out"
    out.mkdir(exist_ok=True)
    cmd = [
        sys.executable, str(SCRIPT),
        "--phenotype-csv", str(pheno_path),
        "--covariates-csv", str(cov_path),
        "--vcf-samples", str(samples_path),
        "--out-dir", str(out),
        "--min-cases", "50",
        "--n-phenotypes", "0",
        *extra,
    ]
    r = subprocess.run(cmd, capture_output=True, text=True)
    assert r.returncode == 0, r.stdout + r.stderr
    stats = pd.read_csv(out / "selected_phenotype_stats.tsv", sep="\t")
    audit = pd.read_csv(out / "phenotype_selection_audit.tsv", sep="\t")
    chosen = (out / "selected_phenotypes.txt").read_text().split()
    return stats, audit, chosen


def test_sex_specific_is_flagged_by_default_and_not_restricted():
    with tempfile.TemporaryDirectory() as td:
        stats, _audit, chosen = _run(Path(td))
    row = stats.set_index("phenotype").loc["CM_185"]
    assert bool(row["sex_specific"]) and row["sex_specific_basis"] == "observed"
    assert pd.isna(row["sex_restricted_to"]) or row["sex_restricted_to"] == ""
    # flag mode leaves the phenotype alone: females are still controls.
    assert row["n_controls"] > N / 2
    assert "CM_185" in chosen


def test_sex_specific_restrict_drops_the_other_sex_to_missing():
    with tempfile.TemporaryDirectory() as td:
        stats, _audit, _chosen = _run(Path(td), "--sex-specific", "restrict")
    row = stats.set_index("phenotype").loc["CM_185"]
    assert row["sex_restricted_to"] == "Male"
    # Controls are now male non-cases only, and the females became missing.
    assert row["n_cases"] == 110
    assert row["n_controls"] + row["n_cases"] == N / 2
    assert row["n_missing"] == N / 2
    other = stats.set_index("phenotype").loc["CM_250"]
    assert other["n_missing"] == 0, "a non-sex-specific phecode is untouched"


def test_definitions_table_cannot_weaken_an_observed_sex_restriction():
    """A "Both" label may not override a case set that is entirely one sex.

    This reverses an earlier decision (the table used to win outright). It had to
    change: PhecodeX labels CA_142.2 "Abnormal mammogram" as Both, yet all 360
    cases in this cohort are female, so sex stayed in the covariate set, perfectly
    separated the cases, and SAIGE step 1 diverged while still writing a .rda
    (theta 2.0e13, intercept -1.7e15, every fitted value 2.22e-16). Every step-2
    shard for that phenotype then aborted on a gamma_q_inv overflow. A label can
    still make the restriction stronger, never weaker.
    """
    with tempfile.TemporaryDirectory() as td:
        tmp = Path(td)
        defs = tmp / "phecode_defs.csv"
        defs.write_text("phecode,sex\nCM_185,Both\nCM_761,Both\nCM_761.1,Both\nCM_250,Both\n")
        stats, _audit, _chosen = _run(tmp, "--phecode-sex-tsv", str(defs),
                                      "--sex-specific", "restrict")
    stats = stats.set_index("phenotype")

    # CM_185's cases are 110/110 male, so the data wins over the Both label.
    row = stats.loc["CM_185"]
    assert bool(row["sex_specific"])
    assert row["sex_specific_basis"] == "observed (definitions said Both)"
    assert row["sex_restricted_to"] == "Male"
    assert row["n_missing"] == N / 2

    # A genuinely both-sex phecode is still left alone, which is what the label
    # is for: it is consulted, and it does not restrict.
    other = stats.loc["CM_250"]
    assert not bool(other["sex_specific"])
    assert other["sex_specific_basis"] == "definitions"
    assert other["n_missing"] == 0


def test_definitions_table_in_the_wrong_dialect_is_refused(tmp_path):
    """A classic-phecode table against PhecodeX columns joins on nothing.

    Left to itself every call falls back to the observed case sex ratio and the
    run looks successful, so this has to fail loudly instead.
    """
    defs = tmp_path / "phecode_defs.csv"
    defs.write_text("phecode,sex\n185,Male\n761,Both\n")
    # _run asserts a zero exit, so the refusal surfaces as the captured output.
    with pytest.raises(AssertionError, match="matched none"):
        _run(tmp_path, "--phecode-sex-tsv", str(defs), "--sex-specific", "restrict")


def test_definitions_table_reads_latin1(tmp_path):
    """The published PhecodeX table is Latin-1, not UTF-8."""
    defs = tmp_path / "phecode_defs.csv"
    defs.write_bytes(
        "phecode,phecode_string,sex\nCM_185,M\u00e9ni\u00e8re,Both\nCM_761,x,Both\n"
        .encode("latin-1")
    )
    mapping = read_phecode_sex(defs)
    assert mapping["CM_185"] == "Both"


def test_nesting_leaves_drops_parents_and_audit_records_why():
    with tempfile.TemporaryDirectory() as td:
        _stats, audit, chosen = _run(Path(td), "--nesting", "leaves")
    assert "CM_761_1" in chosen and "CM_761" not in chosen
    assert audit.set_index("phenotype").loc["CM_761", "excluded_reason"] == "nesting=leaves"


def test_nesting_roots_keeps_only_integer_codes():
    with tempfile.TemporaryDirectory() as td:
        _stats, _audit, chosen = _run(Path(td), "--nesting", "roots")
    assert "CM_761" in chosen and "CM_761_1" not in chosen


def test_r2_pruning_removes_the_duplicate_and_tiers_are_prefixes():
    with tempfile.TemporaryDirectory() as td:
        tmp = Path(td)
        stats, audit, chosen = _run(tmp, "--max-r2", "0.5", "--tier-sizes", "2")
        tier = (tmp / "out" / "selected_phenotypes.top2.txt").read_text().split()
    kept = {c for c in chosen if c.startswith("CM_250")}
    assert len(kept) == 1, f"one of the duplicated pair should survive, got {kept}"
    dropped = audit[audit["excluded_reason"].astype(str).str.startswith("r2>=")].set_index(
        "phenotype"
    )
    # Which member of an exact-duplicate pair survives is set by the ranking
    # tie-break, so assert on the pair rather than on a name.
    pair = {"CM_250", "CM_250_2"}
    gone = (pair - kept).pop()
    assert dropped.loc[gone, "r2"] > 0.99
    assert dropped.loc[gone, "r2_with"] == (pair - {gone}).pop()
    # The parent/child pair is correlated too (child is a subset of parent), so
    # pruning removes one of those as well — that is the point of the filter.
    assert "CM_761" in dropped.index or "CM_761_1" in dropped.index
    assert tier == chosen[:2]
    assert list(stats["safe_name"]) == chosen


def test_audit_covers_every_candidate():
    with tempfile.TemporaryDirectory() as td:
        _stats, audit, chosen = _run(Path(td), "--min-cases", "130")
    assert len(audit) == 5, "one row per candidate phenotype column"
    assert set(audit[audit["selected"]]["phenotype"]) == set(chosen)
    excluded = audit[~audit["selected"]]
    assert (excluded["excluded_reason"] == "n_cases<130").all()


if __name__ == "__main__":
    for name, fn in sorted(globals().items()):
        if name.startswith("test_") and callable(fn):
            fn()
            print(f"  {name} ok")
    print("ok")
