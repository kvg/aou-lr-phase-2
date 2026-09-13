#!/usr/bin/env python3
"""Summarize pb-CpG-tools bedMethyl pileups for the Phase 2 methylation section.

Subcommands:
  summarize     one sample (combined ± hap1/hap2 beds) → JSON + one-row TSV
  extract-chrom one bed, one chromosome (optional start/end), coverage filter → site TSV
  merge         concatenate per-sample TSVs, join covariates, write paper numbers
  concordance   Primrose vs Jasmine site-mean comparison from chrom-site dumps
  plot-s3       hexbin of those site means + haplotype occupancy vs coverage
  inventory     fetch aou2_v1_phased_bams (firecloud) and count bed / stats URIs
  pull-stats    copy per-sample stats_tsv files from that table into a directory
  pull-dumps    copy PbCpgChromSites `sites` files for Primrose/Jasmine concordance
  pick-concordance-samples  25 Primrose + 25 Jasmine IDs for that WDL

Bed format (pb-CpG-tools v3 model mode): chrom, begin, end, mod_score, type, cov, ...
See https://github.com/PacificBiosciences/pb-CpG-tools
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import math
import os
import re
import subprocess
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from io import StringIO
from pathlib import Path
from typing import Any, Iterable, Iterator, Mapping, Sequence, TextIO

DROP_CHROMS = frozenset({"chrM", "chrMT", "MT"})
COV_THRESHOLDS = (1, 5, 10, 15, 20)
DEFAULT_ENTITY_TYPE = "aou2_v1_phased_bams"
DEFAULT_NAMESPACE = "allofus-drc-wgs-LR-prodData"
DEFAULT_WORKSPACE = "AoU_DRC_LongReads_PhaseTwo_Storage"
DEFAULT_ID_COLUMN = "entity:aou2_v1_phased_bams_id"
CALLER_PRIMROSE_PREFIX = "primrose"
CALLER_JASMINE_PREFIX = "jasmine"
CONCORDANCE_DELTA = 10.0
CONCORDANCE_N_PER_GROUP = 25
CONCORDANCE_ENTITY_TYPE = "pbcpg_concordance"
CONCORDANCE_ID_COLUMN = "entity:pbcpg_concordance_id"
STATS_URI_CANDIDATES = (
    "stats_tsv",
    "PbCpgSampleStats.stats_tsv",
    "pbcpg_stats_tsv",
)
SITES_URI_CANDIDATES = (
    "sites",
    "PbCpgChromSites.sites",
    "chr22_sites",
    "chrom_sites",
    "pbcpg_sites",
)


def entity_id_column(entity_type: str) -> str:
    return f"entity:{entity_type}_id"


def open_text(path: Path) -> TextIO:
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return path.open("rt", encoding="utf-8", errors="replace")


def parse_header(handle: TextIO) -> tuple[dict[str, str], list[str] | None, str | None]:
    """Read '#' lines. Return metadata, optional column names, first data line."""
    meta: dict[str, str] = {}
    columns: list[str] | None = None
    leftover: str | None = None
    for raw in handle:
        if raw.startswith("##"):
            body = raw[2:].strip()
            if "=" in body:
                key, value = body.split("=", 1)
                meta[key.strip()] = value.strip().strip('"')
            continue
        if raw.startswith("#"):
            body = raw[1:].strip()
            columns = re.split(r"\s+", body)
            continue
        leftover = raw
        break
    return meta, columns, leftover


def resolve_columns(columns: list[str] | None) -> dict[str, int]:
    """Map column name → index. Fall back to pb-CpG-tools v3 model layout."""
    if columns:
        lowered = {name.lower(): i for i, name in enumerate(columns)}
        chrom = lowered.get("chrom", lowered.get("#chrom", 0))
        begin = lowered.get("begin", lowered.get("start", 1))
        cov = lowered.get("cov", lowered.get("coverage", 5))
        score = lowered.get("mod_score", lowered.get("score", 3))
        return {"chrom": chrom, "begin": begin, "cov": cov, "score": score}
    return {"chrom": 0, "begin": 1, "cov": 5, "score": 3}


def iter_bed_rows(
    path: Path,
    *,
    chrom_keep: set[str] | None = None,
    begin_min: int | None = None,
    begin_max: int | None = None,
) -> tuple[dict[str, str], Iterator[tuple[str, int, int, float]]]:
    handle = open_text(path)
    meta, columns, leftover = parse_header(handle)
    idx = resolve_columns(columns)
    max_i = max(idx.values())

    def rows() -> Iterator[tuple[str, int, int, float]]:
        try:
            lines: Iterable[str]
            if leftover is None:
                lines = handle
            else:
                def _chain() -> Iterator[str]:
                    yield leftover
                    yield from handle

                lines = _chain()
            for raw in lines:
                if not raw or raw.startswith("#"):
                    continue
                parts = raw.rstrip("\n").split("\t")
                if len(parts) == 1:
                    parts = re.split(r"\s+", raw.strip())
                if len(parts) <= max_i:
                    continue
                chrom = parts[idx["chrom"]].strip()
                if chrom_keep is not None and chrom not in chrom_keep:
                    continue
                if chrom in DROP_CHROMS:
                    continue
                try:
                    begin = int(float(parts[idx["begin"]]))
                    cov = int(float(parts[idx["cov"]]))
                    score = float(parts[idx["score"]])
                except ValueError:
                    continue
                if begin_min is not None and begin < begin_min:
                    continue
                if begin_max is not None and begin >= begin_max:
                    continue
                yield chrom, begin, cov, score
        finally:
            handle.close()

    return meta, rows()


def empty_bed_stats() -> dict[str, Any]:
    out: dict[str, Any] = {
        "n_sites": 0,
        "sum_cov": 0,
        "sum_score": 0.0,
        "n_autosome": 0,
        "n_chrX": 0,
        "n_chrY": 0,
        "mean_cov": 0.0,
        "mean_score": 0.0,
    }
    for thresh in COV_THRESHOLDS:
        out[f"n_cov_ge{thresh}"] = 0
    return out


def summarize_bed(path: Path | None, *, chrom_keep: set[str] | None = None) -> dict[str, Any]:
    stats = empty_bed_stats()
    stats["path"] = "" if path is None else str(path)
    stats["present"] = bool(path) and Path(path).is_file()
    if not stats["present"]:
        return stats
    meta, rows = iter_bed_rows(Path(path), chrom_keep=chrom_keep)
    stats["pb_cpg_version"] = meta.get("pb-cpg-tools-version", "")
    stats["pileup_mode"] = meta.get("pileup-mode", "")
    stats["modsites_mode"] = meta.get("modsites-mode", "")
    stats["min_coverage_header"] = meta.get("min-coverage", "")
    stats["basemod_source"] = meta.get("basemod-source", "")
    for chrom, _begin, cov, score in rows:
        stats["n_sites"] += 1
        stats["sum_cov"] += cov
        stats["sum_score"] += score
        if chrom.startswith("chr") and chrom[3:].isdigit() and 1 <= int(chrom[3:]) <= 22:
            stats["n_autosome"] += 1
        elif chrom == "chrX":
            stats["n_chrX"] += 1
        elif chrom == "chrY":
            stats["n_chrY"] += 1
        for thresh in COV_THRESHOLDS:
            if cov >= thresh:
                stats[f"n_cov_ge{thresh}"] += 1
    n = stats["n_sites"]
    stats["mean_cov"] = (stats["sum_cov"] / n) if n else 0.0
    stats["mean_score"] = (stats["sum_score"] / n) if n else 0.0
    return stats


def prefix_stats(prefix: str, stats: Mapping[str, Any]) -> dict[str, Any]:
    skip = {
        "path",
        "present",
        "pb_cpg_version",
        "pileup_mode",
        "modsites_mode",
        "min_coverage_header",
        "basemod_source",
        "sum_cov",
        "sum_score",
    }
    out: dict[str, Any] = {}
    for key, value in stats.items():
        if key in skip:
            continue
        out[f"{prefix}_{key}"] = value
    out[f"{prefix}_present"] = bool(stats.get("present"))
    out[f"{prefix}_path"] = str(stats.get("path") or "")
    return out


def combine_sample_stats(
    sample_id: str,
    combined: Mapping[str, Any],
    hap1: Mapping[str, Any],
    hap2: Mapping[str, Any],
) -> dict[str, Any]:
    n_combined = int(combined.get("n_sites") or 0)
    n_hap1 = int(hap1.get("n_sites") or 0)
    n_hap2 = int(hap2.get("n_sites") or 0)
    n_combined_ge5 = int(combined.get("n_cov_ge5") or 0)
    n_hap1_ge5 = int(hap1.get("n_cov_ge5") or 0)
    n_hap2_ge5 = int(hap2.get("n_cov_ge5") or 0)
    row = {
        "sample_id": sample_id,
        "pb_cpg_version": combined.get("pb_cpg_version") or hap1.get("pb_cpg_version") or "",
        "pileup_mode": combined.get("pileup_mode") or "",
        "modsites_mode": combined.get("modsites_mode") or "",
        "min_coverage_header": combined.get("min_coverage_header") or "",
        "basemod_source": combined.get("basemod_source") or "",
        "haplotype_tracks": bool(hap1.get("present") and hap2.get("present")),
        "haplotype_resolved": bool(n_hap1_ge5 > 0 and n_hap2_ge5 > 0),
        "frac_hap_slots": ((n_hap1 + n_hap2) / (2 * n_combined)) if n_combined else 0.0,
        "frac_hap_slots_ge5": (
            (n_hap1_ge5 + n_hap2_ge5) / (2 * n_combined_ge5) if n_combined_ge5 else 0.0
        ),
    }
    row.update(prefix_stats("combined", combined))
    row.update(prefix_stats("hap1", hap1))
    row.update(prefix_stats("hap2", hap2))
    return row


def write_json(path: Path, payload: Mapping[str, Any]) -> None:
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")


def write_row_tsv(path: Path, row: Mapping[str, Any]) -> None:
    keys = list(row.keys())
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=keys, delimiter="\t")
        writer.writeheader()
        writer.writerow({k: row[k] for k in keys})


def cmd_summarize(args: argparse.Namespace) -> None:
    sample_id = args.sample_id
    combined = summarize_bed(Path(args.combined) if args.combined else None)
    hap1 = summarize_bed(Path(args.hap1) if args.hap1 else None)
    hap2 = summarize_bed(Path(args.hap2) if args.hap2 else None)
    row = combine_sample_stats(sample_id, combined, hap1, hap2)
    out_prefix = Path(args.out_prefix)
    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    write_json(Path(str(out_prefix) + ".json"), row)
    write_row_tsv(Path(str(out_prefix) + ".tsv"), row)
    print(json.dumps({"sample_id": sample_id, "n_sites": row["combined_n_sites"]}, sort_keys=True))


def cmd_extract_chrom(args: argparse.Namespace) -> None:
    chrom = args.chrom.strip()
    min_cov = int(args.min_cov)
    begin_min = int(args.start) if args.start is not None else None
    begin_max = int(args.end) if args.end is not None else None
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    opener = gzip.open if str(out).endswith(".gz") else open
    n_in = 0
    n_out = 0
    sample_id = str(args.sample_id or "")
    with opener(out, "wt", encoding="utf-8") as handle:
        if sample_id:
            handle.write(f"##sample_id={sample_id}\n")
        handle.write("chrom\tbegin\tcov\tmod_score\n")
        _meta, rows = iter_bed_rows(
            Path(args.bed),
            chrom_keep={chrom},
            begin_min=begin_min,
            begin_max=begin_max,
        )
        for chrom_i, begin, cov, score in rows:
            n_in += 1
            if cov < min_cov:
                continue
            handle.write(f"{chrom_i}\t{begin}\t{cov}\t{score:.6g}\n")
            n_out += 1
    print(
        json.dumps(
            {
                "chrom": chrom,
                "start": begin_min,
                "end": begin_max,
                "n_in": n_in,
                "n_out": n_out,
                "out": str(out),
            }
        )
    )


def _read_tables(paths: Sequence[Path]) -> "Any":
    import pandas as pd

    frames = []
    for path in paths:
        frames.append(pd.read_csv(path, sep="\t", dtype={"sample_id": str}))
    if not frames:
        raise FileNotFoundError("No per-sample stats TSVs to merge")
    return pd.concat(frames, ignore_index=True)


def load_stats_dir(stats_dir: Path) -> "Any":
    skip = {
        "pbcpg_sample_stats.tsv",
        "methylation_manuscript_numbers.tsv",
        "labels.tsv",
    }
    paths = [
        path
        for path in sorted(stats_dir.rglob("*.tsv"))
        if path.name not in skip and not path.name.endswith(".table.tsv")
    ]
    return _read_tables(paths)


def _median(series: "Any") -> float:
    import pandas as pd

    s = pd.to_numeric(series, errors="coerce").dropna()
    if s.empty:
        return float("nan")
    return float(s.median())


def _q1(series: "Any") -> float:
    import pandas as pd

    s = pd.to_numeric(series, errors="coerce").dropna()
    if s.empty:
        return float("nan")
    return float(s.quantile(0.25))


def _q3(series: "Any") -> float:
    import pandas as pd

    s = pd.to_numeric(series, errors="coerce").dropna()
    if s.empty:
        return float("nan")
    return float(s.quantile(0.75))


def manuscript_numbers(df: "Any") -> dict[str, Any]:
    n = len(df)
    n_hap_files = int(df["haplotype_tracks"].fillna(False).astype(bool).sum()) if "haplotype_tracks" in df else 0
    n_hap_sites = int(df["haplotype_resolved"].fillna(False).astype(bool).sum()) if "haplotype_resolved" in df else 0
    numbers = {
        "n_samples": n,
        "n_with_combined": int(df["combined_present"].fillna(False).astype(bool).sum())
        if "combined_present" in df
        else n,
        "n_haplotype_tracks": n_hap_files,
        "n_haplotype_resolved": n_hap_sites,
        "median_n_cpg_ge5": _median(df["combined_n_cov_ge5"]),
        "q1_n_cpg_ge5": _q1(df["combined_n_cov_ge5"]),
        "q3_n_cpg_ge5": _q3(df["combined_n_cov_ge5"]),
        "median_n_cpg_ge10": _median(df["combined_n_cov_ge10"]),
        "median_mean_cov": _median(df["combined_mean_cov"]),
        "median_frac_hap_ge5": _median(df["frac_hap_slots_ge5"]),
    }
    if "pb_meth_caller" in df.columns:
        caller = df["pb_meth_caller"].fillna("NA").astype(str)
        numbers["n_primrose"] = int(caller.str.startswith(CALLER_PRIMROSE_PREFIX).sum())
        numbers["n_jasmine"] = int(caller.str.startswith(CALLER_JASMINE_PREFIX).sum())
        numbers["n_mixed_caller"] = int(caller.eq("mixed").sum())
    if "coverage" in df.columns:
        numbers["median_hifi_coverage"] = _median(df["coverage"])
    return numbers


def join_covariates(stats: "Any", covariates_path: Path | None) -> "Any":
    if covariates_path is None or not covariates_path.is_file():
        return stats
    import pandas as pd

    cov = pd.read_csv(covariates_path, dtype={"research_id": str}, low_memory=False)
    if "research_id" not in cov.columns:
        raise ValueError(f"covariates missing research_id: {covariates_path}")
    keep = [
        c
        for c in (
            "research_id",
            "coverage",
            "pb_meth_caller",
            "platform",
            "final_releasable_v9",
            "technology",
            "has_ONT",
            "has_PacBio",
            "lr_phase",
        )
        if c in cov.columns
    ]
    merged = stats.merge(cov[keep], how="left", left_on="sample_id", right_on="research_id")
    return merged


def cmd_merge(args: argparse.Namespace) -> None:
    import pandas as pd

    stats_dir = Path(args.stats_dir)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    df = load_stats_dir(stats_dir)
    cov_path = Path(args.covariates) if args.covariates else None
    df = join_covariates(df, cov_path)
    if args.discovery_only and "final_releasable_v9" in df.columns:
        mask = (
            df["final_releasable_v9"]
            .fillna("false")
            .astype(str)
            .str.strip()
            .str.lower()
            .isin({"true", "1", "yes", "t"})
        )
        if "technology" in df.columns:
            mask &= df["technology"].astype(str).eq("PacBio")
        df = df.loc[mask].copy()
    numbers = manuscript_numbers(df)
    sample_path = out_dir / "pbcpg_sample_stats.tsv"
    numbers_path = out_dir / "methylation_manuscript_numbers.tsv"
    json_path = out_dir / "methylation_manuscript_numbers.json"
    df.to_csv(sample_path, sep="\t", index=False)
    pd.DataFrame([numbers]).to_csv(numbers_path, sep="\t", index=False)
    write_json(json_path, numbers)
    print(json.dumps(numbers, sort_keys=True, default=str))


def pearson_r(xs: Sequence[float], ys: Sequence[float]) -> float:
    n = len(xs)
    if n < 3:
        return float("nan")
    mean_x = sum(xs) / n
    mean_y = sum(ys) / n
    var_x = sum((x - mean_x) ** 2 for x in xs)
    var_y = sum((y - mean_y) ** 2 for y in ys)
    if var_x <= 0 or var_y <= 0:
        return float("nan")
    cov = sum((x - mean_x) * (y - mean_y) for x, y in zip(xs, ys))
    return cov / math.sqrt(var_x * var_y)


def load_site_dump(path: Path) -> tuple[str, dict[tuple[str, int], tuple[int, float]]]:
    sites: dict[tuple[str, int], tuple[int, float]] = {}
    sample_id = path.name.split(".")[0]
    with open_text(path) as handle:
        for raw in handle:
            if raw.startswith("##"):
                body = raw[2:].strip()
                if body.startswith("sample_id="):
                    sample_id = body.split("=", 1)[1].strip()
                continue
            if raw.startswith("chrom") or raw.startswith("#"):
                continue
            parts = raw.rstrip("\n").split("\t")
            if len(parts) < 4:
                continue
            try:
                begin = int(parts[1])
                cov = int(float(parts[2]))
                score = float(parts[3])
            except ValueError:
                continue
            sites[(parts[0], begin)] = (cov, score)
    return sample_id, sites


def paired_concordance(bed_a: Path, bed_b: Path, *, min_cov: int) -> dict[str, Any]:
    meta_a, rows_a = iter_bed_rows(bed_a)
    lookup = {(chrom, begin): (cov, score) for chrom, begin, cov, score in rows_a}
    n_both = 0
    n_conc = 0
    xs: list[float] = []
    ys: list[float] = []
    _meta_b, rows_b = iter_bed_rows(bed_b)
    for chrom, begin, cov_b, score_b in rows_b:
        hit = lookup.get((chrom, begin))
        if hit is None:
            continue
        cov_a, score_a = hit
        if cov_a < min_cov or cov_b < min_cov:
            continue
        n_both += 1
        xs.append(score_a)
        ys.append(score_b)
        if abs(score_a - score_b) < CONCORDANCE_DELTA:
            n_conc += 1
    frac = (n_conc / n_both) if n_both else float("nan")
    return {
        "n_sites": n_both,
        "concordance_frac": frac,
        "concordance_r": pearson_r(xs, ys),
        "min_cov": min_cov,
        "delta": CONCORDANCE_DELTA,
        "basemod_source_a": meta_a.get("basemod-source", ""),
    }


def group_mean_site_pairs(
    dumps: Mapping[str, Path],
    labels: Mapping[str, str],
    *,
    min_cov: int,
    min_n_per_group: int,
) -> tuple[list[float], list[float], dict[str, int]]:
    primrose: dict[tuple[str, int], list[float]] = defaultdict(list)
    jasmine: dict[tuple[str, int], list[float]] = defaultdict(list)
    n_primrose = 0
    n_jasmine = 0
    for sample_id, path in dumps.items():
        caller = labels.get(sample_id, "")
        bucket = None
        if caller.startswith(CALLER_PRIMROSE_PREFIX):
            bucket = primrose
            n_primrose += 1
        elif caller.startswith(CALLER_JASMINE_PREFIX):
            bucket = jasmine
            n_jasmine += 1
        if bucket is None:
            continue
        for key, (cov, score) in load_site_dump(path)[1].items():
            if cov >= min_cov:
                bucket[key].append(score)
    xs: list[float] = []
    ys: list[float] = []
    keys = set(primrose) & set(jasmine)
    for key in keys:
        p_vals = primrose[key]
        j_vals = jasmine[key]
        if len(p_vals) < min_n_per_group or len(j_vals) < min_n_per_group:
            continue
        xs.append(sum(p_vals) / len(p_vals))
        ys.append(sum(j_vals) / len(j_vals))
    counts = {
        "n_primrose_dumps": n_primrose,
        "n_jasmine_dumps": n_jasmine,
        "min_cov": min_cov,
        "min_n_per_group": min_n_per_group,
    }
    return xs, ys, counts


def concordance_from_pairs(
    xs: Sequence[float],
    ys: Sequence[float],
    *,
    extra: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    n_sites = len(xs)
    n_conc = sum(1 for x, y in zip(xs, ys) if abs(x - y) < CONCORDANCE_DELTA)
    frac = (n_conc / n_sites) if n_sites else float("nan")
    result: dict[str, Any] = {
        "n_sites": n_sites,
        "concordance_frac": frac,
        "concordance_r": pearson_r(xs, ys),
        "delta": CONCORDANCE_DELTA,
    }
    if extra:
        result.update(extra)
    return result


def group_mean_concordance(
    dumps: Mapping[str, Path],
    labels: Mapping[str, str],
    *,
    min_cov: int,
    min_n_per_group: int,
) -> dict[str, Any]:
    xs, ys, counts = group_mean_site_pairs(
        dumps, labels, min_cov=min_cov, min_n_per_group=min_n_per_group
    )
    return concordance_from_pairs(xs, ys, extra=counts)


def write_pairs_tsv(path: Path, xs: Sequence[float], ys: Sequence[float]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "wt", encoding="utf-8") as handle:
        handle.write("mean_primrose\tmean_jasmine\n")
        for x, y in zip(xs, ys):
            handle.write(f"{x:.6g}\t{y:.6g}\n")


def load_pairs_tsv(path: Path) -> tuple[list[float], list[float]]:
    xs: list[float] = []
    ys: list[float] = []
    with open_text(path) as handle:
        header = handle.readline()
        if not header:
            return xs, ys
        for raw in handle:
            parts = raw.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            try:
                xs.append(float(parts[0]))
                ys.append(float(parts[1]))
            except ValueError:
                continue
    return xs, ys


def load_concordance_labels(path: Path) -> dict[str, str]:
    labels: dict[str, str] = {}
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            sample = str(row.get("sample_id") or row.get("research_id") or "").strip()
            caller = str(row.get("pb_meth_caller") or row.get("caller") or "").strip()
            if sample:
                labels[sample] = caller
    return labels


def load_dump_paths(dump_dir: Path) -> dict[str, Path]:
    dumps: dict[str, Path] = {}
    for path in sorted(dump_dir.glob("*")):
        if path.suffix not in {".tsv", ".gz"} and not path.name.endswith(".tsv.gz"):
            continue
        if not path.is_file():
            continue
        sample_id, _sites = load_site_dump(path)
        dumps[sample_id] = path
    return dumps


def cmd_concordance(args: argparse.Namespace) -> None:
    xs: list[float] = []
    ys: list[float] = []
    if args.bed_a and args.bed_b:
        result = paired_concordance(Path(args.bed_a), Path(args.bed_b), min_cov=args.min_cov)
    else:
        labels = load_concordance_labels(Path(args.labels))
        dumps = load_dump_paths(Path(args.dumps_dir))
        xs, ys, counts = group_mean_site_pairs(
            dumps, labels, min_cov=args.min_cov, min_n_per_group=args.min_n_per_group
        )
        result = concordance_from_pairs(xs, ys, extra=counts)
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    write_json(out, result)
    pairs_out = Path(args.pairs_out) if args.pairs_out else None
    if pairs_out is None and xs:
        pairs_out = out.with_suffix(".pairs.tsv.gz")
    if pairs_out is not None and xs:
        write_pairs_tsv(pairs_out, xs, ys)
        result["pairs_out"] = str(pairs_out)
    print(json.dumps(result, sort_keys=True))


def _plot_concordance_hexbin(ax, xs: Sequence[float], ys: Sequence[float], stats: Mapping[str, Any]) -> None:
    import numpy as np

    x = np.asarray(xs, dtype=float)
    y = np.asarray(ys, dtype=float)
    vmax = float(max(np.nanmax(x), np.nanmax(y), 1.0)) if x.size else 100.0
    limit = 100.0 if vmax > 1.5 else 1.0
    ax.hexbin(
        x,
        y,
        gridsize=70,
        bins="log",
        cmap="Blues",
        mincnt=1,
        extent=(0, limit, 0, limit),
        linewidths=0,
    )
    ax.plot([0, limit], [0, limit], ls="--", color="#888888", lw=0.8, zorder=3)
    ax.set_xlim(0, limit)
    ax.set_ylim(0, limit)
    ax.set_aspect("equal", adjustable="box")
    unit = "%" if limit > 1.5 else ""
    ax.set_xlabel(f"Primrose site-mean 5mC{unit}")
    ax.set_ylabel(f"Jasmine site-mean 5mC{unit}")
    ax.set_title("Caller concordance (site means, not paired samples)", loc="left", pad=2)
    n = int(stats.get("n_sites") or len(xs))
    frac = stats.get("concordance_frac")
    r = stats.get("concordance_r")
    frac_s = f"{100 * float(frac):.1f}%" if frac is not None and frac == frac else "NA"
    r_s = f"{float(r):.3f}" if r is not None and r == r else "NA"
    ax.text(
        0.04 * limit,
        0.92 * limit,
        f"$r = {r_s}$  ·  {frac_s} with $|\\Delta| < 10$ pp  ·  {n:,} sites",
        fontsize=8,
        color="#333333",
    )


def _plot_haplotype_coverage(ax, stats_path: Path) -> dict[str, Any]:
    import numpy as np
    import pandas as pd

    df = pd.read_csv(stats_path, sep="\t", low_memory=False)
    x_col = "coverage" if "coverage" in df.columns else "combined_mean_cov"
    y_col = "frac_hap_slots_ge5"
    if y_col not in df.columns or x_col not in df.columns:
        raise ValueError(f"{stats_path} needs {x_col} and {y_col}")
    x = pd.to_numeric(df[x_col], errors="coerce")
    y = pd.to_numeric(df[y_col], errors="coerce")
    mask = x.notna() & y.notna()
    x = x.loc[mask].to_numpy(dtype=float)
    y = y.loc[mask].to_numpy(dtype=float)
    ax.hexbin(x, y, gridsize=50, bins="log", cmap="Blues", mincnt=1, linewidths=0)
    median_y = float(np.median(y)) if y.size else float("nan")
    ax.axhline(median_y, color="#888888", ls=":", lw=0.8)
    ax.set_xlabel("HiFi coverage (×)" if x_col == "coverage" else "Mean CpG coverage (×)")
    ax.set_ylabel("Haplotype-slot fraction ≥5×")
    ax.set_title("Haplotype-assigned 5mC versus coverage", loc="left", pad=2)
    ax.set_ylim(0, 1.02)
    ax.text(
        0.04,
        0.08,
        f"median {100 * median_y:.1f}%  ·  n = {int(y.size):,}",
        transform=ax.transAxes,
        fontsize=8,
        color="#333333",
    )
    return {"n": int(y.size), "median_frac_hap_ge5": median_y}


def cmd_plot_s3(args: argparse.Namespace) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    xs: list[float] = []
    ys: list[float] = []
    conc: dict[str, Any] = {}
    if args.pairs:
        xs, ys = load_pairs_tsv(Path(args.pairs))
        conc = concordance_from_pairs(xs, ys)
    elif args.dumps_dir and args.labels:
        labels = load_concordance_labels(Path(args.labels))
        dumps = load_dump_paths(Path(args.dumps_dir))
        xs, ys, counts = group_mean_site_pairs(
            dumps,
            labels,
            min_cov=args.min_cov,
            min_n_per_group=args.min_n_per_group,
        )
        conc = concordance_from_pairs(xs, ys, extra=counts)
        if args.write_pairs:
            write_pairs_tsv(Path(args.write_pairs), xs, ys)
    else:
        raise SystemExit("plot-s3 requires --pairs or --dumps-dir/--labels")

    out_prefix = Path(args.out_prefix)
    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(1, 2, figsize=(10.2, 4.4))
    _plot_concordance_hexbin(axes[0], xs, ys, conc)
    hap: dict[str, Any] = {}
    if args.stats:
        hap = _plot_haplotype_coverage(axes[1], Path(args.stats))
    else:
        axes[1].text(0.5, 0.5, "No --stats TSV", ha="center", va="center")
        axes[1].set_axis_off()
    fig.tight_layout()
    png = Path(str(out_prefix) + ".png")
    pdf = Path(str(out_prefix) + ".pdf")
    fig.savefig(png, dpi=220, bbox_inches="tight")
    fig.savefig(pdf, bbox_inches="tight")
    plt.close(fig)
    summary = {
        "png": str(png),
        "pdf": str(pdf),
        "n_sites": conc.get("n_sites"),
        "concordance_frac": conc.get("concordance_frac"),
        "concordance_r": conc.get("concordance_r"),
        **hap,
    }
    print(json.dumps(summary, sort_keys=True, default=str))


def gsutil_prefix() -> list[str]:
    cmd = ["gsutil"]
    project = os.environ.get("GOOGLE_PROJECT", "").strip()
    if project:
        cmd.extend(["-u", project])
    return cmd


def attr_to_str(value: Any) -> str:
    """Flatten a Terra / firecloud attribute to a string URI or name."""
    if value is None:
        return ""
    if isinstance(value, dict):
        if "entityName" in value:
            return str(value.get("entityName") or "").strip()
        for key in ("path", "value", "url"):
            inner = value.get(key)
            if inner:
                return attr_to_str(inner)
        items = value.get("items")
        if isinstance(items, list) and items:
            return attr_to_str(items[0])
        return ""
    return str(value).strip()


def _entities_to_rows(entities: list[dict[str, Any]], *, id_column: str) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for ent in entities:
        entity_id = str(ent.get("name") or "").strip()
        if not entity_id:
            continue
        attrs = ent.get("attributes") or {}
        row = {id_column: entity_id}
        for key, value in attrs.items():
            row[str(key)] = attr_to_str(value)
        rows.append(row)
    return rows


def _read_entity_tsv(text: str, *, id_column: str) -> "Any":
    import pandas as pd

    df = pd.read_csv(StringIO(text), sep="\t", dtype=str)
    if df.empty:
        raise RuntimeError("Terra entity TSV was empty")
    if id_column not in df.columns:
        candidates = [c for c in df.columns if str(c).startswith("entity:") or str(c).endswith("_id")]
        if not candidates:
            raise ValueError(f"No entity id column in Terra TSV; columns={list(df.columns)}")
        df = df.rename(columns={candidates[0]: id_column})
    return df


def _fetch_entities_tsv(fapi: Any, namespace: str, workspace: str, entity_type: str) -> Any:
    """Download the data table as TSV (Terra UI export). Custom entity types need flexible."""
    try:
        return fapi.get_entities_tsv(namespace, workspace, entity_type, model="flexible")
    except TypeError:
        return fapi.get_entities_tsv(namespace, workspace, entity_type)


def fetch_phased_bams_table(
    *,
    tsv: Path | str | None = None,
    from_firecloud: bool = False,
    namespace: str = DEFAULT_NAMESPACE,
    workspace: str = DEFAULT_WORKSPACE,
    entity_type: str = DEFAULT_ENTITY_TYPE,
    id_column: str = DEFAULT_ID_COLUMN,
) -> "Any":
    """Load `aou2_v1_phased_bams` from a TSV export or the Firecloud/Terra API."""
    import pandas as pd

    if tsv:
        df = pd.read_csv(tsv, sep="\t", dtype=str)
        if id_column not in df.columns:
            candidates = [c for c in df.columns if c.startswith("entity:") or c.endswith("_id")]
            if not candidates:
                raise ValueError(f"No entity id column in {tsv}; columns={list(df.columns)}")
            df = df.rename(columns={candidates[0]: id_column})
        return df
    if not from_firecloud:
        raise ValueError("fetch_phased_bams_table requires tsv= or from_firecloud=True")
    from firecloud import api as fapi

    where = f"{namespace}/{workspace}/{entity_type}"
    tsv_error = ""
    if hasattr(fapi, "get_entities_tsv"):
        resp = _fetch_entities_tsv(fapi, namespace, workspace, entity_type)
        if resp.status_code == 200 and str(getattr(resp, "text", "")).strip():
            return _read_entity_tsv(resp.text, id_column=id_column)
        tsv_error = f"get_entities_tsv {resp.status_code}: {str(getattr(resp, 'text', ''))[:300]}"

    resp = fapi.get_entities(namespace, workspace, entity_type)
    if resp.status_code != 200:
        detail = tsv_error + ("; " if tsv_error else "") + resp.text[:500]
        raise RuntimeError(f"firecloud get_entities failed ({resp.status_code}) for {where}: {detail}")
    entities = resp.json()
    if not isinstance(entities, list):
        raise RuntimeError(f"Unexpected firecloud payload type: {type(entities)}")
    rows = _entities_to_rows(entities, id_column=id_column)
    if not rows:
        raise RuntimeError(f"No entities returned for {where}")
    return pd.DataFrame(rows, dtype=str)


def inventory_from_frame(df: "Any", *, id_column: str) -> dict[str, Any]:
    def nonempty(col: str) -> int:
        if col not in df.columns:
            return 0
        series = df[col].astype(str).str.strip()
        return int(((series != "") & (series.str.lower() != "nan") & series.str.startswith("gs://")).sum())

    n = len(df)
    n_combined = nonempty("combined_bed")
    n_hap1 = nonempty("hap1_bed")
    n_hap2 = nonempty("hap2_bed")
    stats_col = stats_uri_column(df)
    sites_col = sites_uri_column(df)
    if "hap1_bed" in df.columns and "hap2_bed" in df.columns:
        hap1 = df["hap1_bed"].astype(str).str.startswith("gs://")
        hap2 = df["hap2_bed"].astype(str).str.startswith("gs://")
        n_pair = int((hap1 & hap2).sum())
    else:
        n_pair = 0
    return {
        "n_rows": n,
        "n_combined_bed": n_combined,
        "n_hap1_bed": n_hap1,
        "n_hap2_bed": n_hap2,
        "n_haplotype_pair": n_pair,
        "n_stats_tsv": nonempty(stats_col) if stats_col else 0,
        "stats_uri_column": stats_col or "",
        "n_sites": nonempty(sites_col) if sites_col else 0,
        "sites_uri_column": sites_col or "",
        "id_column": id_column,
    }


def gs_uri_column(
    df: "Any",
    candidates: Sequence[str],
    *,
    extra_suffixes: Sequence[str] = (),
) -> str | None:
    columns = list(df.columns)
    ordered: list[str] = []
    for name in candidates:
        if name in columns:
            ordered.append(name)
    for suffix in extra_suffixes:
        ordered.extend(c for c in columns if c not in ordered and str(c).endswith(suffix))
    for col in ordered:
        series = df[col].astype(str).str.strip()
        if int(series.str.startswith("gs://").sum()) > 0:
            return col
    return ordered[0] if ordered else None


def stats_uri_column(df: "Any") -> str | None:
    return gs_uri_column(df, STATS_URI_CANDIDATES, extra_suffixes=("stats_tsv",))


def sites_uri_column(df: "Any") -> str | None:
    return gs_uri_column(df, SITES_URI_CANDIDATES, extra_suffixes=("sites",))


def count_site_dumps(dumps_dir: Path) -> int:
    if not dumps_dir.is_dir():
        return 0
    return len(
        [
            path
            for path in dumps_dir.iterdir()
            if path.is_file()
            and (path.suffix in {".gz", ".tsv"} or path.name.endswith(".tsv.gz"))
        ]
    )


def caller_group(value: Any) -> str:
    text = str(value or "").strip().lower()
    if text.startswith(CALLER_PRIMROSE_PREFIX):
        return "primrose"
    if text.startswith(CALLER_JASMINE_PREFIX):
        return "jasmine"
    return ""


def _copy_gs_jobs(
    jobs_spec: list[tuple[str, str, Path]],
    *,
    jobs: int,
    force: bool,
) -> tuple[int, int]:
    def one(sample_id: str, uri: str, dest: Path) -> str:
        if dest.is_file() and dest.stat().st_size > 0 and not force:
            return "skip"
        dest.parent.mkdir(parents=True, exist_ok=True)
        subprocess.check_call(gsutil_prefix() + ["cp", uri, str(dest)])
        return "ok"

    workers = max(1, int(jobs))
    if workers == 1:
        results = [one(*job) for job in jobs_spec]
    else:
        results = []
        with ThreadPoolExecutor(max_workers=workers) as pool:
            futs = {pool.submit(one, *job): job[0] for job in jobs_spec}
            for fut in as_completed(futs):
                results.append(fut.result())
    n_ok = sum(1 for status in results if status != "skip")
    n_skip = sum(1 for status in results if status == "skip")
    return n_ok, n_skip


def pull_stats_from_table(
    df: "Any",
    stats_dir: Path,
    *,
    id_column: str = DEFAULT_ID_COLUMN,
    jobs: int = 8,
    force: bool = False,
) -> dict[str, Any]:
    """Download per-sample stats_tsv URIs to `{sample_id}.tsv`."""
    stats_dir.mkdir(parents=True, exist_ok=True)
    col = stats_uri_column(df)
    if not col:
        raise ValueError(
            "No stats_tsv column on the data table. Submit PbCpgSampleStats.wdl "
            "with outputs written back to aou2_v1_phased_bams, or set METH_STATS_DIR "
            "to a directory of already-downloaded TSVs."
        )
    if id_column not in df.columns:
        raise ValueError(f"Missing {id_column}; columns={list(df.columns)}")
    jobs_spec: list[tuple[str, str, Path]] = []
    for _, row in df.iterrows():
        sample_id = str(row[id_column]).strip()
        uri = str(row[col]).strip()
        if not sample_id or not uri.startswith("gs://"):
            continue
        dest = stats_dir / f"{sample_id}.tsv"
        jobs_spec.append((sample_id, uri, dest))
    if not jobs_spec:
        raise ValueError(f"Column {col} has no gs:// stats_tsv URIs ({len(df)} table rows)")
    n_ok, n_skip = _copy_gs_jobs(jobs_spec, jobs=jobs, force=force)
    return {
        "stats_uri_column": col,
        "n_table_rows": int(len(df)),
        "n_with_uri": len(jobs_spec),
        "n_downloaded": n_ok,
        "n_skipped": n_skip,
        "stats_dir": str(stats_dir),
    }


def annotate_callers(
    df: "Any",
    *,
    id_column: str,
    covariates: Path | str | None = None,
    stats: Path | str | None = None,
) -> "Any":
    import pandas as pd

    out = df.copy()
    out["_sample_id"] = out[id_column].astype(str)
    if covariates and Path(covariates).is_file():
        cov = pd.read_csv(covariates, dtype={"research_id": str}, low_memory=False)
        keep = [c for c in ("research_id", "pb_meth_caller", "coverage") if c in cov.columns]
        out = out.merge(cov[keep], how="left", left_on="_sample_id", right_on="research_id")
    if stats:
        stats_path = Path(stats)
        if stats_path.is_dir():
            stat_df = load_stats_dir(stats_path)
        elif stats_path.is_file():
            stat_df = pd.read_csv(stats_path, sep="\t", dtype={"sample_id": str}, low_memory=False)
        else:
            stat_df = None
        if stat_df is not None and "sample_id" in stat_df.columns:
            cols = [c for c in ("sample_id", "combined_mean_cov") if c in stat_df.columns]
            out = out.merge(stat_df[cols], how="left", left_on="_sample_id", right_on="sample_id")
    if "pb_meth_caller" in out.columns:
        out["caller_group"] = out["pb_meth_caller"].map(caller_group)
    else:
        out["caller_group"] = ""
    return out


def pick_concordance_samples(
    df: "Any",
    *,
    id_column: str = DEFAULT_ID_COLUMN,
    covariates: Path | str | None = None,
    stats: Path | str | None = None,
    n_per_group: int = CONCORDANCE_N_PER_GROUP,
    seed: int = 1,
    require_sites: bool = False,
) -> "Any":
    """Choose n Primrose + n Jasmine samples (highest mean CpG coverage when known)."""
    import pandas as pd

    labeled = annotate_callers(df, id_column=id_column, covariates=covariates, stats=stats)
    sites_col = sites_uri_column(labeled)
    if require_sites:
        if not sites_col:
            return labeled.iloc[0:0].copy()
        labeled = labeled.loc[labeled[sites_col].astype(str).str.startswith("gs://")].copy()
    labeled = labeled.loc[labeled["caller_group"].isin({"primrose", "jasmine"})].copy()
    if "combined_bed" in labeled.columns:
        labeled = labeled.loc[labeled["combined_bed"].astype(str).str.startswith("gs://")].copy()
    if labeled.empty:
        return labeled
    if "combined_mean_cov" in labeled.columns:
        labeled["_rank"] = pd.to_numeric(labeled["combined_mean_cov"], errors="coerce")
        labeled = labeled.sort_values(["caller_group", "_rank"], ascending=[True, False])
    else:
        labeled = labeled.sample(frac=1.0, random_state=seed).sort_values("caller_group")
    picked = labeled.groupby("caller_group", sort=False).head(int(n_per_group))
    return picked.drop(columns=["_rank"], errors="ignore").reset_index(drop=True)


def write_concordance_entity_tsv(df: "Any", path: Path, *, id_column: str) -> Path:
    keep = df[[id_column]].copy()
    if "combined_bed" in df.columns:
        keep["combined_bed"] = df["combined_bed"]
    keep = keep.rename(columns={id_column: CONCORDANCE_ID_COLUMN})
    path.parent.mkdir(parents=True, exist_ok=True)
    keep.to_csv(path, sep="\t", index=False)
    return path


def pull_dumps_from_table(
    df: "Any",
    dumps_dir: Path,
    *,
    id_column: str = DEFAULT_ID_COLUMN,
    covariates: Path | str | None = None,
    stats: Path | str | None = None,
    n_per_group: int = CONCORDANCE_N_PER_GROUP,
    jobs: int = 8,
    force: bool = False,
) -> dict[str, Any]:
    """Download chr-site dumps for a Primrose/Jasmine subsample."""
    dumps_dir.mkdir(parents=True, exist_ok=True)
    col = sites_uri_column(df)
    if not col:
        raise ValueError(
            "No sites column on the data table. Submit PbCpgChromSites.wdl on a "
            "Primrose+Jasmine subsample and write the `sites` output back to the entity."
        )
    picked = pick_concordance_samples(
        df,
        id_column=id_column,
        covariates=covariates,
        stats=stats,
        n_per_group=n_per_group,
        require_sites=True,
    )
    if picked.empty:
        labeled = annotate_callers(df, id_column=id_column, covariates=covariates, stats=stats)
        picked = labeled.loc[labeled[col].astype(str).str.startswith("gs://")].head(
            int(n_per_group) * 2
        )
    jobs_spec: list[tuple[str, str, Path]] = []
    for _, row in picked.iterrows():
        sample_id = str(row[id_column]).strip()
        uri = str(row[col]).strip()
        if not sample_id or not uri.startswith("gs://"):
            continue
        dest = dumps_dir / f"{sample_id}.sites.tsv.gz"
        jobs_spec.append((sample_id, uri, dest))
    if not jobs_spec:
        raise ValueError(f"Column {col} has no gs:// sites URIs to download")
    n_ok, n_skip = _copy_gs_jobs(jobs_spec, jobs=jobs, force=force)
    n_primrose = int((picked["caller_group"] == "primrose").sum()) if "caller_group" in picked.columns else 0
    n_jasmine = int((picked["caller_group"] == "jasmine").sum()) if "caller_group" in picked.columns else 0
    return {
        "sites_uri_column": col,
        "n_table_rows": int(len(df)),
        "n_with_uri": int(df[col].astype(str).str.startswith("gs://").sum()),
        "n_downloaded": n_ok,
        "n_skipped": n_skip,
        "n_primrose": n_primrose,
        "n_jasmine": n_jasmine,
        "dumps_dir": str(dumps_dir),
    }


def cmd_inventory(args: argparse.Namespace) -> None:
    id_column = args.id_column
    df = fetch_phased_bams_table(
        tsv=args.tsv or None,
        from_firecloud=bool(args.from_firecloud),
        namespace=args.namespace,
        workspace=args.workspace,
        entity_type=args.entity_type,
        id_column=id_column,
    )
    if id_column not in df.columns:
        candidates = [c for c in df.columns if c.startswith("entity:") or c.endswith("_id")]
        id_column = candidates[0] if candidates else id_column
    result = inventory_from_frame(df, id_column=id_column)
    if args.out:
        out = Path(args.out)
        out.parent.mkdir(parents=True, exist_ok=True)
        write_json(out, result)
        table_path = out.with_name(out.stem + ".table.tsv")
        df.to_csv(table_path, sep="\t", index=False)
        result["table_tsv"] = str(table_path)
        write_json(out, result)
    print(json.dumps(result, sort_keys=True))


def cmd_pull_stats(args: argparse.Namespace) -> None:
    df = fetch_phased_bams_table(
        tsv=args.tsv or None,
        from_firecloud=bool(args.from_firecloud),
        namespace=args.namespace,
        workspace=args.workspace,
        entity_type=args.entity_type,
        id_column=args.id_column,
    )
    result = pull_stats_from_table(
        df,
        Path(args.stats_dir),
        id_column=args.id_column,
        jobs=args.jobs,
        force=args.force,
    )
    print(json.dumps(result, sort_keys=True))


def cmd_pull_dumps(args: argparse.Namespace) -> None:
    id_column = args.id_column
    if args.entity_type != DEFAULT_ENTITY_TYPE and id_column == DEFAULT_ID_COLUMN:
        id_column = entity_id_column(args.entity_type)
    df = fetch_phased_bams_table(
        tsv=args.tsv or None,
        from_firecloud=bool(args.from_firecloud),
        namespace=args.namespace,
        workspace=args.workspace,
        entity_type=args.entity_type,
        id_column=id_column,
    )
    result = pull_dumps_from_table(
        df,
        Path(args.dumps_dir),
        id_column=id_column,
        covariates=args.covariates or None,
        stats=args.stats or None,
        n_per_group=args.n_per_group,
        jobs=args.jobs,
        force=args.force,
    )
    print(json.dumps(result, sort_keys=True))


def cmd_pick_concordance_samples(args: argparse.Namespace) -> None:
    df = fetch_phased_bams_table(
        tsv=args.tsv or None,
        from_firecloud=bool(args.from_firecloud),
        namespace=args.namespace,
        workspace=args.workspace,
        entity_type=args.entity_type,
        id_column=args.id_column,
    )
    picked = pick_concordance_samples(
        df,
        id_column=args.id_column,
        covariates=args.covariates or None,
        stats=args.stats or None,
        n_per_group=args.n_per_group,
        seed=args.seed,
    )
    if picked.empty:
        raise RuntimeError(
            "No Primrose/Jasmine samples with combined_bed. Join covariates "
            "(`pb_meth_caller`) and confirm combined_bed URIs."
        )
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    picked.to_csv(out, sep="\t", index=False)
    entity_path = out.with_name(out.stem + ".entities.tsv")
    write_concordance_entity_tsv(picked, entity_path, id_column=args.id_column)
    summary = {
        "n_picked": int(len(picked)),
        "n_primrose": int((picked["caller_group"] == "primrose").sum()),
        "n_jasmine": int((picked["caller_group"] == "jasmine").sum()),
        "samples_tsv": str(out),
        "entities_tsv": str(entity_path),
    }
    print(json.dumps(summary, sort_keys=True))


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="cmd", required=True)

    summarize = sub.add_parser("summarize", help="Per-sample combined/hap1/hap2 stats")
    summarize.add_argument("--sample-id", required=True)
    summarize.add_argument("--combined", default="")
    summarize.add_argument("--hap1", default="")
    summarize.add_argument("--hap2", default="")
    summarize.add_argument("--out-prefix", required=True)
    summarize.set_defaults(func=cmd_summarize)

    extract = sub.add_parser("extract-chrom", help="Write cov/score rows for one chromosome")
    extract.add_argument("--bed", required=True)
    extract.add_argument("--chrom", required=True)
    extract.add_argument("--start", type=int, default=None, help="0-based begin >= start")
    extract.add_argument("--end", type=int, default=None, help="0-based begin < end")
    extract.add_argument("--min-cov", type=int, default=10)
    extract.add_argument("--sample-id", default="")
    extract.add_argument("--out", required=True)
    extract.set_defaults(func=cmd_extract_chrom)

    merge = sub.add_parser("merge", help="Merge per-sample TSVs into manuscript numbers")
    merge.add_argument("--stats-dir", required=True)
    merge.add_argument("--out-dir", required=True)
    merge.add_argument("--covariates", default="")
    merge.add_argument("--discovery-only", action="store_true")
    merge.set_defaults(func=cmd_merge)

    conc = sub.add_parser("concordance", help="Primrose vs Jasmine site concordance")
    conc.add_argument("--bed-a", default="")
    conc.add_argument("--bed-b", default="")
    conc.add_argument("--dumps-dir", default="")
    conc.add_argument("--labels", default="")
    conc.add_argument("--min-cov", type=int, default=10)
    conc.add_argument("--min-n-per-group", type=int, default=5)
    conc.add_argument("--out", required=True)
    conc.add_argument("--pairs-out", default="", help="Site-mean pairs TSV.gz for plotting")
    conc.set_defaults(func=cmd_concordance)

    plot = sub.add_parser("plot-s3", help="Hexbin site-mean concordance and haplotype occupancy")
    plot.add_argument("--pairs", default="", help="mean_primrose/mean_jasmine TSV from concordance")
    plot.add_argument("--dumps-dir", default="")
    plot.add_argument("--labels", default="")
    plot.add_argument("--stats", default="", help="pbcpg_sample_stats.tsv for haplotype vs coverage")
    plot.add_argument("--min-cov", type=int, default=10)
    plot.add_argument("--min-n-per-group", type=int, default=5)
    plot.add_argument("--write-pairs", default="")
    plot.add_argument("--out-prefix", required=True)
    plot.set_defaults(func=cmd_plot_s3)

    inv = sub.add_parser("inventory", help="Count methylation URIs in the Terra table")
    inv.add_argument("--tsv", default="")
    inv.add_argument("--from-firecloud", action="store_true")
    inv.add_argument("--namespace", default=os.environ.get("METH_TERRA_NAMESPACE", DEFAULT_NAMESPACE))
    inv.add_argument("--workspace", default=os.environ.get("METH_TERRA_WORKSPACE", DEFAULT_WORKSPACE))
    inv.add_argument("--entity-type", default=os.environ.get("METH_ENTITY_TYPE", DEFAULT_ENTITY_TYPE))
    inv.add_argument("--id-column", default=DEFAULT_ID_COLUMN)
    inv.add_argument("--out", default="")
    inv.set_defaults(func=cmd_inventory)

    pull = sub.add_parser("pull-stats", help="Download stats_tsv files from the Terra table")
    pull.add_argument("--tsv", default="")
    pull.add_argument("--from-firecloud", action="store_true")
    pull.add_argument("--namespace", default=os.environ.get("METH_TERRA_NAMESPACE", DEFAULT_NAMESPACE))
    pull.add_argument("--workspace", default=os.environ.get("METH_TERRA_WORKSPACE", DEFAULT_WORKSPACE))
    pull.add_argument("--entity-type", default=os.environ.get("METH_ENTITY_TYPE", DEFAULT_ENTITY_TYPE))
    pull.add_argument("--id-column", default=DEFAULT_ID_COLUMN)
    pull.add_argument("--stats-dir", required=True)
    pull.add_argument("--jobs", type=int, default=8)
    pull.add_argument("--force", action="store_true")
    pull.set_defaults(func=cmd_pull_stats)

    dumps = sub.add_parser("pull-dumps", help="Download PbCpgChromSites sites files from the Terra table")
    dumps.add_argument("--tsv", default="")
    dumps.add_argument("--from-firecloud", action="store_true")
    dumps.add_argument("--namespace", default=os.environ.get("METH_TERRA_NAMESPACE", DEFAULT_NAMESPACE))
    dumps.add_argument("--workspace", default=os.environ.get("METH_TERRA_WORKSPACE", DEFAULT_WORKSPACE))
    dumps.add_argument("--entity-type", default=os.environ.get("METH_CONCORDANCE_ENTITY", DEFAULT_ENTITY_TYPE))
    dumps.add_argument("--id-column", default=DEFAULT_ID_COLUMN)
    dumps.add_argument("--dumps-dir", required=True)
    dumps.add_argument("--covariates", default="")
    dumps.add_argument("--stats", default="")
    dumps.add_argument("--n-per-group", type=int, default=CONCORDANCE_N_PER_GROUP)
    dumps.add_argument("--jobs", type=int, default=8)
    dumps.add_argument("--force", action="store_true")
    dumps.set_defaults(func=cmd_pull_dumps)

    pick = sub.add_parser("pick-concordance-samples", help="Choose Primrose+Jasmine rows for PbCpgChromSites")
    pick.add_argument("--tsv", default="")
    pick.add_argument("--from-firecloud", action="store_true")
    pick.add_argument("--namespace", default=os.environ.get("METH_TERRA_NAMESPACE", DEFAULT_NAMESPACE))
    pick.add_argument("--workspace", default=os.environ.get("METH_TERRA_WORKSPACE", DEFAULT_WORKSPACE))
    pick.add_argument("--entity-type", default=os.environ.get("METH_ENTITY_TYPE", DEFAULT_ENTITY_TYPE))
    pick.add_argument("--id-column", default=DEFAULT_ID_COLUMN)
    pick.add_argument("--covariates", default="")
    pick.add_argument("--stats", default="")
    pick.add_argument("--n-per-group", type=int, default=CONCORDANCE_N_PER_GROUP)
    pick.add_argument("--seed", type=int, default=1)
    pick.add_argument("--out", required=True)
    pick.set_defaults(func=cmd_pick_concordance_samples)
    return parser


def main(argv: Sequence[str] | None = None) -> None:
    parser = build_parser()
    args = parser.parse_args(argv)
    if args.cmd in {"inventory", "pull-stats", "pull-dumps", "pick-concordance-samples"} and not args.tsv and not args.from_firecloud:
        parser.error(f"{args.cmd} requires --tsv or --from-firecloud")
    if args.cmd == "concordance":
        paired = bool(args.bed_a and args.bed_b)
        grouped = bool(args.dumps_dir and args.labels)
        if paired == grouped:
            parser.error("concordance requires --bed-a/--bed-b or --dumps-dir/--labels")
    if args.cmd == "plot-s3":
        has_pairs = bool(args.pairs)
        has_dumps = bool(args.dumps_dir and args.labels)
        if has_pairs == has_dumps:
            parser.error("plot-s3 requires --pairs or --dumps-dir/--labels")
    args.func(args)


if __name__ == "__main__":
    main()
