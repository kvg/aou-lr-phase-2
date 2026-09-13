#!/usr/bin/env python3
"""Split VCF sample IDs into per-population keep lists from the covariates table.

FLARE estimates one generations-since-admixture T for the whole gt VCF. AoU is
not one admixed population; this script emits one sample list per continental
group so each FLARE run can have its own T.

HPRC/HGSVC3/GIAB controls (HG*/NA*) overlap the LAI reference panel and are
dropped by default.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import sys
from collections import defaultdict
from pathlib import Path
from typing import Iterable, Optional, TextIO


MISSING = frozenset({"", "NA", "NAN", "NONE", "NULL", ".", "NAN"})
TRUE_VALUES = frozenset({"true", "t", "1", "yes", "y"})


def open_text(path: Path) -> TextIO:
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt", newline="")
    return path.open(newline="")


def norm_id(value: str) -> str:
    return str(value).strip()


def norm_pop(value: str) -> Optional[str]:
    text = str(value).strip()
    if text.upper() in MISSING:
        return None
    return text.upper()


def parse_bool(value: str) -> bool:
    return str(value).strip().lower() in TRUE_VALUES


def looks_like_control_id(sample_id: str) -> bool:
    token = str(sample_id).strip().upper()
    return token.startswith("HG") or token.startswith("NA")


def read_vcf_samples(path: Path) -> list[str]:
    samples: list[str] = []
    seen: set[str] = set()
    with path.open() as fh:
        for line in fh:
            sid = norm_id(line)
            if not sid or sid.startswith("#"):
                continue
            if sid in seen:
                raise SystemExit(f"duplicate VCF sample ID: {sid}")
            seen.add(sid)
            samples.append(sid)
    if not samples:
        raise SystemExit(f"no sample IDs in {path}")
    return samples


def parse_csv_list(raw: str) -> list[str]:
    if not raw or not raw.strip():
        return []
    return [norm_pop(part) or "" for part in raw.split(",") if part.strip()]


def load_covariate_rows(
    covariates: Path,
    id_column: str,
    pop_column: str,
    control_column: str,
) -> dict[str, tuple[Optional[str], bool]]:
    """research_id → (population or None, is_reference_control)."""
    with open_text(covariates) as fh:
        reader = csv.DictReader(fh)
        if reader.fieldnames is None:
            raise SystemExit(f"{covariates}: empty covariates file")
        fields = [f.lstrip("\ufeff") for f in reader.fieldnames]
        reader.fieldnames = fields
        if id_column not in fields:
            raise SystemExit(
                f"{covariates}: missing ID column {id_column!r}; have {fields[:12]}"
            )
        if pop_column not in fields:
            raise SystemExit(
                f"{covariates}: missing population column {pop_column!r}; have {fields[:12]}"
            )
        has_control = control_column in fields
        out: dict[str, tuple[Optional[str], bool]] = {}
        for row in reader:
            rid = norm_id(row.get(id_column) or "")
            if not rid:
                continue
            if rid in out:
                raise SystemExit(f"{covariates}: duplicate {id_column} {rid}")
            flagged = parse_bool(row.get(control_column) or "") if has_control else False
            is_control = flagged or looks_like_control_id(rid)
            out[rid] = (norm_pop(row.get(pop_column) or ""), is_control)
        return out


def split_samples(
    vcf_samples: list[str],
    cov_rows: dict[str, tuple[Optional[str], bool]],
    *,
    include_pops: Optional[Iterable[str]] = None,
    drop_pops: Optional[Iterable[str]] = None,
    min_samples: int = 50,
    unlabeled_pop: str = "",
    exclude_controls: bool = True,
) -> dict:
    include = {p for p in (norm_pop(p) or "" for p in (include_pops or [])) if p}
    drop = {p for p in (norm_pop(p) or "" for p in (drop_pops or [])) if p}
    unlabeled = (norm_pop(unlabeled_pop) or "") if unlabeled_pop else ""

    by_pop: dict[str, list[str]] = defaultdict(list)
    unmatched_vcf: list[str] = []
    dropped: list[tuple[str, str, str]] = []
    excluded_controls: list[str] = []

    for sid in vcf_samples:
        pop, flagged = cov_rows.get(sid, (None, False))
        if pop is not None:
            pop = norm_pop(pop)
        is_control = flagged or looks_like_control_id(sid)
        if exclude_controls and is_control:
            excluded_controls.append(sid)
            dropped.append((sid, pop or "", "control"))
            continue
        if pop is None and unlabeled:
            pop = unlabeled
        if pop is None:
            unmatched_vcf.append(sid)
            continue
        if pop in drop:
            dropped.append((sid, pop, "drop_pop"))
            continue
        if include and pop not in include:
            dropped.append((sid, pop, "not_included"))
            continue
        by_pop[pop].append(sid)

    kept: dict[str, list[str]] = {}
    skipped_small: dict[str, int] = {}
    for pop, samples in sorted(by_pop.items()):
        if len(samples) < min_samples:
            skipped_small[pop] = len(samples)
            continue
        kept[pop] = samples

    n_vcf = len(vcf_samples)
    n_kept = sum(len(v) for v in kept.values())
    return {
        "n_vcf_samples": n_vcf,
        "n_kept": n_kept,
        "frac_kept": (n_kept / n_vcf) if n_vcf else 0.0,
        "n_unmatched_vcf": len(unmatched_vcf),
        "frac_unmatched_vcf": (len(unmatched_vcf) / n_vcf) if n_vcf else 0.0,
        "n_excluded_controls": len(excluded_controls),
        "kept": kept,
        "unmatched_vcf": unmatched_vcf,
        "dropped": dropped,
        "excluded_controls": excluded_controls,
        "skipped_small": skipped_small,
        "include_pops": sorted(include),
        "drop_pops": sorted(drop),
        "min_samples": min_samples,
        "unlabeled_pop": unlabeled or None,
        "exclude_controls": exclude_controls,
    }


def write_outputs(result: dict, out_dir: Path) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    keep_dir = out_dir / "keeps"
    keep_dir.mkdir(exist_ok=True)

    populations: list[str] = []
    manifest_rows: list[list[str]] = []
    for pop, samples in result["kept"].items():
        keep_path = keep_dir / f"{pop}.samples.txt"
        keep_path.write_text("".join(f"{s}\n" for s in samples))
        populations.append(pop)
        manifest_rows.append([pop, str(len(samples)), str(keep_path.name)])

    (out_dir / "populations.txt").write_text("".join(f"{p}\n" for p in populations))
    with (out_dir / "manifest.tsv").open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["population", "n_samples", "keep_list"])
        w.writerows(manifest_rows)

    (out_dir / "unmatched_vcf_samples.txt").write_text(
        "".join(f"{s}\n" for s in result["unmatched_vcf"])
    )
    (out_dir / "excluded_controls.txt").write_text(
        "".join(f"{s}\n" for s in result["excluded_controls"])
    )
    with (out_dir / "dropped_samples.tsv").open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["sample", "population", "reason"])
        w.writerows(result["dropped"])
        for pop, n in result["skipped_small"].items():
            w.writerow([f"# skipped_{pop}_n={n}", pop, "too_small"])

    summary = {
        "n_vcf_samples": result["n_vcf_samples"],
        "n_kept": result["n_kept"],
        "frac_kept": result["frac_kept"],
        "n_unmatched_vcf": result["n_unmatched_vcf"],
        "frac_unmatched_vcf": result["frac_unmatched_vcf"],
        "n_excluded_controls": result["n_excluded_controls"],
        "n_populations": len(populations),
        "populations": {
            pop: len(samples) for pop, samples in result["kept"].items()
        },
        "skipped_small": result["skipped_small"],
        "include_pops": result["include_pops"],
        "drop_pops": result["drop_pops"],
        "min_samples": result["min_samples"],
        "unlabeled_pop": result["unlabeled_pop"],
        "exclude_controls": result["exclude_controls"],
    }
    (out_dir / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")


def main(argv: Optional[list[str]] = None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--covariates", required=True, type=Path)
    p.add_argument("--vcf-samples", required=True, type=Path)
    p.add_argument("--out-dir", required=True, type=Path)
    p.add_argument("--id-column", default="research_id")
    p.add_argument("--pop-column", default="population")
    p.add_argument("--control-column", default="is_reference_control")
    p.add_argument(
        "--include-pops",
        default="",
        help="Comma-separated populations to keep (default: all that pass min-samples)",
    )
    p.add_argument(
        "--drop-pops",
        default="",
        help="Comma-separated populations to omit (e.g. OTH,MID)",
    )
    p.add_argument("--min-samples", type=int, default=50)
    p.add_argument(
        "--unlabeled-pop",
        default="",
        help="Assign VCF samples missing a covariate label to this population",
    )
    p.add_argument(
        "--exclude-controls",
        dest="exclude_controls",
        action="store_true",
        help="Drop HPRC/HGSVC3/GIAB controls (default)",
    )
    p.add_argument(
        "--keep-controls",
        dest="exclude_controls",
        action="store_false",
        help="Leave reference-control samples in the FLARE target",
    )
    p.set_defaults(exclude_controls=True)
    p.add_argument(
        "--max-unmatched-vcf-frac",
        type=float,
        default=0.02,
        help="Fail if this fraction of VCF samples have no usable population",
    )
    args = p.parse_args(argv)

    if args.min_samples < 1:
        raise SystemExit("--min-samples must be >= 1")

    vcf_samples = read_vcf_samples(args.vcf_samples)
    cov_rows = load_covariate_rows(
        args.covariates, args.id_column, args.pop_column, args.control_column
    )
    result = split_samples(
        vcf_samples,
        cov_rows,
        include_pops=parse_csv_list(args.include_pops),
        drop_pops=parse_csv_list(args.drop_pops),
        min_samples=args.min_samples,
        unlabeled_pop=args.unlabeled_pop,
        exclude_controls=args.exclude_controls,
    )
    write_outputs(result, args.out_dir)

    if not result["kept"]:
        raise SystemExit(
            "no populations left after filters "
            f"(min_samples={args.min_samples}, skipped={result['skipped_small']})"
        )
    unmatched_frac = result["frac_unmatched_vcf"]
    if unmatched_frac > args.max_unmatched_vcf_frac:
        raise SystemExit(
            f"unmatched VCF samples {unmatched_frac:.3f} > "
            f"--max-unmatched-vcf-frac {args.max_unmatched_vcf_frac} "
            f"({result['n_unmatched_vcf']}/{result['n_vcf_samples']})"
        )

    pops = ", ".join(f"{k}={len(v)}" for k, v in result["kept"].items())
    print(
        f"kept {result['n_kept']}/{result['n_vcf_samples']} samples in "
        f"{len(result['kept'])} populations: {pops}; "
        f"excluded_controls={result['n_excluded_controls']}",
        file=sys.stderr,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
