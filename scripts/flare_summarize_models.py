#!/usr/bin/env python3
"""Parse per-population FLARE .model files and concatenate .global.anc.gz."""

from __future__ import annotations

import argparse
import csv
import gzip
import sys
from pathlib import Path
from typing import Optional, TextIO

# Allow ``python scripts/flare_summarize_models.py`` without PYTHONPATH.
_SCRIPTS = Path(__file__).resolve().parent
if str(_SCRIPTS) not in sys.path:
    sys.path.insert(0, str(_SCRIPTS))

from flare_model import Model, parse_model, pop_of_model_path  # noqa: E402

# Back-compat re-exports for tests / callers.
pop_from_model_path = pop_of_model_path


def open_text(path: Path) -> TextIO:
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")
    return path.open()


def write_models_tsv(models: list[tuple[str, Model]], dest: Path) -> None:
    """Write population, T, prop_<anc>, and mean mu_<anc> (miscopy) columns."""
    if not models:
        raise SystemExit("no models to summarize")
    ancestries = models[0][1].ancestries
    dest.parent.mkdir(parents=True, exist_ok=True)
    with dest.open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(
            [
                "population",
                "t_gen",
                *[f"prop_{a}" for a in ancestries],
                *[f"mu_{a}" for a in ancestries],
            ]
        )
        for pop, parsed in models:
            if parsed.ancestries != ancestries:
                raise SystemExit(
                    f"{pop}: ancestry list {parsed.ancestries} != {ancestries}"
                )
            mean_mu = parsed.mean_miscopy()
            w.writerow(
                [
                    pop,
                    f"{parsed.t_gen:.6g}",
                    *[f"{m:.8g}" for m in parsed.props],
                    *[f"{m:.8g}" for m in mean_mu],
                ]
            )


def merge_global_anc(paths: list[Path], dest: Path) -> int:
    header: Optional[str] = None
    n = 0
    dest.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(dest, "wt") as out:
        for path in paths:
            with open_text(path) as fh:
                first = True
                for line in fh:
                    if first:
                        if header is None:
                            header = line
                            out.write(line if line.endswith("\n") else line + "\n")
                        elif line.strip() != header.strip():
                            raise SystemExit(
                                f"{path}: global.anc header {line.strip()!r} "
                                f"!= {header.strip()!r}"
                            )
                        first = False
                        continue
                    if not line.strip():
                        continue
                    out.write(line if line.endswith("\n") else line + "\n")
                    n += 1
    if header is None:
        raise SystemExit("no global.anc headers")
    return n


def main(argv: Optional[list[str]] = None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--models", nargs="+", type=Path, required=True)
    p.add_argument("--pops", nargs="*", default=[])
    p.add_argument("--global-anc", nargs="*", type=Path, default=[])
    p.add_argument("--out-models-tsv", type=Path, required=True)
    p.add_argument("--out-global-anc", type=Path)
    args = p.parse_args(argv)

    if args.pops and len(args.pops) != len(args.models):
        raise SystemExit("--pops must match --models 1:1")
    # Prefer explicit shard pops; filename is a fallback only.
    pops = list(args.pops) if args.pops else [pop_of_model_path(path) for path in args.models]
    parsed = [(pop, parse_model(path)) for pop, path in zip(pops, args.models)]
    write_models_tsv(parsed, args.out_models_tsv)
    n_global = 0
    if args.out_global_anc:
        if not args.global_anc:
            raise SystemExit("--out-global-anc requires --global-anc")
        n_global = merge_global_anc(args.global_anc, args.out_global_anc)
    print(
        "wrote",
        args.out_models_tsv,
        "T=" + ", ".join(f"{pop}:{block.t_gen:.2f}" for pop, block in parsed),
        f"global_n={n_global}",
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
