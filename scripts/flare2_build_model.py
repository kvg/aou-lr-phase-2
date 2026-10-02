#!/usr/bin/env python3
"""Build a FLARE2 model (clustering step) with stable ancestry labels.

FLARE2 (Browning, Temple & Browning 2025) is three steps on FLARE 0.6:

1. ``flare panel-probs=true`` writes a ``.panels`` file
2. upstream ``create_model_file.py nanc panels out`` clusters copying
   probabilities with a Gaussian mixture and writes an initial ``.model``
3. ``flare model=<that file> update-p=true`` estimates ancestry

This wraps step 2. Upstream names clusters ``anc_0 .. anc_{n-1}`` in GMM
order, and the GMM has no fixed seed, so two runs (two population shards, two
chromosomes) can label the same ancestry differently. FELIX needs one label
set for the whole cohort and genome, so we:

* seed numpy before running the upstream script (deterministic GMM),
* rename each cluster by its dominant reference panel (``anc0_afr`` …),
* reorder ancestries canonically (by dominant panel, then weight) and permute
  every per-ancestry block to match,
* gate on the upstream lag-1 autocorrelation (spurious ``nanc`` check).

Train the model once (pooled cohort subsample, one chromosome) and apply it to
every population shard and chromosome via ``FlareByPopulation.template_model``.

Example::

    python3 scripts/flare2_build_model.py \\
      --create-model-script /opt/flare2/create_model_file.py \\
      --panels train.panels --nanc 5 --seed 12345 \\
      --min-autocorr 0.25 --out-prefix chr20.flare2_nanc5
"""

from __future__ import annotations

import argparse
import json
import re
import subprocess
import sys
from pathlib import Path
from typing import Any, Optional, Sequence

sys.path.insert(0, str(Path(__file__).resolve().parent))
from flare_model import Model, parse_model, write_model  # noqa: E402

_MIN_AUTOCORR_RE = re.compile(r"#\s*min auto-correlation is\s+(\S+)")
_MAX_CROSSCORR_RE = re.compile(r"#\s*max cross-correlation is\s+(\S+)")


_BOOTSTRAP = (
    "import runpy, sys\n"
    "import numpy as np\n"
    "np.random.seed(int(sys.argv[1]))\n"
    "sys.argv = sys.argv[2:]\n"
    "runpy.run_path(sys.argv[0], run_name='__main__')\n"
)


def run_upstream(script: Path, nanc: int, panels: Path, raw_prefix: str, seed: int) -> None:
    """Run upstream create_model_file.py in a child process with a seeded numpy RNG.

    A child process (not in-process ``runpy``) because upstream never closes its
    output handles; they only flush at interpreter exit.
    """
    subprocess.run(
        [sys.executable, "-c", _BOOTSTRAP, str(seed), str(script), str(nanc), str(panels), raw_prefix],
        check=True,
    )


def read_correlations(model_path: Path) -> dict[str, Optional[float]]:
    out: dict[str, Optional[float]] = {"min_autocorr": None, "max_crosscorr": None}
    for line in model_path.read_text().splitlines():
        m = _MIN_AUTOCORR_RE.match(line.strip())
        if m:
            out["min_autocorr"] = float(m.group(1))
        m = _MAX_CROSSCORR_RE.match(line.strip())
        if m:
            out["max_crosscorr"] = float(m.group(1))
    return out


def canonical_labels(model: Model) -> tuple[list[int], list[str], list[dict[str, Any]]]:
    """Return (order, new_names, label_rows) for a cluster model.

    ``order[k]`` is the old ancestry index placed at new index ``k``.
    """
    info = []
    for i, row in enumerate(model.panel_weights):
        j = max(range(len(row)), key=lambda c: row[c])
        info.append((model.panels[j].lower(), -row[j], i, row[j]))
    info.sort()
    order = [old for _, _, old, _ in info]
    names = [f"anc{k}_{info[k][0]}" for k in range(len(info))]
    rows = []
    for k, (panel, _, old, weight) in enumerate(info):
        rows.append(
            {
                "index": k,
                "name": names[k],
                "upstream_label": model.ancestries[old],
                "dominant_panel": panel,
                "dominant_weight": round(weight, 6),
                "panel_weights": {
                    model.panels[j]: round(model.panel_weights[old][j], 6)
                    for j in range(len(model.panels))
                },
            }
        )
    return order, names, rows


def relabel(model: Model, order: Sequence[int], names: Sequence[str]) -> Model:
    out = Model(
        ancestries=list(names),
        panels=list(model.panels),
        t_gen=model.t_gen,
        props=[model.props[i] for i in order],
        panel_weights=[list(model.panel_weights[i]) for i in order],
        miscopy=[list(model.miscopy[i]) for i in order],
        ibd_rates=[model.ibd_rates[i] for i in order],
    )
    out.validate()
    return out


def main(argv: Optional[list[str]] = None) -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--create-model-script", required=True, type=Path,
                   help="Upstream browning-lab/flare create_model_file.py")
    p.add_argument("--panels", required=True, type=Path, help="FLARE panel-probs=true .panels file")
    p.add_argument("--nanc", required=True, type=int)
    p.add_argument("--seed", type=int, default=12345)
    p.add_argument("--min-autocorr", type=float, default=0.25,
                   help="Fail when upstream min lag-1 autocorrelation is below this (<=0 disables)")
    p.add_argument("--out-prefix", required=True)
    args = p.parse_args(argv)

    if args.nanc < 1:
        raise SystemExit("--nanc must be >= 1")
    raw_prefix = f"{args.out_prefix}.upstream"
    run_upstream(args.create_model_script, args.nanc, args.panels, raw_prefix, args.seed)
    raw_model = Path(f"{raw_prefix}.model")
    if not raw_model.is_file():
        raise SystemExit(f"upstream script did not write {raw_model}")

    corr = read_correlations(raw_model)
    model = parse_model(raw_model)
    order, names, label_rows = canonical_labels(model)
    out_model = Path(f"{args.out_prefix}.model")
    write_model(relabel(model, order, names), out_model, preserve_raw=False)

    min_ac = corr["min_autocorr"]
    passed = True
    reason = ""
    if args.min_autocorr > 0:
        if min_ac is None:
            passed, reason = False, "upstream model has no min auto-correlation line"
        elif min_ac < args.min_autocorr:
            passed, reason = False, (
                f"min lag-1 autocorrelation {min_ac:.3f} < {args.min_autocorr} "
                f"(nanc={args.nanc} likely spurious)"
            )
    dominant = [r["dominant_panel"] for r in label_rows]
    summary = {
        "nanc": args.nanc,
        "seed": args.seed,
        "panels_file": str(args.panels),
        "model": str(out_model),
        "upstream_model": str(raw_model),
        **corr,
        "min_autocorr_gate": args.min_autocorr,
        "passed": passed,
        "reason": reason,
        "duplicate_dominant_panels": sorted({d for d in dominant if dominant.count(d) > 1}),
        "ancestries": label_rows,
    }
    Path(f"{args.out_prefix}.summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    with Path(f"{args.out_prefix}.labels.tsv").open("w") as fh:
        fh.write("index\tname\tupstream_label\tdominant_panel\tdominant_weight\n")
        for r in label_rows:
            fh.write(
                f"{r['index']}\t{r['name']}\t{r['upstream_label']}\t"
                f"{r['dominant_panel']}\t{r['dominant_weight']}\n"
            )
    print(json.dumps({k: summary[k] for k in ("nanc", "min_autocorr", "max_crosscorr", "passed")}),
          file=sys.stderr)
    if not passed:
        raise SystemExit(f"FLARE2 model gate failed: {reason}")
    print("wrote", out_model)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
