#!/usr/bin/env python3
"""Build the inputs JSON for felix/wdl/FelixPilotStep2Resume.wdl from an earlier FelixPilot run.

The resume workflow restarts at Step2 using outputs the earlier run already produced, so a Step2
fix can be tested without redoing the GRM build. Point --run-root at the workflow directory of
that run (the one holding call-MakeGRM/, call-Pack/ and call-Null/), e.g.

  gs://fc-secure-.../submissions/<submission-id>/FelixPilot/<workflow-id>

Needs gsutil on PATH (the Terra VM has it). It lists the three call directories, picks the
output files, and fails rather than guessing if a file is missing or appears more than once
(for example a retried shard with attempt-2/ or cacheCopy/ copies).

  python3 felix/scripts/make_resume_inputs.py \\
      --run-root gs://.../FelixPilot/<workflow-id> \\
      --phenotypes-file gs://.../tractor_mix_pilot/selected_phenotypes.txt \\
      --base-config felix/configs/felix.pilot.inputs.limited.json \\
      --limit 2 --out resume.inputs.json

--limit N keeps the first N phenotypes, for a quick test. Shard i of call-Null belongs to line i
of the phenotypes file, which is how FelixPilot.wdl scatters.
"""

from __future__ import annotations

import argparse
import json
import posixpath
import re
import subprocess
import sys
from pathlib import Path

PREFIX = "FelixPilotStep2Resume"
# Inputs copied from the FelixPilot config when present (same names in both workflows).
COPIED = ["pheno_cov", "run_felix_step2_script", "summarize_script", "chrom", "num_ancs", "min_mac",
          "pvalcutoff_of_haplotype", "docker"]


def gsutil_ls(gsutil: str, prefix: str) -> list[str]:
    proc = subprocess.run([gsutil, "ls", prefix.rstrip("/") + "/**"], capture_output=True, text=True)
    if proc.returncode != 0:
        raise SystemExit(f"gsutil ls {prefix}/** failed: {proc.stderr.strip() or proc.stdout.strip()}")
    return [
        ln.strip() for ln in proc.stdout.splitlines()
        if ln.strip() and not ln.strip().endswith(":") and not ln.strip().endswith("/")
    ]


def read_lines(gsutil: str, path: str) -> list[str]:
    if path.startswith("gs://"):
        proc = subprocess.run([gsutil, "cat", path], capture_output=True, text=True)
        if proc.returncode != 0:
            raise SystemExit(f"gsutil cat {path} failed: {proc.stderr.strip()}")
        text = proc.stdout
    else:
        text = Path(path).read_text()
    return [ln.strip() for ln in text.splitlines() if ln.strip()]


def pick(urls: list[str], name: str, where: str) -> str:
    hits = [u for u in urls if posixpath.basename(u) == name]
    if len(hits) != 1:
        detail = "".join(f"\n    {u}" for u in hits) if hits else ""
        raise SystemExit(f"{where}: expected exactly one {name}, found {len(hits)}{detail}")
    return hits[0]


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--run-root", required=True, help="gs://.../submissions/<id>/FelixPilot/<workflow-id>")
    ap.add_argument("--phenotypes-file", required=True, help="the selected_phenotypes.txt the run used (local or gs://)")
    ap.add_argument("--base-config", required=True, help="FelixPilot inputs JSON to copy pheno_cov, scripts, etc. from")
    ap.add_argument("--limit", type=int, default=0, help="keep only the first N phenotypes (0 = all)")
    ap.add_argument("--out", required=True)
    ap.add_argument("--gsutil", default="gsutil")
    args = ap.parse_args()

    root = args.run_root.rstrip("/")
    phenotypes = read_lines(args.gsutil, args.phenotypes_file)
    if args.limit > 0:
        phenotypes = phenotypes[: args.limit]
    if not phenotypes:
        raise SystemExit("no phenotypes")

    pack = gsutil_ls(args.gsutil, f"{root}/call-Pack")
    grm = gsutil_ls(args.gsutil, f"{root}/call-MakeGRM")
    null = gsutil_ls(args.gsutil, f"{root}/call-Null")

    by_shard: dict[int, list[str]] = {}
    for u in null:
        m = re.search(r"/call-Null/shard-(\d+)/", u)
        if m:
            by_shard.setdefault(int(m.group(1)), []).append(u)

    rdas, vrs, sus = [], [], []
    for i, p in enumerate(phenotypes):
        if i not in by_shard:
            raise SystemExit(f"call-Null/shard-{i}/ has no files, so {p} has no null model in this run")
        where = f"call-Null/shard-{i} ({p})"
        rdas.append(pick(by_shard[i], f"{p}.null.rda", where))
        vrs.append(pick(by_shard[i], f"{p}.varianceRatio.txt", where))
        sus.append(pick(by_shard[i], f"{p}.samples.txt", where))

    base = json.loads(Path(args.base_config).read_text())
    inputs = {
        f"{PREFIX}.packed_tar": pick(pack, "felixla_packed.tar.gz", "call-Pack"),
        f"{PREFIX}.sparse_grm_mtx": pick(grm, "saige_sparseGRM.mtx", "call-MakeGRM"),
        f"{PREFIX}.sparse_grm_sample_ids": pick(grm, "saige_sparseGRM.sampleIDs.txt", "call-MakeGRM"),
        f"{PREFIX}.phenotypes": phenotypes,
        f"{PREFIX}.null_rdas": rdas,
        f"{PREFIX}.variance_ratios": vrs,
        f"{PREFIX}.samples_used": sus,
    }
    for key in COPIED:
        src = f"FelixPilot.{key}"
        if src in base:
            inputs[f"{PREFIX}.{key}"] = base[src]
    for required in ("pheno_cov", "run_felix_step2_script", "summarize_script"):
        if f"{PREFIX}.{required}" not in inputs:
            raise SystemExit(f"{args.base_config} has no FelixPilot.{required}")

    Path(args.out).write_text(json.dumps(inputs, indent=2) + "\n")
    print(f"wrote {args.out}: {len(phenotypes)} phenotype(s), first {phenotypes[0]} -> shard-0, last {phenotypes[-1]} -> shard-{len(phenotypes) - 1}", file=sys.stderr)
    print(f"  packed_tar  {inputs[f'{PREFIX}.packed_tar']}", file=sys.stderr)


if __name__ == "__main__":
    main()
