#!/usr/bin/env python3
"""
Integrate Phase callset partitions into one site table.

Rules:
  - Concatenate main (+ optional bnd/large)
  - If a site key appears in both main and large, keep large and mark/drop main
  - Site key: chrom, pos, ref, alt (fallback: id)
"""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from sv_site_utils import iter_sites, open_text, site_row, write_site_header  # noqa: E402


def site_key(s: dict) -> tuple:
    return (s["chrom"], int(s["pos"]), s["ref"], s["alt"])


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--main", required=True)
    p.add_argument("--bnd")
    p.add_argument("--large")
    p.add_argument("--out", required=True)
    p.add_argument("--manifest", required=True)
    args = p.parse_args()

    large_keys: set[tuple] = set()
    n_large = 0
    if args.large:
        for s in iter_sites(args.large):
            large_keys.add(site_key(s))
            n_large += 1

    n_main = 0
    n_bnd = 0
    n_suppressed = 0
    n_output = 0

    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile("w", delete=False, encoding="utf-8") as tmp:
        tmp_path = tmp.name
        if args.large or args.bnd:
            for s in iter_sites(args.main):
                n_main += 1
                s["source_vcf"] = "main"
                if site_key(s) in large_keys:
                    n_suppressed += 1
                    continue
                s["suppressed_main_duplicate"] = False
                tmp.write(site_row(s))
                n_output += 1
        else:
            for s in iter_sites(args.main):
                n_main += 1
                s["source_vcf"] = "main"
                s["suppressed_main_duplicate"] = False
                tmp.write(site_row(s))
                n_output += 1
        if args.bnd:
            for s in iter_sites(args.bnd):
                n_bnd += 1
                s["source_vcf"] = "bnd"
                s["suppressed_main_duplicate"] = False
                tmp.write(site_row(s))
                n_output += 1
        if args.large:
            for s in iter_sites(args.large):
                s["source_vcf"] = "large"
                s["suppressed_main_duplicate"] = False
                tmp.write(site_row(s))
                n_output += 1

    sorted_path = tmp_path + ".sorted"
    env = os.environ.copy()
    env["LC_ALL"] = "C"
    subprocess.check_call(
        [
            "sort",
            "-s",
            "-t",
            "\t",
            "-k1,1",
            "-k2,2n",
            "-k10,10",
            "-k4,4",
            "-o",
            sorted_path,
            tmp_path,
        ],
        env=env,
    )
    os.unlink(tmp_path)

    with open_text(args.out, "wt") as out, open(sorted_path, "rt", encoding="utf-8") as body:
        write_site_header(out)
        for line in body:
            out.write(line)
    os.unlink(sorted_path)

    manifest = {
        "n_main_input": n_main,
        "n_bnd_input": n_bnd,
        "n_large_input": n_large,
        "n_main_suppressed_as_large_duplicate": n_suppressed,
        "n_output_sites": n_output,
        "overlap_policy": "large_over_main",
    }
    Path(args.manifest).write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
