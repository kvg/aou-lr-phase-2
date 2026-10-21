#!/usr/bin/env python3
"""Build the Terra entity table ``TrgtHaplotypePhase`` runs on.

The workflow needs four per-sample paths: the haplotagged BAM and its index,
and the TRGT VCF and its index. They live in two different places — the BAMs in
a Terra table (``aou2_v1_phased_bams``), the TRGT calls in
``tractor_mix/resources/trgt_table.txt`` — and neither is a superset of the
other. This joins them and emits one row per sample that has **both**.

    make_trgt_phase_table.py --bams-table aou2_v1_phased_bams.tsv \\
        --trgt-table tractor_mix/resources/trgt_table.txt \\
        --out trgt_phase.tsv

Everything else the workflow takes is the same for every sample — the phased
per-chromosome callsets are multi-sample files, and the reference, catalog,
thresholds and image are run-level — so they belong in the method config, not
in this table. A column that is identical in 12,000 rows is a column that can
disagree with itself.

Rows are emitted only for the intersection, and the three populations are
reported separately: a sample with a BAM but no TRGT call cannot be phased, and
a TRGT call with no BAM has no reads to phase it from. Neither is an error, but
both are a count you want to see rather than discover as a failed shard.

One row is one SAMPLE, not one person. ``trgt_table.txt`` holds 12,421 rows
over 12,393 distinct person ids, and the 28 repeats are not resubmissions: they
are people with a genuine sample from each of two sequencing sites, genotyped
once per sample. TRGT was submitted one batch per site, which is why the
submission id in each URI identifies the site:

    terra-d45cc5ca  BI (10,023)      terra-e7fd2863  BCM (762)
    terra-f6a059a8  HA (1,234)       terra-f6671367  UW (402)

That mapping is derived from ``--sample-list`` by majority vote, not hardcoded,
so adding a batch does not require editing this file. Keying by sample makes
each of those 28 pairs two rows rather than a collision to resolve, which is the
point: deduplication to one genome per person is an analysis decision and stays
in FELIX's ``analysis_samples.txt``, and the pairs then serve as replicate
controls on repeat-call reproducibility.

Two different ids therefore appear in the output. The **entity id** is
``<person>_<SITE>``, matching the ``sample-hifi-hg38-all-cohorts`` convention.
The **sample column** is the bare person id, because that is the name
``bcftools view -s`` has to find in the joint phased callset — which carries one
column per person, not per sample. For a cross-site pair both rows therefore
subset the same phased genotypes, so the comparison between them isolates read
and genotyping variation with the phase frame held fixed.

Two TRGT URIs landing on the same ``<person>_<SITE>`` key is a hard error
rather than a policy choice: the site mapping is meant to make that impossible,
so it happening means an assumption here is wrong and picking one silently would
hide it.

Nothing here checks that the objects exist. It is a join over path strings.
"""

from __future__ import annotations

import argparse
import collections
import gzip
import re
import sys
from pathlib import Path
from typing import Optional, TextIO

OUT_COLUMNS = ("sample", "bam", "bam_index", "trgt_vcf", "trgt_vcf_index")
SUBMISSION_RE = re.compile(r"/(terra-[0-9a-f]+)/")


def open_text(path: str, mode: str = "rt") -> TextIO:
    if str(path).endswith(".gz"):
        return gzip.open(path, mode)
    return open(path, mode)


def read_trgt_table(path: str) -> list[tuple[str, str]]:
    """[(person_id, VCF URI)] in file order, duplicates preserved.

    Duplicates are not resolved here. A person id appearing twice means two
    samples from two sites, and which row is which is decided by the site
    mapping, not by file order.
    """
    rows: list[tuple[str, str]] = []
    seen: set[tuple[str, str]] = set()
    with open_text(path) as fh:
        for n, line in enumerate(fh, start=1):
            line = line.rstrip("\n")
            if not line:
                continue
            fields = line.split("\t")
            if len(fields) < 2:
                raise SystemExit(f"{path}:{n}: expected 'id<TAB>uri', got {line!r}")
            key = (fields[0].strip(), fields[1].strip())
            if key not in seen:
                seen.add(key)
                rows.append(key)
    return rows


def read_sample_list(path: str) -> dict[str, list[str]]:
    """person_id -> sorted sites, from a list of ``<person>_<SITE>`` ids.

    Accepts a Terra set-membership export (two columns, ids in the second) or a
    bare one-id-per-line list.
    """
    sites: dict[str, set[str]] = {}
    with open_text(path) as fh:
        for line in fh:
            line = line.rstrip("\n").strip()
            if not line:
                continue
            fields = [f.strip() for f in line.split("\t")]
            ident = fields[-1]
            if ident.startswith(("entity:", "membership:")) or "_" not in ident:
                continue
            # rsplit, not split: the SITE is the suffix, so a person id that
            # itself contains an underscore must still parse correctly.
            person, site = ident.rsplit("_", 1)
            sites.setdefault(person, set()).add(site)
    if not sites:
        raise SystemExit(f"{path}: no '<person>_<SITE>' ids found")
    return {p: sorted(s) for p, s in sites.items()}


def submission_sites(rows: list[tuple[str, str]],
                     sites: dict[str, list[str]],
                     min_purity: float) -> tuple[dict[str, str], list[str]]:
    """Map each TRGT submission to a sequencing site by majority vote.

    Only people with exactly one site vote, since an ambiguous person is what the
    mapping is for. Purity is reported and enforced: TRGT was batched per site,
    so a submission whose voters disagree means that is no longer true and the
    mapping must not be guessed.
    """
    votes: dict[str, collections.Counter] = collections.defaultdict(collections.Counter)
    for person, uri in rows:
        m = SUBMISSION_RE.search(uri)
        if m is None:
            continue
        s = sites.get(person)
        if s is not None and len(s) == 1:
            votes[m.group(1)][s[0]] += 1
    mapping, notes = {}, []
    for subm, counter in sorted(votes.items()):
        site, n = counter.most_common(1)[0]
        total = sum(counter.values())
        purity = n / total if total else 0.0
        mapping[subm] = site
        notes.append(f"  {subm}  ->  {site:<4} {n:,}/{total:,} unambiguous voters"
                     f" ({purity:.4f})")
        if purity < min_purity:
            raise SystemExit(
                f"submission {subm} does not map cleanly to one site: "
                f"{dict(counter)}. TRGT was batched per site; if that has changed, "
                f"pass --min-site-purity to override deliberately.")
    return mapping, notes


def key_by_sample(rows: list[tuple[str, str]],
                  sites: dict[str, list[str]],
                  mapping: dict[str, str]) -> tuple[dict[str, tuple[str, str]], dict[str, int]]:
    """``<person>_<SITE>`` -> (person_id, URI), plus a tally of what could not be keyed.

    Resolved per person, not per row, so a row whose submission has no site can
    still be placed by elimination: if a person's other rows account for all but
    one of their sites, the remaining row takes the remaining site. That matters
    because ``submission_sites`` votes only with unambiguous people, so a batch
    consisting solely of cross-site people would otherwise map to nothing.

    A person with one site and one row takes it outright; the submission is only
    consulted when there is something to disambiguate.
    """
    by_person: dict[str, list[str]] = collections.defaultdict(list)
    for person, uri in rows:
        by_person[person].append(uri)

    out: dict[str, tuple[str, str]] = {}
    tally: collections.Counter = collections.Counter()
    for person, uris in by_person.items():
        person_sites = sites.get(person)
        if not person_sites:
            tally["person_not_in_sample_list"] += len(uris)
            continue
        if len(uris) == 1 and len(person_sites) == 1:
            out[f"{person}_{person_sites[0]}"] = (person, uris[0])
            continue

        assigned: dict[str, str] = {}
        leftover: list[str] = []
        taken: set[str] = set()
        for uri in uris:
            m = SUBMISSION_RE.search(uri)
            site = mapping.get(m.group(1)) if m else None
            if site is not None and site not in person_sites:
                tally["submission_site_not_a_site_of_person"] += 1
                continue
            if site is not None and site not in taken:
                assigned[uri] = site
                taken.add(site)
            else:
                leftover.append(uri)

        remaining = [s for s in person_sites if s not in taken]
        if len(leftover) == len(remaining):
            for uri, site in zip(leftover, remaining):
                assigned[uri] = site
        elif leftover:
            tally["could_not_disambiguate"] += len(leftover)

        for uri, site in assigned.items():
            key = f"{person}_{site}"
            if key in out:
                tally["key_collision"] += 1
                continue
            out[key] = (person, uri)
    return out, tally


def read_bams_table(path: str, bam_column: str, bai_column: str) -> dict[str, dict[str, str]]:
    """entity id -> {bam, bam_index}, from a Terra table export.

    The id is the first column, whose header Terra writes as ``entity:<type>_id``.
    Rows whose BAM or index cell is empty are dropped and counted by the caller:
    an empty cell in Terra is a real state, not a missing file.
    """
    rows: dict[str, dict[str, str]] = {}
    with open_text(path) as fh:
        header: Optional[list[str]] = None
        for n, line in enumerate(fh, start=1):
            line = line.rstrip("\n")
            if not line:
                continue
            fields = line.split("\t")
            if header is None:
                header = [f.strip() for f in fields]
                if not header[0].startswith("entity:"):
                    raise SystemExit(
                        f"{path}: first column header is {header[0]!r}, expected 'entity:<type>_id'. "
                        "Export the table from Terra rather than hand-building it.")
                for col in (bam_column, bai_column):
                    if col not in header:
                        raise SystemExit(f"{path}: no column {col!r}. Columns: {', '.join(header)}")
                continue
            row = dict(zip(header, fields))
            sample = fields[0].strip()
            if not sample:
                continue
            rows[sample] = {
                "bam": row.get(bam_column, "").strip(),
                "bam_index": row.get(bai_column, "").strip(),
            }
    if header is None:
        raise SystemExit(f"{path}: empty")
    return rows


def run(args: argparse.Namespace) -> int:
    bams = read_bams_table(args.bams_table, args.bam_column, args.bai_column)
    rows = read_trgt_table(args.trgt_table)
    sites = read_sample_list(args.sample_list)
    mapping, notes = submission_sites(rows, sites, args.min_site_purity)
    trgt, tally = key_by_sample(rows, sites, mapping)

    if tally.get("key_collision"):
        raise SystemExit(
            f"{tally['key_collision']} TRGT rows collided on a <person>_<SITE> key. "
            "The site mapping is supposed to make that impossible, so this means an "
            "assumption in this script is wrong; it is not resolved silently.")

    incomplete = sorted(s for s, r in bams.items() if not r["bam"] or not r["bam_index"])
    usable_bams = {s: r for s, r in bams.items() if r["bam"] and r["bam_index"]}

    both = set(usable_bams) & set(trgt)
    bam_only = sorted(set(usable_bams) - set(trgt))
    trgt_only = sorted(set(trgt) - set(usable_bams))

    # Emit in the BAM table's order so the new table lines up with the old one.
    order = [s for s in bams if s in both]
    with open(args.out, "wt") as fh:
        fh.write(f"entity:{args.entity_type}_id\t" + "\t".join(OUT_COLUMNS) + "\n")
        for key in order:
            person, uri = trgt[key]
            fh.write("\t".join((
                key, person, usable_bams[key]["bam"], usable_bams[key]["bam_index"],
                uri, uri + args.index_suffix)) + "\n")

    print("TRGT submission -> sequencing site (majority vote, not hardcoded):")
    print("\n".join(notes))
    print(f"\nsample list:                    {len(sites):,} people,"
          f" {sum(len(v) for v in sites.values()):,} samples")
    multi = {p_: v for p_, v in sites.items() if len(v) > 1}
    print(f"  people at more than one site:  {len(multi):,}"
          + (f"  (e.g. {', '.join(sorted(multi)[:3])})" if multi else ""))
    print(f"TRGT rows:                      {len(rows):,}")
    print(f"  keyed to a sample:            {len(trgt):,}")
    for reason, n in sorted(tally.items()):
        print(f"  dropped, {reason}: {n:,}")
    print(f"BAM table rows:                 {len(bams):,}")
    print(f"  with both BAM and index:      {len(usable_bams):,}")
    print(f"  missing a BAM or index cell:  {len(incomplete):,}"
          + (f"  (e.g. {', '.join(incomplete[:5])})" if incomplete else ""))
    print()
    print(f"ROWS WRITTEN (both sources):    {len(order):,}")
    print(f"  BAM but no TRGT call:         {len(bam_only):,}"
          + (f"  (e.g. {', '.join(bam_only[:5])})" if bam_only else ""))
    print(f"  TRGT call but no usable BAM:  {len(trgt_only):,}"
          + (f"  (e.g. {', '.join(trgt_only[:5])})" if trgt_only else ""))

    written_people = {trgt[k][0] for k in order}
    pairs = sorted(p_ for p_ in multi
                   if sum(1 for k in order if trgt[k][0] == p_) > 1)
    print(f"\npeople represented:             {len(written_people):,}")
    print(f"  with both site samples kept:  {len(pairs):,}"
          + (f"  (e.g. {', '.join(pairs[:3])})" if pairs else ""))
    if pairs:
        print("  -> these are replicate pairs, not duplicates; one genome per person is "
              "FELIX's analysis_samples.txt to enforce, not this table's.", file=sys.stderr)
    print(f"\nwrote {args.out}")
    if bam_only or trgt_only:
        print("Samples outside the intersection are absent from this table by construction, so "
              "they cannot be submitted by accident.", file=sys.stderr)
    return 0


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--bams-table", required=True,
                   help="Terra export of the haplotagged-BAM table (aou2_v1_phased_bams)")
    p.add_argument("--trgt-table", required=True,
                   help="person_id<TAB>TRGT VCF URI, no header")
    p.add_argument("--entity-type", default="trgt_phase",
                   help="Entity type for the new table, used in the 'entity:<type>_id' header")
    p.add_argument("--bam-column", default="gcs_bam_path")
    p.add_argument("--bai-column", default="gcs_bai_path")
    p.add_argument("--index-suffix", default=".tbi")
    p.add_argument(
        "--sample-list", required=True,
        help="List of '<person>_<SITE>' sample ids -- a Terra set-membership export "
             "(e.g. sample-hifi-hg38-all-cohorts) or one id per line. Supplies the "
             "site for each person and the submission->site vote.")
    p.add_argument(
        "--min-site-purity", type=float, default=0.95,
        help="Floor on the fraction of a TRGT submission's unambiguous people that "
             "agree on one site (default 0.95). TRGT was batched per site; a low "
             "value means that is no longer true.")
    p.add_argument("--out", required=True)
    return p


def main(argv: Optional[list[str]] = None) -> int:
    return run(build_parser().parse_args(argv))


if __name__ == "__main__":
    raise SystemExit(main())
