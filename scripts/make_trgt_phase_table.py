#!/usr/bin/env python3
"""Build the Terra entity table ``TrgtHaplotypePhase`` runs on.

The workflow needs four per-sample paths: the haplotagged BAM and its index,
and the TRGT VCF and its index. They live in two different places — the BAMs in
a Terra table (``aou2_v1_phased_bams``), the TRGT calls in
``tractor_mix/resources/trgt_table.txt`` — and neither is a superset of the
other. This joins them and emits one row per sample that has **both**.

    make_trgt_phase_table.py --bams-table aou2_v1_phased_bams.tsv \\
        --trgt-table tractor_mix/resources/trgt_table.txt \\
        --sample-list sample-hifi-hg38-all-cohorts_set_membership.tsv \\
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

One row is one PERSON, because there is one haplotagged BAM per person. The
``aou2_v1_phased_bams`` export carries 12,261 numeric participant ids (exactly
the releasable cohort) plus 293 non-numeric 1000 Genomes and GIAB control ids,
one row each, and all 28 cross-site people appear once. The choice of which
sample represents a person was therefore already made upstream, when the
haplotagged BAMs were produced; this script does not get to revisit it.

``trgt_table.txt`` nonetheless holds 12,421 rows over 12,393 people. The 28
repeats are not resubmissions: they are people with a genuine sample from each
of two sequencing sites, genotyped once per sample. TRGT was submitted one batch
per site, so the submission id in each URI identifies the site:

    terra-d45cc5ca  BI (10,023)      terra-e7fd2863  BCM (762)
    terra-f6a059a8  HA (1,234)       terra-f6671367  UW (402)

That mapping is derived from ``--sample-list`` by majority vote rather than
hardcoded, and is used only to break those 28 ties in a stated order
(``--prefer-site``) instead of by file position. The 28 pairs are BI+UW 19,
BI+HA 3, BI+BCM 3 and HA+UW 3, so a BI-first order resolves all of them.

The stakes on that tie-break are low and worth being explicit about: the
workflow RE-GENOTYPES from the BAM, and the TRGT VCF is used only to emit a
comparison slice of the production call. So the choice affects what the run is
compared against, not what it computes. It is made deterministic because an
arbitrary one would be invisible, not because it changes a result.

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


def choose_uri(rows: list[tuple[str, str]],
               sites: dict[str, list[str]],
               mapping: dict[str, str],
               prefer: list[str]) -> tuple[dict[str, str], dict[str, str], dict[str, int]]:
    """person_id -> one VCF URI, plus the site chosen for each resolved tie.

    A person with one URI keeps it. A person with several is resolved by site
    preference: the URI whose submission maps to the earliest site in ``prefer``
    wins. Anything unresolvable is tallied and dropped rather than guessed.
    """
    by_person: dict[str, list[str]] = collections.defaultdict(list)
    for person, uri in rows:
        by_person[person].append(uri)

    chosen: dict[str, str] = {}
    chosen_site: dict[str, str] = {}
    tally: collections.Counter = collections.Counter()
    rank = {s: i for i, s in enumerate(prefer)}
    for person, uris in by_person.items():
        if len(uris) == 1:
            chosen[person] = uris[0]
            continue
        cands = []
        for uri in uris:
            m = SUBMISSION_RE.search(uri)
            site = mapping.get(m.group(1)) if m else None
            if site is None:
                continue
            cands.append((rank.get(site, len(rank)), site, uri))
        if not cands:
            # No site for any candidate. Dropping the person would be
            # disproportionate -- they have a BAM and a usable TRGT call, and the
            # tie only decides which production call the run is compared against.
            # So fall back to a deterministic order over the URIs themselves,
            # which is reproducible rather than dependent on file position, and
            # say so loudly in the tally.
            chosen[person] = sorted(uris)[0]
            tally["tie_resolved_by_uri_order"] += 1
            continue
        cands.sort()
        if len(cands) < len(uris):
            tally["tie_partially_sited"] += 1
        _, site, uri = cands[0]
        chosen[person] = uri
        chosen_site[person] = site
        tally["tie_resolved_by_site"] += 1
    return chosen, chosen_site, tally


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

    by_person: dict[str, list[str]] = collections.defaultdict(list)
    for person, uri in rows:
        by_person[person].append(uri)
    ties = sorted(p_ for p_, u in by_person.items() if len(u) > 1)

    prefer = [s.strip() for s in args.prefer_site.split(",") if s.strip()]
    notes: list[str] = []
    if args.sample_list:
        sites = read_sample_list(args.sample_list)
        mapping, notes = submission_sites(rows, sites, args.min_site_purity)
    elif ties:
        raise SystemExit(
            f"{len(ties)} people have more than one TRGT VCF (e.g. {', '.join(ties[:3])}) "
            "and --sample-list was not given, so the site of each cannot be determined. "
            "Pass --sample-list so the tie is broken by a stated site order rather than "
            "by file position.")
    else:
        sites, mapping = {}, {}
    trgt, chosen_site, tally = choose_uri(rows, sites, mapping, prefer)

    participants = {s for s in bams if s.isdigit()}
    controls = sorted(set(bams) - participants)
    usable = {s: r for s, r in bams.items()
              if s in participants and r["bam"] and r["bam_index"]}
    incomplete = sorted(s for s in participants if not bams[s]["bam"] or not bams[s]["bam_index"])

    both = set(usable) & set(trgt)
    bam_only = sorted(set(usable) - set(trgt))
    trgt_only = sorted(set(trgt) - set(usable))

    order = [s for s in bams if s in both]
    with open(args.out, "wt") as fh:
        fh.write(f"entity:{args.entity_type}_id\t" + "\t".join(OUT_COLUMNS) + "\n")
        for person in order:
            fh.write("\t".join((
                person, person, usable[person]["bam"], usable[person]["bam_index"],
                trgt[person], trgt[person] + args.index_suffix)) + "\n")

    if notes:
        print("TRGT submission -> sequencing site (majority vote, not hardcoded):")
        print("\n".join(notes))
        print()
    print(f"BAM table rows:                 {len(bams):,}")
    print(f"  participant ids (numeric):    {len(participants):,}")
    print(f"  control ids (1000G/GIAB):     {len(controls):,}"
          + (f"  (e.g. {', '.join(controls[:3])})" if controls else ""))
    print(f"  missing a BAM or index cell:  {len(incomplete):,}"
          + (f"  (e.g. {', '.join(incomplete[:5])})" if incomplete else ""))
    print(f"TRGT rows:                      {len(rows):,}")
    print(f"  distinct people:              {len(by_person):,}")
    print(f"  people with >1 TRGT VCF:      {len(ties):,}")
    for reason, n in sorted(tally.items()):
        print(f"    {reason}: {n:,}")
    if chosen_site:
        picked = collections.Counter(chosen_site.values())
        print(f"    site chosen for ties:       "
              f"{', '.join(f'{s} {n}' for s, n in picked.most_common())}"
              f"   (preference {' > '.join(prefer)})")
    print()
    print(f"ROWS WRITTEN (both sources):    {len(order):,}")
    print(f"  BAM but no TRGT call:         {len(bam_only):,}"
          + (f"  (e.g. {', '.join(bam_only[:5])})" if bam_only else ""))
    print(f"  TRGT call but no BAM:         {len(trgt_only):,}"
          + (f"  (e.g. {', '.join(trgt_only[:5])})" if trgt_only else ""))
    print(f"\nwrote {args.out}")
    if trgt_only:
        print("A TRGT call with no haplotagged BAM has no reads to phase from, so those "
              "people are absent by construction and cannot be submitted by accident.",
              file=sys.stderr)
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
        "--sample-list", default=None,
        help="List of '<person>_<SITE>' sample ids -- a Terra set-membership export "
             "(e.g. sample-hifi-hg38-all-cohorts) or one id per line. Supplies the "
             "site for each person and the submission->site vote. Required only "
             "when some person has more than one TRGT VCF.")
    p.add_argument(
        "--prefer-site", default="BI,HA,BCM,UW",
        help="Order in which sites win a tie when a person has a TRGT call from "
             "more than one (default BI,HA,BCM,UW). Affects only which production "
             "call the run is compared against, since re-genotyping is done from "
             "the BAM.")
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
