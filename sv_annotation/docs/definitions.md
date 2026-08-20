# Manuscript definitions (locked)

These definitions drive both the Phase 1 vs Phase 2 site-count table and the
Ebert-style cumulative discovery figures.

## Site counts

- Counts are **total unique sites** in the cohort callset(s), not median SVs per genome.
- Phase 2 mid-pass vs high-pass is **not** used for site totals; all Phase 2 samples
  entered one combined callset.
- Phase 1 and Phase 2 discovery curves are plotted as **separate panels** (not overlaid).

## Callset partitions

Callsets are **pre-partitioned by size** at the source (not split in this pipeline):

| Partition | Phase 1 | Phase 2 | Size / type |
|-----------|---------|---------|-------------|
| `main` | yes | yes | Resolved DEL/DUP/INS/INV, **20 bp ≤ \|SVLEN\| ≤ 10 kb** |
| `bnd` | no | yes | Breakends |
| `large` | no | yes | Resolved events **> 10 kb** (ultralong partition) |

If a site appears in both `main` and `large`, count it under **`large` only**.

Because `main` is bounded at 10 kb, **CADD-SV on `main` is effectively scoring events under 10 kb**
(DEL/DUP/INS/INV only; CADD-SV itself requires \|SVLEN\| ≥ 50 bp). The `large` and `bnd`
partitions are not CADD-scored.

## Classes reported

- DEL, DUP, INS, INV from `main` (and size-eligible large events as their native types)
- Breakends from `bnd` (reported separately; no SVLEN dual cutoffs)
- Large events from `large` (>10 kb partition)
- Do **not** report CPX or MCNV

## Size cutoffs (main resolved classes)

- `ge50`: \|SVLEN\| ≥ 50 bp (canonical SV)
- `ge20`: \|SVLEN\| ≥ 20 bp (broader LR sensitivity)
- Indels <20 bp belong in the SNV/indel table block, not here

## Repetitive region (union)

A site is **repetitive** if it overlaps the union of these GRCh38 tracks
(configurable via `configs/repetitive_beds.json`):

1. UCSC `rmsk` (RepeatMasker)
2. UCSC `simpleRepeat`
3. UCSC `genomicSuperDups` (segmental duplications)

Non-repetitive = complement of that union. Underlying per-track hit flags are
retained for QC even though figures use the binary `REGION_CLASS`.

## Frequency strata (Ebert et al. 2021)

Computed within the plotted cohort after the documented no-call policy:

| Label | Definition |
|-------|------------|
| `singleton` | AC = 1 |
| `polymorphic` | AC ≥ 2 and AF < 0.5 |
| `major` | AF ≥ 0.5 and not present in every sample |
| `shared` | present (alt carrier) in all samples |

## CADD-SV pathogenicity bins

CADD-SV PHRED scores are generated in-pipeline. For discovery facets we use:

| Bin | PHRED |
|-----|-------|
| `low` | < 10 |
| `mid` | ≥ 10 and < 20 |
| `high` | ≥ 20 |
| `unscored` | missing (`bnd`; `large`; unsupported class; or main sites with \|SVLEN\| < 50 bp) |

There is no universal pathogenicity threshold for CADD-SV; these bins are for
stratified discovery plots, not hard filters.

## Discovery sample order

Samples are ordered by ancestry label (non-African groups first, then African),
matching the Ebert-style layout. Exact order key is written to the run manifest.

Default ancestry sort key (configurable):

1. Non-African: `eas`, `amr`, `eur`, `sas`, `oth` / unknown
2. African: `afr`

## No-call policy (AF / carriers)

- Default: **no-calls do not contribute alleles** (`bcftools +fill-tags` / query
  using called genotypes only). Missing GT is not converted to 0/0.
- A site is a carrier for discovery if it has at least one alt allele called
  (GT with alt). Hom-ref and no-call are non-carriers.
