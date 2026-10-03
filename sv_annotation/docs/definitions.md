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

## Sequence-context strata (site counts)

Deletions and insertions in Table 2 report the **percent of sites** in each
exclusive class (US / RM / SD / SR), plus CMRG, separately for |SVLEN|
≥50 bp and |SVLEN| ≥20 bp (`deletions_context_pct_ge50` /
`deletions_context_pct_ge20`, and the same for insertions). Within each
row, US, RM, SD, and SR partition that size bin and sum to 100. CMRG is
independent, so the five displayed percentages do not sum to 100.

The label is the breakpoint rule from Xuefang Zhao's
`annotate_genomic_context.sh` (gatk-sv branch `xz_fixes_3`; Zhao et al.,
AJHG 2021). `region_class` is one of `US`, `RM`, `SD`, `SR`.

| Label | Definition |
|-------|------------|
| `US` | Neither breakpoint (and, for a long CNV, not the body) meets a repeat class |
| `RM` | RepeatMasker, unless SD or SR also applies |
| `SD` | Segmental duplication, unless SR also applies |
| `SR` | Simple repeat |
| `CMRG` | Reference span overlaps GIAB challenging medically relevant gene coordinates |

Either breakpoint can set the class. Breakpoints are the 0-based BED start
(VCF `POS` − 1) and the BED end (VCF `END`). Priority is SR, then SD, then RM.

`DEL`, `DUP`, and `CNV` with reference span **> 5 kb** ignore breakpoints. A
class applies only when the merged track covers more than half of the body,
with the same priority; otherwise the site is `US`. Insertions always use
breakpoints, including when `|SVLEN|` is large. If `end` was filled from
`SVLEN` (`end − pos = |SVLEN|`), both insertion breakpoints are the anchor
at `POS`.

`hit_rmsk`, `hit_simpleRepeat`, and `hit_genomicSuperDups` are the evidence
flags before the priority collapse, so more than one can be true.
`hit_cmrg` is span overlap and is not part of `region_class`. A site can be
`US` and CMRG, or `SR` and CMRG.

Discovery strata use this four-way `region_class`, not a repetitive/unique
union. Tracks are the merged GRCh38 BEDs from that script
(`configs/repetitive_beds.json`).

## Context tracks

Stage with `notebooks/terra/sv_01_stage_repeat_tracks.ipynb`:

1. `hg38.RM.sorted.merged.bed.gz` (RepeatMasker)
2. `hg38.SR.sorted.merged.bed.gz` (simple repeats)
3. `hg38.SD.sorted.merged.bed.gz` (segmental duplications)

These are the merged beds from gatk-sv `benchmark_scripts/input` on
`xz_fixes_3`, not the raw UCSC table dumps. CMRG stays the GIAB v1.00 gene
BED and is not one of the three context tracks.

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
