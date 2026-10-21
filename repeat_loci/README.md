# Repeat loci (Figure 2C)

Panel C counts TRGT allele lengths at TRExplorer v1.0.1 loci. `TrgtPlvi` ranks
every locus by PLVI; the three examples are chosen from that ranking once the
allele counts exist. `TrgtAlleleCounts` counts haplotypes at each distinct
length. Run it once on the Phase 1 TRGT list and once on the Phase 2 list,
with the same catalog. `panel_c_repeat_alleles.py` turns the two allele tables
and the PLVI table into the rarefaction curve and the example histograms.

`RepeatLocusDosage` sums integrated SNV/indel/SV records. That is not the
Panel C input. `TrgtLocusConcordance` on chr20 is already done.

`TrgtHaplotypePhase` is unrelated to Panel C: it is the prototype for re-basing
the FELIX repeat-dosage arm on TRGT lengths. See
[TRGT haplotype phasing](#trgt-haplotype-phasing) below.

| Step | Workflow | Script | Output |
|---|---|---|---|
| Allele tables | `wdl/TrgtAlleleCounts.wdl` (Phase 1 and Phase 2) | `trgt_allele_counts.py` | `<label>.alleles.tsv.gz` (locus × dosage × haplotype count) and `<label>.summary.json` |
| Example ranking | `wdl/TrgtPlvi.wdl` (Phase 2) | `trgt_plvi.py` | PLVI per TRGT locus, ranked within motif-length groups |
| Figure inputs | `panel_c_repeat_alleles.py` | same | Rarefaction curves, example-locus alleles, summary |

No output carries sample ids.

## TRGT allele counts

A haplotype's dosage is `AL - (END - POS + 1)`: the TRGT allele length minus
the reference span, in bp. Repeat units are that change divided by the first
motif's length. Two sequences of the same length are one allele, including a
length that is not a whole number of copies. A missing genotype counts as
missing, not as reference. A catalog locus absent from a VCF counts as two
missing haplotypes. A hemizygous call counts as one haplotype.

Both phases use the catalog BED from `trgt_plvi.py catalog` (one TRExplorer
VCF). Phase 2 reuses the VCF list `TrgtPlvi` is already running. Phase 1 is
`trgt.phase1.txt`: `sample_id` and the GCS path, tab-separated, no header,
same shape as `trgt_table.txt`. Example inputs are in `configs/`. Shards of
50 VCFs are localized.
The merge writes the allele table above (`locus_id` is the TRID).

Docker: `us-central1-docker.pkg.dev/broad-dsp-lrma/aou-lr/aou-sv-annotation:0.1.6`.

## Integrated dosage (not used for Panel C)

Every length-changing record inside a TRExplorer v1.0.1 locus (± `flank` bp)
contributes its signed length change to that haplotype's locus dosage, so
split or partially merged records at one repeat are summed into a single
allele. The same dosage is what a repeat-aware association test would use
later (`--out-vcf`, FORMAT `GT[:AN1:AN2]`, INFO `RU_DOSAGE`). The chr20
concordance workflow is `wdl/TrgtLocusConcordance.wdl`.

> **Association input, not the figure.** `aggregate_repeat_loci.py --out-vcf
> --with-ancestry` is what the FELIX repeat-dosage branch tests. The coding
> rules below are association conventions — see
> [`felix/REPEAT_DOSAGE.md`](../felix/REPEAT_DOSAGE.md) §6.1 before changing
> any of them.

## Rules

- **Sign.** Each ALT is signed by its own length: `len(ALT) - len(REF)` for
  sequence-resolved alleles, else `SVLEN` / `END` / `SVTYPE` for symbolic ones.
  INV and BND add 0.
- **Assignment.** A record goes to the catalog locus it overlaps most, within
  `flank` bp. Records with a change above `max_bp` are ignored.
- **Ambiguity.** A sample with two or more unphased heterozygous length
  changes at one locus is set missing (both haplotypes), because the changes
  cannot be assigned to haplotypes; `hap_unphased_ambiguous` counts them. One
  unphased het plus any homozygous records is resolved. The integrated inputs
  are unphased; `--ignore-phase` applies the same rule to a phased callset.
- **Size split.** Each input has a role: `small` keeps alleles with
  |change| < `split_bp` (50), `sv` keeps the rest, `all` keeps everything.
  An event called by both DeepVariant and the SV caller at the same size
  therefore counts once. Phase 1 SVs start at 50 bp and Phase 2 `v3_main` at
  20 bp, so 50 bp puts 20–49 bp events on the DeepVariant side in both phases.
  `alleles_out_of_band` counts the alleles zeroed this way.
- **Filters.** `apply_filters` (default `PASS,.`) keeps records with those
  FILTER values, as in the Figure 2A length spectrum.
- **Possible duplicates.** The same event can still be called at different
  sizes, e.g. 45 bp by DeepVariant and 52 bp by the SV caller. A haplotype
  carrying both a < 50 bp and a ≥ 50 bp change in the same direction is
  counted in `hap_possible_duplicate`, and `TrgtLocusConcordance` measures how
  often the sum overshoots TRGT.

## Inputs

| Phase | Small (`small`) | SV (`sv`) |
|---|---|---|
| 2 | DeepVariant + GLnexus shards, `GL_INTERVAL_set` in `AoU_DRC_LongReads_PhaseTwo_Storage` (`chrom_vcfs`) | `integrated_sv_callset/v3_main.bcf` (`genome_vcfs`) |
| 1 | `phase1/Jointcall_dvcf.g.vcf.bgz.QualFT40.vcf.gz` (`genome_vcfs`) | `scratch/kvg/phase1/concat_annotated.sens_07.bcf` (`genome_vcfs`) |

- `chrom_vcfs` are localized per task. `genome_vcfs` are cut per chromosome
  to length-changing records in `SliceChrom` (one streamed read; GCS tokens
  expire after about an hour, so whole-chromosome aggregation does not stream).
- The aggregator keeps samples present in every input; the notebook checks
  that each phase's two callsets share sample ids.
- Catalog BED: `trgt_plvi.py catalog --vcf <one TRGT VCF>` (chrom, start,
  end, TRID, motifs).

Example configs are in `configs/`. Docker:
`us-central1-docker.pkg.dev/broad-dsp-lrma/aou-lr/aou-sv-annotation:0.1.6`
(bcftools, python3, numpy).

## Tests

```bash
python -m pytest scripts/test_aggregate_repeat_loci.py scripts/test_trgt_plvi.py \
  scripts/test_trgt_allele_counts.py scripts/test_compare_trgt_locus_dosage.py \
  scripts/test_panel_c_repeat_alleles.py
```

## TRGT haplotype phasing

Prototype for the repeat-dosage re-base, not part of Figure 2C. The locus-level
`ru` arm infers a haplotype's repeat length by summing the integrated callset's
records at a catalog locus, which double-counts an event described once as a
GLnexus indel and once as a merged SV record: 99.84% exact where the repeat does
not change length, 8.10% exact where it expands by 50 bp or more. TRGT measures
the repeat directly. The obstacle is that the production TRGT VCFs are raw TRGT
output — `/`-separated genotypes, no `PS` — so no allele can be attributed to
the haplotype carrying a FLARE ancestry call. HiPhase was run on the small and
SV BCFs only, never on the TRGT VCF, so there is no phased sibling.

`wdl/TrgtHaplotypePhase.wdl` resolves it at read level for one sample over one
window: `bcftools` pulls the sample out of the ligated phased callset (this
defines the haplotype frame — slot 0 is the haplotype `AN1` describes),
`samtools` reads the window out of the haplotagged BAM and strips HiPhase's
`HP`/`PS` (they are in HiPhase's block frame, not SHAPEIT4's), `whatshap
haplotag` re-tags against the phased callset, `trgt` re-genotypes from the
tagged reads, and `phase_trgt_by_haplotag.py` joins read → haplotype against
read → allele. A single HiFi molecule carries both the repeat allele and the
flanking heterozygous SNVs that fix its haplotype, so the join is a measurement.

| Step | Workflow | Script | Output |
|---|---|---|---|
| Haplotype-resolved TRGT | `wdl/TrgtHaplotypePhase.wdl` | `phase_trgt_by_haplotag.py` | `<sample>.hap_alleles.tsv.gz` (locus × haplotype × length), `<sample>.trgt.phased.vcf.gz`, `<sample>.phase_summary.json` |
| Verification | `notebooks/terra/repeat_01_trgt_haplotype_phase.ipynb` | same | tag census, resolution tally, orientation agreement by length stratum, exposure |

Image: `repeat_loci/docker/Dockerfile` (whatshap + TRGT + samtools/bcftools),
built with `repeat_loci/build_docker.sh` (Cloud Build by default, `--local` for
a laptop build) and pushed as `aou-trgt-phase:0.1.0`. Nothing from the repo is
baked in — the WDL takes `phase_trgt_by_haplotag.py` as a File input from the
workspace bucket, as the other repeat_loci workflows do — so the build context
is just `repeat_loci/`. Example inputs:
`configs/trgt_haplotype_phase.inputs.json.example`.

`TrgtHaplotypePhase.docker` is pinned by **digest**, not tag, because a
re-pushed tag would silently change what ran:

    aou-trgt-phase@sha256:4a6d36d095aa7021431ce65bdea84780be5de58e23683869c01d2b58f890ac88

built 2026-10-20 from Cloud Build `c32d8ec2-81ed-4609-9ec5-29ec8a389137`, with
`trgt` 1.5.1 (the design needs ≥ 0.7.0 for `HP` during genotyping and ≥ 0.5.0
for read-to-allele assignments in the spanning BAM), `whatshap` 2.3 and
samtools/bcftools 1.21. After any rebuild, re-pin with the digest
`build_docker.sh` prints on exit.

### The `trgt_phase` table

The workflow needs four per-sample paths, and they live in two places: the
haplotagged BAM and its index in the `aou2_v1_phased_bams` Terra table, the TRGT
VCF in `tractor_mix/resources/trgt_table.txt`. Neither is a superset of the
other, so `trgt_phase` is a purpose-built table holding one row per sample that
has **both**:

    entity:trgt_phase_id  sample  bam  bam_index  trgt_vcf  trgt_vcf_index

Build it with:

    python3 scripts/make_trgt_phase_table.py \
      --bams-table aou2_v1_phased_bams.tsv \
      --trgt-table tractor_mix/resources/trgt_table.txt \
      --out trgt_phase.tsv

Nothing else the workflow takes belongs in the table. The phased per-chromosome
callsets are multi-sample files, identical for every row; the reference,
catalog, thresholds and image are run-level. A column that is the same in 12,000
rows is a column that can disagree with itself, so those stay in the config:

| Workflow input | Source |
|---|---|
| `sample_id` | `${this.sample}` |
| `bam` / `bam_index` | `${this.bam}` / `${this.bam_index}` |
| `production_trgt_vcf` / `_index` | `${this.trgt_vcf}` / `${this.trgt_vcf_index}` |
| `phased_vcf` / `_index`, `region`, everything else | literal in the config |

`phased_vcf` is per-chromosome and `region` must lie on that chromosome. Nothing
checks this: a chr20 window against a chr21 phased VCF fails at the first
`bcftools view`, so change the two together.

#### Two repeat catalogs, and they are not interchangeable

This trips you once. There are two files describing TR_Explorer 1.01 loci:

| File | Shape | Used by |
|---|---|---|
| TRGT repeat definitions | 4 fields: `chrom start end ID=<trid>;MOTIFS=<m>;STRUC=<s>` | `trgt --repeats`, i.e. `TrgtHaplotypePhase.trgt_repeats_bed` |
| `refs/repeat_loci/trexplorer_v1.0.1.catalog.bed.gz` | 5 fields: `chrom start end TRID motifs` | `aggregate_repeat_loci.py`, `trgt_allele_counts.py`, `trgt_plvi.py`, i.e. every other workflow here |

The 5-field file is **derived from** a TRGT VCF by `trgt_plvi.py catalog`, and in
deriving it `STRUC` was dropped. `STRUC` cannot be reconstructed from a motif
list at a compound locus — chr20:10,013,431 in the prototype window carries two
TRIDs and two motifs — so the derived catalog cannot be converted back. Feeding
it to `trgt --repeats` produces one `Locus Processing: Error` per locus
("Expected 4 fields ... found 5") and no genotypes.

`TrgtHaplotypePhase` therefore takes the TRGT-format definitions, and refuses a
file whose first row is not 4 fields containing `ID=`, naming the derived
catalog as the likely confusion.

Rebuild the definitions from any production TRGT **output** VCF — TRGT echoes
the definition it was given into INFO, so `TRID`, `MOTIFS` and `STRUC` are all
there (`STRUC=<TR22>` for a simple locus, `STRUC=<VC2>` for a compound one):

    CATALOG=gs://.../refs/repeat_loci/trexplorer_v1.0.1.catalog.bed.gz
    N=$(gsutil cat "$CATALOG" | gunzip | wc -l)

    python3 scripts/trgt_defs_from_vcf.py \
      --vcf gs://.../TRGT/TR_Explorer_1.01/1000234_trgt.TR_Explorer_1.01.vcf.gz \
      --expect-loci "$N" \
      --out trexplorer_v1.0.1.trgt_definitions.bed.gz

    gsutil cp trexplorer_v1.0.1.trgt_definitions.bed.gz \
      gs://.../refs/repeat_loci/

This is by construction the definition set those calls were genotyped against,
which is what keeps a re-genotyped call locus-for-locus comparable with the
delivered one — the notebook's re-genotyping check depends on it.

`--expect-loci` is the guard that matters: **one sample's VCF holds only the
loci TRGT emitted for that sample**, so a short file is possible. If the count
comes up short, either take the union over several samples or find the original
TRExplorer release BED; do not quietly genotype against a subset.

Two things the generator reports rather than hides:

- **The three populations.** Rows are written only for the intersection; samples
  with a BAM but no TRGT call, and TRGT calls with no usable BAM, are counted
  separately. Neither is an error — a sample with no TRGT call cannot be phased,
  and a TRGT call with no reads has nothing to phase it from — but both are a
  number you want up front rather than as a failed shard. Samples outside the
  intersection are absent from the table by construction and so cannot be
  submitted by accident.
- **Duplicate TRGT calls.** `trgt_table.txt` has **12,421 rows but 12,393
  distinct person ids**: 28 samples were genotyped under two different Terra
  submissions and both copies are listed (`trgt_ids.txt` and `trgt_files.txt`
  carry the same duplication; the rows split 10,023 / 1,234 / 762 / 402 across
  `terra-d45cc5ca`, `terra-f6a059a8`, `terra-e7fd2863`, `terra-f6671367`). One
  column holds one URI, so `--on-duplicate` makes the choice explicit: the
  default refuses and names them, `first` or `last` resolves. A duplicate
  outside the intersection never blocks the build.

Count TRGT samples as 12,393, not 12,421. Coverage still exceeds the 12,261
releasable participants, but by 132 rather than 160.

**Scaling past one window.** The workflow is single-sample, single-window by
construction, and `region` must carry explicit coordinates — a bare contig like
`chr20` leaves the task's `START` parsing with the contig name. Running a cohort
genome-wide therefore means either one submission per chromosome with explicit
coordinates, or a WDL change to scatter over a chromosome list and take
`Array[File] phased_vcfs`. Either way it is a new method snapshot, so decide it
before the next upload rather than after.

**Submit the first run with `allele_tag` empty.** The spanning BAM's
read-to-allele tag is TRGT-version-dependent, so it is never defaulted: an
empty `allele_tag` runs the join in `--census-only` mode, which writes the SAM
tag census and stops. Read the tag off `phase_summary`, then re-submit with it
set. A guessed tag would produce an empty join rather than an error.

Orientation only changes a per-ancestry dosage where a locus is
length-heterozygous **and** the sample's two haplotypes carry different local
ancestry — a same-ancestry swap leaves `DS{k}` unchanged, because it sums over
one bin. Section 5 of the notebook measures that product. It sizes the
genome-wide job, and a small value argues for a pooled length test plus a
per-ancestry carrier test instead of a per-ancestry length test.

Nothing here touches `FelixGenome.wdl` or any cached result.
