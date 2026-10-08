# Repeat loci (Figure 2C)

Panel C counts TRGT allele lengths at TRExplorer v1.0.1 loci. `TrgtPlvi` ranks
every locus by PLVI; the three examples are chosen from that ranking once the
allele counts exist. `TrgtAlleleCounts` counts haplotypes at each distinct
length. Run it once on the Phase 1 TRGT list and once on the Phase 2 list,
with the same catalog. `panel_c_repeat_alleles.py` turns the two allele tables
and the PLVI table into the rarefaction curve and the example histograms.

`RepeatLocusDosage` sums integrated SNV/indel/SV records. That is not the
Panel C input. `TrgtLocusConcordance` on chr20 is already done.

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
