# Repeat-locus dosage (Figure 2C)

Per-locus repeat alleles from the SNV/indel and SV callsets, without TRGT
calls. Every length-changing record inside a TRExplorer v1.0.1 locus
(± `flank` bp) contributes its signed length change to that haplotype's locus
dosage, so split or partially merged records at one repeat are summed into a
single allele. The same dosage is what a repeat-aware association test would
use later (`--out-vcf`, FORMAT `GT[:AN1:AN2]`, INFO `RU_DOSAGE`).

Driven by [`notebooks/terra/fig2_02_panel_c_repeat_loci.ipynb`](../notebooks/terra/fig2_02_panel_c_repeat_loci.ipynb).

> **Also the association input.** `aggregate_repeat_loci.py --out-vcf
> --with-ancestry` is what the FELIX repeat-dosage branch tests, so the loci in
> Figure 2C and the loci in the association scan are the same objects. The
> coding rules below (sign, period, ambiguity, size split) are therefore
> association conventions too — see
> [`felix/REPEAT_DOSAGE.md`](../felix/REPEAT_DOSAGE.md) §6.1 before changing any
> of them.

| Step | Workflow | Script | Output |
|---|---|---|---|
| 1. Double-count check | `wdl/TrgtLocusConcordance.wdl` | `aggregate_repeat_loci.py`, `compare_trgt_locus_dosage.py` | Integrated vs TRGT locus dosage on one chromosome for a participant subset |
| 2. Allele tables | `wdl/RepeatLocusDosage.wdl` (run once for Phase 1, once for Phase 2) | `aggregate_repeat_loci.py` | `<label>.alleles.tsv.gz` (locus × dosage × haplotype count) and `<label>.summary.json` |
| 3. Example loci | `wdl/TrgtPlvi.wdl` | `trgt_plvi.py` | PLVI per TRGT locus, ranked within motif-length groups |
| Figure inputs | notebook | `panel_c_repeat_alleles.py` | Rarefaction curves, example-locus alleles, summary |

No output carries sample ids.

## Rules

- **Sign.** Each ALT is signed by its own length: `len(ALT) - len(REF)` for
  sequence-resolved alleles, else `SVLEN` / `END` / `SVTYPE` for symbolic ones.
  INV and BND add 0.
- **Assignment.** A record goes to the catalog locus it overlaps most, within
  `flank` bp. Records with a change above `max_bp` are ignored.
- **Ambiguity.** A sample with two or more unphased heterozygous length
  changes at one locus is set missing (both haplotypes), because the changes
  cannot be assigned to haplotypes; `hap_unphased_ambiguous` counts them. One
  unphased het plus any homozygous records is resolved. All Figure 2C inputs
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
  scripts/test_compare_trgt_locus_dosage.py scripts/test_panel_c_repeat_alleles.py
```
