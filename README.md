# AoU-LR Phase 2 analyses

Reproducible notebooks, CLIs, and Cromwell/WDL workflows for the All of Us
long-read Phase 2 manuscript. Companion LaTeX:
[`kvg/aou-lr-phase-2-manuscript`](https://github.com/kvg/aou-lr-phase-2-manuscript).

Most analyses run on **Terra**. Workbench-only notebooks (Locityper,
ExpansionHunter) live under `notebooks/rw/`. Person-level data are All of Us
controlled-tier; this repository does not publish callsets or covariates.

## Layout

| Path | What |
|---|---|
| [`notebooks/terra/`](notebooks/terra/) | Terra analysis notebooks (primary entrypoints) |
| [`notebooks/rw/`](notebooks/rw/) | Verily Workbench notebooks (`wb workflow` CLI) |
| [`scripts/`](scripts/) | CLIs used by notebooks and WDLs |
| [`locityper/`](locityper/) | Locityper targeted genotyping (VWB; Isaac / EichlerLab stream WDL) |
| [`expansion_hunter/`](expansion_hunter/) | ExpansionHunter STR genotyping (VWB; minicram + bw2 fork) |
| [`tractor_mix/`](tractor_mix/) | Tractor-Mix / SAIGE association (WDL, Docker, configs) |
| [`felix/`](felix/) | FELIX LAI GWAS (WDL, Docker, configs) |
| [`sv_annotation/`](sv_annotation/) | SV site table, CADD-SV, discovery plots |
| [`snv_stats/`](snv_stats/) | DeepVariant + GLnexus SNV/indel `bcftools stats` WDL |
| [`methylation_stats/`](methylation_stats/) | pb-CpG-tools bedMethyl coverage / haplotype / concordance |
| [`flare/`](flare/) | Per-population FLARE (fix T≈120 from cohort-wide EM) |
| [`bam_to_contig/`](bam_to_contig/) | Haplotig sequence at GRCh38 loci (Julie / EichlerLab) |

Package-level runbooks (image tags, input JSON, calibration gates):

- [Tractor-Mix / SAIGE](tractor_mix/README.md)
- [FELIX](felix/README.md)
- [SV annotation](sv_annotation/README.md)
- [SNV/indel stats (bcftools WDL)](snv_stats/README.md)
- [Methylation maps (pb-CpG-tools)](methylation_stats/README.md)
- [FLARE by population](flare/README.md)
- [Propagate annotations](propagate_annotations/README.md)
- [Haplotig locus extraction](bam_to_contig/README.md)
- [Locityper (Verily Workbench)](locityper/README.md)
- [ExpansionHunter (Verily Workbench)](expansion_hunter/README.md)
- [Terra notebook bootstrap](notebooks/terra/README.md)
- [Verily Workbench notebooks](notebooks/rw/README.md)

## How to run a Terra notebook

1. Upload `notebooks/terra/` and `scripts/` to the workspace as siblings
   (Jupyter cwd may be the workspace root or the notebooks folder).
2. Stage CLIs to the bucket (do this whenever scripts change):

   ```bash
   gsutil -m rsync -r scripts/ "$WORKSPACE_BUCKET/scripts/"
   ```

   Association WDLs also need `./scripts/stage_tractor_scripts.sh` and/or
   `./scripts/stage_felix_scripts.sh` (see those package READMEs).
3. Open the notebook. The first code cell rsyncs `$WORKSPACE_BUCKET/scripts/`
   onto the VM **before** importing (Terra’s persistent `edit/scripts/` is
   often stale). Set `TERRA_SYNC_SCRIPTS=false` only if you are iterating on
   local VM copies.
4. Many Hail / bcftools notebooks are **dry-run by default**. Set
   `SNV_RUN_PIPELINE=true`, `PCA_RUN_PIPELINE=true`, or `METH_RUN_PIPELINE=true`
   (see each intro cell) before a real run.

Person-level covariates used by Table 1 and association inputs:

- Local checkout: `tractor_mix/covariates.source_rebuilt.csv.gz`
- Terra: `$WORKSPACE_BUCKET/covariates/covariates.source_rebuilt.csv.gz`
- Override: `$AOU_COVARIATES`

## Manuscript map

Numbers in the paper should be regenerated from these outputs, then copied
into the LaTeX. Do not edit table TSVs by hand.

| Manuscript item | Entrypoint | Output |
|---|---|---|
| Table 1 cohort (`tab:cohort`) | `tractor_04_table1_cohort_summary.ipynb` | `summaries/manuscript/table1_cohort_summary.{tsv,md}` |
| Table S1 pedigrees (`tab:pedigrees`) | same notebook | `summaries/manuscript/tableS_phase2_pedigrees.{tsv,md}` |
| Table 2 SNV/indel rows (`tab:callset`) | **`snv_stats/wdl/BcftoolsGlnexusStats.wdl`**, then `snv_00_merge_glnexus_stats.ipynb` | `snv_indel_site_counts.{tsv,md}`, `snv_indel_sample_qc.tsv` |
| Table 2 assembly rows | QUAST downloads (no notebook yet) | `scripts/download_phase1_hifiasm_quast.sh`, `scripts/download_phase2_quast.sh` |
| Table 2 SV rows + discovery plots | `sv_03_manuscript_stats.ipynb` after `AnnotateSvCallset.wdl` | `sv_annotation/outputs/` (or Terra `sv_outputs/`) |
| Genome-wide methylation maps | **`methylation_stats/wdl/PbCpgSampleStats.wdl`** on `aou2_v1_phased_bams`, then `meth_00_merge_pbcpg_stats.ipynb`. Concordance: `PbCpgChromSites.wdl` on a Primrose vs Jasmine subsample. | `methylation_manuscript_numbers.{tsv,json}` |
| Ancestry PCs (`lr_PC*`) | `tractor_05` → `tractor_06` → `tractor_07` | `$WORKSPACE_BUCKET/pca/<run_label>/` |
| Association λGC / QQ / Manhattan | Tractor-Mix / FELIX / SAIGE WDLs, then `tractor_09` | workflow `results_tsvs` + notebook figures |

Per-sample SNV/indel QC (`snv_indel_sample_qc.tsv`) is **all VCF samples**.
Join to covariates locally to restrict to Phase 2 PacBio discovery, mid-pass,
or high-pass.

## Recommended order

Pipelines are independent after covariates exist, except where noted.

### 1. Cohort and Table 1

| Step | Notebook / script |
|---|---|
| Rebuild covariates from CDR / technical sheets | `tractor_00_cov_rebuild_source.ipynb` |
| Covariate atlas (optional QA) | `tractor_03_cov_summarize.ipynb` |
| Table 1 + pedigree supplement | `tractor_04_table1_cohort_summary.ipynb` |

High-pass size is the proxy `HIGH_PASS_N = 1133` in `tractor_04` (all
PacBio+ONT dual-tech discovery samples, plus highest-coverage PacBio-only
fill). There is no mid/high-pass column in the covariates table.

### 2. Long-read PCA (association covariates)

Hail Cloud Environment. Checkpoints are resumable under `PCA_RUN_LABEL`.

| Step | Notebook |
|---|---|
| Global autosomal SNV PCA | `tractor_05_pca_deepvariant_long_read.ipynb` |
| Fill remaining ancestry labels | `tractor_06_fill_lr_ancestry.ipynb` |
| Within-population PCs | `tractor_07_pca_within_population.ipynb` |
| Merge PCs into covariates | `scripts/merge_lr_global_pcs_into_covariates.py`, `scripts/merge_lr_pop_pcs_into_covariates.py` |

`tractor_05` reads `GL_INTERVAL_set` `VCF` shards (DeepVariant + GLnexus).
Do **not** use FastFilter `output_vcf` unless you want that filtered callset.
The PCA MatrixTable is **not** the Table 2 SNV catalog (common SNVs only).

### 3. SNV / indel catalog (Table 2)

Prefer a **per-row WDL** on `GL_INTERVAL_set` (Terra logs notebooks out at 24 h):

| Step | Entrypoint |
|---|---|
| `bcftools stats` on one chrom | `snv_stats/wdl/BcftoolsGlnexusStats.wdl` (`chrom=this.GL_INTERVAL_set_id`, `vcf=this.VCF`) |
| Merge shard `stats` from the Terra table | `snv_00_merge_glnexus_stats.ipynb` (`SNV_RUN_PIPELINE=true`) |

Submit the WDL against nuclear rows of `GL_INTERVAL_set` (skip `chrM`). After
all shards succeed, the `stats` column is written back to the table. Merge
with the notebook, or:

```bash
python3 scripts/snv_bcftools_sample_qc.py --from-firecloud --pull-stats --out-dir summaries/manuscript
```

Details: [`snv_stats/README.md`](snv_stats/README.md).

Per-sample QC is **all VCF samples**. Join `snv_indel_sample_qc.tsv` to
covariates locally for mid/high-pass means.

### 4. Assemblies (Table 2 QUAST)

Requester-pays GCS. Set `GOOGLE_PROJECT` to a billing project.

```bash
export GOOGLE_PROJECT=...
./scripts/download_phase1_hifiasm_quast.sh
./scripts/download_phase2_quast.sh          # JOBS=8 by default
```

### 5. Structural variants

| Step | Notebook / WDL |
|---|---|
| CADD-SV bundle (≥60 GB disk, once) | `sv_00_prepare_caddsv_annotations.ipynb` |
| RepeatMasker / simpleRepeat / segdup BEDs | `sv_01_stage_repeat_tracks.ipynb` |
| Sample ancestry TSVs | `sv_02_stage_sample_ancestry.ipynb` |
| Annotate + count | `sv_annotation/wdl/AnnotateSvCallset.wdl` |
| Site-count table + discovery plots | `sv_03_manuscript_stats.ipynb` |

Details: [`sv_annotation/README.md`](sv_annotation/README.md).

### 6. Methylation maps

Per-sample pb-CpG-tools pileups live on `aou2_v1_phased_bams` (`combined_bed`,
`hap1_bed`, `hap2_bed`). Prefer the WDL (Terra logs notebooks out at 24 h):

| Step | Entrypoint |
|---|---|
| Coverage + haplotype counts, one sample | `methylation_stats/wdl/PbCpgSampleStats.wdl` |
| Merge sample stats | `meth_00_merge_pbcpg_stats.ipynb` (`METH_RUN_PIPELINE=true`) |
| Primrose vs Jasmine chr22 dumps (subsample) | `methylation_stats/wdl/PbCpgChromSites.wdl` |

Details: [`methylation_stats/README.md`](methylation_stats/README.md).

### 7. Association (Tractor-Mix, SAIGE, FELIX)

Shared cohort from `tractor_01_prepare_inputs.ipynb`, then WDLs. GRM QC:
`tractor_08_qc_grm.ipynb`. Genome-wide post-workflow views:
`tractor_09_genome_post_workflow.ipynb`.

FLARE ancestry-switch vs `GQ`/`DP`, and old vs new per-pop FLARE (switch rates / tract lengths):
`flare_01_switch_gq_dp.ipynb`. Re-infer LAI per covariates population (own T):
`flare/wdl/FlareByPopulation.wdl` — start with the 1 Mb chr22 smoke test (AFR+AMR, controls excluded); see [`flare/README.md`](flare/README.md).

Details: [`tractor_mix/README.md`](tractor_mix/README.md),
[`felix/README.md`](felix/README.md),
[`flare/README.md`](flare/README.md),
[`propagate_annotations/README.md`](propagate_annotations/README.md).

### 8. Haplotig locus extraction (Julie)

Per-sample haplotig FASTA at GRCh38 intervals (Jiadong STRs, CYP2D6–7,
CEL–CELP): [`bam_to_contig/wdl/BamToContig.wdl`](bam_to_contig/wdl/BamToContig.wdl)
on the `sample-hifi-hg38-all-cohorts` table. Smoke-test **one row**, then four
cohort submissions. Details: [`bam_to_contig/README.md`](bam_to_contig/README.md).

### 9. Locityper (Verily Workbench)

Isaac’s per-sample genotype WDL, adapted for VWB Cromwell. Nearline CRAMs
are subset once with `print_reads` (not `samtools view` per shard). The
Workbench Workflows GUI is not used. Open
[`notebooks/rw/locityper_00_prep_reference.ipynb`](notebooks/rw/locityper_00_prep_reference.ipynb)
then
[`notebooks/rw/locityper_01_run_stream.ipynb`](notebooks/rw/locityper_01_run_stream.ipynb)
in a VWB Jupyter app and submit with `wb workflow job run`. Details:
[`locityper/README.md`](locityper/README.md).

### 10. ExpansionHunter (Verily Workbench)

Per-sample STR genotyping with the bw2 ExpansionHunter fork. Nearline CRAMs
are subset once with `make_minicram_for_expansion_hunter` (catalog regions +
mates), then EH runs `--analysis-mode optimized-streaming` on the local
minicram. Reuse the assembly38 FASTA from the Locityper prep notebook. Open
[`notebooks/rw/expansion_hunter_01_run.ipynb`](notebooks/rw/expansion_hunter_01_run.ipynb)
in a VWB Jupyter app and submit with `wb workflow job run`. Details:
[`expansion_hunter/README.md`](expansion_hunter/README.md).

## Ad hoc notebooks (not manuscript tables)

`tractor_cov_explore.ipynb`, `scratch_lr_pcs_hgna.ipynb`,
`explore_lr_pop_mahalanobis_qc.ipynb` — exploratory. Do not treat their
outputs as frozen paper numbers.

## Production data tables

Default Terra workspace for joint-call shards and FLARE URIs:
`allofus-drc-wgs-LR-prodData` / `AoU_DRC_LongReads_PhaseTwo_Storage`.

| Table | Used for |
|---|---|
| `GL_INTERVAL_set` (`VCF`, `VCF_idx`, `stats`) | DeepVariant + GLnexus chrom shards (`snv_00`, `tractor_05`) |
| `aou2_v1_phased_bams` (`combined_bed`, `hap1_bed`, `hap2_bed`) | pb-CpG-tools methylation pileups (`meth_00`, `PbCpgSampleStats.wdl`) |
| `aou_lr_chrom` (`model_chr_anc_vcf`) | FLARE LAI VCFs (`scripts/resolve_flare_uris.py`) |
| `sample-hifi-hg38-all-cohorts` (`asm_bam_h*`, `asm_h*`) | Haplotig BAM + FASTA for `BamToContig.wdl` |

```bash
python scripts/resolve_gl_interval_manifest.py --from-firecloud
python scripts/resolve_flare_uris.py --from-firecloud --scan-chrom chr22 --grm-chroms chr1 chr22
```

## Updating this document

When you add a notebook or a paper table/figure, add a row here **and** in
[`notebooks/terra/README.md`](notebooks/terra/README.md). Keep package READMEs
for WDL/Docker details; keep this file as the reviewer-facing map of
entrypoints and order.
