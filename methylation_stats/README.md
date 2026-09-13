# Genome-wide methylation maps (pb-CpG-tools)

Numbers and figures for the main-text **Genome-wide methylation maps**
section. Per-read 5mC (Primrose / Jasmine MM/ML tags) is already in the
haplotagged BAMs. This package summarizes the **pileups** from
`aligned_bam_to_cpg_scores` v3.0.0 (`--pileup-mode model --modsites-mode denovo`).

The resource is meant as (1) a haplotype-resolved methylation **panel of
normals** for disease studies and (2) a later **GWAS / meQTL** input. This
thread filled cohort QC and Fig. S3; it did not run meQTLs or a WGBS
comparison.

## What fills each claim

| Claim | Source |
|---|---|
| Median CpGs ≥5× / ≥10×, median mean coverage | `PbCpgSampleStats.wdl` → `merge` |
| All 12,261 discovery samples haplotype-resolved | hap1/hap2 site counts (`frac_hap_slots_ge5`) |
| Primrose vs Jasmine concordance | `PbCpgChromSites.wdl` → `concordance` (site means) |
| Fig. S3 | `plot-s3`, restyled in the manuscript repo |

Discovery-only merge (n = 12,261 PacBio): median 28,871,895 CpGs ≥5×
(IQR 28,548,008–29,049,209), 25,640,013 ≥10×, mean CpG coverage 16.9×.
Haplotype-slot occupancy ≥5×: all samples, median 75.2%.

Caller mix (`pb_meth_caller`): Jasmine 7,124, Primrose 5,049 (1.4.0 3,832 +
1.2.0 1,217), mixed 88. Mixed rows used different movies, so there are **no
paired Primrose/Jasmine pileups on the same reads**. Concordance is
**between-group site means** (chr22, coverage ≥10×, |Δ| < 10 percentage
points). Locked numbers: 25 + 25 genomes, 613,729 sites, 93.3%, Pearson
r = 0.986.

## Terra table

Workspace: `allofus-drc-wgs-LR-prodData` /
`AoU_DRC_LongReads_PhaseTwo_Storage`. Entity: `aou2_v1_phased_bams`.

Pileup columns: `combined_bed`, `hap1_bed`, `hap2_bed`. Stats WDL writes
`stats_tsv` / `stats_json`. Chrom-site dumps live in `sites` (full table
already has ~12,261 URIs from `PbCpgChromSites`; do not resubmit that WDL
on every row just to plot).

On Terra, `edit/scripts/` is often stale. Stage from git, then let
`meth_00_merge_pbcpg_stats.ipynb` rsync `$WORKSPACE_BUCKET/scripts/`
before import:

```bash
gsutil -m rsync -r scripts/ "$WORKSPACE_BUCKET/scripts/"
```

Flags: `METH_RUN_PIPELINE=true` (pull stats + merge),
`METH_RUN_CONCORDANCE=true` (25+25 `sites` + concordance),
`METH_RUN_PLOT=true` (hexbins). Firecloud table fetch does not need the
pipeline flag.

## WDL

[`wdl/PbCpgSampleStats.wdl`](wdl/PbCpgSampleStats.wdl) — one sample:

| Input | Data table |
|---|---|
| `sample_id` | `this.aou2_v1_phased_bams_id` |
| `combined_bed` | `this.combined_bed` |
| `hap1_bed` | `this.hap1_bed` |
| `hap2_bed` | `this.hap2_bed` |
| `pbcpg_stats_py` | `$WORKSPACE_BUCKET/scripts/pbcpg_stats.py` |

Image: `aou-sv-annotation:0.1.6`. Example:
[`configs/pbcpg_sample_stats.inputs.json.example`](configs/pbcpg_sample_stats.inputs.json.example).

```bash
python3 scripts/pbcpg_stats.py inventory --from-firecloud \
  --out summaries/manuscript/pbcpg_table_inventory.json

python3 scripts/pbcpg_stats.py pull-stats --from-firecloud \
  --stats-dir methylation_stats/outputs/shards --jobs 8

python3 scripts/pbcpg_stats.py merge \
  --stats-dir methylation_stats/outputs/shards \
  --out-dir summaries/manuscript \
  --covariates tractor_mix/covariates.source_rebuilt.csv.gz \
  --discovery-only
```

Writes `methylation_manuscript_numbers.{tsv,json}` and
`pbcpg_sample_stats.tsv` (includes `coverage`, `pb_meth_caller`, `has_ONT`
when covariates join).

## Concordance

`pull-dumps` subsamples 25 Primrose + 25 Jasmine `sites` URIs from
`aou2_v1_phased_bams` (highest `combined_mean_cov` when stats exist; mixed
skipped). A separate `pbcpg_concordance` table is optional.

```bash
python3 scripts/pbcpg_stats.py concordance \
  --dumps-dir methylation_stats/outputs/chr22_dumps \
  --labels summaries/manuscript/pbcpg_concordance_labels.tsv \
  --out summaries/manuscript/methylation_concordance.json
```

Also writes `methylation_concordance.pairs.tsv.gz` (`mean_primrose`,
`mean_jasmine`) for plotting. Scores in dumps are 0–100.

## Fig. S3 (current)

Two real panels in
`aou-lr-phase-2-manuscript/figures/figS3_methylation.{png,pdf}`:

- **A** hexbin of Primrose vs Jasmine **site means** (not paired samples)
- **B** haplotype-slot fraction ≥5× vs HiFi `coverage` (n = 12,261)

Renderer: `scripts/mockup_figures/figS3_methylation.py` (Science style).
Inputs (gitignored) live at
`aou-lr-phase-2-manuscript/data/figS3/`:

- `methylation_concordance.pairs.tsv.gz`
- `pbcpg_sample_stats.tsv`
- `methylation_concordance.json`

```bash
# Terra / analysis repo, default matplotlib:
python3 scripts/pbcpg_stats.py plot-s3 \
  --pairs summaries/manuscript/methylation_concordance.pairs.tsv.gz \
  --stats summaries/manuscript/pbcpg_sample_stats.tsv \
  --out-prefix summaries/manuscript/figS3_methylation

# Manuscript restyle:
.venv/bin/python scripts/render_mockup_figures.py figS3_methylation
```

Panel B vertical lines are the **Table 1 means** (16.3× / 33.2×) of the
reconstructed mid-pass / high-pass groups, not a coverage cutoff.

## Mid-pass vs high-pass

There is **no** mid/high-pass column in covariates. Table 1
(`notebooks/terra/tractor_04_table1_cohort_summary.ipynb`) rebuilds it:

1. Universe: `final_releasable_v9` and `technology == PacBio` (n = 12,261).
2. High-pass (n = 1,133): all `has_ONT` dual-tech (n = 876), then the
   257 PacBio-only samples with highest `coverage`.
3. Mid-pass (n = 11,128): the rest (no ONT).

`GC` (BI/HA ~15×, BCM/UW ~30×) is design intent, not the classifier.

## Later (not done)

Resource intent: panel of normals + GWAS/meQTL. Do not expand S3 into a
QTL figure.

| Panel | Why | How |
|---|---|---|
| Imprint hap1/hap2 tracks (IGF2/H19 or SNRPN) | Shows the resource is haplotype-true | One high-pass sample; `extract-chrom` on hap1 and hap2 beds. IGF2/H19 hg38 window: `chr11 --start 1990000 --end 2160000` |
| `combined_mean_score` by `pb_meth_caller` (and `platform`) | Callers still need to be a GWAS/meQTL covariate; S3A is site-mean concordance only | Box/violin from `pbcpg_sample_stats.tsv` |

```bash
python3 scripts/pbcpg_stats.py extract-chrom \
  --bed <hap1.bed.gz> --chrom chr11 --start 1990000 --end 2160000 \
  --min-cov 5 --out igf2.hap1.tsv.gz
```

Repeat for hap2. Pick a high-pass sample from the Table 1 rule above
(ONT overlap, or top-coverage PacBio-only).

## Do not (from this thread)

- Treat S3A as overlapping-sample / paired-read concordance.
- Draw fake short-read bisulfite bars (no AoU WGBS). Sequence-class
  ascertainment needs RepeatMasker/segdup intersects if revived.
- Run `PbCpgChromSites` on all 12k again; `sites` is already on the table.
- Put meQTL / SV–methylation colocalization here (Fig. S5 / QTL section).
- Commit `data/figS3/` or `pbcpg_sample_stats.tsv` (research IDs).

## Tests

```bash
python3 -m pytest methylation_stats/tests -q
```
