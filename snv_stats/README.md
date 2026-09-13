# DeepVariant + GLnexus SNV/indel stats

Preferred path for manuscript Table 2 small-variant rows and per-sample
Ti/Tv / het/hom: Terra data-table launch of `bcftools stats`, not Hail.

## WDL (recommended)

[`wdl/BcftoolsGlnexusStats.wdl`](wdl/BcftoolsGlnexusStats.wdl) processes **one
chromosome**. Run it on the `GL_INTERVAL_set` table so Terra submits one
workflow per row:

| Input | Data table |
|---|---|
| `chrom` | `this.GL_INTERVAL_set_id` |
| `vcf` | `this.VCF` |

Skip `chrM`. Keep chr1–22, X, Y.

Image: `aou-sv-annotation:0.1.6` (bcftools). Example root entity JSON:
[`configs/bcftools_glnexus.inputs.json.example`](configs/bcftools_glnexus.inputs.json.example).

Optional knobs: `cpu`, `memory_gb`, `preemptible` (set `0` if preemptibles
die mid-scan), `disk_gb_floor` / `disk_gb_multiplier` if a shard OOMs on disk.

Each row writes `stats` (`bcftools stats -s - -f PASS,.`). GLnexus usually
leaves FILTER as `.` (unfiltered), so `-f PASS` alone yields 0 records.

After the nuclear
chromosomes succeed, merge from the Terra table:

```bash
python3 scripts/snv_bcftools_sample_qc.py --from-firecloud --pull-stats --out-dir summaries/manuscript
```

Or open `notebooks/terra/snv_00_merge_glnexus_stats.ipynb` and set
`SNV_RUN_PIPELINE=true`. That produces `snv_indel_site_counts.{tsv,md}` and
`snv_indel_sample_qc.tsv` (all VCF samples; join to covariates later).

Catalog indel bins `<20 bp` and `<50 bp` come from the bcftools IDD section.
Per-sample PSC `n_indel` is not length-binned.
