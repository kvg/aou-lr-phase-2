# propagate-annotations

Three Terra entry points:

1. **`AnnotateFlareGqDp.wdl`** — FLARE backbone, region-indexed smoke test.
   Copies `FORMAT/GQ` and `FORMAT/DP` from the DeepVariant/GLnexus VCF onto
   matching FLARE sites. Keeps FLARE `GT`/`AN1`/`AN2`. This is the path for
   ancestry-switch vs genotype-quality QC on a 1 Mb interval.
2. **`PropagateFlareAncestry.wdl`** — **target** backbone (phased SNV/indel,
   SV, or joint). Copies covering FLARE `AN1`/`AN2` onto every target site,
   including interstitial sites that FLARE dropped. Keeps target `GT`.
   This is the path toward a `GT:AN1:AN2` VCF for FELIXla / extract-tracts.
3. **`PropagateAnnotations.wdl`** — full-chromosome **union** of SNV/indel +
   SV + FLARE into one annotated BCF (no `GT`; sparse FORMAT). Use this when
   you want every site from all three callsets, not a FLARE-only VCF.

Memory stays O(samples × FORMAT fields × open streams), not O(sites).

The union merger is designed for FLARE flip / bad-genotype QC: sites keep
`GQ`/`DP`/`AD` from the matching callset (SNV/indel or SV), and overlapping
SNV sites get `AN1`/`AN2`/`ANP*` from FLARE. The union path does not
propagate genotypes; field names are not renamed.

## Performance (v0.3)

| Optimization | Effect |
|--------------|--------|
| **Sparse FORMAT** | Each site emits only tags from matched sources (no `SV_*` dots on SNV sites) |
| **BCF output** | Streams VCF text into `bcftools view -Ob` (much smaller than `.vcf.gz`) |
| **Rayon `--threads`** | Parallel sample-column assembly within each site |
| **WDL region shards** | `n_shards` (default 8) scatter by **equal variant count** (SNV density; FLARE/SV fallback), then concat |
| **Uncompressed shard inputs** | `bcftools view -Ov` into the Rust parser (avoids gunzip) |
| **Progress every 100k sites** | sites/s + elapsed minutes on stderr |

**WDL prep:** whole-genome SV BCFs are subset with `bcftools view -r/-t
{chromosome}`; per-chrom SNV/indel and FLARE are copied and indexed only (no
full-file `bcftools sort` — that was multi-hour on wide callsets).

## Inputs (at least one required)

| Flag | Description |
|------|-------------|
| `--snv-indel-vcf` | Joint SNV/indel callset (`.vcf` / `.vcf.gz`) |
| `--sv-vcf` | SV partition (convert BCF → VCF first, or let WDL do it) |
| `--flare-vcf` | FLARE `*.anc.vcf.gz` |
| `--output` | `.bcf` (preferred), `.vcf`, or `.vcf.gz` |
| `--threads N` | Rayon workers (0 = all CPUs) |
| `--dense-format` | Legacy: emit every FORMAT tag with `.` fillers |
| `--region CHR:START-END` | Optional filter (WDL shards pass this) |

**Output sites** = sorted union of all variant keys from the provided inputs.
Samples = union of sample columns (SNV order first, then extras from SV/FLARE).

## Algorithm

Each input is a one-record stream. Take the minimum `(CHROM, POS, REF, ALT)`,
merge FORMAT/INFO from matching streams, write one row, advance. Unsorted
inputs fail with a clear error.

## Output FORMAT fields

| Source | Fields (original names) |
|--------|-------------------------|
| SNV/indel | `GQ`, `DP`, `AD` |
| SV | `GQ`, `DP`, `AD` |
| FLARE | `AN1`, `AN2`, `ANP1`, `ANP2` (if present) |

No `GT`. Sparse mode emits only tags from sources that matched the site.
Selected INFO tags are copied with original names (`AF`, `SVTYPE`, …).

INFO `SRC` lists which sources matched each row.

## Build & test

```bash
cd propagate_annotations
cargo test
cargo build --release
```

## Docker

```bash
./build_docker.sh              # Cloud Build → Artifact Registry
./build_docker.sh --local      # local build
```

Image: `us-central1-docker.pkg.dev/broad-dsp-lrma/aou-lr/aou-propagate-annotations:0.3.1`

## Sanity-check stats

Every run writes `{prefix}.propagate_stats.json` and `.propagate_stats.tsv`
(WDL merges shard stats in `ConcatShards`):

| Metric | Meaning |
|--------|---------|
| `union_sites` / `output_sites` | Must be equal |
| `{source}_source_sites` | Sites in that input VCF |
| `{source}_union_sites_with_source` | Union rows annotated by that source |
| `snv_flare_overlap_rate` | Fraction of SNV/indel sites that also have FLARE |

**Fail-fast:** WDL default `min_match_rate_flare = 0.95` (checked after concat).

## Terra / WDL — region GQ/DP/RNC onto FLARE (`AnnotateFlareGqDp`)

Workflow: [`wdl/AnnotateFlareGqDp.wdl`](wdl/AnnotateFlareGqDp.wdl)

Stage the annotator, then submit on a 1 Mb (or full-contig) interval. Both
VCFs must be tabix/csi-indexed; the task uses `bcftools view -r` then a
linear streaming merge (`scripts/annotate_flare_gq_dp.py`).

On Terra (Pipelines API) the VCF and index inputs are
`localization_optional`
([Cromwell docs](https://cromwell.readthedocs.io/en/latest/optimizations/FileLocalization/)):
Cromwell leaves them as `gs://` URIs and htslib range-reads the region.
Set `gcs_project` to the workspace Google project so requester-pays AoU
buckets work (`GCS_REQUESTER_PAYS_PROJECT`). The same command accepts local
paths if Cromwell does localize (miniwdl). Disk is
`ceil((size(flare)+size(dv)) × multiplier) + floor` (defaults `1.0` / `80`).
The annotator streams into `bcftools view -Oz`; DeepVariant is restricted to
FLARE positions so a GLnexus chrom shard is not copied in full.

```bash
gsutil -m rsync -r scripts/ "$WORKSPACE_BUCKET/scripts/"
```

| WDL input | Typical source |
|-----------|----------------|
| `flare_vcf` / `flare_vcf_index` | `aou_lr_chrom.model_chr_anc_vcf` (+ `.tbi`) |
| `deepvariant_vcf` / `deepvariant_vcf_index` | `GL_INTERVAL_set.VCF` / `VCF_idx` |
| `region` | `chr22:16000000-17000000` (or `chr22` for the whole contig) |
| `annotate_script` | `gs://$WORKSPACE_BUCKET/scripts/annotate_flare_gq_dp.py` |
| `gcs_project` | Terra workspace Google project (requester-pays) |
| `min_match_rate` | warn (default) if FLARE–DeepVariant site overlap is below this (default `0.95`); set `fail_on_low_match=true` to fail after writing outputs |

Outputs: `{prefix}.flare.gq_dp.vcf.gz` (+ `.tbi`) and match-rate stats.
Default `format_tags` is `GQ,DP,RNC`. GLnexus `RNC` is two characters per
sample (`..` = called; `I` = incomplete gVCF). That VCF is a valid input to
`scripts/flare_switch_qc.py` / `notebooks/terra/flare_01_switch_gq_dp.ipynb`.

Example: `configs/annotate_flare_gq_dp.inputs.json.example`

Local smoke (no Terra):

```bash
python3 scripts/annotate_flare_gq_dp.py \
  --flare flare.region.vcf.gz \
  --deepvariant dv.region.vcf.gz \
  --output annotated.vcf \
  --stats-json stats.json
python3 -m pytest scripts/test_annotate_flare_gq_dp.py -q
```

## Terra / WDL — FLARE ancestry onto interstitial target sites (`PropagateFlareAncestry`)

Workflow: [`wdl/PropagateFlareAncestry.wdl`](wdl/PropagateFlareAncestry.wdl)

FLARE only emits LAI markers (sites that overlap the reference panel and pass
MAF/MAC). FELIXla and `extract-tracts-flare` need `AN1`/`AN2` on the **target**
sites you want to test — SNVs/indels between markers, and SVs. This workflow
keeps the target VCF as the backbone (phased `GT` stays) and copies covering
FLARE ancestry in one streaming pass.

Interval rule matches FELIXla:

- first FLARE site on a contig covers `1..current_lai_pos`
- later FLARE sites cover `(previous_lai_pos, current_lai_pos]`
- the last FLARE state extends through the last target site

The WDL left-pads the FLARE extract to `CHR:1-END` so a 1 Mb target window
still sees the marker before `START`. Use `aou-flare-bcftools:1.24` (htslib
GCS retry). Default `min_sample_overlap` is `0` because per-population FLARE
drops EUR/EAS/SAS and reference controls.

```bash
gsutil -m rsync -r scripts/ "$WORKSPACE_BUCKET/scripts/"
```

| WDL input | Typical source |
|-----------|----------------|
| `flare_vcf` / `flare_vcf_index` | `{prefix}.anc.vcf.gz` from `FlareByPopulation` (or `model_chr_anc_vcf`) |
| `target_vcf` / `target_vcf_index` | phased chrom VCF (`gt` / joint SNV+SV) |
| `region` | `chr22:26897597-27897597` (or `chr22` for the whole contig) |
| `propagate_script` | `gs://$WORKSPACE_BUCKET/scripts/propagate_flare_ancestry.py` |
| `gcs_project` | Terra workspace Google project (requester-pays) |
| `min_coverage_rate` | warn (default) if target sites with ancestry are below this (default `0.95`); set `fail_on_low_coverage=true` to fail after writing outputs |

Outputs: `{prefix}.ancestry.vcf.gz` (+ `.tbi`) and coverage stats.
`INFO/FLARE_POS` is the covering marker. Samples absent from FLARE get `AN1=.,AN2=.`.

Example: `configs/propagate_flare_ancestry.inputs.json.example`

Local smoke (no Terra):

```bash
python3 scripts/propagate_flare_ancestry.py \
  --flare flare.region.vcf.gz \
  --target phased.region.vcf.gz \
  --output annotated.vcf \
  --stats-json stats.json
python3 -m pytest scripts/test_propagate_flare_ancestry.py -q
```

## Terra / WDL — full-chrom union (`PropagateAnnotations`)

Workflow: `wdl/PropagateAnnotations.wdl`

| WDL input | Typical column |
|-----------|----------------|
| `chromosome` | `aou_lr_chrom_id` (e.g. `chr22`) |
| `snv_indel_vcf` | `vcf_gz` |
| `sv_vcf` | integrated SV partition (whole-genome BCF) |
| `sv_vcf_index` | optional `.csi` (fast `-r`; without it WDL streams via `-t`) |
| `flare_vcf` | `model_chr_anc_vcf` |
| `sv_is_whole_genome` | `true` when `sv_vcf` spans all chromosomes |
| `n_shards` | equal-variant parallelism (default 8; balances acrocentric chroms) |
| `disk_gb_floor` / `disk_gb_multiplier` | PrepareChrom + Concat disk = `ceil(size(inputs)×mult)+floor` |
| `shard_disk_gb_floor` / `shard_disk_gb_multiplier` | PropagateShard disk from chrom VCF sizes (same formula) |
| `prefix` | e.g. `aou_lr_phase2_v1.chr22` |

Outputs: `{prefix}.annotated.bcf` (+ `.csi`) and merged stats.

Example: `configs/example.inputs.json`

## Notes

- Inputs must be sorted by `(CHROM, POS)` (standard VCF/TBI order). Several
  ALTs at one POS need not be REF/ALT-sorted — the merger sorts that locally.
- SV-only sites omit FLARE FORMAT tags (sparse mode).
- SNV-only sites omit ancestry tags unless FLARE also matches.
- Multi-allelic sites must use the same REF/ALT representation across inputs.
- Whole-genome SV BCFs: pass `sv_vcf_index` (`.csi`) for fast chrom extract;
  if omitted, PrepareChrom uses `bcftools view -t` (no index, streams the file).
