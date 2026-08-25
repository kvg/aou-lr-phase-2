# propagate-annotations

Union SNV/indel, SV, and FLARE VCFs into one **fully annotated per-chromosome
callset** via a **sorted streaming multi-way merge**, then (on Terra) **scatter
across genomic shards** and concat to BCF.

Memory stays O(samples × FORMAT fields × open streams), not O(sites).

Designed for FLARE flip / bad-genotype QC: sites keep `GQ`/`DP`/`AD` from the
matching callset (SNV/indel or SV), and overlapping SNV sites get
`AN1`/`AN2`/`ANP*` from FLARE. Genotypes are not propagated; field names are
not renamed.

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

## Terra / WDL

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
