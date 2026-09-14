# Haplotig locus extraction (`bam_to_contig`)

Per-sample Terra workflow: map a GRCh38 `BED` through haplotig-vs-reference
BAMs and `faidx` those intervals from the haplotype FASTAs. One workflow row
= one sample (hap1 + hap2 concatenated).

Julie’s four cohort jobs (Jiadong UTR/CDS STRs, CYP2D6–7, CEL–CELP) all use
this WDL against `sample-hifi-hg38-all-cohorts` (~12,425 rows). Assemblies
and BAMs are GRCh38 (`AlignDiploidAssemblyToRef`), so the region BEDs are
GRCh38 — not CHM13.

Upstream scripts: [EichlerLab/AoU_WDL `bam_to_contig`](https://github.com/EichlerLab/AoU_WDL/tree/main/bam_to_contig)
([`scripts/bam_to_contig/SOURCE.md`](../scripts/bam_to_contig/SOURCE.md)).
The WDL here is Julie’s **per-sample** Terra adaptation.

Laptop `gsutil` cannot reach AoU Research Program buckets. Do every staging
and submit step on a **Terra Workbench notebook VM** (pmi-ops login).

## What you need on GitHub first

These files are not on `main` until you commit and push. `00_sync_repo`
clones GitHub, so **push before you sync**.

## 1. Sync the repo onto Terra

In the workspace, run [`notebooks/terra/00_sync_repo.ipynb`](../notebooks/terra/00_sync_repo.ipynb)
(or `python3 scripts/terra_sync_repo.py --ref main` on the VM). That will:

- rsync helper scripts → `$WORKSPACE_BUCKET/scripts/bam_to_contig/`
- copy the WDL → `$WORKSPACE_BUCKET/wdl/bam_to_contig/wdl/BamToContig.wdl`
- rsync region BEDs → `$WORKSPACE_BUCKET/bam_to_contig/regions/`

Check the **WDL manual import checklist** at the end of the notebook. If
`BamToContig.wdl` is listed (new/changed), Terra → **Workflows** → import /
replace using that GCS URI. AoU Methods-repo Create is usually 403.

## 2. Import the sample table (once)

Data table TSV: [`configs/sample-hifi-hg38-all-cohorts.tsv`](configs/sample-hifi-hg38-all-cohorts.tsv)
(12,425 entities: `{research_id}_{cohort}`).

In Terra: **Data → Import Data → Data Table** → that TSV. Skip if
`sample-hifi-hg38-all-cohorts` already exists with `asm_bam_h1` / `asm_h1`
columns.

Do **not** add this table to the default `00_sync_repo` upsert list — it is
12k rows.

## 3. Smoke test one sample

Before 12k rows: Workflows → `bam_to_contig` (or whatever you named the
import) → **Run**.

- Root entity type: `sample-hifi-hg38-all-cohorts`
- Choose **1 row**
- Inputs: [`configs/bam_to_contig.inputs.json.example`](configs/bam_to_contig.inputs.json.example)
  with `gs://WORKSPACE_BUCKET` replaced by the real `WORKSPACE_BUCKET`
- `flank_bp` / `regions_bed`: use Job 3 (`cyp2d6_7`, flank 30000) — one
  interval, cheapest check

Confirm `final_out` (`all_contigs_concat.fa`) has `>sample_hap1:` /
`>sample_hap2:` records. Empty FASTA usually means no haplotig fully covered
the locus (the Python filter requires a single alignment spanning the BED
interval before flank).

## 4. Launch Julie’s four jobs

Same WDL, same sample table, **one submission per job**. Only `regions_bed`
and `flank_bp` change. Select **all rows**.

| Job | `flank_bp` | `regions_bed` (under `$WORKSPACE_BUCKET/bam_to_contig/regions/`) | Julie’s cost note |
|---|---:|---|---|
| 1. Jiadong triplet UTR | 150 | `jiadong_variable_strs_utr.bed` (2559 intervals) | ~$400 |
| 2. Jiadong triplet CDS | 150 | `jiadong_variable_strs_cds.bed` (2309 intervals) | same as Job 1 |
| 3. CYP2D6–CYP2D7 | 30000 | `CYP2D6_through_7_GRCh38.bed` (`chr22:42126499-42144483`) | ~$120 |
| 4. CEL–CELP | 30000 | `CEL_through_CELP_GRCh38.bed` (`chr9:133061981-133087091`) | ~$120 |

`flank_bp` is the Python `locus_buffer` (fetch + extract padding). The shell
wrapper’s `bedtools merge` does not add a second pad.

Shared File inputs (constant across jobs):

```
bam_to_contig_bash_script = ${WORKSPACE_BUCKET}/scripts/bam_to_contig/bam_to_contig.sh
make_contig_bed_script    = ${WORKSPACE_BUCKET}/scripts/bam_to_contig/make_contig_bed.py
concat_contigs_script     = ${WORKSPACE_BUCKET}/scripts/bam_to_contig/concat_contigs.py
run_faidx_script          = ${WORKSPACE_BUCKET}/scripts/bam_to_contig/run_samtools_faidx_for_sample.sh
```

Data-table inputs:

| WDL input | Table column |
|---|---|
| `sample_name` | `this.sample-hifi-hg38-all-cohorts_id` |
| `bam_hap1` / `bai_hap1` | `this.asm_bam_h1` / `this.asm_bai_h1` |
| `bam_hap2` / `bai_hap2` | `this.asm_bam_h2` / `this.asm_bai_h2` |
| `haplotig_fasta_hap1` / `_hap2` | `this.asm_h1` / `this.asm_h2` |

Image: `us.gcr.io/broad-dsp-lrma/lr-talon:5.0` (pysam, bedtools, samtools).
Do not swap it without those tools.

## 5. Concatenate and send Julie the FASTAs

Each row writes one `all_contigs_concat.fa` (hap1 + hap2). After a
submission finishes, on the Terra VM:

```bash
gsutil cat \
  "$WORKSPACE_BUCKET/submissions/<submission_id>/bam_to_contig/*/call-ConcatContigs/**/all_contigs_concat.fa" \
  > jiadong_utr.all_contigs_concat.fa
```

Repeat per submission. Copy the four FASTAs somewhere Julie can read (workspace
bucket folder is enough) and send her the URIs.

Failed rows will be missing from that glob; check the submission UI and
resubmit failures only.

## WDL

[`wdl/BamToContig.wdl`](wdl/BamToContig.wdl) — workflow name `bam_to_contig`.

| Step | What |
|---|---|
| `RunBamToContig` (scatter hap1/hap2) | `make_contig_bed.py` walks CIGAR over the BED ± flank; `bedtools merge -d 100`; emit contig `chr:start-end` |
| `RunFaidx` | `samtools faidx` those intervals from the haplotig FASTA (gunzip if needed; reverse-complement names keep `/rc`) |
| `ConcatContigs` | concatenate the two haplotype FASTAs; headers become `{sample}_{hap}:{seq}` |

Haplotypes with no single haplotig covering the full BED interval are
dropped by the Python filter (empty that haplotype’s FASTA).
