# Per-population FLARE

The popout [`flare.wdl`](https://github.com/broadinstitute/popout/blob/main/workflows/flare/wdl/flare.wdl)
ran all 12 347 Phase 2 samples as one FLARE `gt` VCF. EM then estimated a
single generations-since-admixture **T ≈ 120** (started at `gen=10`) and
study-wide µ equal to cohort composition (EUR 41%, AFR 24%, …). That T
predicts ~0.8 Mb mean tracts; it is not African American / Latino admixture.

This workflow splits the same phased chrom VCF by covariates `population`,
drops HPRC/HGSVC3/GIAB **reference controls** (they overlap the LAI panel),
runs FLARE in each group with **`update-p=true`**, and merges `AN1`/`AN2`.

## WDL

[`wdl/FlareByPopulation.wdl`](wdl/FlareByPopulation.wdl)

| Step | What |
|---|---|
| `SplitSamples` | `bcftools query -l` + `scripts/flare_split_samples.py` |
| `CutRefRegion` | if `region` is set: one tabix slice of the LAI ref (shared by every FLARE shard). Retries on curl/HTTP flakes. |
| `SubsetOnePopulation` | **scatter** one shard per keep-list: `bcftools view -S` (and `-r` when `region` is set). Streams `gs://` (`localization_optional`). bcftools/htslib **1.24** reconnects mid-stream on curl 92 ([htslib#1987](https://github.com/samtools/htslib/pull/1987)). Shell + Cromwell `maxRetries` are a backstop. `subset_preemptible` default 5 covers GCE reclaiming a ~90 min VM. Seven pops ⇒ seven parallel GCS reads. |
| `flare_task` | FLARE 0.6 docker from popout; own EM / T / µ per group (same scatter) |
| `MergePopulationFlare` | `bcftools merge -m none` + concatenate `.global.anc.gz` + `models.tsv` |

FLARE forbids missing GT, so we subset samples rather than convert no-calls
to `./.`. Resource scaling follows popout (heap ≈ 2.2× gt GB), with a **small
bucket** (4 CPU / 16 GB) when the subset gt is &lt; 1 GB so 1 Mb tests do not
request a 48 GB VM. **`probs=true` raises floors** (≈72g `-Xmx` / 96 GB VM
even for &lt;1 GB subsets) because full-chrom ANP fwd/bwd OOMs the small
bucket — compressed size underestimates peak heap. Override with
`flare_memory_override` **and** `flare_xmx_gb_override` together if needed.
Magicwand/W&B is omitted.

Image for split/subset/merge: `aou-flare-bcftools:1.24` (bcftools/htslib 1.24).
Do **not** reuse `aou-sv-annotation:0.1.6` here — that image is bcftools 1.20
and has no GCS HTTP retry. FLARE image:
`us-docker.pkg.dev/broad-dsde-methods/popout/flare:latest`.

Build / push once (same Artifact Registry as the other AoU-LR images):

```bash
cd flare
./build_docker.sh                 # Cloud Build → aou-flare-bcftools:1.24
```

## 1 Mb smoke test (do this first)

Colleague chr1 was ~16 h because FLARE indexed a full chromosome of target
**and** reference. Pass `region` so the LAI ref is sliced once and each
population shard tabix-cuts the same interval from gt (`view -S -r`).

The example JSON uses the same interval as switch QC:

`chr22:26897597-27897597` (~1 Mb, ~5.7k FLARE sites in the original run)

AFR+AMR only, controls excluded. Expect minutes, not hours. Then read
`*.models.tsv` (`t_gen`). Clear `region` only after T looks sane.

## Inputs

Example: [`configs/chr22.inputs.json.example`](configs/chr22.inputs.json.example).

Stage scripts first (from a Terra notebook — see
[`notebooks/terra/00_sync_repo.ipynb`](../notebooks/terra/00_sync_repo.ipynb);
laptop `gsutil` cannot reach AoU Research Program buckets):

```bash
# on the Terra VM after 00_sync_repo, or:
python3 scripts/terra_sync_repo.py --ref main --upsert-tables flare_lai_exp
# equivalent bucket staging only:
gsutil -m rsync -r scripts/ "$WORKSPACE_BUCKET/scripts/"
```

| Input | Notes |
|---|---|
| `covariates` | `covariates.source_rebuilt.csv.gz` (`research_id`, `population`, `is_reference_control`) |
| `gt_vcf` / `gt_vcf_index` | Phased per-chrom target (same file the original FLARE used) |
| `ref_vcf` / `ref_vcf_index` | gnomAD LAI (or other) for **this** chrom; index required for `region` |
| `ref_panel` | `sample <tab> panel` |
| `map_file` | PLINK cM map; **chromosome column must match the VCF** (`chr22`, not `chrchr22`) |
| `split_script` / `summarize_script` / `flare_model_script` | `$WORKSPACE_BUCKET/scripts/flare_split_samples.py`, `flare_summarize_models.py`, and `flare_model.py` |
| `region` | `chr22:start-end` smoke test (empty = full contig, streamed) |
| `subset_preemptible` | Default **5**. Full-chrom subset is ~90 min on a preemptible VM. |
| `max_retries` | Default **5**. Backstop if htslib’s mid-stream retry still gives up. curl 92 is a task failure, not a preemption. |
| `include_pops` | Optional comma list, e.g. `AFR,AMR` |
| `drop_pops` | Optional, e.g. `OTH,MID` |
| `exclude_controls` | Default **true**. Drops `is_reference_control` and `HG*`/`NA*` IDs |
| `min_samples` | Default 50; smaller groups are skipped |
| `update_p` | Default **true** (popout left this at FLARE’s false) |
| `em` | **true** = estimate T/µ (`gen` is EM start only). **false** = pin T at `gen` / `gen_by_pop` (no `.model` required) |
| `gen` | Default T (or EM start). Literature pin: AFR ~8, AMR ~12 |
| `gen_by_pop` | Optional overrides, e.g. `AFR:8,AMR:12`. Blank = use `gen` for every shard |
| `gen_by_pop_allow_default` | Default **false**. When `gen_by_pop` is set, every shard pop must appear unless this is true |
| `props_by_pop` | Optional admixture proportions, e.g. `AFR:AFR=0.80,EUR=0.20` (case-insensitive ancestry labels) |
| `inherit_props` / `inherit_panel_weights` / `inherit_mu` | When a source `.model` is present, keep or reset each block independently (`mu` = miscopy matrix) |
| `template_model` | Required when resetting/overriding blocks without a matched `pop_models` entry |
| `allow_unrepresented_pops` | Default **true** (temporary). Warn — do not fail — when a shard pop has no same-named panel ancestry (MID today). Revisit with a MID panel / remap later. |
| `probs` | Default **false** (no `ANP`). Set **true** only when posterior dosages are needed |
| `pop_models` | Empty for explore / pin-T. Apply: basenames must encode pop as `*.<POP>.model` or Cromwell `*.<POP>.out.model` / `.in.model`. Rewritten under `in_model/` so outputs never collide |
| `seed` | Fixed `12345` ⇒ identical outputs show **determinism**, not independent corroboration |

`population` is uppercase AFR/AMR/EAS/EUR/SAS/OTH/MID (legacy / soft labels).
VCF sample IDs must match `research_id`. Split on hard `ancestry_pred`
(lowercase in covariates; uppercased to shard names; no AoU `oth` after Rule C).
Method config: `pop_column = ancestry_pred`. **MID policy (for now):** keep MID
shards with `allow_unrepresented_pops=true` (WDL default); panel lacks a MID
ancestry, so validation warns only. A dedicated MID panel or remap may replace
this later.
## Terra experiment table (10 Mb method grid)

Use a dedicated table so each row is one LAI recipe on the same window.
Do not overload `aou_lr_chrom` with this.

1. Import [`configs/lai_exp.tsv`](configs/lai_exp.tsv) (chr22 paths already filled
   from `aou_lr_chrom`). Creates table `flare_lai_exp`.
2. Set column types: `em` / `update_p` / `exclude_controls` / `inherit_*` /
   `probs` / `allow_unrepresented_pops` / `gen_by_pop_allow_default` → boolean;
   `gen` / `min_maf` → number; `min_mac` / `seed` → number;
   VCF/map/`ref_panel`/`template_model`/`exclude_regions`/`include_sites` → file;
   `pop_models` → array of file (`[]` on most rows).
3. Re-import [`wdl/FlareByPopulation.wdl`](wdl/FlareByPopulation.wdl).
4. Root entity: `flare_lai_exp`. Config:
   [`configs/lai_exp.inputs.json.example`](configs/lai_exp.inputs.json.example)
   (`this.region`, `this.em`, `this.gen`, `this.gen_by_pop`, `this.props_by_pop`,
   `this.inherit_*`, `this.ref_panel`, `this.pop_models`, …).
5. Launch selected rows with `pop_column=ancestry_pred` (no OTH shard). MID is
   allowed via `allow_unrepresented_pops=true`. Re-import the updated
   `lai_exp.tsv` so `gen_by_pop` / `props_by_pop` no longer list OTH.
6. Score with [`notebooks/terra/flare_02_lai_exp_compare.ipynb`](../notebooks/terra/flare_02_lai_exp_compare.ipynb)
   (Firecloud-fetches `flare_lai_exp`, switch-QCs finished rows).

Fixed window (legacy grid): `chr22:26897597-36897597` (~10 Mb). Chr20 rows use
full contig (`region` empty) with paths from `aou_lr_chrom`
(`aou_lr_phase2_v1.chr20` gt/ref/map). Context mask is staged at
`gs://fc-secure-8f7d6a20-04ce-40d7-8c88-aececeac3e09/refs/flare/flare_context_mask.rmsk_sr_sd.bed.gz`
(via `flare_03_stage_site_filters.ipynb`). HiFi↔srWGS concordance `include_sites`
ladder step is deferred (no `chr20_filter_concordance` row for now). Call-QC
`include_sites` (GQ/DP/RNC) is Part 7 — stage with `flare_build_call_qc_sites.py`
then fill `chr22_filter_call_qc` / `chr20_filter_call_qc`. FLARE2 rows need the
WDL two-pass mode (Part 6) before launch.

| Row id | Intent |
|---|---|
| `em_baseline_afr_amr` | Current per-pop EM |
| `pin_gen8_afr_amr` / `pin_gen12_afr_amr` | Single pinned T |
| `pin_gen_by_pop` | AFR:8 and AMR:12 together |
| `em_strict_mac` / `pin_gen8_strict_mac` | Stricter `min_maf`/`min_mac` |
| `pin_afr_only_gen8` / `pin_amr_only_gen12` | One population only |
| `em_all_pops` | Negative control |
| `pin_all_pops_gen_by_pop` | All pops; pin AFR:8 AMR:12 others:8 (flat µ) |
| `pin_all_pops_em_mu` | **Discard prior conclusions** (old `out_model` glob picked rewritten input). Re-run only after model I/O fix |
| `flare2_amr_nanc2` / `flare2_amr_nanc3` | **Blocked until FLARE2 WDL**; AMR clustering pilot |
| `flare2_mid_nanc2` | **Blocked until FLARE2 WDL**; MID poorly-matched panel pilot |
| `flare2_afr_nanc2` | **Blocked until FLARE2 WDL**; AFR control (expect little gain) |
| `chr20_flat_props_pin_t` | Full chr20; FLARE defaults + pinned T; `probs=true` |
| `chr20_em_props_defaults` | Full chr20; `props_by_pop` from EM; default weights/µ |
| `chr20_em_props_inherit_weights` | Full chr20; EM props + inherit panel weights |
| `chr20_em_props_inherit_weights_mu` | Full chr20; + inherit miscopy µ |
| `chr20_filter_*` | Site-restriction ladder (baseline → biallelic → context mask) |
| `chr22_filter_call_qc` / `chr20_filter_call_qc` | Ladder + GQ/DP/RNC `include_sites` (fill path after staging) |

`pop_models` is empty (`[]`) except when attaching prior EM models. One method
config covers every row; no per-row JSON override.

Prior `pin_all_pops_em_mu` switch-QC conclusions are invalid (defect: input
model globbed as `out_model`) — re-interpret only after the model I/O fix.
MID is allowed via `allow_unrepresented_pops=true`. HiFi–srWGS concordance → µ
measurement is deferred.

For a HQ-thinned `gt`, change that row’s `gt_vcf` / `gt_vcf_index` only.

## What “sane T” looks like

| Group | T that would look sane |
|---|---|
| AFR (African American) | ~6–15 |
| AMR (Latino) | ~8–20 |
| EUR / EAS / SAS if you add them later | small (a few generations) if they are mostly unadmixed |

A 1 Mb window makes T noisier than a chromosome, but 120 vs 12 is obvious.
If AFR T is still ~100 on the slice, splitting the cohort was not enough.

Do **not** reuse high-T EM `.model` files from noisy chr1 runs for genome-wide apply.

## Terra chrom table (train once, apply)

You do **not** need a dummy extra row. **chr22 is the training row.** Other autosomes are apply rows. Other rows must be **`em=false`**, not `em=true` — otherwise each chrom re-estimates \(T\).

Terra will not automatically point chr1 at chr22’s outputs. After the train submission succeeds, copy the `per_pop_model` URIs into columns (or workspace attributes) and reference those on every apply row.

Suggested `aou_lr_chrom` columns:

| Column | chr22 (train) | chr1, chr3, … (apply) |
|---|---|---|
| `lai_em` | true | false |
| `lai_model_AFR` | (empty) | gs://…/`*.AFR.model` from the chr22 job |
| `lai_model_AMR` | (empty) | gs://…/`*.AMR.model` from the chr22 job |
| `gt` / `ref` / `map` | this chrom | this chrom |

Workflow inputs:

- `em` = `this.lai_em`
- `pop_models` = `[this.lai_model_AFR, this.lai_model_AMR]` on apply rows; leave unset on chr22
- Same `include_pops=AFR,AMR`, `drop_pops=OTH,MID`, `seed`, `min_maf`, `min_mac`
- Pin `flare_cpu_override` to the same value on train and apply (FLARE results depend on `nthreads`)

If a model file is present for a population, the WDL forces `em=false` for that pop even if `lai_em` is true.

Apply example: [`configs/chrom.apply.inputs.json.example`](configs/chrom.apply.inputs.json.example).

## Outputs

| File | Contents |
|---|---|
| `{prefix}.anc.vcf.gz` | Merged samples, same FLARE FORMAT (`GT`, `AN1`, `AN2`) |
| `{prefix}.global.anc.gz` | Concatenated per-sample global ancestry |
| `{prefix}.models.tsv` | One row per population: `t_gen`, `prop_*` (admixture), `mu_*` (mean miscopy) |
| `excluded_controls` | HG/NA (and flagged) samples left out of the target |
| `per_pop_*` | Unmerged FLARE VCF / `.out.model` / `.in.model` / `.log` |

Samples dropped for missing labels, controls, or `min_samples` are **absent**
from the merged VCF. See `unmatched_vcf_samples` and `sample_manifest`.

## Evaluate old vs new

[`notebooks/terra/flare_01_switch_gq_dp.ipynb`](../notebooks/terra/flare_01_switch_gq_dp.ipynb)
re-scans the old annotated backbone and this merged anc VCF on **shared
samples** (so dropped controls / EUR-EAS do not confound the rate).

Set `NEW_FLARE_VCF` to `{prefix}.anc.vcf.gz` and `NEW_FLARE_MODELS` to
`{prefix}.models.tsv`. Optional `COMPARE_REGION` (the 1 Mb smoke-test window
is `chr22:26897597-27897597`). Full chr1: leave the region empty.

The fair 1 Mb metric is **switches per hap per Mb** (≈ T/100 at 1 cM/Mb).
Tract-length ECDFs wait for a chromosome-scale VCF.

Stage `scripts/flare_switch_qc.py` before re-running the notebook:

```bash
gsutil -m rsync -r scripts/ "$WORKSPACE_BUCKET/scripts/"
```

## Ancestry onto interstitial sites

FLARE only writes LAI markers. To put `AN1`/`AN2` on the rest of a phased
chrom VCF (SNVs between markers, SVs), use
[`propagate_annotations/wdl/PropagateFlareAncestry.wdl`](../propagate_annotations/wdl/PropagateFlareAncestry.wdl).
Same 1 Mb window works: `chr22:26897597-27897597`.

## Tests

```bash
python3 scripts/test_flare_model.py
python3 scripts/test_flare_split_samples.py
python3 scripts/test_flare_switch_qc.py
python3 scripts/test_flare_site_stats.py
python3 scripts/test_flare_build_call_qc_sites.py
python3 scripts/test_flare_lai_association_metrics.py
```

Synthetic full-format model fixtures live under
`scripts/testdata/flare_models/` (replace with GCS downloads when credentials
allow). Expected switch rate uses `(T/100) * (1 - Σ p²)`; see
`scripts/flare_switch_qc.py`.

Site-restriction addendum (Parts 4–5), FLARE2 (Part 6), call-QC GQ/DP/RNC
(Part 7), and **association-facing recipe selection (Part 8)** — see
[`flare_lai_fix_instructions.md`](flare_lai_fix_instructions.md).

**Scoring rule (Part 8):** switch / tract metrics are diagnostics only.
Recipes must clear allele–ancestry concordance **and** pedigree Mendelian
gates; among survivors, Tractor null-λ (closest to 1) decides the winner
(`scripts/flare_score_allele_ancestry.py`,
`scripts/flare_score_mendelian_lai.py`,
`scripts/flare_lai_null_lambda.py`; notebook `flare_02` Part 8).

Null-λ phenotypes are **directly simulated**
(`L = β·q + γ_GC + u_family + ε`, thresholded to an anchor case rate) — not
shuffled. Run `flare_lai_null_lambda.py preflight` first: (1) check that fixed
global ancestry `q` is not absorbed by the null PC set (drop PCs for this
pilot if multivariate R² is too high); (2) require a non-trivial anchor
phenotype log-OR before locking β_mid. Score at lo/mid/hi β; trust a winner
only if ranking is stable across magnitudes.

Stage filtered-site helpers (`flare_site_stats.py`,
`flare_build_indel_flanks.py`, `flare_build_call_qc_sites.py`) and panel
builder (`flare_build_af_panel.py`) with the other scripts.
Diagnostic per-Mb / per-marker rates still use
`flare_lai_exp.enrich_compare_row`.
