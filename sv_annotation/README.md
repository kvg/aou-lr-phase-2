# SV annotation pipeline (AoU-LR Phase 1 / Phase 2)

Manuscript-first pipeline that:

1. Normalizes Phase 1 (main) and Phase 2 (main + breakends + large-event) SV callsets
2. Builds a unified site table with AF/carriers, repetitive-union context, and CADD-SV bins
3. Emits the Phase 1 vs Phase 2 **total site** comparison table
4. Emits ancestry-ordered Ebert-style cumulative discovery tables (optionally stratified)

See [docs/definitions.md](docs/definitions.md) for locked counting rules.

## Layout

```
notebooks/terra/sv_00_prepare_caddsv_annotations.ipynb
notebooks/terra/sv_01_stage_repeat_tracks.ipynb
notebooks/terra/sv_02_stage_sample_ancestry.ipynb
notebooks/terra/sv_03_manuscript_stats.ipynb
scripts/                         # shared with Tractor-Mix; baked into the Docker image
sv_annotation/
  docs/definitions.md
  docs/data_dictionary.md
  schemas/
  configs/
  wdl/AnnotateSvCallset.wdl
  tests/
```

## Local smoke test

```bash
cd sv_annotation
python3 -m pytest tests/ -q
```

## Terra / Cromwell

1. Build / push the Docker image (same Artifact Registry as Tractor-Mix):

```bash
cd sv_annotation
./build_docker.sh                 # Cloud Build → aou-sv-annotation:0.1.6
# or: ./build_docker.sh --local --push
```

2. Build the CADD-SV annotation bundle once on a machine/Workbench with ≥60 GB
   free (Docker Desktop local disks are usually too small). Prefer the Workbench
   notebook [`../notebooks/terra/sv_00_prepare_caddsv_annotations.ipynb`](../notebooks/terra/sv_00_prepare_caddsv_annotations.ipynb),
   which uploads to `$WORKSPACE_BUCKET/refs/caddsv/`. Alternatively run
   `../scripts/prepare_caddsv_annotations.sh` on a large VM.
3. Stage RepeatMasker / simple-repeat / segdup BEDs with
   [`../notebooks/terra/sv_01_stage_repeat_tracks.ipynb`](../notebooks/terra/sv_01_stage_repeat_tracks.ipynb)
   (`$WORKSPACE_BUCKET/refs/grch38/`).
4. Stage sample ancestry TSVs with
   [`../notebooks/terra/sv_02_stage_sample_ancestry.ipynb`](../notebooks/terra/sv_02_stage_sample_ancestry.ipynb)
   (`$WORKSPACE_BUCKET/metadata/`). Prefer setting `SV_PHASE1_MAIN_VCF` /
   `SV_PHASE2_MAIN_VCF` so IDs match the callsets.
5. Import `wdl/AnnotateSvCallset.wdl` into a Terra Method.
6. Fill `configs/phase2.inputs.json.example` (and phase1) with real GCS paths:
   main / bnd / large VCF (or BCF) callsets, repeat BEDs, CADD-SV tar,
   ancestry TSV, and the docker image.
   The container runs `bcftools view`, so **BCF is fine**.
7. Run Phase 2 with `main_vcf`, `bnd_vcf`, `large_vcf`; Phase 1 with `main_vcf` only.
   CADD-SV runs and attaches scores during this single workflow submission.
8. Merge phase count TSVs with `../scripts/merge_manuscript_counts.py`.
9. Plot discovery curves with `../scripts/plot_discovery.py` or
   `../notebooks/terra/sv_03_manuscript_stats.ipynb`.

## Core script chain

```text
extract_sites_from_vcf.py
  -> fill_af_and_carriers.py
  -> annotate_repeat_context.py
  -> prepare_caddsv_input.py -> run_caddsv.py -> attach_caddsv_scores.py
  -> integrate_partitions.py   # large-over-main precedence
  -> manuscript_site_counts.py
  -> ebert_discovery.py
  -> plot_discovery.py
```

Optional landscape producers (TR / remap / MEI / external AF) format TSVs for later
`bcftools annotate` transfer; they are not required for the main table or Ebert plots.

Repeat-unit dosage (FELIX Score B path) is a separate VCF-INFO annotator:

```bash
python3 ../scripts/annotate_repeat_units.py \
  --vcf JOINT.vcf.gz \
  --simple-repeat-bed simpleRepeat.bed.gz \
  --out JOINT.ru.vcf.gz
```
