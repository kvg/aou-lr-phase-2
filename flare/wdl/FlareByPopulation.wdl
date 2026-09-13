version 1.0

# Per-population FLARE on one chromosome, then merge sample columns back.
#
# The existing popout flare.wdl runs all 12k AoU samples as one "admixed
# population". FLARE then estimates a single T (generations since admixture)
# for everyone; EM on this cohort produced T≈120 and short tracts. This
# workflow splits the gt VCF by covariates `population`, runs FLARE with its
# own EM (+ update-p) in each group, and bcftools-merges AN1/AN2.
#
# Based on https://github.com/broadinstitute/popout/blob/main/workflows/flare/wdl/flare.wdl
# (resource scaling, FLARE docker, command-line flags). Magicwand/W&B omitted.
#
#   covariates     = covariates.source_rebuilt.csv.gz  (research_id, population)
#   gt_vcf         = phased per-chrom target VCF (no missing GT)
#   ref_vcf        = gnomAD LAI (or other) reference VCF for this chrom
#   ref_panel      = sample <tab> panel
#   map_file       = PLINK cM map; chromosome IDs must match the VCF (chr1 not chrchr1)
#   region         = optional chr:start-end smoke test (subset gt AND ref)
#   split_script   = gs://$WORKSPACE_BUCKET/scripts/flare_split_samples.py
#   summarize_script = gs://$WORKSPACE_BUCKET/scripts/flare_summarize_models.py
#   flare_model_script = gs://$WORKSPACE_BUCKET/scripts/flare_model.py
#   flare_site_stats_script = gs://$WORKSPACE_BUCKET/scripts/flare_site_stats.py
#
# Target gt_vcf is a phased VCF with no missing genotypes (AoU Phase 2
# ligation). A call-rate filter would be a no-op here; SubsetOnePopulation
# verifies F_MISSING==0 on the shard and skips implementing one.
#
# Scatter is per keep-list: each shard subsets that population from the gt VCF
# (streaming gs://) and runs FLARE. Full-chrom stream is ~90 min. bcftools
# 1.24 / htslib 1.23+ reconnects mid-transfer on GCS HTTP/2 GOAWAY (curl 92;
# PR 1987). Shell + Cromwell retries remain as a backstop. When region is set,
# the reference is sliced once and shared. Cutting only gt still leaves FLARE
# loading a full-chrom ref.
#
# Fixed seed=12345 ⇒ identical outputs show determinism, not independent
# corroboration across "replicate" submissions.

workflow FlareByPopulation {
  input {
    File covariates
    File gt_vcf
    File gt_vcf_index
    File ref_vcf
    File ref_vcf_index
    File ref_panel
    File map_file
    File split_script
    File summarize_script
    File flare_model_script
    File flare_site_stats_script
    String output_prefix

    String id_column = "research_id"
    String pop_column = "ancestry_pred"
    String include_pops = ""
    String drop_pops = ""
    Int min_samples = 50
    String unlabeled_pop = ""
    Float max_unmatched_vcf_frac = 0.02
    Boolean exclude_controls = true
    String region = ""

    # Train / explore: em=true estimates T (gen is only the EM start).
    # Pin T for association-grade LAI: em=false, pop_models empty, set gen
    # (and optional gen_by_pop). Apply genome-wide: em=false + pop_models from
    # a prior run (basenames must end in .<POP>.model).
    # With pop_models, FLARE always runs em=false. Model blocks are assembled by
    # flare_model.py (T always rewritten; props / panel weights / miscopy µ
    # controlled independently via props_by_pop + inherit_*).
    Boolean em = true
    Array[File] pop_models = []
    File? template_model
    # Per-pop admixture proportions, e.g. "AFR:AFR=0.80,EUR=0.20;EAS:EAS=0.97"
    String props_by_pop = ""
    Boolean inherit_props = true
    Boolean inherit_panel_weights = true
    Boolean inherit_mu = true
    Boolean update_p = true
    # false ⇒ no ANP posterior dosages. Diagnostic runs that need global props
    # from AN1/AN2 lengths do not need probs; set true only when ANP is required.
    Boolean probs = false
    Boolean array = false
    Float min_maf = 0.005
    Int min_mac = 50
    Float gen = 10.0
    # Optional per-population T overrides, e.g. "AFR:8,AMR:12". Blank = use gen
    # for every shard. Also applied when rewriting T inside a matched pop_model.
    String gen_by_pop = ""
    # If true, missing gen_by_pop keys fall back to gen; if false, every shard
    # pop must appear when gen_by_pop is non-empty.
    Boolean gen_by_pop_allow_default = false
    # If false, fail when a shard pop has no same-named panel in ref_panel
    # (e.g. MID without a matching ancestry). If true, warn only.
    # Default true for now: MID is kept as a shard without a MID panel label;
    # revisit with a dedicated panel / remap later.
    Boolean allow_unrepresented_pops = true
    Int seed = 12345

    # Part 4 site restriction — applied identically to target and reference.
    Boolean biallelic_snvs_only = false
    File? include_sites
    File? exclude_regions
    Int indel_flank_bp = 0
    File indel_flank_script
    # Fail if post-filter density < this (sites per Mb of region span).
    Float min_sites_per_mb = 1000.0
    # Part 4.2: target-side filters after sample subset (0 = disabled).
    Float target_min_maf = 0.0
    Float target_hwe_pval = 0.0

    String gcs_project = ""
    String bcftools_docker = "us-central1-docker.pkg.dev/broad-dsp-lrma/aou-lr/aou-flare-bcftools:1.24"
    String flare_docker = "us-docker.pkg.dev/broad-dsde-methods/popout/flare:latest"

    Int split_cpu = 2
    Int split_memory_gb = 4
    Int subset_cpu = 8
    Int subset_memory_gb = 16
    Int subset_disk_gb_floor = 80
    Float subset_disk_gb_multiplier = 0.0
    # Full-chrom stream ~90 min. htslib 1.23+ retries curl 92 mid-stream;
    # subset_preemptible covers GCE reclaiming the VM. maxRetries is a backstop.
    Int subset_stream_attempts = 3
    Int max_retries = 5
    Int preemptible = 1
    Int subset_preemptible = 5
    Int flare_preemptible = 0

    Int? flare_cpu_override
    # Runtime memory string (e.g. "96 GB"). Does not change -Xmx by itself —
    # set flare_xmx_gb_override together, or rely on probs-aware auto sizing.
    String? flare_memory_override
    Int? flare_xmx_gb_override
    Int? flare_disk_gb_override
    String flare_disk_type = "HDD"
  }

  call SplitSamples {
    input:
      covariates = covariates,
      gt_vcf = gt_vcf,
      gt_vcf_index = gt_vcf_index,
      split_script = split_script,
      id_column = id_column,
      pop_column = pop_column,
      include_pops = include_pops,
      drop_pops = drop_pops,
      min_samples = min_samples,
      unlabeled_pop = unlabeled_pop,
      max_unmatched_vcf_frac = max_unmatched_vcf_frac,
      exclude_controls = exclude_controls,
      gcs_project = gcs_project,
      docker = bcftools_docker,
      cpu = split_cpu,
      memory_gb = split_memory_gb,
      preemptible = preemptible
  }

  call ValidateFlareInputs {
    input:
      keep_lists = SplitSamples.keep_lists,
      ref_panel = ref_panel,
      gen_by_pop = gen_by_pop,
      gen_by_pop_allow_default = gen_by_pop_allow_default,
      allow_unrepresented_pops = allow_unrepresented_pops,
      docker = bcftools_docker,
      preemptible = preemptible
  }

  # Merge user exclude BED with optional indel ± flank BED from the target so
  # ref and gt share one effective exclude set (Part 4.1).
  call PrepExcludeRegions {
    input:
      exclude_regions = exclude_regions,
      indel_flank_bp = indel_flank_bp,
      gt_vcf = gt_vcf,
      gt_vcf_index = gt_vcf_index,
      region = region,
      indel_flank_script = indel_flank_script,
      gcs_project = gcs_project,
      docker = bcftools_docker,
      cpu = subset_cpu,
      memory_gb = subset_memory_gb,
      disk_gb_floor = subset_disk_gb_floor,
      preemptible = preemptible,
      max_retries = max_retries
  }

  call HashSiteFilters {
    input:
      biallelic_snvs_only = biallelic_snvs_only,
      include_sites = include_sites,
      effective_exclude = PrepExcludeRegions.effective_exclude,
      docker = bcftools_docker,
      preemptible = preemptible
  }

  # Prep ref whenever a region is set OR any site filter is active so gt and
  # ref always see the same site treatment.
  Boolean prep_ref = (region != "")
    || biallelic_snvs_only
    || defined(include_sites)
    || (PrepExcludeRegions.has_exclude)
    || (indel_flank_bp > 0)

  if (prep_ref) {
    call CutRefRegion {
      input:
        ref_vcf = ref_vcf,
        ref_vcf_index = ref_vcf_index,
        region = region,
        biallelic_snvs_only = biallelic_snvs_only,
        include_sites = include_sites,
        exclude_regions = PrepExcludeRegions.effective_exclude,
        site_filter_digest = HashSiteFilters.digest,
        min_sites_per_mb = min_sites_per_mb,
        gcs_project = gcs_project,
        docker = bcftools_docker,
        cpu = subset_cpu,
        memory_gb = subset_memory_gb,
        disk_gb_floor = subset_disk_gb_floor,
        stream_attempts = subset_stream_attempts,
        preemptible = preemptible,
        max_retries = max_retries
    }
  }
  File flare_ref_vcf = select_first([CutRefRegion.vcf, ref_vcf])

  scatter (keep_list in SplitSamples.keep_lists) {
    call SubsetOnePopulation {
      input:
        gt_vcf = gt_vcf,
        gt_vcf_index = gt_vcf_index,
        keep_list = keep_list,
        region = region,
        biallelic_snvs_only = biallelic_snvs_only,
        include_sites = include_sites,
        exclude_regions = PrepExcludeRegions.effective_exclude,
        site_filter_digest = HashSiteFilters.digest,
        min_sites_per_mb = min_sites_per_mb,
        target_min_maf = target_min_maf,
        target_hwe_pval = target_hwe_pval,
        gcs_project = gcs_project,
        docker = bcftools_docker,
        cpu = subset_cpu,
        memory_gb = subset_memory_gb,
        disk_gb_floor = subset_disk_gb_floor,
        disk_gb_multiplier = subset_disk_gb_multiplier,
        stream_attempts = subset_stream_attempts,
        preemptible = subset_preemptible,
        max_retries = max_retries
    }
    call flare_task {
      input:
        ref_vcf = flare_ref_vcf,
        ref_panel = ref_panel,
        gt_vcf = SubsetOnePopulation.subset_vcf,
        map_file = map_file,
        output_prefix = output_prefix,
        em = em,
        pop_models = pop_models,
        template_model = template_model,
        props_by_pop = props_by_pop,
        inherit_props = inherit_props,
        inherit_panel_weights = inherit_panel_weights,
        inherit_mu = inherit_mu,
        flare_model_script = flare_model_script,
        flare_site_stats_script = flare_site_stats_script,
        update_p = update_p,
        probs = probs,
        array = array,
        min_maf = min_maf,
        min_mac = min_mac,
        gen = gen,
        gen_by_pop = gen_by_pop,
        seed = seed,
        validation_ok = ValidateFlareInputs.ok,
        cpu_override = flare_cpu_override,
        memory_override = flare_memory_override,
        xmx_gb_override = flare_xmx_gb_override,
        disk_size_gb_override = flare_disk_gb_override,
        disk_type = flare_disk_type,
        preemptible = flare_preemptible,
        docker_image = flare_docker
    }
  }

  call ConcatSubsetCounts {
    input:
      tables = SubsetOnePopulation.counts_tsv,
      docker = bcftools_docker,
      preemptible = preemptible
  }

  call MergeSiteStats {
    input:
      tables = flare_task.site_stats,
      site_classes = SubsetOnePopulation.site_classes,
      site_stats_script = flare_site_stats_script,
      output_prefix = output_prefix,
      docker = bcftools_docker,
      preemptible = preemptible
  }

  call MergePopulationFlare {
    input:
      anc_vcfs = flare_task.anc_vcf,
      global_ancs = flare_task.global_anc,
      models = flare_task.out_model,
      pops = flare_task.population,
      summarize_script = summarize_script,
      flare_model_script = flare_model_script,
      output_prefix = output_prefix,
      docker = bcftools_docker,
      cpu = subset_cpu,
      preemptible = preemptible
  }

  output {
    File anc_vcf = MergePopulationFlare.anc_vcf
    File anc_vcf_index = MergePopulationFlare.anc_vcf_index
    File global_anc = MergePopulationFlare.global_anc
    File models_tsv = MergePopulationFlare.models_tsv
    File site_stats_tsv = MergeSiteStats.site_stats
    Array[File] per_pop_anc_vcf = flare_task.anc_vcf
    Array[File] per_pop_model = flare_task.out_model
    Array[File] per_pop_in_model = flare_task.in_model
    Array[File] per_pop_log = flare_task.log
    Array[File] per_pop_site_stats = flare_task.site_stats
    File sample_manifest = SplitSamples.manifest
    File sample_summary_json = SplitSamples.summary_json
    File unmatched_vcf_samples = SplitSamples.unmatched_vcf_samples
    File excluded_controls = SplitSamples.excluded_controls
    File subset_counts = ConcatSubsetCounts.subset_counts
    File validation_log = ValidateFlareInputs.log
    File site_filter_digest = HashSiteFilters.digest_file
    File effective_exclude = PrepExcludeRegions.effective_exclude
  }

  meta {
    description: "FLARE per covariates population, then merge ancestry VCFs."
    allowNestedInputs: true
  }

  parameter_meta {
    probs: "false ⇒ no ANP posterior dosages. Diagnostic runs that only need AN1/AN2 (and optional global props from allele counts) can leave this false; set true only when ANP is required. probs=true raises flare_task heap floors (full-chrom OOM otherwise)."
    seed: "Fixed seed (default 12345) makes identical inputs produce identical outputs — that shows determinism, not independent corroboration."
    flare_model_script: "scripts/flare_model.py — parse/write .model blocks and assemble per-shard inputs under in_model/."
    flare_memory_override: "Optional Cromwell memory string (e.g. \"96 GB\"). Pair with flare_xmx_gb_override so -Xmx tracks the VM."
    flare_xmx_gb_override: "Optional -Xmx in GB for flare.jar. Defaults to probs-aware auto sizing from subset gt size."
  }
}

task SplitSamples {
  input {
    File covariates
    File gt_vcf
    File gt_vcf_index
    File split_script
    String id_column
    String pop_column
    String include_pops
    String drop_pops
    Int min_samples
    String unlabeled_pop
    Float max_unmatched_vcf_frac
    Boolean exclude_controls
    String gcs_project
    String docker
    Int cpu
    Int memory_gb
    Int preemptible
  }

  parameter_meta {
    gt_vcf: {
      description: "Phased target VCF. Header-only read on Terra when streamed.",
      localization_optional: true
    }
    gt_vcf_index: {
      description: "Tabix/CSI sibling of gt_vcf.",
      localization_optional: true
    }
  }

  command <<<
    set -euo pipefail
    python3 - <<'PY'
import os
from pathlib import Path

def export_token():
    uri = "~{gt_vcf}"
    if not uri.startswith("gs://"):
        return
    token = os.environ.get("GCS_OAUTH_TOKEN", "")
    if not token:
        try:
            import subprocess
            token = subprocess.check_output(
                ["gcloud", "auth", "application-default", "print-access-token"],
                text=True,
                stderr=subprocess.DEVNULL,
            ).strip()
        except Exception:
            token = ""
    if not token:
        import json, urllib.request
        for url in (
            "http://metadata.google.internal/computeMetadata/v1/instance/service-accounts/default/token",
            "http://169.254.169.254/computeMetadata/v1/instance/service-accounts/default/token",
        ):
            req = urllib.request.Request(url, headers={"Metadata-Flavor": "Google"})
            try:
                with urllib.request.urlopen(req, timeout=5) as resp:
                    token = json.load(resp)["access_token"]
                    break
            except Exception:
                continue
    if not token:
        raise SystemExit("error: bcftools needs GCS_OAUTH_TOKEN to read ~{gt_vcf}")
    Path("gcs.env").write_text(f"export GCS_OAUTH_TOKEN='{token}'\n")

project = "~{gcs_project}".strip() or os.environ.get("GOOGLE_PROJECT", "")
lines = []
if project:
    lines.append(f"export GCS_REQUESTER_PAYS_PROJECT={project}")
Path("gcs.env").write_text("")
export_token()
if lines:
    Path("gcs.env").write_text("\n".join(lines) + "\n" + Path("gcs.env").read_text())
PY
    # shellcheck disable=SC1091
    if [[ -s gcs.env ]]; then
      source gcs.env
    fi

    bcftools query -l "~{gt_vcf}" > vcf_samples.txt
    python3 "~{split_script}" \
      --covariates "~{covariates}" \
      --vcf-samples vcf_samples.txt \
      --out-dir split \
      --id-column "~{id_column}" \
      --pop-column "~{pop_column}" \
      --include-pops "~{include_pops}" \
      --drop-pops "~{drop_pops}" \
      --min-samples ~{min_samples} \
      --unlabeled-pop "~{unlabeled_pop}" \
      --max-unmatched-vcf-frac ~{max_unmatched_vcf_frac} \
      ~{true="--exclude-controls" false="--keep-controls" exclude_controls}
    test -s split/populations.txt
    ls -l split/keeps/*.samples.txt
  >>>

  output {
    Array[File] keep_lists = glob("split/keeps/*.samples.txt")
    File populations = "split/populations.txt"
    File manifest = "split/manifest.tsv"
    File summary_json = "split/summary.json"
    File unmatched_vcf_samples = "split/unmatched_vcf_samples.txt"
    File excluded_controls = "split/excluded_controls.txt"
    File dropped_samples = "split/dropped_samples.tsv"
    File vcf_samples = "vcf_samples.txt"
  }

  runtime {
    docker: docker
    cpu: cpu
    memory: memory_gb + " GB"
    disks: "local-disk 20 HDD"
    preemptible: preemptible
  }
}

# Merge user exclude BED with optional indel ± flank intervals from the target.
task PrepExcludeRegions {
  input {
    File? exclude_regions
    Int indel_flank_bp
    File gt_vcf
    File gt_vcf_index
    String region
    File indel_flank_script
    String gcs_project
    String docker
    Int cpu
    Int memory_gb
    Int disk_gb_floor
    Int preemptible
    Int max_retries
  }

  Boolean has_user_exclude = defined(exclude_regions)
  File user_exclude_or_dummy = select_first([exclude_regions, indel_flank_script])

  parameter_meta {
    gt_vcf: {
      description: "Phased target VCF for indel flank discovery.",
      localization_optional: true
    }
    gt_vcf_index: {
      description: "Tabix/CSI sibling of gt_vcf.",
      localization_optional: true
    }
  }

  command <<<
    set -euo pipefail
    export HTS_RETRY_MAX="${HTS_RETRY_MAX:-8}"
    export HTS_RETRY_DELAY="${HTS_RETRY_DELAY:-500}"
    export HTS_RETRY_MAX_DELAY="${HTS_RETRY_MAX_DELAY:-60000}"
    cat > gcs_auth.py <<'PY'
import json, os, subprocess, sys, urllib.request
from pathlib import Path

def fresh_token():
    for cmd in (
        ["gcloud", "auth", "print-access-token"],
        ["gcloud", "auth", "application-default", "print-access-token"],
    ):
        try:
            tok = subprocess.check_output(cmd, text=True, stderr=subprocess.DEVNULL).strip()
            if tok:
                return tok
        except Exception:
            pass
    for url in (
        "http://metadata.google.internal/computeMetadata/v1/instance/service-accounts/default/token",
        "http://169.254.169.254/computeMetadata/v1/instance/service-accounts/default/token",
    ):
        req = urllib.request.Request(url, headers={"Metadata-Flavor": "Google"})
        try:
            with urllib.request.urlopen(req, timeout=5) as resp:
                tok = json.load(resp).get("access_token") or ""
                if tok:
                    return tok
        except Exception:
            continue
    return ""

need_gcs = any(u.startswith("gs://") for u in sys.argv[1:])
token = fresh_token() if need_gcs else ""
if need_gcs and not token:
    raise SystemExit("error: bcftools needs a GCS OAuth token to read gs:// VCFs")
project = "~{gcs_project}".strip() or os.environ.get("GOOGLE_PROJECT", "")
parts = []
if project:
    parts.append(f"export GCS_REQUESTER_PAYS_PROJECT={project}")
if token:
    parts.append(f"export GCS_OAUTH_TOKEN='{token}'")
Path("gcs.env").write_text(("\n".join(parts) + "\n") if parts else "")
PY
    GT="~{gt_vcf}"
    GT_IDX="~{gt_vcf_index}"
    refresh_gcs_auth() {
      python3 gcs_auth.py "${GT}"
      if [[ -s gcs.env ]]; then
        # shellcheck disable=SC1091
        source gcs.env
      fi
      echo "$(date -u +%Y-%m-%dT%H:%M:%SZ) refreshed GCS auth" >&2
    }
    refresh_gcs_auth

    python3 -c "import pathlib,sys; pathlib.Path('region.txt').write_text(sys.argv[1].strip().strip(chr(34)).strip(chr(39)) + chr(10))" "~{region}"
    REGION="$(cat region.txt)"

    : > effective_exclude.bed
    if [[ "~{has_user_exclude}" == "true" && -s "~{user_exclude_or_dummy}" ]]; then
      cat "~{user_exclude_or_dummy}" >> effective_exclude.bed
      echo "merged user exclude_regions ($(wc -l < effective_exclude.bed) lines so far)" >&2
    fi

    if [[ ~{indel_flank_bp} -gt 0 ]]; then
      refresh_gcs_auth
      python3 "~{indel_flank_script}" \
        --vcf "${GT}##idx##${GT_IDX}" \
        --region "${REGION}" \
        --flank ~{indel_flank_bp} \
        --out indel_flanks.bed
      if [[ -s indel_flanks.bed ]]; then
        cat indel_flanks.bed >> effective_exclude.bed
        echo "merged indel flanks (±~{indel_flank_bp} bp)" >&2
      fi
    fi

    if [[ -s effective_exclude.bed ]]; then
      echo "true" > has_exclude.txt
      echo "effective_exclude.bed lines=$(wc -l < effective_exclude.bed)" >&2
    else
      echo "false" > has_exclude.txt
      : > effective_exclude.bed
      echo "effective_exclude.bed empty" >&2
    fi
    ls -lh effective_exclude.bed has_exclude.txt >&2
  >>>

  output {
    File effective_exclude = "effective_exclude.bed"
    Boolean has_exclude = read_boolean("has_exclude.txt")
  }

  runtime {
    docker: docker
    cpu: cpu
    memory: memory_gb + " GB"
    disks: "local-disk " + disk_gb_floor + " HDD"
    preemptible: preemptible
    maxRetries: max_retries
  }
}

# Digest gate so CutRefRegion and SubsetOnePopulation share identical site filters.
task HashSiteFilters {
  input {
    Boolean biallelic_snvs_only
    File? include_sites
    File effective_exclude
    String docker
    Int preemptible
  }

  Boolean has_include = defined(include_sites)
  File include_or_dummy = select_first([include_sites, effective_exclude])

  command <<<
    set -euo pipefail
    python3 - <<'PY'
import hashlib
from pathlib import Path

biallelic = "~{biallelic_snvs_only}".lower() == "true"
has_include = "~{has_include}" == "true"
include_path = Path("~{include_or_dummy}")
exclude_path = Path("~{effective_exclude}")

if has_include:
    include_hash = hashlib.sha256(include_path.read_bytes()).hexdigest()
else:
    include_hash = hashlib.sha256(b"NONE").hexdigest()
exclude_hash = hashlib.sha256(exclude_path.read_bytes()).hexdigest()
digest_input = f"{str(biallelic).lower()}:{include_hash}:{exclude_hash}"
digest = hashlib.sha256(digest_input.encode()).hexdigest()
Path("digest.txt").write_text(digest + "\n")
Path("digest_meta.txt").write_text(
    f"biallelic={biallelic}\ninclude={include_hash}\nexclude={exclude_hash}\n"
)
print("site_filter_digest", digest, flush=True)
PY
    cat digest.txt digest_meta.txt >&2
  >>>

  output {
    String digest = read_string("digest.txt")
    File digest_file = "digest.txt"
  }

  runtime {
    docker: docker
    cpu: 1
    memory: "2 GB"
    disks: "local-disk 10 HDD"
    preemptible: preemptible
  }
}

# Shared GCS auth + contig aliasing for streamed bcftools. A Cromwell
# GCS_OAUTH_TOKEN dies after ~1 h; mint a fresh one (gcloud/metadata) instead
# of reusing the env token. One shard = one bcftools open, so refresh at start.
task CutRefRegion {
  input {
    File ref_vcf
    File ref_vcf_index
    String region
    Boolean biallelic_snvs_only = false
    File? include_sites
    File exclude_regions
    String site_filter_digest
    Float min_sites_per_mb = 1000.0
    String gcs_project
    String docker
    Int cpu
    Int memory_gb
    Int disk_gb_floor
    Int stream_attempts
    Int preemptible
    Int max_retries
  }

  Boolean has_include = defined(include_sites)
  File include_or_dummy = select_first([include_sites, exclude_regions])

  parameter_meta {
    ref_vcf: {
      description: "LAI reference VCF; tabix-sliced, not fully localized.",
      localization_optional: true
    }
    ref_vcf_index: {
      description: "Tabix/CSI sibling of ref_vcf.",
      localization_optional: true
    }
  }

  command <<<
    set -euo pipefail
    export HTS_RETRY_MAX="${HTS_RETRY_MAX:-8}"
    export HTS_RETRY_DELAY="${HTS_RETRY_DELAY:-500}"
    export HTS_RETRY_MAX_DELAY="${HTS_RETRY_MAX_DELAY:-60000}"
    echo "$(date -u +%Y-%m-%dT%H:%M:%SZ) $(bcftools --version 2>/dev/null | head -n 1) HTS_RETRY_MAX=${HTS_RETRY_MAX}" >&2
    cat > gcs_auth.py <<'PY'
import json, os, subprocess, sys, urllib.request
from pathlib import Path

def fresh_token():
    for cmd in (
        ["gcloud", "auth", "print-access-token"],
        ["gcloud", "auth", "application-default", "print-access-token"],
    ):
        try:
            tok = subprocess.check_output(cmd, text=True, stderr=subprocess.DEVNULL).strip()
            if tok:
                return tok
        except Exception:
            pass
    for url in (
        "http://metadata.google.internal/computeMetadata/v1/instance/service-accounts/default/token",
        "http://169.254.169.254/computeMetadata/v1/instance/service-accounts/default/token",
    ):
        req = urllib.request.Request(url, headers={"Metadata-Flavor": "Google"})
        try:
            with urllib.request.urlopen(req, timeout=5) as resp:
                tok = json.load(resp).get("access_token") or ""
                if tok:
                    return tok
        except Exception:
            continue
    return ""

need_gcs = any(u.startswith("gs://") for u in sys.argv[1:])
token = fresh_token() if need_gcs else ""
if need_gcs and not token:
    raise SystemExit("error: bcftools needs a GCS OAuth token to read gs:// VCFs")
project = "~{gcs_project}".strip() or os.environ.get("GOOGLE_PROJECT", "")
parts = []
if project:
    parts.append(f"export GCS_REQUESTER_PAYS_PROJECT={project}")
if token:
    parts.append(f"export GCS_OAUTH_TOKEN='{token}'")
Path("gcs.env").write_text(("\n".join(parts) + "\n") if parts else "")
PY
    SRC="~{ref_vcf}"
    IDX="~{ref_vcf_index}"
    refresh_gcs_auth() {
      python3 gcs_auth.py "${SRC}"
      if [[ -s gcs.env ]]; then
        # shellcheck disable=SC1091
        source gcs.env
      fi
      echo "$(date -u +%Y-%m-%dT%H:%M:%SZ) refreshed GCS auth" >&2
    }
    bcftools_view_retry() {
      local out="$1"
      shift
      local attempt=1
      local max=~{stream_attempts}
      local delay=20
      if [[ "${max}" -lt 1 ]]; then
        max=1
      fi
      while [[ "${attempt}" -le "${max}" ]]; do
        refresh_gcs_auth
        echo "$(date -u +%Y-%m-%dT%H:%M:%SZ) bcftools view attempt ${attempt}/${max} -> ${out}" >&2
        rm -f "${out}"
        if bcftools view -o "${out}" "$@"; then
          return 0
        fi
        echo "$(date -u +%Y-%m-%dT%H:%M:%SZ) bcftools view failed (attempt ${attempt}/${max})" >&2
        if [[ "${attempt}" -eq "${max}" ]]; then
          return 1
        fi
        sleep "${delay}"
        delay=$((delay * 2))
        attempt=$((attempt + 1))
      done
    }
    refresh_gcs_auth

    python3 - <<'PY'
import hashlib
from pathlib import Path

biallelic = "~{biallelic_snvs_only}".lower() == "true"
has_include = "~{has_include}" == "true"
include_path = Path("~{include_or_dummy}")
exclude_path = Path("~{exclude_regions}")
if has_include:
    include_hash = hashlib.sha256(include_path.read_bytes()).hexdigest()
else:
    include_hash = hashlib.sha256(b"NONE").hexdigest()
exclude_hash = hashlib.sha256(exclude_path.read_bytes()).hexdigest()
digest_input = f"{str(biallelic).lower()}:{include_hash}:{exclude_hash}"
computed = hashlib.sha256(digest_input.encode()).hexdigest()
expected = "~{site_filter_digest}".strip()
if computed != expected:
    raise SystemExit(f"site filter digest mismatch: {computed} != {expected}")
print("site filter digest ok", computed, flush=True)
PY

    python3 -c "import pathlib,sys; pathlib.Path('region.txt').write_text(sys.argv[1].strip().strip(chr(34)).strip(chr(39)) + chr(10))" "~{region}"
    REGION="$(cat region.txt)"

    SITE_FILTER_ARGS=()
    if [[ "~{has_include}" == "true" ]]; then
      SITE_FILTER_ARGS+=(-T "~{include_or_dummy}")
    fi
    if [[ -s "~{exclude_regions}" ]]; then
      SITE_FILTER_ARGS+=(-T "^~{exclude_regions}")
    fi
    BIALLELIC_ARGS=()
    if [[ "~{biallelic_snvs_only}" == "true" ]]; then
      BIALLELIC_ARGS=(-m2 -M2 -v snps)
    fi
    site_filters_active=false
    if [[ "~{has_include}" == "true" ]] || [[ -s "~{exclude_regions}" ]] || [[ "~{biallelic_snvs_only}" == "true" ]]; then
      site_filters_active=true
    fi

    REF_REGION=""
    if [[ -n "${REGION}" ]]; then
      bcftools view -h "${SRC}##idx##${IDX}" > ref.header
      python3 - "${REGION}" ref.header ref.region <<'PY'
import re, sys
from pathlib import Path

region, header_path, out_path = sys.argv[1], sys.argv[2], sys.argv[3]
header = Path(header_path).read_text(errors="replace")
contigs = []
for line in header.splitlines():
    if line.startswith("##contig=<"):
        m = re.search(r"ID=([^,>]+)", line)
        if m:
            contigs.append(m.group(1))
print(f"{Path(header_path).name} contigs: {', '.join(contigs[:20])}{' …' if len(contigs) > 20 else ''}", file=sys.stderr)
if not contigs:
    raise SystemExit(f"{header_path}: no ##contig=<ID=...> in header")
chrom, rest = (region.split(":", 1) + [""])[:2]
aliases = [chrom]
if chrom.startswith("chr"):
    aliases.append(chrom[3:])
else:
    aliases.append("chr" + chrom)
matched = next((a for a in aliases if a in set(contigs)), None)
if matched is None:
    raise SystemExit(
        f"region chrom {chrom!r} is not in this VCF. Contigs: {contigs}. "
        "Use the gnomAD ref for the SAME chromosome as region."
    )
resolved = matched if not rest else f"{matched}:{rest}"
if resolved != region:
    print(f"rewrote region {region} -> {resolved}", file=sys.stderr)
Path(out_path).write_text(resolved + "\n")
print(resolved)
PY
      REF_REGION="$(cat ref.region)"
    fi

    VIEW_SRC="${SRC}##idx##${IDX}"
    n_pre=0
    if [[ -n "${REF_REGION}" ]]; then
      echo "Extract ref region ${REF_REGION} (pre site filters) from ${SRC}" >&2
      bcftools_view_retry region_only.vcf.gz -r "${REF_REGION}" -Oz --threads ~{cpu} \
        "${SRC}##idx##${IDX}"
      bcftools index -t region_only.vcf.gz
      n_pre=$(bcftools index -n region_only.vcf.gz 2>/dev/null || echo 0)
      VIEW_SRC="region_only.vcf.gz"
    elif [[ "${site_filters_active}" == "true" ]]; then
      n_pre=$(bcftools index -n "${SRC}##idx##${IDX}" 2>/dev/null || echo 0)
      if [[ "${n_pre}" == "0" ]]; then
        n_pre=$(bcftools view -H --threads 1 "${SRC}##idx##${IDX}" | wc -l | tr -d ' ')
      fi
    fi

    echo "Apply ref site filters region=${REF_REGION:-FULL} include=~{has_include} exclude_nonempty=$([[ -s '~{exclude_regions}' ]] && echo yes || echo no) biallelic=~{biallelic_snvs_only}" >&2
    bcftools_view_retry region_ref.vcf.gz \
      ${REF_REGION:+-r "${REF_REGION}"} \
      "${SITE_FILTER_ARGS[@]}" \
      "${BIALLELIC_ARGS[@]}" \
      -Oz --threads ~{cpu} \
      "${VIEW_SRC}"
    bcftools index -t region_ref.vcf.gz
    n_post=$(bcftools index -n region_ref.vcf.gz 2>/dev/null || echo 0)
    if [[ "${site_filters_active}" != "true" && -z "${REF_REGION}" ]]; then
      n_pre="${n_post}"
    elif [[ "${site_filters_active}" != "true" && -n "${REF_REGION}" ]]; then
      n_pre="${n_post}"
    fi
    echo "ref sites pre=${n_pre} post=${n_post}" >&2
    echo -e "ref\t${n_pre}\t${n_post}" > ref_counts.tsv
    cat ref_counts.tsv >&2

    if [[ "${n_post}" == "0" ]]; then
      echo "error: no sites in ref after filters" >&2
      exit 1
    fi

    python3 - "${REF_REGION}" ref.header ~{min_sites_per_mb} ref_counts.tsv <<'PY'
import re, sys
from pathlib import Path

region = sys.argv[1].strip()
header_path = Path(sys.argv[2])
min_density = float(sys.argv[3])
counts = Path(sys.argv[4]).read_text().strip().split()
if len(counts) < 3:
    raise SystemExit("bad ref_counts.tsv")
n_post = int(counts[2])

def span_mb(reg, header_text):
    reg = reg.strip()
    if not reg:
        return None
    chrom, rest = (reg.split(":", 1) + [""])[:2]
    if rest and "-" in rest:
        start_s, end_s = rest.split("-", 1)
        return (int(end_s) - int(start_s) + 1) / 1e6
    for line in header_text.splitlines():
        if not line.startswith("##contig=<"):
            continue
        m_id = re.search(r"ID=([^,>]+)", line)
        if not m_id or m_id.group(1) != chrom:
            continue
        m_len = re.search(r"length=([0-9]+)", line)
        if m_len:
            return int(m_len.group(1)) / 1e6
    return None

header_text = header_path.read_text(errors="replace") if header_path.is_file() else ""
span = span_mb(region, header_text)
if span is None or span <= 0:
    print("WARNING: skip ref density check (empty region or unknown contig span)", file=sys.stderr)
    raise SystemExit(0)
density = n_post / span
print(f"ref density={density:.1f} sites/Mb span={span:.3f}Mb min={min_density}", file=sys.stderr)
if density < min_density:
    raise SystemExit(
        f"ref post-filter density {density:.1f} sites/Mb < min {min_density} "
        f"({n_post} sites / {span:.3f} Mb)"
    )
PY

    ls -lh region_ref.vcf.gz* ref_counts.tsv >&2
  >>>

  output {
    File vcf = "region_ref.vcf.gz"
    File index = "region_ref.vcf.gz.tbi"
    File counts_tsv = "ref_counts.tsv"
  }

  runtime {
    docker: docker
    cpu: cpu
    memory: memory_gb + " GB"
    disks: "local-disk " + disk_gb_floor + " HDD"
    preemptible: preemptible
    maxRetries: max_retries
  }
}

task SubsetOnePopulation {
  input {
    File gt_vcf
    File gt_vcf_index
    File keep_list
    String region
    Boolean biallelic_snvs_only = false
    File? include_sites
    File exclude_regions
    String site_filter_digest
    Float min_sites_per_mb = 1000.0
    Float target_min_maf = 0.0
    Float target_hwe_pval = 0.0
    String gcs_project
    String docker
    Int cpu
    Int memory_gb
    Int disk_gb_floor
    Float disk_gb_multiplier
    Int stream_attempts
    Int preemptible
    Int max_retries
  }

  Boolean has_include = defined(include_sites)
  File include_or_dummy = select_first([include_sites, exclude_regions])

  # flare_task reads the population from the gt VCF basename (AFR.vcf.gz → AFR).
  String population = basename(keep_list, ".samples.txt")
  Int disk_gb = ceil(size(gt_vcf, "GB") * disk_gb_multiplier) + disk_gb_floor

  parameter_meta {
    gt_vcf: {
      description: "Phased target VCF streamed from GCS on Terra PAPI.",
      localization_optional: true
    }
    gt_vcf_index: {
      description: "Tabix/CSI sibling of gt_vcf.",
      localization_optional: true
    }
  }

  command <<<
    set -euo pipefail
    export HTS_RETRY_MAX="${HTS_RETRY_MAX:-8}"
    export HTS_RETRY_DELAY="${HTS_RETRY_DELAY:-500}"
    export HTS_RETRY_MAX_DELAY="${HTS_RETRY_MAX_DELAY:-60000}"
    echo "$(date -u +%Y-%m-%dT%H:%M:%SZ) $(bcftools --version 2>/dev/null | head -n 1) HTS_RETRY_MAX=${HTS_RETRY_MAX}" >&2
    cat > gcs_auth.py <<'PY'
import json, os, subprocess, sys, urllib.request
from pathlib import Path

def fresh_token():
    for cmd in (
        ["gcloud", "auth", "print-access-token"],
        ["gcloud", "auth", "application-default", "print-access-token"],
    ):
        try:
            tok = subprocess.check_output(cmd, text=True, stderr=subprocess.DEVNULL).strip()
            if tok:
                return tok
        except Exception:
            pass
    for url in (
        "http://metadata.google.internal/computeMetadata/v1/instance/service-accounts/default/token",
        "http://169.254.169.254/computeMetadata/v1/instance/service-accounts/default/token",
    ):
        req = urllib.request.Request(url, headers={"Metadata-Flavor": "Google"})
        try:
            with urllib.request.urlopen(req, timeout=5) as resp:
                tok = json.load(resp).get("access_token") or ""
                if tok:
                    return tok
        except Exception:
            continue
    return ""

need_gcs = any(u.startswith("gs://") for u in sys.argv[1:])
token = fresh_token() if need_gcs else ""
if need_gcs and not token:
    raise SystemExit("error: bcftools needs a GCS OAuth token to read gs:// VCFs")
project = "~{gcs_project}".strip() or os.environ.get("GOOGLE_PROJECT", "")
parts = []
if project:
    parts.append(f"export GCS_REQUESTER_PAYS_PROJECT={project}")
if token:
    parts.append(f"export GCS_OAUTH_TOKEN='{token}'")
Path("gcs.env").write_text(("\n".join(parts) + "\n") if parts else "")
PY
    SRC="~{gt_vcf}"
    IDX="~{gt_vcf_index}"
    KEEP="~{keep_list}"
    pop="~{population}"
    refresh_gcs_auth() {
      python3 gcs_auth.py "${SRC}"
      if [[ -s gcs.env ]]; then
        # shellcheck disable=SC1091
        source gcs.env
      fi
      echo "$(date -u +%Y-%m-%dT%H:%M:%SZ) refreshed GCS auth" >&2
    }
    bcftools_view_retry() {
      local out="$1"
      shift
      local attempt=1
      local max=~{stream_attempts}
      local delay=20
      if [[ "${max}" -lt 1 ]]; then
        max=1
      fi
      while [[ "${attempt}" -le "${max}" ]]; do
        refresh_gcs_auth
        echo "$(date -u +%Y-%m-%dT%H:%M:%SZ) bcftools view attempt ${attempt}/${max} -> ${out}" >&2
        rm -f "${out}"
        if bcftools view -o "${out}" "$@"; then
          return 0
        fi
        echo "$(date -u +%Y-%m-%dT%H:%M:%SZ) bcftools view failed (attempt ${attempt}/${max})" >&2
        if [[ "${attempt}" -eq "${max}" ]]; then
          return 1
        fi
        sleep "${delay}"
        delay=$((delay * 2))
        attempt=$((attempt + 1))
      done
    }
    refresh_gcs_auth

    python3 - <<'PY'
import hashlib
from pathlib import Path

biallelic = "~{biallelic_snvs_only}".lower() == "true"
has_include = "~{has_include}" == "true"
include_path = Path("~{include_or_dummy}")
exclude_path = Path("~{exclude_regions}")
if has_include:
    include_hash = hashlib.sha256(include_path.read_bytes()).hexdigest()
else:
    include_hash = hashlib.sha256(b"NONE").hexdigest()
exclude_hash = hashlib.sha256(exclude_path.read_bytes()).hexdigest()
digest_input = f"{str(biallelic).lower()}:{include_hash}:{exclude_hash}"
computed = hashlib.sha256(digest_input.encode()).hexdigest()
expected = "~{site_filter_digest}".strip()
if computed != expected:
    raise SystemExit(f"site filter digest mismatch: {computed} != {expected}")
print("site filter digest ok", computed, flush=True)
PY

    python3 -c "import pathlib,sys; pathlib.Path('region.txt').write_text(sys.argv[1].strip().strip(chr(34)).strip(chr(39)) + chr(10))" "~{region}"
    REGION="$(cat region.txt)"

    SITE_FILTER_ARGS=()
    if [[ "~{has_include}" == "true" ]]; then
      SITE_FILTER_ARGS+=(-T "~{include_or_dummy}")
    fi
    if [[ -s "~{exclude_regions}" ]]; then
      SITE_FILTER_ARGS+=(-T "^~{exclude_regions}")
    fi
    BIALLELIC_ARGS=()
    if [[ "~{biallelic_snvs_only}" == "true" ]]; then
      BIALLELIC_ARGS=(-m2 -M2 -v snps)
    fi
    site_filters_active=false
    if [[ "~{has_include}" == "true" ]] || [[ -s "~{exclude_regions}" ]] || [[ "~{biallelic_snvs_only}" == "true" ]]; then
      site_filters_active=true
    fi

    bcftools view -h "${SRC}##idx##${IDX}" > gt.header
    GT_REGION=""
    if [[ -n "${REGION}" ]]; then
      python3 - "${REGION}" gt.header gt.region <<'PY'
import re, sys
from pathlib import Path

region, header_path, out_path = sys.argv[1], sys.argv[2], sys.argv[3]
header = Path(header_path).read_text(errors="replace")
contigs = []
for line in header.splitlines():
    if line.startswith("##contig=<"):
        m = re.search(r"ID=([^,>]+)", line)
        if m:
            contigs.append(m.group(1))
print(f"{Path(header_path).name} contigs: {', '.join(contigs[:20])}{' …' if len(contigs) > 20 else ''}", file=sys.stderr)
if not contigs:
    raise SystemExit(f"{header_path}: no ##contig=<ID=...> in header")
chrom, rest = (region.split(":", 1) + [""])[:2]
aliases = [chrom]
if chrom.startswith("chr"):
    aliases.append(chrom[3:])
else:
    aliases.append("chr" + chrom)
matched = next((a for a in aliases if a in set(contigs)), None)
if matched is None:
    raise SystemExit(
        f"region chrom {chrom!r} is not in this VCF. Contigs: {contigs}. "
        "Use the phased gt for the SAME chromosome as region "
        "(chr22 window needs aou_lr_phase2_v1.chr22.vcf.gz, not chr1)."
    )
resolved = matched if not rest else f"{matched}:{rest}"
if resolved != region:
    print(f"rewrote region {region} -> {resolved}", file=sys.stderr)
Path(out_path).write_text(resolved + "\n")
print(resolved)
PY
      GT_REGION="$(cat gt.region)"
    fi

    n_pre=0
    if [[ -n "${GT_REGION}" ]]; then
      echo "Count ${pop} region ${GT_REGION} pre site filters" >&2
      bcftools_view_retry region_only.vcf.gz -S "${KEEP}" -r "${GT_REGION}" -Oz --threads ~{cpu} \
        "${SRC}##idx##${IDX}"
      bcftools index -t region_only.vcf.gz
      n_pre=$(bcftools index -n region_only.vcf.gz 2>/dev/null || echo 0)
    elif [[ "${site_filters_active}" == "true" ]]; then
      n_pre=$(bcftools index -n "${SRC}##idx##${IDX}" 2>/dev/null || echo 0)
      if [[ "${n_pre}" == "0" ]]; then
        n_pre=$(bcftools view -H -S "${KEEP}" --threads 1 "${SRC}##idx##${IDX}" | wc -l | tr -d ' ')
      fi
    fi

    echo "Subset ${pop} region=${GT_REGION:-FULL} samples=$(wc -l < "${KEEP}")" >&2
    bcftools_view_retry "${pop}.prefilter.vcf.gz" \
      -S "${KEEP}" \
      ${GT_REGION:+-r "${GT_REGION}"} \
      "${SITE_FILTER_ARGS[@]}" \
      "${BIALLELIC_ARGS[@]}" \
      -Oz --threads ~{cpu} \
      "${SRC}##idx##${IDX}"
    bcftools index -t "${pop}.prefilter.vcf.gz"
    n_post=$(bcftools index -n "${pop}.prefilter.vcf.gz" 2>/dev/null || echo 0)
    if [[ "${site_filters_active}" != "true" ]]; then
      n_pre="${n_post}"
    fi
    n_s=$(bcftools query -l "${pop}.prefilter.vcf.gz" | wc -l | tr -d ' ')

    echo "Call-rate check on ${pop} (phased complete GT assumed; no call-rate filter applied)" >&2
    bcftools +fill-tags "${pop}.prefilter.vcf.gz" -Oz -o "${pop}.cr_check.vcf.gz" -- -t F_MISSING
    bcftools index -t "${pop}.cr_check.vcf.gz"
    n_missing_sites=$(bcftools view -H -i 'F_MISSING>0' "${pop}.cr_check.vcf.gz" | wc -l | tr -d ' ')
    if [[ "${n_missing_sites}" -gt 0 ]]; then
      echo "WARNING: ${pop} has ${n_missing_sites} sites with F_MISSING>0 (call-rate filter skipped)" >&2
    else
      echo "${pop}: F_MISSING==0 on all sites; call-rate filter would be a no-op" >&2
    fi

    n_removed_maf=0
    n_removed_hwe=0
    CURRENT="${pop}.prefilter.vcf.gz"
    if awk "BEGIN{exit !(~{target_min_maf} > 0)}"; then
      bcftools +fill-tags "${CURRENT}" -Oz -o "${pop}.maf_tag.vcf.gz" -- -t MAF
      bcftools index -t "${pop}.maf_tag.vcf.gz"
      n_before=$(bcftools index -n "${pop}.maf_tag.vcf.gz" 2>/dev/null || echo 0)
      bcftools view -e 'MAF<~{target_min_maf}' "${pop}.maf_tag.vcf.gz" -Oz -o "${pop}.maf_filt.vcf.gz"
      bcftools index -t "${pop}.maf_filt.vcf.gz"
      n_after=$(bcftools index -n "${pop}.maf_filt.vcf.gz" 2>/dev/null || echo 0)
      n_removed_maf=$((n_before - n_after))
      CURRENT="${pop}.maf_filt.vcf.gz"
      echo "${pop} target_min_maf=~{target_min_maf} removed=${n_removed_maf}" >&2
    fi
    if awk "BEGIN{exit !(~{target_hwe_pval} > 0)}"; then
      bcftools +fill-tags "${CURRENT}" -Oz -o "${pop}.hwe_tag.vcf.gz" -- -t HWE
      bcftools index -t "${pop}.hwe_tag.vcf.gz"
      n_before=$(bcftools index -n "${pop}.hwe_tag.vcf.gz" 2>/dev/null || echo 0)
      bcftools view -e 'HWE<~{target_hwe_pval}' "${pop}.hwe_tag.vcf.gz" -Oz -o "${pop}.hwe_filt.vcf.gz"
      bcftools index -t "${pop}.hwe_filt.vcf.gz"
      n_after=$(bcftools index -n "${pop}.hwe_filt.vcf.gz" 2>/dev/null || echo 0)
      n_removed_hwe=$((n_before - n_after))
      CURRENT="${pop}.hwe_filt.vcf.gz"
      echo "${pop} target_hwe_pval=~{target_hwe_pval} removed=${n_removed_hwe}" >&2
    fi

    cp -f "${CURRENT}" "${pop}.vcf.gz"
    if [[ -f "${CURRENT}.tbi" ]]; then
      cp -f "${CURRENT}.tbi" "${pop}.vcf.gz.tbi"
    else
      bcftools index -t "${pop}.vcf.gz"
    fi
    n_final=$(bcftools index -n "${pop}.vcf.gz" 2>/dev/null || echo 0)

    # Part 4.0 variant-class breakdown on the post-filter target shard.
    n_bial_snv=$(bcftools view -H -v snps -m2 -M2 "${pop}.vcf.gz" | wc -l | tr -d ' ')
    n_indels=$(bcftools view -H -v indels "${pop}.vcf.gz" | wc -l | tr -d ' ')
    n_bial_any=$(bcftools view -H -m2 -M2 "${pop}.vcf.gz" | wc -l | tr -d ' ')
    n_multi=$((n_final - n_bial_any))
    if [[ "${n_multi}" -lt 0 ]]; then n_multi=0; fi
    frac_bial="NA"
    if [[ "${n_final}" -gt 0 ]]; then
      frac_bial=$(awk -v a="${n_bial_snv}" -v b="${n_final}" 'BEGIN{printf "%.6g", a/b}')
    fi
    echo -e "population\tn_sites_input\tn_biallelic_snv\tn_indels\tn_multiallelic\tfrac_biallelic_snv" > "${pop}.site_classes.tsv"
    echo -e "${pop}\t${n_final}\t${n_bial_snv}\t${n_indels}\t${n_multi}\t${frac_bial}" >> "${pop}.site_classes.tsv"
    cat "${pop}.site_classes.tsv" >&2

    echo -e "${pop}\t${n_s}\t${n_pre}\t${n_final}\t${n_removed_maf}\t${n_removed_hwe}" > counts.tsv
    cat counts.tsv >&2

    if [[ "${n_final}" == "0" ]]; then
      echo "error: no sites in ${pop} subset after site filters" >&2
      exit 1
    fi

    python3 - "${GT_REGION}" gt.header ~{min_sites_per_mb} counts.tsv <<'PY'
import re, sys
from pathlib import Path

region = sys.argv[1].strip()
header_path = Path(sys.argv[2])
min_density = float(sys.argv[3])
parts = Path(sys.argv[4]).read_text().strip().split("\t")
if len(parts) < 4:
    raise SystemExit("bad counts.tsv")
n_post = int(parts[3])

def span_mb(reg, header_text):
    reg = reg.strip()
    if not reg:
        return None
    chrom, rest = (reg.split(":", 1) + [""])[:2]
    if rest and "-" in rest:
        start_s, end_s = rest.split("-", 1)
        return (int(end_s) - int(start_s) + 1) / 1e6
    for line in header_text.splitlines():
        if not line.startswith("##contig=<"):
            continue
        m_id = re.search(r"ID=([^,>]+)", line)
        if not m_id or m_id.group(1) != chrom:
            continue
        m_len = re.search(r"length=([0-9]+)", line)
        if m_len:
            return int(m_len.group(1)) / 1e6
    return None

header_text = header_path.read_text(errors="replace") if header_path.is_file() else ""
span = span_mb(region, header_text)
if span is None or span <= 0:
    print("WARNING: skip subset density check (empty region or unknown contig span)", file=sys.stderr)
    raise SystemExit(0)
density = n_post / span
print(f"{parts[0]} density={density:.1f} sites/Mb span={span:.3f}Mb min={min_density}", file=sys.stderr)
if density < min_density:
    raise SystemExit(
        f"{parts[0]} post-filter density {density:.1f} sites/Mb < min {min_density} "
        f"({n_post} sites / {span:.3f} Mb)"
    )
PY

    ls -lh "${pop}.vcf.gz"* counts.tsv >&2
  >>>

  output {
    File subset_vcf = "~{population}.vcf.gz"
    File subset_index = "~{population}.vcf.gz.tbi"
    File counts_tsv = "counts.tsv"
    File site_classes = "~{population}.site_classes.tsv"
  }

  runtime {
    docker: docker
    cpu: cpu
    memory: memory_gb + " GB"
    disks: "local-disk " + disk_gb + " HDD"
    preemptible: preemptible
    maxRetries: max_retries
  }
}

task ConcatSubsetCounts {
  input {
    Array[File] tables
    String docker
    Int preemptible
  }

  command <<<
    set -euo pipefail
    echo -e "population\tn_samples\tn_sites_pre_filter\tn_sites_post_filter\tn_removed_target_maf\tn_removed_target_hwe" > subset_counts.tsv
    if [[ ~{length(tables)} -gt 0 ]]; then
      for f in ~{sep=" " tables}; do
        tail -n +2 "${f}" >> subset_counts.tsv
      done
    fi
    cat subset_counts.tsv >&2
  >>>

  output {
    File subset_counts = "subset_counts.tsv"
  }

  runtime {
    docker: docker
    cpu: 1
    memory: "2 GB"
    disks: "local-disk 10 HDD"
    preemptible: preemptible
  }
}

task ValidateFlareInputs {
  input {
    Array[File] keep_lists
    File ref_panel
    String gen_by_pop
    Boolean gen_by_pop_allow_default
    Boolean allow_unrepresented_pops
    String docker
    Int preemptible
  }

  command <<<
    set -euo pipefail
    python3 - <<'PY'
from pathlib import Path
import sys

keep_lists = """~{sep="\n" keep_lists}""".strip().splitlines()
pops = []
for path in keep_lists:
    path = path.strip()
    if not path:
        continue
    # keep list basename is <POP>.keep.txt or similar; prefer first line of
    # sibling naming: files from SplitSamples are named {POP}.samples.txt
    name = Path(path).name
    pop = name.split(".")[0]
    pops.append(pop)

panel_pops = set()
for line in Path("~{ref_panel}").read_text().splitlines():
    line = line.strip()
    if not line or line.startswith("#"):
        continue
    bits = line.split()
    if len(bits) >= 2:
        panel_pops.add(bits[1])

# Panel labels are often lowercase (eas/afr); shard pops are uppercase.
panel_upper = {p.upper() for p in panel_pops}
missing = sorted({p for p in pops if p.upper() not in panel_upper})
allow = "~{allow_unrepresented_pops}".lower() == "true"
if missing:
    msg = (
        f"shard pops with no same-named panel ancestry: {missing}; "
        f"panel labels={sorted(panel_pops)}"
    )
    if allow:
        print("WARNING:", msg, file=sys.stderr)
    else:
        raise SystemExit(msg)

gen_raw = """~{gen_by_pop}""".strip()
allow_default = "~{gen_by_pop_allow_default}".lower() == "true"
if gen_raw:
    overrides = {}
    for part in gen_raw.split(","):
        part = part.strip()
        if not part:
            continue
        if ":" not in part:
            raise SystemExit(f"gen_by_pop entry {part!r} must look like POP:T")
        key, val = part.split(":", 1)
        overrides[key.strip()] = float(val.strip())
    missing_gen = sorted(p for p in pops if p not in overrides)
    if missing_gen and not allow_default:
        raise SystemExit(
            f"gen_by_pop missing shard pops {missing_gen}; "
            f"set gen_by_pop_allow_default=true to fall back to gen"
        )
    if missing_gen:
        print(
            f"WARNING: gen_by_pop missing {missing_gen}; using workflow gen",
            file=sys.stderr,
        )

Path("ok.txt").write_text("ok\n")
Path("validation.log").write_text(
    f"pops={pops}\npanel={sorted(panel_pops)}\n"
    f"allow_unrepresented_pops={allow}\ngen_by_pop={gen_raw!r}\n"
)
print("ValidateFlareInputs ok", pops, file=sys.stderr)
PY
  >>>

  output {
    String ok = read_string("ok.txt")
    File log = "validation.log"
  }

  runtime {
    docker: docker
    cpu: 1
    memory: "2 GB"
    disks: "local-disk 10 HDD"
    preemptible: preemptible
  }
}

task flare_task {
  input {
    File ref_vcf
    File ref_panel
    File gt_vcf
    File map_file
    String output_prefix

    Boolean em = true
    Array[File] pop_models = []
    File? template_model
    String props_by_pop = ""
    Boolean inherit_props = true
    Boolean inherit_panel_weights = true
    Boolean inherit_mu = true
    File flare_model_script
    File flare_site_stats_script
    Boolean update_p = true
    Boolean probs = false
    Boolean array = false
    Float min_maf = 0.005
    Int min_mac = 50
    Float gen = 10.0
    String gen_by_pop = ""
    Int seed = 12345
    String validation_ok = "ok"

    Int? cpu_override
    String? memory_override
    Int? xmx_gb_override
    Int? disk_size_gb_override
    String disk_type = "HDD"
    Int preemptible = 0
    String docker_image = "us-docker.pkg.dev/broad-dsde-methods/popout/flare:latest"
  }

  Float ref_gb = size(ref_vcf, "GB")
  Float gt_gb = size(gt_vcf, "GB")

  # Compressed gt GB underestimates full-chrom HMM + ANP (probs) peak heap.
  # Stack OOM was in EstimatedGlobalAncProportions during fwd/bwd with probs.
  Int bucket_cpu = if probs then (
                     if gt_gb < 1.0 then 8
                     else if gt_gb > 100.0 then 32
                     else if gt_gb > 50.0 then 24
                     else if gt_gb > 20.0 then 16
                     else 12
                   ) else (
                     if gt_gb < 1.0 then 4
                     else if gt_gb > 100.0 then 32
                     else if gt_gb > 50.0 then 24
                     else if gt_gb > 20.0 then 16
                     else 8
                   )
  Int bucket_mem_gb = if probs then (
                        if gt_gb < 1.0 then 96
                        else if gt_gb > 100.0 then 256
                        else if gt_gb > 50.0 then 192
                        else if gt_gb > 20.0 then 128
                        else 96
                      ) else (
                        if gt_gb < 1.0 then 16
                        else if gt_gb > 100.0 then 192
                        else if gt_gb > 50.0 then 128
                        else if gt_gb > 20.0 then 64
                        else 48
                      )
  Int bucket_xmx_gb = if probs then (
                        if gt_gb < 1.0 then 72
                        else if gt_gb > 100.0 then 192
                        else if gt_gb > 50.0 then 144
                        else if gt_gb > 20.0 then 96
                        else 72
                      ) else (
                        if gt_gb < 1.0 then 12
                        else if gt_gb > 100.0 then 144
                        else if gt_gb > 50.0 then 96
                        else if gt_gb > 20.0 then 48
                        else 36
                      )

  Float heap_probs_mult = if probs then 3.0 else 1.0
  Float predicted_heap_gb = gt_gb * 2.2 * heap_probs_mult
  Int sized_xmx_gb = ceil(predicted_heap_gb * 1.3)
  Int sized_mem_gb = ceil(sized_xmx_gb * 1.34)
  Int sized_cpu = ceil(sized_mem_gb * 1.0 / 8.0)

  Int auto_cpu = if sized_cpu > bucket_cpu then sized_cpu else bucket_cpu
  Int auto_xmx_gb = if sized_xmx_gb > bucket_xmx_gb then sized_xmx_gb else bucket_xmx_gb
  Int auto_mem_int = if sized_mem_gb > bucket_mem_gb then sized_mem_gb else bucket_mem_gb
  String auto_memory = "~{auto_mem_int} GB"

  Float output_multiplier = if probs then 4.0 else 2.5
  Int auto_disk = ceil(ref_gb + gt_gb * output_multiplier) + (if gt_gb < 1.0 then 20 else 50)

  Int cpu = select_first([cpu_override, auto_cpu])
  String memory = select_first([memory_override, auto_memory])
  Int xmx_gb = select_first([xmx_gb_override, auto_xmx_gb])
  Int disk_size_gb = select_first([disk_size_gb_override, auto_disk])

  Boolean has_template = defined(template_model)
  # select_first so optional File is always a concrete path in the command block.
  File template_model_or_dummy = select_first([template_model, flare_model_script])

  command <<<
    set -euo pipefail
    echo "validation_ok=~{validation_ok}" >&2
    POP="$(basename "~{gt_vcf}" .vcf.gz)"
    PREFIX="~{output_prefix}.${POP}"
    echo "FLARE population=${POP} prefix=${PREFIX}" >&2
    # Fixed seed=12345 ⇒ identical outputs show determinism, not independent
    # corroboration across replicate submissions.
    NUM_REF_SAMPLES=$(awk 'NF && $1 !~ /^#/' ~{ref_panel} | wc -l | tr -d ' ')
    NUM_REF_POPS=$(awk 'NF && $1 !~ /^#/ {print $2}' ~{ref_panel} | sort -u | wc -l | tr -d ' ')
    echo "ref panel: ${NUM_REF_SAMPLES} samples across ${NUM_REF_POPS} populations" >&2

    mkdir -p in_model
    IN_MODEL="in_model/${PREFIX}.model"
    SENTINEL="in_model/NONE"
    cat > model_paths.txt <<'END_MODELS'
~{sep="\n" pop_models}
END_MODELS

    POP_MODEL="$(
      python3 - "${POP}" model_paths.txt "~{flare_model_script}" <<'PY'
import sys
from pathlib import Path
sys.path.insert(0, str(Path(sys.argv[3]).resolve().parent))
from flare_model import match_pop_model

pop = sys.argv[1]
paths = [ln.strip() for ln in Path(sys.argv[2]).read_text().splitlines() if ln.strip()]
try:
    got = match_pop_model(pop, paths, require=bool(paths))
except ValueError as exc:
    raise SystemExit(exc)
if got is not None:
    print(got)
PY
    )" || exit 1

    POP_MODEL_ARGS=()
    if [[ -n "${POP_MODEL}" ]]; then
      POP_MODEL_ARGS=(--pop-model "${POP_MODEL}")
    fi
    TEMPLATE_ARGS=()
    if [[ "~{has_template}" == "true" ]]; then
      TEMPLATE_ARGS=(--template-model "~{template_model_or_dummy}")
    fi

    python3 "~{flare_model_script}" build-input \
      --pop "${POP}" \
      --gen "~{gen}" \
      --gen-by-pop "~{gen_by_pop}" \
      --props-by-pop "~{props_by_pop}" \
      ~{true="--inherit-props" false="--no-inherit-props" inherit_props} \
      ~{true="--inherit-panel-weights" false="--no-inherit-panel-weights" inherit_panel_weights} \
      ~{true="--inherit-mu" false="--no-inherit-mu" inherit_mu} \
      "${POP_MODEL_ARGS[@]}" \
      "${TEMPLATE_ARGS[@]}" \
      --out "${IN_MODEL}"

    EM_FLAG="~{em}"
    MODEL=""
    if [[ -s "${IN_MODEL}" ]]; then
      MODEL="$(realpath "${IN_MODEL}")"
      # Matching / rewritten model always pins parameters.
      EM_FLAG="false"
    else
      # Empty build-input output ⇒ omit model=; honor workflow em.
      printf 'NONE\n' > "${SENTINEL}"
      IN_MODEL="${SENTINEL}"
    fi
    GEN="$(
      python3 - "${POP}" "~{gen}" "~{gen_by_pop}" <<'PY'
import sys
pop, gen_s, raw = sys.argv[1], sys.argv[2], sys.argv[3].strip()
gen = float(gen_s)
if raw:
    for part in raw.split(","):
        part = part.strip()
        if not part:
            continue
        key, val = part.split(":", 1)
        if key.strip() == pop:
            gen = float(val.strip())
            break
print(f"{gen:g}")
PY
    )"

    MODEL_ARG=()
    if [[ -n "${MODEL}" ]]; then
      MODEL_ARG=("model=${MODEL}")
    fi
    echo "population=${POP} em=${EM_FLAG} gen=${GEN} model=${MODEL:-'(none)'} in_model=${IN_MODEL}" >&2
    if [[ -s "${IN_MODEL}" && "${IN_MODEL}" != "${SENTINEL}" ]]; then
      echo "----- input model -----" >&2
      cat "${IN_MODEL}" >&2
      echo "----- end input model -----" >&2
    fi

    echo "flare_task resources: cpu=~{cpu} memory=~{memory} xmx=~{xmx_gb}g gt_gb=~{gt_gb} probs=~{probs}" >&2
    java -Xmx~{xmx_gb}g -jar /opt/flare/flare.jar \
      ref=~{ref_vcf} \
      ref-panel=~{ref_panel} \
      gt=~{gt_vcf} \
      map=~{map_file} \
      out="${PREFIX}" \
      em="${EM_FLAG}" \
      "${MODEL_ARG[@]}" \
      update-p=~{update_p} \
      probs=~{probs} \
      array=~{array} \
      nthreads=~{cpu} \
      min-maf=~{min_maf} \
      min-mac=~{min_mac} \
      gen="${GEN}" \
      seed=~{seed}

    # Glob only FLARE outputs under cwd — never in_model/.
    python3 - "${PREFIX}" <<'PY'
import sys
from pathlib import Path
prefix = sys.argv[1]
cwd = Path(".")
cands = []
for p in cwd.glob(f"{prefix}.*.model"):
    cands.append(p)
exact = Path(f"{prefix}.model")
if exact.is_file():
    cands.append(exact)
# Unique, exclude in_model/
uniq = []
seen = set()
for p in cands:
    if "in_model" in p.parts:
        continue
    key = str(p.resolve())
    if key in seen:
        continue
    seen.add(key)
    uniq.append(p)
if len(uniq) != 1:
    raise SystemExit(
        f"expected exactly one output model for {prefix}; got {[str(u) for u in uniq]}"
    )
out = Path(f"{prefix}.out.model")
out.write_bytes(uniq[0].read_bytes())
print(uniq[0], "->", out, file=sys.stderr)
PY
    cp -L "${IN_MODEL}" "${PREFIX}.in.model"
    ls -lh "${PREFIX}".* in_model/ || true
    printf '%s\n' "${POP}" > pop.txt

    python3 "~{flare_site_stats_script}" summarize \
      --pop "${POP}" \
      --log "${PREFIX}.log" \
      --out "${PREFIX}.site_stats.tsv"
    # VCF class breakdown runs in SubsetOnePopulation (bcftools docker). FLARE
    # image may lack bcftools; retained marker count comes from the log here.
    ls -lh "${PREFIX}.site_stats.tsv" >&2
  >>>

  output {
    String population = read_string("pop.txt")
    File anc_vcf = glob("*.anc.vcf.gz")[0]
    File global_anc = glob("*.global.anc.gz")[0]
    File out_model = glob("*.out.model")[0]
    File in_model = glob("*.in.model")[0]
    File log = glob("*.log")[0]
    File site_stats = glob("*.site_stats.tsv")[0]
  }

  runtime {
    docker: docker_image
    cpu: cpu
    memory: memory
    disks: "local-disk ~{disk_size_gb} ~{disk_type}"
    preemptible: preemptible
  }
}

task MergeSiteStats {
  input {
    Array[File] tables
    Array[File] site_classes
    File site_stats_script
    String output_prefix
    String docker
    Int preemptible
  }

  command <<<
    set -euo pipefail
    python3 "~{site_stats_script}" merge \
      --tables ~{sep=" " tables} \
      --out "~{output_prefix}.site_stats.log_only.tsv"
    python3 - <<'PY'
import csv
from pathlib import Path

log_rows = {}
with open("~{output_prefix}.site_stats.log_only.tsv") as fh:
    for row in csv.DictReader(fh, delimiter="\t"):
        log_rows[row["population"]] = row

class_paths = """~{sep="\n" site_classes}""".strip().splitlines()
class_rows = {}
for path in class_paths:
    path = path.strip()
    if not path:
        continue
    with open(path) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            class_rows[row["population"]] = row

fields = [
    "population",
    "n_sites_input",
    "n_sites_flare_retained",
    "frac_biallelic_snv",
    "n_biallelic_snv",
    "n_indels",
    "n_multiallelic",
]
pops = sorted(set(log_rows) | set(class_rows))
out = Path("~{output_prefix}.site_stats.tsv")
with out.open("w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
    w.writeheader()
    for pop in pops:
        lr = log_rows.get(pop, {})
        cr = class_rows.get(pop, {})
        w.writerow(
            {
                "population": pop,
                "n_sites_input": cr.get("n_sites_input") or lr.get("n_sites_input") or "",
                "n_sites_flare_retained": lr.get("n_sites_flare_retained") or "",
                "frac_biallelic_snv": cr.get("frac_biallelic_snv") or lr.get("frac_biallelic_snv") or "",
                "n_biallelic_snv": cr.get("n_biallelic_snv") or "",
                "n_indels": cr.get("n_indels") or "",
                "n_multiallelic": cr.get("n_multiallelic") or "",
            }
        )
print("wrote", out)
PY
    cat "~{output_prefix}.site_stats.tsv" >&2
  >>>

  output {
    File site_stats = "~{output_prefix}.site_stats.tsv"
  }

  runtime {
    docker: docker
    cpu: 1
    memory: "2 GB"
    disks: "local-disk 10 HDD"
    preemptible: preemptible
  }
}

task MergePopulationFlare {
  input {
    Array[File] anc_vcfs
    Array[File] global_ancs
    Array[File] models
    Array[String] pops
    File summarize_script
    # Co-localize next to summarize_script so ``import flare_model`` works
    # (Cromwell places each File input in its own directory).
    File flare_model_script
    String output_prefix
    String docker
    Int cpu
    Int preemptible
  }

  Int n_vcf = length(anc_vcfs)
  Int disk_gb = 80 + ceil(size(anc_vcfs, "GB") * 3.0)

  command <<<
    set -euo pipefail
    mkdir -p _scripts
    cp -L "~{summarize_script}" _scripts/flare_summarize_models.py
    cp -L "~{flare_model_script}" _scripts/flare_model.py
    python3 _scripts/flare_summarize_models.py \
      --models ~{sep=" " models} \
      --pops ~{sep=" " pops} \
      --global-anc ~{sep=" " global_ancs} \
      --out-models-tsv "~{output_prefix}.models.tsv" \
      --out-global-anc "~{output_prefix}.global.anc.gz"
    cat "~{output_prefix}.models.tsv" >&2

    mkdir -p indexed
    i=0
    for vcf in ~{sep=" " anc_vcfs}; do
      i=$((i + 1))
      out="indexed/p${i}.anc.vcf.gz"
      cp -L "${vcf}" "${out}"
      bcftools index -t "${out}"
    done

    python3 - <<'PY'
from pathlib import Path
import subprocess
paths = sorted(Path("indexed").glob("*.anc.vcf.gz"))
counts = []
for p in paths:
    n = subprocess.check_output(["bcftools", "index", "-n", str(p)], text=True).strip()
    counts.append((p.name, n))
    print(p.name, n, flush=True)
uniq = {n for _, n in counts}
if len(uniq) != 1:
    raise SystemExit(f"FLARE marker counts differ across populations: {counts}")
PY

    shopt -s nullglob
    VCFS=(indexed/*.anc.vcf.gz)
    if [[ ${#VCFS[@]} -eq 0 ]]; then
      echo "no indexed ancestry VCFs" >&2
      exit 1
    fi
    if [[ "~{n_vcf}" -eq 1 ]]; then
      cp -L "${VCFS[0]}" "~{output_prefix}.anc.vcf.gz"
      if [[ -f "${VCFS[0]}.tbi" ]]; then
        cp -L "${VCFS[0]}.tbi" "~{output_prefix}.anc.vcf.gz.tbi"
      else
        bcftools index -t "~{output_prefix}.anc.vcf.gz"
      fi
    else
      bcftools merge -m none --threads ~{cpu} -Oz -o "~{output_prefix}.anc.vcf.gz" "${VCFS[@]}"
      bcftools index -t "~{output_prefix}.anc.vcf.gz"
    fi
    echo -n "merged samples: "
    bcftools query -l "~{output_prefix}.anc.vcf.gz" | wc -l
    echo -n "merged sites: "
    bcftools index -n "~{output_prefix}.anc.vcf.gz"
  >>>

  output {
    File anc_vcf = "~{output_prefix}.anc.vcf.gz"
    File anc_vcf_index = "~{output_prefix}.anc.vcf.gz.tbi"
    File global_anc = "~{output_prefix}.global.anc.gz"
    File models_tsv = "~{output_prefix}.models.tsv"
  }

  runtime {
    docker: docker
    cpu: cpu
    memory: "16 GB"
    disks: "local-disk " + disk_gb + " HDD"
    preemptible: preemptible
  }
}
