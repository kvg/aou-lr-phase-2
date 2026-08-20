version 1.0

# Self-contained Terra Methods workflow for AoU-LR SV annotation + manuscript tables.
# The companion Docker image contains all project scripts, bcftools, and CADD-SV.

workflow AnnotateSvCallset {
  input {
    String phase  # phase1 | phase2
    File main_vcf          # VCF, VCF.gz, or BCF
    File? bnd_vcf
    File? large_vcf        # large / ultralong partition

    File rmsk_bed
    File simple_repeat_bed
    File segdup_bed

    File sample_ancestry_tsv
    File caddsv_annotations_tar

    String prefix
    String docker
  }

  call ProcessPartition as MainPart {
    input:
      vcf = main_vcf,
      source_vcf = "main",
      phase = phase,
      prefix = prefix + ".main",
      rmsk_bed = rmsk_bed,
      simple_repeat_bed = simple_repeat_bed,
      segdup_bed = segdup_bed,
      docker = docker,
      cpu = 4,
      memory_gb = 64,
      disk_gb = 500,
      preemptible = 0
  }

  call CaddSvAndAttach as MainCadd {
    input:
      sites = MainPart.sites,
      caddsv_bed = MainPart.caddsv_bed,
      caddsv_annotations_tar = caddsv_annotations_tar,
      prefix = prefix + ".main",
      docker = docker
  }

  if (defined(bnd_vcf)) {
    call ProcessPartition as BndPart {
      input:
        vcf = select_first([bnd_vcf]),
        source_vcf = "bnd",
        phase = phase,
        prefix = prefix + ".bnd",
        rmsk_bed = rmsk_bed,
        simple_repeat_bed = simple_repeat_bed,
        segdup_bed = segdup_bed,
        docker = docker,
        cpu = 4,
        memory_gb = 16,
        disk_gb = 200,
        preemptible = 0
    }
  }

  if (defined(large_vcf)) {
    call ProcessPartition as LargePart {
      input:
        vcf = select_first([large_vcf]),
        source_vcf = "large",
        phase = phase,
        prefix = prefix + ".large",
        rmsk_bed = rmsk_bed,
        simple_repeat_bed = simple_repeat_bed,
        segdup_bed = segdup_bed,
        docker = docker,
        cpu = 4,
        memory_gb = 32,
        disk_gb = 400,
        preemptible = 0
    }

    # CADD-SV is not run on the large/ultralong partition: ~225k intervals,
    # each >10 kb, make the 236-job annotation DAG too heavy (OOM / timeout).
    # Discovery CADD strata use main only; large sites stay cadd_sv_bin=unscored.
  }

  call IntegrateAndSummarize {
    input:
      main_sites = MainCadd.scored_sites,
      main_carriers = MainPart.carriers,
      bnd_sites = BndPart.sites,
      large_sites = LargePart.sites,
      sample_ancestry_tsv = sample_ancestry_tsv,
      phase = phase,
      prefix = prefix,
      docker = docker,
      memory_gb = 64,
      disk_gb = 400
  }

  output {
    File unified_sites = IntegrateAndSummarize.unified_sites
    File integrate_manifest = IntegrateAndSummarize.integrate_manifest
    File manuscript_counts_tsv = IntegrateAndSummarize.manuscript_counts_tsv
    File manuscript_counts_json = IntegrateAndSummarize.manuscript_counts_json
    File discovery_tsv = IntegrateAndSummarize.discovery_tsv
    File discovery_region_tsv = IntegrateAndSummarize.discovery_region_tsv
    File discovery_cadd_tsv = IntegrateAndSummarize.discovery_cadd_tsv
    File main_caddsv_bed = MainPart.caddsv_bed
    File main_carriers = MainPart.carriers
  }
}

task ProcessPartition {
  input {
    File vcf
    String source_vcf
    String phase
    String prefix
    File rmsk_bed
    File simple_repeat_bed
    File segdup_bed
    String docker
    Int cpu = 4
    Int memory_gb = 32
    Int disk_gb = 400
    Int preemptible = 0
    Int max_retries = 1
  }

  command <<<
    set -euo pipefail
    echo "[$(date -Is)] ProcessPartition ~{source_vcf} start" >&2
    python3 /opt/aou_sv/scripts/extract_sites_from_vcf.py \
      --vcf "~{vcf}" \
      --source-vcf ~{source_vcf} \
      --phase ~{phase} \
      --out ~{prefix}.raw.sites.tsv
    echo "[$(date -Is)] extract_sites done" >&2

    python3 /opt/aou_sv/scripts/fill_af_and_carriers.py \
      --vcf "~{vcf}" \
      --sites ~{prefix}.raw.sites.tsv \
      --out-sites ~{prefix}.af.sites.tsv \
      --out-carriers ~{prefix}.carriers.tsv
    echo "[$(date -Is)] fill_af done" >&2

    python3 /opt/aou_sv/scripts/annotate_repeat_context.py \
      --sites ~{prefix}.af.sites.tsv \
      --rmsk-bed "~{rmsk_bed}" \
      --simple-repeat-bed "~{simple_repeat_bed}" \
      --segdup-bed "~{segdup_bed}" \
      --out ~{prefix}.region.sites.tsv
    echo "[$(date -Is)] annotate_repeat done" >&2

    python3 /opt/aou_sv/scripts/prepare_caddsv_input.py \
      --sites ~{prefix}.region.sites.tsv \
      --out-bed ~{prefix}.caddsv.bed

    cp ~{prefix}.region.sites.tsv ~{prefix}.sites.tsv
    echo "[$(date -Is)] ProcessPartition ~{source_vcf} done" >&2
  >>>

  output {
    File sites = "~{prefix}.sites.tsv"
    File carriers = "~{prefix}.carriers.tsv"
    File caddsv_bed = "~{prefix}.caddsv.bed"
  }

  runtime {
    docker: docker
    cpu: cpu
    memory: memory_gb + " GiB"
    disks: "local-disk " + disk_gb + " SSD"
    preemptible: preemptible
    maxRetries: max_retries
  }
}

task CaddSvAndAttach {
  input {
    File sites
    File caddsv_bed
    File caddsv_annotations_tar
    String prefix
    String docker
    Int cpu = 4
    Int memory_gb = 32
    Int disk_gb = 300
    Int preemptible = 0
  }

  command <<<
    set -euo pipefail

    # CADD-SV retrieve_seq_for_nt.py pysam-fetches the full DUP/INV span
    # (then truncates). Ultralong intervals OOM; DELs only use 96 bp flanks.
    awk -F '\t' '($4 != "DUP" && $4 != "INV") || (($3 - $2) <= 1000000)' \
      "~{caddsv_bed}" > caddsv.capped.bed
    echo "[CaddSvAndAttach] input lines=$(wc -l < "~{caddsv_bed}") capped=$(wc -l < caddsv.capped.bed)" >&2

    python3 /opt/aou_sv/scripts/run_caddsv.py \
      --input-bed caddsv.capped.bed \
      --annotations-tar "~{caddsv_annotations_tar}" \
      --out-scores ~{prefix}.caddsv.scores.tsv \
      --threads ~{cpu} \
      --conda-prefix /caddsv-conda

    python3 /opt/aou_sv/scripts/attach_caddsv_scores.py \
      --sites "~{sites}" \
      --scores ~{prefix}.caddsv.scores.tsv \
      --out ~{prefix}.caddsv.sites.tsv
  >>>

  output {
    File scored_sites = "~{prefix}.caddsv.sites.tsv"
    File scores = "~{prefix}.caddsv.scores.tsv"
  }

  runtime {
    docker: docker
    cpu: cpu
    memory: memory_gb + " GiB"
    disks: "local-disk " + disk_gb + " SSD"
    preemptible: preemptible
  }
}

task IntegrateAndSummarize {
  input {
    File main_sites
    File main_carriers
    File? bnd_sites
    File? large_sites
    File sample_ancestry_tsv
    String phase
    String prefix
    String docker
    Int cpu = 2
    Int memory_gb = 64
    Int disk_gb = 400
    Int preemptible = 0
  }

  command <<<
    set -euo pipefail
    BND_ARG=()
    LARGE_ARG=()
    if [[ -n "~{bnd_sites}" ]]; then BND_ARG=(--bnd "~{bnd_sites}"); fi
    if [[ -n "~{large_sites}" ]]; then LARGE_ARG=(--large "~{large_sites}"); fi

    python3 /opt/aou_sv/scripts/integrate_partitions.py \
      --main "~{main_sites}" \
      "${BND_ARG[@]}" \
      "${LARGE_ARG[@]}" \
      --out ~{prefix}.unified.sites.tsv \
      --manifest ~{prefix}.integrate.manifest.json

    if [[ "~{phase}" == "phase2" ]]; then
      python3 /opt/aou_sv/scripts/manuscript_site_counts.py \
        --phase2-sites ~{prefix}.unified.sites.tsv \
        --out-tsv ~{prefix}.manuscript_counts.tsv \
        --out-json ~{prefix}.manuscript_counts.json
    else
      python3 /opt/aou_sv/scripts/manuscript_site_counts.py \
        --phase1-sites ~{prefix}.unified.sites.tsv \
        --out-tsv ~{prefix}.manuscript_counts.tsv \
        --out-json ~{prefix}.manuscript_counts.json
    fi

    # Discovery: stream main scored sites + carriers in lockstep (same order).
    # Do not use the unified table here — it reorders/drops main∩large rows and
    # would break the 1:1 join with main.carriers.tsv.
    python3 /opt/aou_sv/scripts/ebert_discovery.py \
      --sites "~{main_sites}" \
      --carriers "~{main_carriers}" \
      --sample-ancestry "~{sample_ancestry_tsv}" \
      --strata-list none,region,cadd \
      --include-sources main \
      --out-prefix ~{prefix}
  >>>

  output {
    File unified_sites = "~{prefix}.unified.sites.tsv"
    File integrate_manifest = "~{prefix}.integrate.manifest.json"
    File manuscript_counts_tsv = "~{prefix}.manuscript_counts.tsv"
    File manuscript_counts_json = "~{prefix}.manuscript_counts.json"
    File discovery_tsv = "~{prefix}.discovery.tsv"
    File discovery_region_tsv = "~{prefix}.discovery.region.tsv"
    File discovery_cadd_tsv = "~{prefix}.discovery.cadd.tsv"
  }

  runtime {
    docker: docker
    cpu: cpu
    memory: memory_gb + " GiB"
    disks: "local-disk " + disk_gb + " HDD"
    preemptible: preemptible
  }
}
