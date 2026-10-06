version 1.0

# Scan the v3 SV companions once (genome-wide, as in Table 2) and merge them
# with the per-chromosome GLnexus VariantLengthSpectrum outputs.
#
#   main_vcf  = gs://…/v3_main.bcf       resolved DEL / INS with |SVLEN| >= sv_min_bp
#   large_vcf = gs://…/v3_ultralong.bcf  every record
#   bnd_vcf   = gs://…/v3_bnd.bcf        count only (no length)

workflow MergeVariantLengthSpectrum {
  input {
    Array[File] chrom_bins
    Array[File] chrom_summaries

    File main_vcf
    File? large_vcf
    File? bnd_vcf

    File rmsk_bed
    File simple_repeat_bed
    File segdup_bed

    File summarize_py
    File sv_site_utils_py
    File annotate_repeat_context_py
    File merge_py

    Int sv_min_bp = 20

    String docker = "us-central1-docker.pkg.dev/broad-dsp-lrma/aou-lr/aou-sv-annotation:0.1.6"
    Int companion_cpu = 2
    Int companion_memory_gb = 16
    Int companion_disk_gb_floor = 50
    Int preemptible = 1
  }

  call SummarizeCompanions {
    input:
      main_vcf = main_vcf,
      large_vcf = large_vcf,
      bnd_vcf = bnd_vcf,
      rmsk_bed = rmsk_bed,
      simple_repeat_bed = simple_repeat_bed,
      segdup_bed = segdup_bed,
      summarize_py = summarize_py,
      sv_site_utils_py = sv_site_utils_py,
      annotate_repeat_context_py = annotate_repeat_context_py,
      sv_min_bp = sv_min_bp,
      docker = docker,
      cpu = companion_cpu,
      memory_gb = companion_memory_gb,
      disk_gb_floor = companion_disk_gb_floor,
      preemptible = preemptible
  }

  call Merge {
    input:
      bin_files = flatten([chrom_bins, [SummarizeCompanions.bins]]),
      summary_files = flatten([chrom_summaries, [SummarizeCompanions.summary]]),
      merge_py = merge_py,
      docker = docker,
      preemptible = preemptible
  }

  output {
    File companion_bins = SummarizeCompanions.bins
    File companion_summary = SummarizeCompanions.summary
    File bins = Merge.bins
    File summary = Merge.summary
  }

  meta {
    description: "Figure 2 length spectrum: v3 SV companions + per-chromosome GLnexus counts, summed."
    allowNestedInputs: true
  }
}

task SummarizeCompanions {
  input {
    File main_vcf
    File? large_vcf
    File? bnd_vcf
    File rmsk_bed
    File simple_repeat_bed
    File segdup_bed
    File summarize_py
    File sv_site_utils_py
    File annotate_repeat_context_py
    Int sv_min_bp
    String docker
    Int cpu
    Int memory_gb
    Int disk_gb_floor
    Int preemptible
  }

  Int large_disk = if defined(large_vcf) then ceil(size(select_first([large_vcf]), "GB") * 1.2) else 0
  Int bnd_disk = if defined(bnd_vcf) then ceil(size(select_first([bnd_vcf]), "GB") * 1.2) else 0
  Int disk_gb = ceil(size(main_vcf, "GB") * 1.2) + large_disk + bnd_disk + disk_gb_floor

  command <<<
    set -euo pipefail
    mkdir -p "$PWD/py"
    cp "~{summarize_py}" "$PWD/py/summarize_variant_length_spectrum.py"
    cp "~{sv_site_utils_py}" "$PWD/py/sv_site_utils.py"
    cp "~{annotate_repeat_context_py}" "$PWD/py/annotate_repeat_context.py"
    export PYTHONPATH="$PWD/py:${PYTHONPATH:-}"

    EXTRA=()
    if [[ -n "~{large_vcf}" ]]; then
      EXTRA+=(--large-vcf "~{large_vcf}")
    fi
    if [[ -n "~{bnd_vcf}" ]]; then
      EXTRA+=(--bnd-vcf "~{bnd_vcf}")
    fi

    python3 "$PWD/py/summarize_variant_length_spectrum.py" \
      --main-vcf "~{main_vcf}" \
      --rmsk-bed "~{rmsk_bed}" \
      --simple-repeat-bed "~{simple_repeat_bed}" \
      --segdup-bed "~{segdup_bed}" \
      --sv-min-bp ~{sv_min_bp} \
      --out-bins companions.length_spectrum.bins.tsv \
      --out-summary companions.length_spectrum.summary.json \
      "${EXTRA[@]}"
  >>>

  output {
    File bins = "companions.length_spectrum.bins.tsv"
    File summary = "companions.length_spectrum.summary.json"
  }

  runtime {
    docker: docker
    cpu: cpu
    memory: memory_gb + " GiB"
    disks: "local-disk " + disk_gb + " HDD"
    preemptible: preemptible
  }
}

task Merge {
  input {
    Array[File] bin_files
    Array[File] summary_files
    File merge_py
    String docker
    Int preemptible
  }

  command <<<
    set -euo pipefail
    MERGE_PY="~{merge_py}"
    echo "[merge] using ${MERGE_PY}" >&2
    if ! python3 "${MERGE_PY}" --help 2>&1 | grep -q -- '--bins'; then
      echo "[merge] ERROR: merge_py does not accept --bins." >&2
      echo "[merge] Point MergeVariantLengthSpectrum.merge_py at" >&2
      echo "[merge]   scripts/merge_variant_length_spectrum.py" >&2
      echo "[merge] not merge_manuscript_counts.py (phase1/phase2 TSV merge)." >&2
      python3 "${MERGE_PY}" --help >&2 || true
      exit 2
    fi
    python3 "${MERGE_PY}" \
      --bins ~{sep=' ' bin_files} \
      --summaries ~{sep=' ' summary_files} \
      --out-bins fig2_length_spectrum.bins.tsv \
      --out-summary fig2_length_spectrum.summary.json
  >>>

  output {
    File bins = "fig2_length_spectrum.bins.tsv"
    File summary = "fig2_length_spectrum.summary.json"
  }

  runtime {
    docker: docker
    cpu: 1
    memory: "4 GiB"
    disks: "local-disk 20 HDD"
    preemptible: preemptible
  }
}
