version 1.0

# Genome-wide FELIX for Terra Methods Repository (self-contained; no imports).
#
# Shape:
#   CheckChromVcfPairs (length-only; do not pass Array[File] to the check task)
#   BuildSaigePlinkAndSparseGRM once
#   FitFelixNull once per phenotype (fit_felix_null.R → null_export for RU scorer)
#   PackFelixla once per chromosome (joint_vcf OR phase_vcf + flare_vcf arrays)
#   RunFelixStep2 once per (phenotype × chromosome)
#   Concat → one merged TSV per phenotype
#   SummarizeFelixResults → phenotype-named copies, λGC / QQ / Manhattan
#
# Optional repeat-dosage branch:
#   ru_engine="felix" (default, needs repeat_catalog_bed):
#     AggregateRepeatLoci → WriteRuAdmixedVcf → RunFelixStep2RuVcf
#     One test per catalog repeat locus, REF-relative units, SPA + carrier QC.
#     ru_biallelic_arm (default true) runs the same loci a second time with
#     coding="carrier" (one unit per non-reference haplotype), and CompareRuArms
#     builds the 2x2 — pooled/local-ancestry x biallelic/repeat-dosage — from
#     the two merged TSVs. The pooled column is FELIX's own ancestry-collapsed
#     pass (p.value_ancALL = P_hom_admixed), not a separate run, so all four
#     arms share one Step 1 null, GRM, sample set and locus set.
#   ru_engine="rust" (needs simple_repeat_bed), comparison only:
#     AnnotateRuTest → ExtractRuDosage (dosage / collapse / split) → ScoreRuDosage
#   (tractor-mix-score --mode felix) per (phenotype × chromosome × encoding)
#
# Terra chrom-set launch: parallel Arrays
#   chroms      <- explicit chr names
#   flare_vcfs  <- FLARE anc VCFs (PackFelixla)
#   phase_vcfs  <- joint phased VCFs (PackFelixla + RU annotator)
# or joint_vcfs filling both phase and FLARE when GT:AN1:AN2 is in one file.

task CheckChromVcfPairs {
  input {
    Int n_chroms
    Int n_joint_vcfs
    Int n_phase_vcfs
    Int n_flare_vcfs
    String docker = "us-central1-docker.pkg.dev/broad-dsp-lrma/aou-lr/felix-pilot:0.2.0"
  }

  command <<<
    set -euo pipefail
    N_CHROM=~{n_chroms}
    N_JOINT=~{n_joint_vcfs}
    N_PHASE=~{n_phase_vcfs}
    N_FLARE=~{n_flare_vcfs}

    if [[ "${N_CHROM}" -lt 1 ]]; then
      echo "Need at least one chromosome" >&2
      exit 1
    fi

    # Pack reads phase_vcfs + flare_vcfs when given, else joint_vcfs; the repeat
    # branch reads joint_vcfs when given, else phase_vcfs. Any given array must
    # match chroms, and phase_vcfs and flare_vcfs come as a pair.
    if [[ "${N_JOINT}" -gt 0 && "${N_JOINT}" -ne "${N_CHROM}" ]]; then
      echo "joint_vcfs length (${N_JOINT}) != chroms length (${N_CHROM})" >&2
      exit 1
    fi
    if [[ "${N_PHASE}" -gt 0 || "${N_FLARE}" -gt 0 ]]; then
      if [[ "${N_PHASE}" -ne "${N_CHROM}" || "${N_FLARE}" -ne "${N_CHROM}" ]]; then
        echo "phase_vcfs (${N_PHASE}) and flare_vcfs (${N_FLARE}) must both match chroms (${N_CHROM})" >&2
        exit 1
      fi
    elif [[ "${N_JOINT}" -eq 0 ]]; then
      echo "Supply joint_vcfs, or both phase_vcfs and flare_vcfs" >&2
      exit 1
    fi
    echo "OK: ${N_CHROM} chrom; joint=${N_JOINT} phase=${N_PHASE} flare=${N_FLARE}"
  >>>

  output {
    Int n_chroms_out = n_chroms
  }

  runtime {
    docker: docker
    cpu: 1
    memory: "1 GB"
    disks: "local-disk 10 HDD"
    preemptible: 3
  }
}

task BuildSaigePlinkAndSparseGRM {
  input {
    Array[File] grm_vcfs
    File analysis_samples
    File build_script
    File make_plink_keep_script
    Float relatedness_cutoff = 0.05
    Int num_random_markers = 2000
    String docker = "us-central1-docker.pkg.dev/broad-dsp-lrma/aou-lr/felix-pilot:0.2.0"
    Int cpu = 8
    Int memory_gb = 32
    Int disk_gb = 300
    Int preemptible = 0
  }

  command <<<
    set -euo pipefail
    chmod +x "~{build_script}"
    export MAKE_PLINK_KEEP_PY="~{make_plink_keep_script}"
    "~{build_script}" \
      --analysis-samples "~{analysis_samples}" \
      --relatedness-cutoff ~{relatedness_cutoff} \
      --n-threads ~{cpu} \
      --num-random-markers ~{num_random_markers} \
      --out-prefix saige_grm/plink \
      -- \
      ~{sep=" " grm_vcfs}

    cp -L saige_grm/sparseGRM.mtx saige_sparseGRM.mtx
    cp -L saige_grm/sparseGRM.sampleIDs.txt saige_sparseGRM.sampleIDs.txt
    cp saige_grm/plink.bed saige_plink.bed
    cp saige_grm/plink.bim saige_plink.bim
    cp saige_grm/plink.fam saige_plink.fam
    cp saige_grm/saige_sparse_grm.paths.tsv .
  >>>

  output {
    File plink_bed = "saige_plink.bed"
    File plink_bim = "saige_plink.bim"
    File plink_fam = "saige_plink.fam"
    File sparse_grm_mtx = "saige_sparseGRM.mtx"
    File sparse_grm_sample_ids = "saige_sparseGRM.sampleIDs.txt"
    File paths_tsv = "saige_sparse_grm.paths.tsv"
  }

  runtime {
    docker: docker
    cpu: cpu
    memory: memory_gb + " GB"
    disks: "local-disk " + disk_gb + " HDD"
    preemptible: preemptible
  }
}

# Guard for precomputed_null_*: the arrays are indexed against
# read_lines(selected_phenotypes), so a length mismatch would silently pair a
# phenotype with another phenotype's null model. Counts only — passing the
# arrays themselves would localize every null file just to count them.
task CheckNullArrays {
  input {
    Int n_phenotypes
    Int n_null_rdas
    Int n_variance_ratios
    Int n_samples_used
    Boolean precomputed
    String docker = "us-central1-docker.pkg.dev/broad-dsp-lrma/aou-lr/felix-pilot:0.2.0"
  }

  command <<<
    set -euo pipefail
    N=~{n_phenotypes}
    for pair in "null_rdas:~{n_null_rdas}" "variance_ratios:~{n_variance_ratios}" \
                "samples_used:~{n_samples_used}"; do
      name="${pair%%:*}"; count="${pair##*:}"
      if [[ "${count}" -ne "${N}" ]]; then
        echo "precomputed ${name}: ${count} entries for ${N} phenotypes" >&2
        exit 1
      fi
    done
    if [[ "~{precomputed}" == "true" ]]; then
      echo "Using ${N} precomputed null fits; FitFelixNull skipped."
    else
      echo "Fitting ${N} null models in this run."
    fi
    printf '%s\n' "${N}" > n_phenotypes.txt
  >>>

  output {
    Int n_phenotypes_out = read_int("n_phenotypes.txt")
  }

  runtime {
    docker: docker
    cpu: 1
    memory: "2 GB"
    disks: "local-disk 10 HDD"
    preemptible: 3
  }
}

task FitFelixNull {
  input {
    File pheno_cov
    String phenotype
    File covariate_columns
    File analysis_samples
    File plink_bed
    File plink_bim
    File plink_fam
    File sparse_grm_mtx
    File sparse_grm_sample_ids
    File fit_null_script
    # Optional override for export_felix_null.R. Unset, fit_felix_null.R runs the copy
    # baked into the image (/opt/felix_scripts), because Cromwell localizes only declared
    # inputs. Set it to a staged copy to pick up script fixes without an image rebuild.
    File? export_null_script
    String trait_type = "binary"
    String docker = "us-central1-docker.pkg.dev/broad-dsp-lrma/aou-lr/felix-pilot:0.2.0"
    Int cpu = 8
    Int memory_gb = 32
    Int disk_gb = 100
    Int preemptible = 0
  }

  command <<<
    set -euo pipefail
    mkdir -p null
    ln -sf "~{plink_bed}" null/plink.bed
    ln -sf "~{plink_bim}" null/plink.bim
    ln -sf "~{plink_fam}" null/plink.fam

    Rscript "~{fit_null_script}" \
      --pheno-cov "~{pheno_cov}" \
      --phenotype "~{phenotype}" \
      --covariates "~{covariate_columns}" \
      --analysis-samples "~{analysis_samples}" \
      --plink-prefix null/plink \
      --sparse-grm "~{sparse_grm_mtx}" \
      --sparse-grm-ids "~{sparse_grm_sample_ids}" \
      --trait-type "~{trait_type}" \
      --step1-r /usr/local/bin/step1_fitNULLGLMM.R \
      ~{"--export-r " + export_null_script} \
      --n-threads ~{cpu} \
      --out-prefix "null/~{phenotype}"

    cp "null/~{phenotype}.rda" "~{phenotype}.null.rda" || \
      cp $(ls -1 null/~{phenotype}*.rda | head -n1) "~{phenotype}.null.rda"
    cp "null/~{phenotype}.varianceRatio.txt" "~{phenotype}.varianceRatio.txt" || \
      cp $(ls -1 null/~{phenotype}*.varianceRatio.txt | head -n1) "~{phenotype}.varianceRatio.txt"
    cp "null/~{phenotype}.pheno.tsv" "~{phenotype}.pheno.tsv"
    cp "null/~{phenotype}.samples.txt" "~{phenotype}.samples.txt"
    cp "null/~{phenotype}.null_meta.tsv" "~{phenotype}.null_meta.tsv"

    NULL_EXPORT="null/~{phenotype}.null_export"
    if [[ ! -d "${NULL_EXPORT}" ]]; then
      echo "missing ${NULL_EXPORT} (fit_felix_null.R should export CSC layout)" >&2
      exit 1
    fi
    tar czf null_export.tar.gz -C null "$(basename "${NULL_EXPORT}")"
  >>>

  output {
    File null_rda = "~{phenotype}.null.rda"
    File variance_ratio = "~{phenotype}.varianceRatio.txt"
    File pheno_used = "~{phenotype}.pheno.tsv"
    File samples_used = "~{phenotype}.samples.txt"
    File null_meta = "~{phenotype}.null_meta.tsv"
    File null_export_tar = "null_export.tar.gz"
  }

  runtime {
    docker: docker
    cpu: cpu
    memory: memory_gb + " GB"
    disks: "local-disk " + disk_gb + " HDD"
    preemptible: preemptible
  }
}

task PackFelixla {
  input {
    File phase_vcf
    File flare_vcf
    File analysis_samples
    String chrom
    Int num_ancs = 5
    String docker = "us-central1-docker.pkg.dev/broad-dsp-lrma/aou-lr/felix-pilot:0.2.0"
    Int cpu = 1
    Int memory_gb = 6
    Int disk_gb_floor = 100
    Float disk_gb_multiplier = 3.5
    Int preemptible = 0
  }

  Int disk_gb = ceil((size(phase_vcf, "GB") + size(flare_vcf, "GB")) * disk_gb_multiplier) + disk_gb_floor

  command <<<
    set -euo pipefail
    mkdir -p pack_vcf packed
    CHROM=$(python3 -c 'import sys; print(sys.argv[1].strip().strip(chr(34)).strip(chr(39)))' "~{chrom}")
    PREFIX="packed/${CHROM}"

    fix_vcf() {
      local in_vcf="$1"
      local out_vcf="$2"
      local hdr
      hdr=$(mktemp)
      bcftools view -h "${in_vcf}" \
        | awk '/^##FORMAT=<ID=GT,/{ if (gt++) next } { print }' \
        > "${hdr}"
      bcftools reheader -h "${hdr}" -o "${out_vcf}" "${in_vcf}"
      bcftools index -t "${out_vcf}" || bcftools index -c "${out_vcf}"
      rm -f "${hdr}"
    }

    if [[ "~{phase_vcf}" == "~{flare_vcf}" ]]; then
      fix_vcf "~{phase_vcf}" pack_vcf/joint.vcf.gz
      PHASE_IN=pack_vcf/joint.vcf.gz
      FLARE_IN=pack_vcf/joint.vcf.gz
    else
      fix_vcf "~{phase_vcf}" pack_vcf/phase.vcf.gz
      fix_vcf "~{flare_vcf}" pack_vcf/flare.vcf.gz
      PHASE_IN=pack_vcf/phase.vcf.gz
      FLARE_IN=pack_vcf/flare.vcf.gz
    fi

    felixla \
      --phase-vcf "${PHASE_IN}" \
      --flare-vcf "${FLARE_IN}" \
      --n-ancestries ~{num_ancs} \
      --keep "~{analysis_samples}" \
      --make-felixla \
      --out "${PREFIX}"

    if [[ ! -f "${PREFIX}.meta" || ! -f "${PREFIX}.samples" ]]; then
      echo "felixla did not write ${PREFIX}.meta / .samples" >&2
      ls -lh packed/ || true
      exit 1
    fi

    printf '%s\n' "${CHROM}" > packed/chrom.txt
    tar czf felixla_packed.tar.gz -C packed .
  >>>

  output {
    File packed_tar = "felixla_packed.tar.gz"
    String chrom_id = read_string("packed/chrom.txt")
  }

  runtime {
    docker: docker
    cpu: cpu
    memory: memory_gb + " GB"
    disks: "local-disk " + disk_gb + " HDD"
    preemptible: preemptible
  }
}

task RunFelixStep2 {
  input {
    File packed_tar
    String chrom
    String phenotype
    File null_rda
    File variance_ratio
    File samples_used
    File sparse_grm_mtx
    File sparse_grm_sample_ids
    File run_step2_script
    Int num_ancs = 5
    Int min_mac = 50
    Float pvalcutoff_of_haplotype = 0.05
    String docker = "us-central1-docker.pkg.dev/broad-dsp-lrma/aou-lr/felix-pilot:0.2.0"
    Int cpu = 1
    Int memory_gb = 6
    Int disk_gb_floor = 50
    Float disk_gb_multiplier = 4.0
    Int preemptible = 1
  }

  Int disk_gb = ceil(size(packed_tar, "GB") * disk_gb_multiplier) + disk_gb_floor

  command <<<
    set -euo pipefail
    mkdir -p packed step2
    CHROM=$(python3 -c 'import sys; print(sys.argv[1].strip().strip(chr(34)).strip(chr(39)))' "~{chrom}")
    PHENO=$(python3 -c 'import sys; print(sys.argv[1].strip().strip(chr(34)).strip(chr(39)))' "~{phenotype}")
    tar xzf "~{packed_tar}" -C packed
    PREFIX="packed/${CHROM}"

    ln -sf "~{null_rda}" "step2/${PHENO}.rda"
    ln -sf "~{variance_ratio}" "step2/${PHENO}.varianceRatio.txt"

    Rscript "~{run_step2_script}" \
      --felixla-prefix "${PREFIX}" \
      --chrom "${CHROM}" \
      --null-prefix "step2/${PHENO}" \
      --sample-file "~{samples_used}" \
      --sparse-grm "~{sparse_grm_mtx}" \
      --sparse-grm-ids "~{sparse_grm_sample_ids}" \
      --min-mac ~{min_mac} \
      --n-ancestries ~{num_ancs} \
      --pvalcutoff-of-haplotype ~{pvalcutoff_of_haplotype} \
      --n-threads 1 \
      --out-tsv "${PHENO}.felix.tsv" \
      --out-raw "${PHENO}.felix.raw.txt"

    cp "${PHENO}.felix.tsv" results.felix.tsv
    printf '%s\n' "${PHENO}.felix.tsv" > results.output_name.txt
  >>>

  output {
    File results_tsv = "results.felix.tsv"
    File results_name = "results.output_name.txt"
  }

  runtime {
    docker: docker
    cpu: cpu
    memory: memory_gb + " GB"
    disks: "local-disk " + disk_gb + " HDD"
    preemptible: preemptible
  }
}

task ConcatPhenotypeScores {
  input {
    String phenotype
    Array[File] shard_tsvs
    String docker = "us-central1-docker.pkg.dev/broad-dsp-lrma/aou-lr/felix-pilot:0.2.0"
    Int cpu = 1
    Int memory_gb = 4
    Int disk_gb_floor = 20
    # Localized shards + the merged output + sort's spill directory, so three
    # copies of the data live on disk at once, not two.
    Float disk_gb_multiplier = 4.0
    Int preemptible = 1
  }

  Int disk_gb = ceil(size(shard_tsvs, "GB") * disk_gb_multiplier) + disk_gb_floor

  command <<<
    set -euo pipefail
    PHENO=$(python3 -c 'import sys; print(sys.argv[1].strip().strip(chr(34)).strip(chr(39)))' "~{phenotype}")
    OUT="${PHENO}.felix.tsv.gz"
    export PHENO OUT
    cat > shards.txt <<'EOF'
~{sep="\n" shard_tsvs}
EOF

    # Header: take the first shard's and fail if any other differs.
    HDR=""
    NSHARD=0
    while IFS= read -r f; do
      [[ -n "${f}" ]] || continue
      NSHARD=$((NSHARD + 1))
      h=$(head -n 1 "${f}")
      if [[ -z "${h}" ]]; then
        echo "empty shard: ${f}" >&2
        exit 1
      fi
      if [[ -z "${HDR}" ]]; then
        HDR="${h}"
      elif [[ "${h}" != "${HDR}" ]]; then
        echo "header mismatch in ${f}" >&2
        echo "  got: ${h}" >&2
        echo "  exp: ${HDR}" >&2
        exit 1
      fi
    done < shards.txt
    if [[ -z "${HDR}" ]]; then
      echo "no shard TSVs to concatenate" >&2
      exit 1
    fi

    # Body: stream every shard through an EXTERNAL sort that spills to disk.
    #
    # This used to be a python block that concatenated to a temp file, then did
    # tmp.read_text().splitlines() and sorted() on the result -- i.e. it held
    # the whole phenotype's genome-wide results in memory as python strings.
    # On the per-record arm that is ~27 GB per phenotype against a 4 GB
    # request, so the OOM killer took it after an hour of useful work ("Killed"
    # with no traceback, then "Required file output merged.felix.tsv does not
    # exist"). coreutils sort is bounded by -S and spills to -T instead.
    #
    # Key matches the old python one: chromosome lexicographically, position
    # numerically. LC_ALL=C keeps the chromosome ordering byte-wise, which is
    # what python's str comparison did.
    # Written GZIPPED and written ONCE. Uncompressed, the per-record arm is
    # ~25.65 GB per phenotype (measured on BI_164) and gzip -6 gives 3.43x on a
    # real result table, so this is ~24 TB -> ~7 TB across 951 phenotypes. It
    # also makes the disk-bound Summarize ~3.4x faster, because that task's cost
    # is bytes moved, not parsing.
    #
    # The `cp "${OUT}" merged.felix.tsv` that used to end this task kept a
    # SECOND full copy in the same call directory -- Cromwell retains the whole
    # execution directory, so both were billed. The phenotype-named form is
    # recorded in merged.output_name.txt, which is all any consumer needed it
    # for; nothing has to carry the bytes twice to preserve the name.
    #
    # Body rows are counted by awk inside the pipeline rather than by wc -l on
    # the result, which would mean decompressing the whole output again.
    SORT_MEM_MB=$(( ~{memory_gb} * 1024 / 2 ))
    TAB=$(printf '\t')
    mkdir -p sorttmp
    {
      printf '%s\n' "${HDR}"
      while IFS= read -r f; do
        [[ -n "${f}" ]] || continue
        tail -n +2 "${f}"
      done < shards.txt \
        | awk 'NF{n++; print} END{printf "%d\n", n+0 > "body_rows.txt"}' \
        | LC_ALL=C sort -t "${TAB}" -k1,1 -k2,2n -S "${SORT_MEM_MB}M" -T sorttmp
    } | gzip -6 > merged.felix.tsv.gz
    rm -rf sorttmp

    NROWS=$(cat body_rows.txt)
    echo "Wrote merged.felix.tsv.gz as ${OUT} (${NROWS} variant rows from ${NSHARD} shards for ${PHENO})"
    printf '%s\n' "${OUT}" > merged.output_name.txt
  >>>

  output {
    File merged_tsv = "merged.felix.tsv.gz"
    File merged_name = "merged.output_name.txt"
  }

  runtime {
    docker: docker
    cpu: cpu
    memory: memory_gb + " GB"
    disks: "local-disk " + disk_gb + " HDD"
    preemptible: preemptible
    maxRetries: 2
  }
}

task SummarizeFelixResults {
  input {
    Array[File] result_tsvs
    Array[String] phenotypes
    File pheno_cov
    File summarize_script
    # Names this arm's outputs so the main and repeat summaries can sit side by
    # side in one bundle: "main", "ru", "ru_biallelic".
    String arm = "main"
    Float p_threshold = 0.00000005
    Int top_n = 50
    String docker = "us-central1-docker.pkg.dev/broad-dsp-lrma/aou-lr/felix-pilot:0.2.0"
    Int cpu = 2
    Int memory_gb = 16
    Int disk_gb_floor = 50
    Float disk_gb_multiplier = 2.0
    Int preemptible = 1
  }

  Int disk_gb = ceil(size(result_tsvs, "GB") * disk_gb_multiplier) + disk_gb_floor

  command <<<
    set -euo pipefail
    mkdir -p summary
    cat > phenotypes.txt <<'EOF'
~{sep="\n" phenotypes}
EOF
    cat > results.txt <<'EOF'
~{sep="\n" result_tsvs}
EOF

    mapfile -t RESULT_ARR < results.txt
    mapfile -t PHENO_ARR < phenotypes.txt
    if [[ "${#RESULT_ARR[@]}" -ne "${#PHENO_ARR[@]}" ]]; then
      echo "results (${#RESULT_ARR[@]}) != phenotypes (${#PHENO_ARR[@]})" >&2
      exit 1
    fi

    python3 "~{summarize_script}" \
      --results "${RESULT_ARR[@]}" \
      --phenotypes "${PHENO_ARR[@]}" \
      --pheno-cov "~{pheno_cov}" \
      --out-dir summary \
      --arm "~{arm}" \
      --top-n ~{top_n} \
      --p-threshold ~{p_threshold}

    ls -lh summary/results_by_phenotype/ || true
    wc -l summary/results_manifest.tsv summary/calibration_summary.tsv

    # Small review bundle: manifests, calibration table, QQ/Manhattan and top
    # hits, but NOT results_by_phenotype/ (full per-phenotype scans, ~400 MB
    # each). Those stay in the bucket for the notebook to read.
    ARM_DIR="~{arm}"
    mkdir -p "bundle/${ARM_DIR}"
    cp summary/results_manifest.tsv summary/calibration_summary.tsv "bundle/${ARM_DIR}/"
    if [[ -d summary/qc ]]; then
      cp -r summary/qc "bundle/${ARM_DIR}/qc"
    fi
    tar czf "~{arm}.summary.tar.gz" -C bundle "${ARM_DIR}"
    du -sh "~{arm}.summary.tar.gz"
  >>>

  output {
    File results_manifest = "summary/results_manifest.tsv"
    File calibration_summary = "summary/calibration_summary.tsv"
    Array[File] qq_plots = glob("summary/qc/*/qq_cct.png")
    Array[File] manhattan_plots = glob("summary/qc/*/manhattan_cct.png")
    Array[File] top_hits_tables = glob("summary/qc/*/top_hits.tsv")
    File summary_tar = "~{arm}.summary.tar.gz"
  }

  runtime {
    docker: docker
    cpu: cpu
    memory: memory_gb + " GB"
    # SSD, not HDD. This task is bytes-moved, not parsing: it localizes the
    # merged table once and reads it twice. pd-standard gives 0.12 MiB/s per GiB
    # and disk_gb = size x 2.0 + 50, so gzipping the input to 7.48 GB also
    # shrinks the disk to 65 GB and DROPS throughput to 7.8 MB/s -- compression
    # alone only takes 1.80 h to 0.82 h. pd-ssd is 0.48 MiB/s per GiB, so the
    # same 65 GB runs at 31.2 MB/s and the task finishes in 0.20 h for the same
    # ~$0.003. Raising disk_gb_floor to 250 on HDD reaches the same place, but
    # that is a task input (so it is hashed) and asks for a 265 GB disk against
    # a shared regional quota, where this asks for 65 GB.
    disks: "local-disk " + disk_gb + " SSD"
    preemptible: preemptible
    maxRetries: 2
  }
}

task AnnotateRuTest {
  input {
    File source_vcf
    String chrom
    File simple_repeat_bed
    File annotate_script
    String docker = "us-central1-docker.pkg.dev/broad-dsp-lrma/aou-lr/felix-pilot:0.2.0"
    Int cpu = 2
    Int memory_gb = 8
    Int disk_gb_floor = 50
    Float disk_gb_multiplier = 2.5
    Int preemptible = 1
  }

  Int disk_gb = ceil(size(source_vcf, "GB") * disk_gb_multiplier) + disk_gb_floor

  command <<<
    set -euo pipefail
    CHROM=$(python3 -c 'import sys; print(sys.argv[1].strip().strip(chr(34)).strip(chr(39)))' "~{chrom}")
    OUT="annotated/annotated.ru.vcf.gz"
    mkdir -p annotated

    python3 "~{annotate_script}" \
      --vcf "~{source_vcf}" \
      --simple-repeat-bed "~{simple_repeat_bed}" \
      --out "${OUT}"

    bcftools index -t "${OUT}" || bcftools index -c "${OUT}"
    printf '%s\n' "${CHROM}" > annotated/chrom.txt
  >>>

  output {
    File annotated_vcf = "annotated/annotated.ru.vcf.gz"
    String chrom_id = read_string("annotated/chrom.txt")
  }

  runtime {
    docker: docker
    cpu: cpu
    memory: memory_gb + " GB"
    disks: "local-disk " + disk_gb + " HDD"
    preemptible: preemptible
  }
}

task ExtractRuDosage {
  input {
    File annotated_vcf
    String chrom
    File analysis_samples
    Int num_ancs
    String docker = "us-central1-docker.pkg.dev/broad-dsp-lrma/aou-lr/felix-pilot:0.2.0"
    Int cpu = 4
    Int memory_gb = 16
    Int disk_gb_floor = 100
    Float disk_gb_multiplier = 3.5
    Int preemptible = 0
  }

  Int disk_gb = ceil(size(annotated_vcf, "GB") * disk_gb_multiplier) + disk_gb_floor

  command <<<
    set -euo pipefail
    mkdir -p extract ordered
    CHROM=$(python3 -c 'import sys; print(sys.argv[1].strip().strip(chr(34)).strip(chr(39)))' "~{chrom}")

    extract-tracts-flare \
      --vcf "~{annotated_vcf}" \
      --num-ancs ~{num_ancs} \
      --output-dir extract \
      --compress-output \
      --samples "~{analysis_samples}"

    python3 - <<'PY'
import re
import shutil
from pathlib import Path

src = Path("extract")
od = Path("ordered")
od.mkdir(parents=True, exist_ok=True)

for p in sorted(src.glob("*.dosage.txt.gz")):
    if ".collapse." in p.name or ".split." in p.name:
        continue
    m = re.search(r"\.anc(\d+)\.dosage\.txt\.gz$", p.name)
    if not m:
        raise SystemExit(f"unexpected dosage name: {p.name}")
    shutil.copy2(p, od / f"anc_{int(m.group(1)):02d}.dosage.txt.gz")
for p in sorted(src.glob("*.collapse.anc*.dosage.txt.gz")):
    m = re.search(r"\.anc(\d+)\.dosage\.txt\.gz$", p.name)
    if m:
        shutil.copy2(p, od / f"collapse_anc_{int(m.group(1)):02d}.dosage.txt.gz")
for p in sorted(src.glob("*.split.anc*.dosage.txt.gz")):
    m = re.search(r"\.anc(\d+)\.dosage\.txt\.gz$", p.name)
    if m:
        shutil.copy2(p, od / f"split_anc_{int(m.group(1)):02d}.dosage.txt.gz")

order = src / "dosage_sample_order.txt"
if not order.exists():
    raise SystemExit("missing extract/dosage_sample_order.txt")
shutil.copy2(order, od / "dosage_sample_order.txt")
PY

    cp "~{analysis_samples}" ordered/analysis_samples.txt
    python3 -c "import pathlib,sys; pathlib.Path('ordered/chrom.txt').write_text(sys.argv[1].strip().strip(chr(34)).strip(chr(39)) + chr(10))" "${CHROM}"
  >>>

  output {
    String chrom_id = read_string("ordered/chrom.txt")
    Array[File] dosage_files = glob("ordered/anc_*.dosage.txt.gz")
    Array[File] collapse_files = glob("ordered/collapse_anc_*.dosage.txt.gz")
    Array[File] split_files = glob("ordered/split_anc_*.dosage.txt.gz")
  }

  runtime {
    docker: docker
    cpu: cpu
    memory: memory_gb + " GB"
    disks: "local-disk " + disk_gb + " HDD"
    preemptible: preemptible
  }
}

task ScoreRuDosage {
  input {
    String phenotype
    String chrom
    String encoding
    File null_export_tar
    Array[File] dosage_files
    Array[File] collapse_files
    Array[File] split_files
    Int score_threads = 8
    Int chunk_size = 2048
    String docker = "us-central1-docker.pkg.dev/broad-dsp-lrma/aou-lr/felix-pilot:0.2.0"
    Int cpu = 8
    Int memory_gb = 32
    Int disk_gb = 200
    Int preemptible = 0
  }

  command <<<
    set -euo pipefail
    CHROM=$(python3 -c 'import sys; print(sys.argv[1].strip().strip(chr(34)).strip(chr(39)))' "~{chrom}")
    PHENO=$(python3 -c 'import sys; print(sys.argv[1].strip().strip(chr(34)).strip(chr(39)))' "~{phenotype}")
    ENC=$(python3 -c 'import sys; print(sys.argv[1].strip().strip(chr(34)).strip(chr(39)))' "~{encoding}")
    OUT="${PHENO}.${CHROM}.ru_${ENC}.felix.tsv"

    tar xzf "~{null_export_tar}"
    NULL_DIR=$(tar tzf "~{null_export_tar}" | head -n1 | cut -d/ -f1)
    if [[ ! -d "${NULL_DIR}" ]]; then
      NULL_DIR=null_export
    fi

    DOSAGE_ARGS=()
    case "${ENC}" in
      dosage)
        for f in ~{sep=" " dosage_files}; do DOSAGE_ARGS+=("$f"); done
        ;;
      collapse)
        for f in ~{sep=" " collapse_files}; do DOSAGE_ARGS+=("$f"); done
        ;;
      split)
        for f in ~{sep=" " split_files}; do DOSAGE_ARGS+=("$f"); done
        ;;
      *)
        echo "encoding must be dosage, collapse, or split; got ${ENC}" >&2
        exit 1
        ;;
    esac

    if [[ "${#DOSAGE_ARGS[@]}" -lt 1 ]]; then
      echo "no dosage files for encoding=${ENC}" >&2
      exit 1
    fi

    tractor-mix-score \
      --null-export "${NULL_DIR}" \
      --dosage-files "${DOSAGE_ARGS[@]}" \
      --out "${OUT}" \
      --mode felix \
      --threads ~{score_threads} \
      --chunk-size ~{chunk_size}

    cp "${OUT}" results.ru.felix.tsv
    printf '%s\n' "${OUT}" > results.output_name.txt
  >>>

  output {
    File results_tsv = "results.ru.felix.tsv"
    File results_name = "results.output_name.txt"
  }

  runtime {
    docker: docker
    cpu: cpu
    memory: memory_gb + " GB"
    disks: "local-disk " + disk_gb + " HDD"
    preemptible: preemptible
  }
}

# One tested unit per catalog repeat locus: every length-changing record inside a
# TRExplorer interval contributes its signed length change to that haplotype's
# locus dosage (scripts/aggregate_repeat_loci.py --out-vcf). This replaces the
# record-level annotator for ru_engine="felix" — see felix/REPEAT_DOSAGE.md §6.
task AggregateRepeatLoci {
  input {
    File source_vcf
    File? source_vcf_index
    String chrom
    File repeat_catalog_bed
    File analysis_samples
    File aggregate_script
    Int num_ancs = 5
    Int flank = 10
    Int max_bp = 100000
    String apply_filters = "PASS,."
    String docker = "us-central1-docker.pkg.dev/broad-dsp-lrma/aou-lr/felix-pilot:0.2.0"
    Int cpu = 1
    Int memory_gb = 6
    Int disk_gb_floor = 100
    Float disk_gb_multiplier = 3.0
    Int preemptible = 0
  }

  Int disk_gb = ceil(size(source_vcf, "GB") * disk_gb_multiplier) + disk_gb_floor

  command <<<
    set -euo pipefail
    CHROM=$(python3 -c 'import sys; print(sys.argv[1].strip().strip(chr(34)).strip(chr(39)))' "~{chrom}")
    mkdir -p loci src

    # The aggregator seeks with bcftools, so the source VCF needs an index next
    # to it. Cromwell localizes the index (when given) to a different directory,
    # so symlink both into one place; index here when none was provided.
    SRC="src/$(basename "~{source_vcf}")"
    ln -s "~{source_vcf}" "${SRC}"
    IDX="~{default='' source_vcf_index}"
    if [[ -n "${IDX}" ]]; then
      ln -s "${IDX}" "src/$(basename "${IDX}")"
    else
      bcftools index -t "${SRC}" || bcftools index -c "${SRC}"
    fi

    # --num-ancs fails loudly if FLARE emitted a code the downstream tests
    # cannot represent, instead of silently folding it into another ancestry.
    python3 "~{aggregate_script}" \
      --vcf "${SRC}" \
      --catalog-bed "~{repeat_catalog_bed}" \
      --chrom "${CHROM}" \
      --samples "~{analysis_samples}" \
      --apply-filters "~{apply_filters}" \
      --flank ~{flank} \
      --max-bp ~{max_bp} \
      --with-ancestry \
      --num-ancs ~{num_ancs} \
      --out-vcf loci/loci.vcf \
      --out-alleles loci/loci.alleles.tsv.gz \
      --out-summary loci/loci.summary.json

    bcftools view -Oz -o loci/loci.vcf.gz loci/loci.vcf
    bcftools index -c loci/loci.vcf.gz
    rm -f loci/loci.vcf
    cat loci/loci.summary.json
  >>>

  output {
    File locus_vcf = "loci/loci.vcf.gz"
    File locus_vcf_index = "loci/loci.vcf.gz.csi"
    File alleles_tsv = "loci/loci.alleles.tsv.gz"
    File summary_json = "loci/loci.summary.json"
  }

  runtime {
    docker: docker
    cpu: cpu
    memory: memory_gb + " GB"
    disks: "local-disk " + disk_gb + " HDD"
    preemptible: preemptible
  }
}

# Repeat-SV dosages for FELIX's admixed VCF path: DS{k} / ANC{k} from RU_TEST
# loci (REF-relative repeat units, scaled per locus; scripts/write_admixed_dosage_vcf.py).
#
# coding = "ru"      repeat units relative to REF (the test arm)
# coding = "carrier" 1 per non-reference repeat haplotype (the biallelic
#                    comparator arm). Same loci, same carriers, same ANC{k} —
#                    only the weighting by repeat count differs, so the two
#                    runs differ in exactly one thing.
task WriteRuAdmixedVcf {
  input {
    File annotated_vcf
    String chrom
    File analysis_samples
    File write_script
    Int num_ancs = 5
    String missing_policy = "ref"
    Float max_hap_missing = 0.1
    String coding = "ru"
    String docker = "us-central1-docker.pkg.dev/broad-dsp-lrma/aou-lr/felix-pilot:0.2.0"
    Int cpu = 1
    Int memory_gb = 6
    Int disk_gb_floor = 50
    Float disk_gb_multiplier = 3.0
    Int preemptible = 0
  }

  Int disk_gb = ceil(size(annotated_vcf, "GB") * disk_gb_multiplier) + disk_gb_floor

  command <<<
    set -euo pipefail
    CHROM=$(python3 -c 'import sys; print(sys.argv[1].strip().strip(chr(34)).strip(chr(39)))' "~{chrom}")
    mkdir -p ru
    python3 "~{write_script}" \
      --vcf "~{annotated_vcf}" \
      --num-ancs ~{num_ancs} \
      --samples "~{analysis_samples}" \
      --missing "~{missing_policy}" \
      --max-hap-missing ~{max_hap_missing} \
      --coding "~{coding}" \
      --out ru/ru.admixed.vcf \
      --stats-json ru/ru.admixed.stats.json
    bcftools view -Oz -o ru/ru.admixed.vcf.gz ru/ru.admixed.vcf
    bcftools index -c ru/ru.admixed.vcf.gz
    rm -f ru/ru.admixed.vcf
    cat ru/ru.admixed.stats.json
    printf '%s\n' "${CHROM}" > ru/chrom.txt
  >>>

  output {
    File admixed_vcf = "ru/ru.admixed.vcf.gz"
    File admixed_vcf_index = "ru/ru.admixed.vcf.gz.csi"
    File stats_json = "ru/ru.admixed.stats.json"
  }

  runtime {
    docker: docker
    cpu: cpu
    memory: memory_gb + " GB"
    disks: "local-disk " + disk_gb + " HDD"
    preemptible: preemptible
  }
}

# FELIX step 2 on the repeat dosage VCF with carrier-count QC
# (FELIX_DOSAGE_QC=carrier; needs the patched image, felix/patches/).
#
# dosage_qc stays "carrier" for the biallelic comparator arm as well: both
# codings mark the same carrier haplotypes, so one threshold keeps the two runs
# on nearly the same locus set (they can still differ where a sample's gain and
# loss cancel under the repeat-unit coding; CompareRuArms reports that as
# n_loci_*_only). "none" restores FELIX's allele-count QC, which is meaningful
# only for the carrier coding — sum(DS) is not an allele count under
# REF-relative repeat units — and widens the gap between the arms' locus sets.
task RunFelixStep2RuVcf {
  input {
    File admixed_vcf
    File admixed_vcf_index
    String chrom
    String phenotype
    String dosage_qc = "carrier"
    File null_rda
    File variance_ratio
    File samples_used
    File sparse_grm_mtx
    File sparse_grm_sample_ids
    File run_step2_script
    Int num_ancs = 5
    Int min_carriers = 20
    Float pvalcutoff_of_haplotype = 0.05
    String docker = "us-central1-docker.pkg.dev/broad-dsp-lrma/aou-lr/felix-pilot:0.2.0"
    Int cpu = 1
    Int memory_gb = 6
    Int disk_gb = 50
    Int preemptible = 3
  }

  command <<<
    set -euo pipefail
    mkdir -p step2
    CHROM=$(python3 -c 'import sys; print(sys.argv[1].strip().strip(chr(34)).strip(chr(39)))' "~{chrom}")
    PHENO=$(python3 -c 'import sys; print(sys.argv[1].strip().strip(chr(34)).strip(chr(39)))' "~{phenotype}")
    ln -sf "~{null_rda}" "step2/${PHENO}.rda"
    ln -sf "~{variance_ratio}" "step2/${PHENO}.varianceRatio.txt"
    ln -sf "~{admixed_vcf}" ru.admixed.vcf.gz
    ln -sf "~{admixed_vcf_index}" ru.admixed.vcf.gz.csi

    # stderr is DISCARDED.
    #
    # FELIX's VCF reader prints one line per multiallelic record
    # (src/VCF.cpp:181, "Warning: skipping multiallelic variant"), and the
    # repeat-dosage VCF is multiallelic BY CONSTRUCTION -- one ALT per distinct
    # length change, ~5.8 per locus. That is ~2.4M lines per phenotype per arm,
    # across two arms and 951 phenotypes. Under the GCP Batch backend the agent
    # ingests runnable stderr into Cloud Logging as batch_task_logs at
    # $0.50/GiB, charged once at ingestion and unrecoverable by deleting.
    #
    # The message is wrong twice over, so nothing diagnostic is lost: VCF.cpp
    # prints it and then does NOT skip -- there is no continue and no return, it
    # falls through to t_alt = alts()[0] -- and the statistic does not come from
    # the ALT at all. m_fmtField is DS here, so t_altCounts accumulates the DS
    # values (VCF.cpp:223), and write_admixed_dosage_vcf.py has already summed
    # every ALT's contribution into one per-haplotype dosage per DS{k}. The
    # multiple ALTs are informational; only the reported ALT string takes
    # alts()[0].
    #
    # This replaces an earlier capture-count-gzip-and-tail-on-failure block.
    # That version had to grep the stream, and a grep that matched nothing
    # exited 1 under `set -o pipefail` and killed the task on exactly the
    # HEALTHY case. `2> /dev/null` has no pattern to miss and no pipeline to
    # break, which is why it is preferred here: the task runs ~42,000 times
    # across both arms and the simplest construct that cannot fail wins. The
    # cost is that a failing shard leaves no stderr anywhere -- diagnose those
    # by re-running the single shard with this redirect removed.
    #
    # --min-mac is a carrier count here (samples with a non-zero dosage).
    Rscript "~{run_step2_script}" \
      --vcf-file ru.admixed.vcf.gz \
      --vcf-index ru.admixed.vcf.gz.csi \
      --dosage-qc "~{dosage_qc}" \
      --chrom "${CHROM}" \
      --null-prefix "step2/${PHENO}" \
      --sample-file "~{samples_used}" \
      --sparse-grm "~{sparse_grm_mtx}" \
      --sparse-grm-ids "~{sparse_grm_sample_ids}" \
      --min-mac ~{min_carriers} \
      --n-ancestries ~{num_ancs} \
      --pvalcutoff-of-haplotype ~{pvalcutoff_of_haplotype} \
      --n-threads 1 \
      --out-tsv "${PHENO}.ru.felix.tsv" \
      --out-raw "${PHENO}.ru.felix.raw.txt" \
      2> /dev/null

    mv "${PHENO}.ru.felix.tsv" results.ru.felix.tsv
  >>>

  output {
    File results_tsv = "results.ru.felix.tsv"
  }

  runtime {
    docker: docker
    cpu: cpu
    memory: memory_gb + " GB"
    disks: "local-disk " + disk_gb + " HDD"
    preemptible: preemptible
  }
}

# 2x2 comparison table for one phenotype: local ancestry x repeat-dosage
# coding, from the two merged result TSVs (felix/scripts/compare_ru_arms.py).
task CompareRuArms {
  input {
    String phenotype
    File dosage_tsv
    File biallelic_tsv
    File compare_script
    Int num_ancs = 5
    Float p_threshold = 0.00000005
    String docker = "us-central1-docker.pkg.dev/broad-dsp-lrma/aou-lr/felix-pilot:0.2.0"
    Int cpu = 1
    # 8 GB OOM-kills this genome-wide. compare_ru_arms.py reads both merged
    # result tables into pandas and inner-joins them, so peak is the two inputs
    # plus the join result, all live at once. Measured on a 492,000-row fixture
    # tiled from the chr22 merged table: 0.20 GB of interpreter and pandas plus
    # 3.31 GB per million shared loci, so 8 GB dies above ~2.35M loci.
    # memory_gb is a task INPUT here, so this re-runs CompareRuArms and nothing
    # else; the call site does not bind it, so the default governs.
    Int memory_gb = 32
    Int disk_gb_floor = 20
    Float disk_gb_multiplier = 3.0
    Int preemptible = 1
  }

  Int disk_gb = ceil((size(dosage_tsv, "GB") + size(biallelic_tsv, "GB")) * disk_gb_multiplier) + disk_gb_floor

  command <<<
    set -euo pipefail
    PHENO=$(python3 -c 'import sys; print(sys.argv[1].strip().strip(chr(34)).strip(chr(39)))' "~{phenotype}")
    mkdir -p arms
    python3 "~{compare_script}" \
      --dosage-tsv "~{dosage_tsv}" \
      --biallelic-tsv "~{biallelic_tsv}" \
      --phenotype "${PHENO}" \
      --num-ancs ~{num_ancs} \
      --p-threshold ~{p_threshold} \
      --out-tsv "arms/${PHENO}.ru_arms.tsv.gz" \
      --out-summary-json "arms/${PHENO}.ru_arms.summary.json" \
      --out-md "arms/${PHENO}.ru_arms.md"

    mv "arms/${PHENO}.ru_arms.tsv.gz" arms_table.tsv.gz
    cp "arms/${PHENO}.ru_arms.summary.json" arms_summary.json
    mv "arms/${PHENO}.ru_arms.md" arms_table.md
  >>>

  output {
    File arm_table = "arms_table.tsv.gz"
    File arm_summary = "arms_summary.json"
    File arm_markdown = "arms_table.md"
  }

  runtime {
    docker: docker
    cpu: cpu
    memory: memory_gb + " GB"
    disks: "local-disk " + disk_gb + " HDD"
    preemptible: preemptible
    maxRetries: 2
  }
}

# One small tar.gz per run: every table, JSON and plot needed to decide whether
# the run worked, with the per-phenotype scan TSVs left in the bucket. Arrays
# are renamed by chromosome on the way in, so the 22 loci.summary.json files do
# not collide.
task BundleRunSummary {
  input {
    Array[File] arm_tars
    Array[String] chroms
    Array[File] locus_summaries = []
    Array[File] admixed_stats = []
    Array[File] biallelic_stats = []
    Array[File] arm_tables = []
    Array[File] arm_summaries = []
    Array[File] arm_markdown = []
    String prefix = "felix_run"
    String docker = "us-central1-docker.pkg.dev/broad-dsp-lrma/aou-lr/felix-pilot:0.2.0"
    Int cpu = 1
    Int memory_gb = 4
    Int disk_gb = 50
    Int preemptible = 1
  }

  command <<<
    set -euo pipefail
    mkdir -p bundle/repeat_qc bundle/arm_compare

    # each arm tar already carries its own top-level directory (main/, ru/, ...)
    ARM_TARS=(~{sep=" " arm_tars})
    for t in ${ARM_TARS[@]+"${ARM_TARS[@]}"}; do
      tar xzf "${t}" -C bundle
    done

    cat > chroms.txt <<'CHROM_EOF'
~{sep="\n" chroms}
CHROM_EOF
    mapfile -t CHROM_ARR < chroms.txt

    # rename per-chromosome QC files so 22 identically named files can coexist
    copy_by_chrom() {
      local dest="$1"; local tag="$2"; shift 2
      local i=0
      for f in ${@+"$@"}; do
        [[ -n "${f}" ]] || continue
        local c="${CHROM_ARR[$i]:-shard${i}}"
        cp "${f}" "${dest}/${c}.${tag}"
        i=$((i + 1))
      done
    }
    LOCUS_SUMMARIES=(~{sep=" " locus_summaries})
    ADMIXED_STATS=(~{sep=" " admixed_stats})
    BIALLELIC_STATS=(~{sep=" " biallelic_stats})
    copy_by_chrom bundle/repeat_qc loci.summary.json ${LOCUS_SUMMARIES[@]+"${LOCUS_SUMMARIES[@]}"}
    copy_by_chrom bundle/repeat_qc ru.admixed.stats.json ${ADMIXED_STATS[@]+"${ADMIXED_STATS[@]}"}
    copy_by_chrom bundle/repeat_qc ru_biallelic.admixed.stats.json ${BIALLELIC_STATS[@]+"${BIALLELIC_STATS[@]}"}

    # arm-comparison outputs are per phenotype; keep their own basenames
    ARM_FILES=(~{sep=" " arm_tables} ~{sep=" " arm_summaries} ~{sep=" " arm_markdown})
    for f in ${ARM_FILES[@]+"${ARM_FILES[@]}"}; do
      cp "${f}" bundle/arm_compare/
    done

    rmdir bundle/repeat_qc bundle/arm_compare 2>/dev/null || true
    find bundle -type f | sort > bundle/MANIFEST.txt
    wc -l < bundle/MANIFEST.txt
    tar czf "~{prefix}.summary.tar.gz" -C bundle .
    du -sh "~{prefix}.summary.tar.gz"
  >>>

  output {
    File summary_tar = "~{prefix}.summary.tar.gz"
  }

  runtime {
    docker: docker
    cpu: cpu
    memory: memory_gb + " GB"
    disks: "local-disk " + disk_gb + " HDD"
    preemptible: preemptible
    maxRetries: 2
  }
}

workflow FelixGenome {
  input {
    Array[String] chroms

  # Parallel per-chrom VCFs (same length / order as chroms).
  # Supply joint_vcfs OR both phase_vcfs and flare_vcfs, or all three. felixla
  # (Pack) assigns each phase record to its covering FLARE marker itself, so the
  # ordinary scan needs only the raw phased callset plus the FLARE anc VCF. When
  # phase_vcfs and flare_vcfs are given they are what Pack reads; the repeat
  # branch reads joint_vcfs (GT:AN1:AN2) when given, else phase_vcfs.
    Array[File]? joint_vcfs
    Array[File]? phase_vcfs
    Array[File]? flare_vcfs
    Int num_ancs = 5

    File analysis_samples
    File pheno_cov
    File selected_phenotypes
    File covariate_columns
    String trait_type = "binary"

    Array[File] grm_vcfs

    File build_saige_grm_script
    File make_plink_keep_script
    File fit_felix_null_script
    File? export_null_script
    File run_felix_step2_script
    File summarize_script

  # Step-1 nulls from felix/wdl/FitFelixNulls.wdl, in the order of
  # selected_phenotypes (that workflow's null_manifest records the order it
  # used; feed both workflows the same phenotype file). Supplying
  # precomputed_null_rdas skips FitFelixNull entirely — the GRM is still built,
  # because Pack and Step 2 need the PLINK set. The other three arrays must be
  # supplied with it; CheckNullArrays fails the run if the lengths disagree.
    Array[File]? precomputed_null_rdas
    Array[File]? precomputed_variance_ratios
    Array[File]? precomputed_samples_used
    Array[File]? precomputed_null_export_tars

  # Optional repeat-dosage branch. Both engines read GT:AN1:AN2, so pass
  # joint_vcfs (PropagateFlareAncestry output), not phase_vcfs alone.
  #
  # "felix" (default): one test per catalog repeat locus. Needs
  #   repeat_catalog_bed (TRExplorer intervals: chrom start end TRID motifs),
  #   aggregate_repeat_loci_script and write_admixed_dosage_script, plus the
  #   patched image for carrier QC. felix/REPEAT_DOSAGE.md.
  # "rust": legacy record-level tractor-mix-score encodings (no SPA, one test
  #   per VCF record), for the encoding comparison only. Needs
  #   simple_repeat_bed and annotate_repeat_units_script.
    String ru_engine = "felix"

    File? repeat_catalog_bed
    File? aggregate_repeat_loci_script
    File? write_admixed_dosage_script

  # Biallelic comparator arm (ru_engine="felix" only). A second pass of the
  # same loci through the same Step 1 null, GRM and samples, with DS{k} coded
  # as the non-reference haplotype count instead of repeat units. Together with
  # FELIX's own ancestry-collapsed test (p.value_ancALL = P_hom_admixed) it
  # fills the 2x2: pooled/local-ancestry x biallelic/repeat-dosage, where
  # pooled+biallelic is the conventional-GWAS baseline. Needs
  # compare_ru_arms_script for the summary table.
    Boolean ru_biallelic_arm = true
    File? compare_ru_arms_script
    Float ru_arm_p_threshold = 0.00000005
    Int ru_min_carriers = 20
    Int ru_flank = 10
    Int ru_max_bp = 100000
    String ru_apply_filters = "PASS,."
    String ru_missing_policy = "ref"
    Float ru_max_hap_missing = 0.1

    File? simple_repeat_bed
    File? annotate_repeat_units_script

    Float relatedness_cutoff = 0.05
    Int num_random_markers = 2000
    Int min_mac = 50
    Float pvalcutoff_of_haplotype = 0.05
    Float summarize_p_threshold = 0.00000005
    Int summarize_top_n = 50

    Int pack_cpu = 1
    Int pack_memory_gb = 6
    Int pack_disk_gb_floor = 100
    Float pack_disk_gb_multiplier = 3.5
    Int step2_disk_gb_floor = 50
    Float step2_disk_gb_multiplier = 4.0
    Int ru_disk_gb_floor = 100
    Float ru_disk_gb_multiplier = 3.5
    Int score_threads = 8
    Int chunk_size = 2048
    Int score_disk_gb = 200

    String docker = "us-central1-docker.pkg.dev/broad-dsp-lrma/aou-lr/felix-pilot:0.2.0"
  }

  Array[String] phenotypes = read_lines(selected_phenotypes)
  Array[String] ru_encodings = ["dosage", "collapse", "split"]

  call CheckChromVcfPairs as Check {
    input:
      n_chroms = length(chroms),
      n_joint_vcfs = length(select_first([joint_vcfs, []])),
      n_phase_vcfs = length(select_first([phase_vcfs, []])),
      n_flare_vcfs = length(select_first([flare_vcfs, []])),
      docker = docker
  }

  call BuildSaigePlinkAndSparseGRM as MakeGRM {
    input:
      grm_vcfs = grm_vcfs,
      analysis_samples = analysis_samples,
      build_script = build_saige_grm_script,
      make_plink_keep_script = make_plink_keep_script,
      relatedness_cutoff = relatedness_cutoff,
      num_random_markers = num_random_markers,
      docker = docker
  }

  # Step-1 nulls: fit here, or take them from a FitFelixNulls run. They depend
  # only on the GRM, the covariates and the phenotype column, so a long
  # phenotype list is worth fitting once in that workflow and reusing.
  if (!defined(precomputed_null_rdas)) {
    scatter (pheno in phenotypes) {
      call FitFelixNull as Null {
        input:
          pheno_cov = pheno_cov,
          phenotype = pheno,
          covariate_columns = covariate_columns,
          analysis_samples = analysis_samples,
          plink_bed = MakeGRM.plink_bed,
          plink_bim = MakeGRM.plink_bim,
          plink_fam = MakeGRM.plink_fam,
          sparse_grm_mtx = MakeGRM.sparse_grm_mtx,
          sparse_grm_sample_ids = MakeGRM.sparse_grm_sample_ids,
          fit_null_script = fit_felix_null_script,
          export_null_script = export_null_script,
          trait_type = trait_type,
          docker = docker
      }
    }
  }

  # Index these, not Null.*, everywhere below. Order is the contract: the
  # arrays must be in the order of read_lines(selected_phenotypes), which is
  # what FitFelixNulls writes and records in its null_manifest.
  Array[File] null_rda_by_pheno =
    select_first([precomputed_null_rdas, Null.null_rda])
  Array[File] variance_ratio_by_pheno =
    select_first([precomputed_variance_ratios, Null.variance_ratio])
  Array[File] samples_used_by_pheno =
    select_first([precomputed_samples_used, Null.samples_used])
  Array[File] null_export_tar_by_pheno =
    select_first([precomputed_null_export_tars, Null.null_export_tar])

  call CheckNullArrays as CheckNulls {
    input:
      n_phenotypes = length(phenotypes),
      n_null_rdas = length(null_rda_by_pheno),
      n_variance_ratios = length(variance_ratio_by_pheno),
      n_samples_used = length(samples_used_by_pheno),
      precomputed = defined(precomputed_null_rdas),
      docker = docker
  }

  scatter (i in range(length(chroms))) {
    String chrom = chroms[i]
    # Cromwell cannot index an optional array; select_first unwraps it
    # (CheckChromVcfPairs has already checked which inputs were given).
    File phase_for_pack = if defined(phase_vcfs) then select_first([phase_vcfs])[i] else select_first([joint_vcfs])[i]
    File flare_for_pack = if defined(flare_vcfs) then select_first([flare_vcfs])[i] else select_first([joint_vcfs])[i]
    File phase_for_ru = if defined(joint_vcfs) then select_first([joint_vcfs])[i] else select_first([phase_vcfs])[i]

    call PackFelixla as Pack {
      input:
        phase_vcf = phase_for_pack,
        flare_vcf = flare_for_pack,
        analysis_samples = analysis_samples,
        chrom = chrom,
        num_ancs = num_ancs,
        docker = docker,
        cpu = pack_cpu,
        memory_gb = pack_memory_gb,
        disk_gb_floor = pack_disk_gb_floor,
        disk_gb_multiplier = pack_disk_gb_multiplier
    }

    # Locus-level dosage (default). One test per catalog repeat locus: the
    # record-level annotator is not used here, because it tests each VCF record
    # separately and drops records whose alleles are not a whole number of
    # repeat units (felix/REPEAT_DOSAGE.md §6.2, felix/eval/record_vs_locus/).
    if (ru_engine == "felix" && defined(repeat_catalog_bed)
        && defined(aggregate_repeat_loci_script) && defined(write_admixed_dosage_script)) {
      call AggregateRepeatLoci as RuLoci {
        input:
          source_vcf = phase_for_ru,
          chrom = chrom,
          repeat_catalog_bed = select_first([repeat_catalog_bed]),
          analysis_samples = analysis_samples,
          aggregate_script = select_first([aggregate_repeat_loci_script]),
          num_ancs = num_ancs,
          flank = ru_flank,
          max_bp = ru_max_bp,
          apply_filters = ru_apply_filters,
          docker = docker,
          disk_gb_floor = ru_disk_gb_floor,
          disk_gb_multiplier = ru_disk_gb_multiplier
      }

      call WriteRuAdmixedVcf as RuVcf {
        input:
          annotated_vcf = RuLoci.locus_vcf,
          chrom = chrom,
          analysis_samples = analysis_samples,
          write_script = select_first([write_admixed_dosage_script]),
          num_ancs = num_ancs,
          missing_policy = ru_missing_policy,
          max_hap_missing = ru_max_hap_missing,
          docker = docker
      }

      scatter (j in range(CheckNulls.n_phenotypes_out)) {
        call RunFelixStep2RuVcf as RuStep2 {
          input:
            admixed_vcf = RuVcf.admixed_vcf,
            admixed_vcf_index = RuVcf.admixed_vcf_index,
            chrom = chrom,
            phenotype = phenotypes[j],
            null_rda = null_rda_by_pheno[j],
            variance_ratio = variance_ratio_by_pheno[j],
            samples_used = samples_used_by_pheno[j],
            sparse_grm_mtx = MakeGRM.sparse_grm_mtx,
            sparse_grm_sample_ids = MakeGRM.sparse_grm_sample_ids,
            run_step2_script = run_felix_step2_script,
            num_ancs = num_ancs,
            min_carriers = ru_min_carriers,
            pvalcutoff_of_haplotype = pvalcutoff_of_haplotype,
            docker = docker
        }
      }

      # Baseline arm: same loci, same null, same carriers; carrier coding.
      if (ru_biallelic_arm) {
        call WriteRuAdmixedVcf as BaVcf {
          input:
            annotated_vcf = RuLoci.locus_vcf,
            chrom = chrom,
            analysis_samples = analysis_samples,
            write_script = select_first([write_admixed_dosage_script]),
            num_ancs = num_ancs,
            missing_policy = ru_missing_policy,
            max_hap_missing = ru_max_hap_missing,
            coding = "carrier",
            docker = docker
        }

        scatter (j in range(CheckNulls.n_phenotypes_out)) {
          call RunFelixStep2RuVcf as BaStep2 {
            input:
              admixed_vcf = BaVcf.admixed_vcf,
              admixed_vcf_index = BaVcf.admixed_vcf_index,
              chrom = chrom,
              phenotype = phenotypes[j],
              null_rda = null_rda_by_pheno[j],
              variance_ratio = variance_ratio_by_pheno[j],
              samples_used = samples_used_by_pheno[j],
              sparse_grm_mtx = MakeGRM.sparse_grm_mtx,
              sparse_grm_sample_ids = MakeGRM.sparse_grm_sample_ids,
              run_step2_script = run_felix_step2_script,
              num_ancs = num_ancs,
              min_carriers = ru_min_carriers,
              pvalcutoff_of_haplotype = pvalcutoff_of_haplotype,
              docker = docker
          }
        }
      }
    }

    # Legacy record-level encodings, kept for the encoding comparison only.
    if (ru_engine == "rust" && defined(simple_repeat_bed)
        && defined(annotate_repeat_units_script)) {
      call AnnotateRuTest as RuAnnot {
        input:
          source_vcf = phase_for_ru,
          chrom = chrom,
          simple_repeat_bed = select_first([simple_repeat_bed]),
          annotate_script = select_first([annotate_repeat_units_script]),
          docker = docker
      }

      call ExtractRuDosage as RuExtract {
        input:
          annotated_vcf = RuAnnot.annotated_vcf,
          chrom = chrom,
          analysis_samples = analysis_samples,
          num_ancs = num_ancs,
          docker = docker,
          disk_gb_floor = ru_disk_gb_floor,
          disk_gb_multiplier = ru_disk_gb_multiplier
      }

      scatter (j in range(CheckNulls.n_phenotypes_out)) {
        scatter (enc in ru_encodings) {
          call ScoreRuDosage as RuScore {
            input:
              phenotype = phenotypes[j],
              chrom = chrom,
              encoding = enc,
              null_export_tar = null_export_tar_by_pheno[j],
              dosage_files = RuExtract.dosage_files,
              collapse_files = RuExtract.collapse_files,
              split_files = RuExtract.split_files,
              score_threads = score_threads,
              chunk_size = chunk_size,
              docker = docker,
              disk_gb = score_disk_gb
          }
        }
      }
    }

    scatter (j in range(CheckNulls.n_phenotypes_out)) {
      call RunFelixStep2 as Step2 {
        input:
          packed_tar = Pack.packed_tar,
          chrom = chrom,
          phenotype = phenotypes[j],
          null_rda = null_rda_by_pheno[j],
          variance_ratio = variance_ratio_by_pheno[j],
          samples_used = samples_used_by_pheno[j],
          sparse_grm_mtx = MakeGRM.sparse_grm_mtx,
          sparse_grm_sample_ids = MakeGRM.sparse_grm_sample_ids,
          run_step2_script = run_felix_step2_script,
          num_ancs = num_ancs,
          min_mac = min_mac,
          pvalcutoff_of_haplotype = pvalcutoff_of_haplotype,
          docker = docker,
          disk_gb_floor = step2_disk_gb_floor,
          disk_gb_multiplier = step2_disk_gb_multiplier
      }
    }
  }

  Array[Array[File]] results_by_pheno = transpose(Step2.results_tsv)
  # Chroms with RU results (felix engine), then per phenotype across chroms.
  Array[Array[File]] ru_by_chrom = select_all(RuStep2.results_tsv)
  Boolean has_ru = length(ru_by_chrom) > 0

  if (has_ru) {
    Array[Array[File]] ru_by_pheno = transpose(ru_by_chrom)
    scatter (i in range(CheckNulls.n_phenotypes_out)) {
      call ConcatPhenotypeScores as RuConcat {
        input:
          phenotype = phenotypes[i] + ".ru",
          shard_tsvs = ru_by_pheno[i],
          docker = docker
      }
    }
  }

  # Biallelic comparator arm, merged the same way, then the 2x2 table.
  Array[Array[File]] ba_by_chrom = select_all(BaStep2.results_tsv)
  Boolean has_ba = length(ba_by_chrom) > 0

  if (has_ba) {
    Array[Array[File]] ba_by_pheno = transpose(ba_by_chrom)
    scatter (i in range(CheckNulls.n_phenotypes_out)) {
      call ConcatPhenotypeScores as BaConcat {
        input:
          phenotype = phenotypes[i] + ".ru_biallelic",
          shard_tsvs = ba_by_pheno[i],
          docker = docker
      }
    }

    if (defined(compare_ru_arms_script)) {
      scatter (i in range(CheckNulls.n_phenotypes_out)) {
        call CompareRuArms as ArmCompare {
          input:
            phenotype = phenotypes[i],
            dosage_tsv = select_first([RuConcat.merged_tsv])[i],
            biallelic_tsv = BaConcat.merged_tsv[i],
            compare_script = select_first([compare_ru_arms_script]),
            num_ancs = num_ancs,
            p_threshold = ru_arm_p_threshold,
            docker = docker
        }
      }
    }
  }

  scatter (i in range(CheckNulls.n_phenotypes_out)) {
    call ConcatPhenotypeScores as Concat {
      input:
        phenotype = phenotypes[i],
        shard_tsvs = results_by_pheno[i],
        docker = docker
    }
  }

  call SummarizeFelixResults as Summarize {
    input:
      result_tsvs = Concat.merged_tsv,
      phenotypes = phenotypes,
      pheno_cov = pheno_cov,
      summarize_script = summarize_script,
      arm = "main",
      p_threshold = summarize_p_threshold,
      top_n = summarize_top_n,
      docker = docker
  }

  # Same summariser on the repeat-dosage results: identical joint P columns, one
  # row per catalog locus instead of per variant, so lambdaGC here is directly
  # comparable with the main arm's from the same null.
  if (has_ru) {
    call SummarizeFelixResults as SummarizeRu {
      input:
        result_tsvs = select_first([RuConcat.merged_tsv]),
        phenotypes = phenotypes,
        pheno_cov = pheno_cov,
        summarize_script = summarize_script,
        arm = "ru",
        p_threshold = summarize_p_threshold,
        top_n = summarize_top_n,
        docker = docker
    }
  }

  if (has_ba) {
    call SummarizeFelixResults as SummarizeBa {
      input:
        result_tsvs = select_first([BaConcat.merged_tsv]),
        phenotypes = phenotypes,
        pheno_cov = pheno_cov,
        summarize_script = summarize_script,
        arm = "ru_biallelic",
        p_threshold = summarize_p_threshold,
        top_n = summarize_top_n,
        docker = docker
    }
  }

  call BundleRunSummary as Bundle {
    input:
      arm_tars = select_all([Summarize.summary_tar, SummarizeRu.summary_tar, SummarizeBa.summary_tar]),
      chroms = chroms,
      locus_summaries = select_all(RuLoci.summary_json),
      admixed_stats = select_all(RuVcf.stats_json),
      biallelic_stats = select_all(BaVcf.stats_json),
      arm_tables = select_first([ArmCompare.arm_table, []]),
      arm_summaries = select_first([ArmCompare.arm_summary, []]),
      arm_markdown = select_first([ArmCompare.arm_markdown, []]),
      docker = docker
  }

  output {
    Int n_chroms = Check.n_chroms_out
    Int n_phenotypes_used = CheckNulls.n_phenotypes_out
    File plink_bed = MakeGRM.plink_bed
    File sparse_grm_mtx = MakeGRM.sparse_grm_mtx
    Array[File] null_rdas = null_rda_by_pheno
    Array[File] null_export_tars = null_export_tar_by_pheno
    File run_summary_tar = Bundle.summary_tar
    File? ru_calibration_summary = SummarizeRu.calibration_summary
    File? ru_results_manifest = SummarizeRu.results_manifest
    Array[File] results_tsvs = Concat.merged_tsv
    Array[File] results_names = Concat.merged_name
    File results_manifest = Summarize.results_manifest
    File calibration_summary = Summarize.calibration_summary
    Array[File] qq_plots = Summarize.qq_plots
    Array[File] manhattan_plots = Summarize.manhattan_plots
    Array[File] top_hits_tables = Summarize.top_hits_tables
    Array[Array[File]] results_tsvs_by_chrom = Step2.results_tsv
    Array[File] packed_tars = Pack.packed_tar
    Array[String] chrom_ids = Pack.chrom_id
    Array[File]? ru_results_tsvs = RuConcat.merged_tsv
    Array[File?] ru_admixed_vcfs = RuVcf.admixed_vcf
    Array[File?] ru_admixed_stats = RuVcf.stats_json
    # Biallelic comparator arm (ru_biallelic_arm) and the 2x2 table per
    # phenotype: arm p-values per locus, lambda / significant-locus counts and
    # the gained/lost contrasts behind the manuscript's benefit claims.
    Array[File]? ru_biallelic_results_tsvs = BaConcat.merged_tsv
    Array[File?] ru_biallelic_admixed_vcfs = BaVcf.admixed_vcf
    Array[File?] ru_biallelic_admixed_stats = BaVcf.stats_json
    Array[File]? ru_arm_tables = ArmCompare.arm_table
    Array[File]? ru_arm_summaries = ArmCompare.arm_summary
    Array[File]? ru_arm_markdown = ArmCompare.arm_markdown
    # Locus/haplotype tallies travel with the results: missingness, possible
    # double counts between the indel and SV calls, ancestry inconsistency.
    Array[File?] ru_locus_summaries = RuLoci.summary_json
    Array[File?] ru_locus_alleles = RuLoci.alleles_tsv
    Array[File?] ru_locus_vcfs = RuLoci.locus_vcf
    # ru_engine="rust" only: per-chrom, per-phenotype, per-encoding score TSVs
    # for the encoding comparison. Not merged into Summarize.
    Array[Array[Array[File]]?] ru_rust_encoding_tsvs = RuScore.results_tsv
  }

  meta {
    description: "Genome-wide FELIX: shared SAIGE mtx GRM + FELIX nulls; per-chr felixla pack/step2; optional locus-level repeat dosages through FELIX step 2 (carrier QC)."
    allowNestedInputs: true
  }
}
