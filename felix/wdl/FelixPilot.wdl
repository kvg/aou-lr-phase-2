version 1.0

# Self-contained FELIX chr22 pilot for Terra Methods Repository.
# Shape: Check → MakeGRM (SAIGE mtx) → FitNull (FELIX step1) → Pack (felixla)
#        → Step2 (admixed SPA) → concat → summarize (P_cct_admixed_c).
#
# Consume lhu1/felix:latest via the felix-pilot image. Do not fork FELIX.

task CheckInputs {
  input {
    String chrom
    Int n_phenotypes
    String docker = "us-central1-docker.pkg.dev/broad-dsp-lrma/aou-lr/felix-pilot:0.1.0"
  }

  command <<<
    set -euo pipefail
    CHROM=$(python3 -c 'import sys; print(sys.argv[1].strip().strip(chr(34)).strip(chr(39)))' "~{chrom}")
    if [[ -z "${CHROM}" ]]; then
      echo "chrom is empty" >&2
      exit 1
    fi
    if [[ "~{n_phenotypes}" -lt 1 ]]; then
      echo "Need at least one phenotype" >&2
      exit 1
    fi
    printf '%s\n' "${CHROM}" > chrom.txt
    echo "OK: chrom=${CHROM} phenotypes=~{n_phenotypes}"
  >>>

  output {
    String chrom_id = read_string("chrom.txt")
    File chrom_txt = "chrom.txt"
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
    String docker = "us-central1-docker.pkg.dev/broad-dsp-lrma/aou-lr/felix-pilot:0.1.0"
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
    String docker = "us-central1-docker.pkg.dev/broad-dsp-lrma/aou-lr/felix-pilot:0.1.0"
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
      --step1-r /usr/local/bin/step1_fitNULLGLMM.R \
      --n-threads ~{cpu} \
      --out-prefix "null/~{phenotype}"

    cp "null/~{phenotype}.rda" "~{phenotype}.null.rda" || \
      cp $(ls -1 null/~{phenotype}*.rda | head -n1) "~{phenotype}.null.rda"
    cp "null/~{phenotype}.varianceRatio.txt" "~{phenotype}.varianceRatio.txt" || \
      cp $(ls -1 null/~{phenotype}*.varianceRatio.txt | head -n1) "~{phenotype}.varianceRatio.txt"
    cp "null/~{phenotype}.pheno.tsv" "~{phenotype}.pheno.tsv"
    cp "null/~{phenotype}.samples.txt" "~{phenotype}.samples.txt"
    cp "null/~{phenotype}.null_meta.tsv" "~{phenotype}.null_meta.tsv"
  >>>

  output {
    File null_rda = "~{phenotype}.null.rda"
    File variance_ratio = "~{phenotype}.varianceRatio.txt"
    File pheno_used = "~{phenotype}.pheno.tsv"
    File samples_used = "~{phenotype}.samples.txt"
    File null_meta = "~{phenotype}.null_meta.tsv"
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
    String docker = "us-central1-docker.pkg.dev/broad-dsp-lrma/aou-lr/felix-pilot:0.1.0"
    Int cpu = 4
    Int memory_gb = 16
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
      # FLARE / joint VCFs often declare ##FORMAT=<ID=GT,...> twice.
      bcftools view -h "${in_vcf}" \
        | awk '/^##FORMAT=<ID=GT,/{ if (gt++) next } { print }' \
        > "${hdr}"
      echo "FORMAT/GT header lines: $(bcftools view -h "${in_vcf}" | grep -c '^##FORMAT=<ID=GT,' || true) -> $(grep -c '^##FORMAT=<ID=GT,' "${hdr}" || true)"
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

    echo "Packing FELIXla chrom=${CHROM} n_ancs=~{num_ancs}"
    echo "phase=${PHASE_IN}"
    echo "flare=${FLARE_IN}"
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
    echo "Packed samples: $(wc -l < "${PREFIX}.samples")"
    ls -lh packed/

    printf '%s\n' "${CHROM}" > packed/chrom.txt
    printf '%s\n' "${PREFIX}" > packed/prefix.txt
    tar czf felixla_packed.tar.gz -C packed .
  >>>

  output {
    File packed_tar = "felixla_packed.tar.gz"
    String chrom_id = read_string("packed/chrom.txt")
    String felixla_prefix = read_string("packed/prefix.txt")
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
    String docker = "us-central1-docker.pkg.dev/broad-dsp-lrma/aou-lr/felix-pilot:0.1.0"
    Int cpu = 8
    Int memory_gb = 16
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
    if [[ ! -f "${PREFIX}.meta" ]]; then
      echo "missing ${PREFIX}.meta after unpack" >&2
      ls -lh packed/ || true
      exit 1
    fi

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
    File results_raw = "~{phenotype}.felix.raw.txt"
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
    String docker = "us-central1-docker.pkg.dev/broad-dsp-lrma/aou-lr/felix-pilot:0.1.0"
    Int cpu = 1
    Int memory_gb = 4
    Int disk_gb_floor = 20
    Float disk_gb_multiplier = 2.5
    Int preemptible = 1
  }

  Int disk_gb = ceil(size(shard_tsvs, "GB") * disk_gb_multiplier) + disk_gb_floor

  command <<<
    set -euo pipefail
    PHENO=$(python3 -c 'import sys; print(sys.argv[1].strip().strip(chr(34)).strip(chr(39)))' "~{phenotype}")
    OUT="${PHENO}.felix.tsv"
    export PHENO OUT
    LIST=shards.txt
    cat > "${LIST}" <<'EOF'
~{sep="\n" shard_tsvs}
EOF

    python3 - <<'PY'
from pathlib import Path
import os

pheno = os.environ["PHENO"]
out = Path(os.environ["OUT"])
shards = [Path(l.strip()) for l in Path("shards.txt").read_text().splitlines() if l.strip()]
if not shards:
    raise SystemExit("no shard TSVs to concatenate")

header = None
nrows = 0
tmp = Path(str(out) + ".unsorted")
with tmp.open("w", encoding="utf-8") as fo:
    for path in shards:
        with path.open(encoding="utf-8") as fi:
            h = fi.readline()
            if not h:
                raise SystemExit(f"empty shard: {path}")
            if header is None:
                header = h
                fo.write(header)
            elif h != header:
                raise SystemExit(
                    f"header mismatch in {path}:\n  got: {h!r}\n  exp: {header!r}"
                )
            for line in fi:
                if line.strip():
                    fo.write(line)
                    nrows += 1

hdr = header.rstrip("\n")
body = tmp.read_text(encoding="utf-8").splitlines()[1:]

def _sort_key(line: str):
    parts = line.split("\t", 2)
    chrom = parts[0] if parts else ""
    pos = parts[1] if len(parts) > 1 else ""
    return (chrom, int(pos) if pos.isdigit() else pos)

body_sorted = sorted(body, key=_sort_key)
with out.open("w", encoding="utf-8") as fo:
    fo.write(hdr + "\n")
    for line in body_sorted:
        fo.write(line + "\n")
tmp.unlink()
print(f"Wrote {out} ({nrows} variant rows from {len(shards)} shards for {pheno})")
PY
    wc -l "${OUT}"
    cp "${OUT}" merged.felix.tsv
    printf '%s\n' "${OUT}" > merged.output_name.txt
  >>>

  output {
    File merged_tsv = "merged.felix.tsv"
    File merged_name = "merged.output_name.txt"
  }

  runtime {
    docker: docker
    cpu: cpu
    memory: memory_gb + " GB"
    disks: "local-disk " + disk_gb + " HDD"
    preemptible: preemptible
  }
}

task SummarizeFelixResults {
  input {
    Array[File] result_tsvs
    Array[String] phenotypes
    File pheno_cov
    File summarize_script
    Float p_threshold = 0.00000005
    Int top_n = 50
    String p_column = "P_cct_admixed_c"
    String named_suffix = ".felix.tsv"
    String report_title = "FELIX chr22 QC summary"
    String docker = "us-central1-docker.pkg.dev/broad-dsp-lrma/aou-lr/felix-pilot:0.1.0"
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
      --top-n ~{top_n} \
      --p-threshold ~{p_threshold} \
      --p-column "~{p_column}" \
      --named-suffix "~{named_suffix}" \
      --report-title "~{report_title}"

    ls -lh summary/results_by_phenotype/ || true
    wc -l summary/results_manifest.tsv summary/calibration_summary.tsv
  >>>

  output {
    File results_manifest = "summary/results_manifest.tsv"
    File calibration_summary = "summary/calibration_summary.tsv"
    File calibration_summary_md = "summary/calibration_summary.md"
    File lambda_gc_wide = "summary/lambda_gc_wide.tsv"
    File phewas_genomewide_hits = "summary/phewas_genomewide_hits.tsv"
    Array[File] results_tsvs_named = glob("summary/results_by_phenotype/*.felix.tsv")
    Array[File] qq_plots = glob("summary/qc/*/qq_joint_acpass.png")
    Array[File] manhattan_plots = glob("summary/qc/*/manhattan_joint.png")
    Array[File] top_hits_tables = glob("summary/qc/*/top_hits.tsv")
  }

  runtime {
    docker: docker
    cpu: cpu
    memory: memory_gb + " GB"
    disks: "local-disk " + disk_gb + " HDD"
    preemptible: preemptible
  }
}

workflow FelixPilot {
  input {
    # Joint GT:AN1:AN2 VCF fills both felixla flags. Override with phase_vcf /
    # flare_vcf when phase and FLARE are still separate files.
    File? joint_vcf
    File? phase_vcf
    File? flare_vcf
    String chrom = "chr22"
    Int num_ancs = 5

    File analysis_samples
    File pheno_cov
    File selected_phenotypes
    File covariate_columns

    Array[File] grm_vcfs

    File build_saige_grm_script
    File make_plink_keep_script
    File fit_felix_null_script
    File run_felix_step2_script
    File summarize_script

    Float relatedness_cutoff = 0.05
    Int num_random_markers = 2000
    Int min_mac = 50
    Float pvalcutoff_of_haplotype = 0.05
    Float summarize_p_threshold = 0.00000005
    Int summarize_top_n = 50

    Int pack_cpu = 4
    Int pack_memory_gb = 16
    Int pack_disk_gb_floor = 100
    Float pack_disk_gb_multiplier = 3.5

    String docker = "us-central1-docker.pkg.dev/broad-dsp-lrma/aou-lr/felix-pilot:0.1.0"
  }

  Array[String] phenotypes = read_lines(selected_phenotypes)
  File phase_for_pack = select_first([phase_vcf, joint_vcf])
  File flare_for_pack = select_first([flare_vcf, joint_vcf])

  call CheckInputs as Check {
    input:
      chrom = chrom,
      n_phenotypes = length(phenotypes),
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

  call PackFelixla as Pack {
    input:
      phase_vcf = phase_for_pack,
      flare_vcf = flare_for_pack,
      analysis_samples = analysis_samples,
      chrom = Check.chrom_id,
      num_ancs = num_ancs,
      docker = docker,
      cpu = pack_cpu,
      memory_gb = pack_memory_gb,
      disk_gb_floor = pack_disk_gb_floor,
      disk_gb_multiplier = pack_disk_gb_multiplier
  }

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
        docker = docker
    }

    call RunFelixStep2 as Step2 {
      input:
        packed_tar = Pack.packed_tar,
        chrom = Check.chrom_id,
        phenotype = pheno,
        null_rda = Null.null_rda,
        variance_ratio = Null.variance_ratio,
        samples_used = Null.samples_used,
        sparse_grm_mtx = MakeGRM.sparse_grm_mtx,
        sparse_grm_sample_ids = MakeGRM.sparse_grm_sample_ids,
        run_step2_script = run_felix_step2_script,
        num_ancs = num_ancs,
        min_mac = min_mac,
        pvalcutoff_of_haplotype = pvalcutoff_of_haplotype,
        docker = docker
    }
  }

  scatter (i in range(length(phenotypes))) {
    call ConcatPhenotypeScores as Concat {
      input:
        phenotype = phenotypes[i],
        shard_tsvs = [Step2.results_tsv[i]],
        docker = docker
    }
  }

  call SummarizeFelixResults as Summarize {
    input:
      result_tsvs = Concat.merged_tsv,
      phenotypes = phenotypes,
      pheno_cov = pheno_cov,
      summarize_script = summarize_script,
      p_threshold = summarize_p_threshold,
      top_n = summarize_top_n,
      docker = docker
  }

  output {
    String chrom_id = Check.chrom_id
    File plink_bed = MakeGRM.plink_bed
    File plink_bim = MakeGRM.plink_bim
    File plink_fam = MakeGRM.plink_fam
    File sparse_grm_mtx = MakeGRM.sparse_grm_mtx
    File sparse_grm_sample_ids = MakeGRM.sparse_grm_sample_ids
    File packed_tar = Pack.packed_tar
    Array[File] null_rdas = Null.null_rda
    Array[File] null_metas = Null.null_meta
    Array[File] results_tsvs = Concat.merged_tsv
    Array[File] results_names = Concat.merged_name
    Array[File] results_tsvs_named = Summarize.results_tsvs_named
    File results_manifest = Summarize.results_manifest
    File calibration_summary = Summarize.calibration_summary
    File calibration_summary_md = Summarize.calibration_summary_md
    File lambda_gc_wide = Summarize.lambda_gc_wide
    File phewas_genomewide_hits = Summarize.phewas_genomewide_hits
    Array[File] qq_plots = Summarize.qq_plots
    Array[File] manhattan_plots = Summarize.manhattan_plots
    Array[File] top_hits_tables = Summarize.top_hits_tables
  }

  meta {
    description: "FELIX chr22 pilot: SAIGE mtx GRM, FELIX step1 null, felixla pack, admixed step2, P_cct_admixed_c summary."
    allowNestedInputs: true
  }
}
