version 1.0

# Fit one FELIX/SAIGE step-1 null model per phenotype, standalone, for reuse by
# every association run that uses the same cohort and covariate arm.
#
# Why this exists: FelixPilot and FelixGenome each fit their own nulls inside
# the run, so a scan over P phenotypes pays P null fits every time it is
# resubmitted, and a 900-phecode list makes that the dominant cost. The null
# depends only on the GRM, the covariates and the phenotype column - nothing
# from the LAI chain, the repeat branch or any chromosome - so it can be built
# once, here, while the rest of the pipeline is still being settled, and fed to
# FelixGenome through its precomputed_null_* inputs.
#
# The task below is copied VERBATIM from felix/wdl/FelixGenome.wdl so the nulls
# are identical in kind to the ones that workflow fits for itself. Do not edit
# it here without editing it there.
#
# Pair with grm/wdl/BuildSparseGRM.wdl: its plink_bed/bim/fam and
# sparse_grm_mtx/sample_ids outputs are this workflow's inputs.
#
# ORDER IS THE CONTRACT. The output arrays are in the order of
# selected_phenotypes, and FelixGenome indexes them against its own
# read_lines(selected_phenotypes). Feed both workflows the SAME file, and check
# null_manifest before using the arrays anywhere else.
#
# CACHE KEYS. A null fit is keyed on pheno_cov, covariate_columns,
# analysis_samples, the GRM files and fit_null_script. Switching covariate arms
# (limited vs full) invalidates every fit, so settle the arm before running
# this over a long phenotype list. fit_felix_null.R drops covariates that have
# no variation on a given phenotype's samples - a sex-restricted phecode makes
# `sex` constant - and records them in null_meta.

workflow FitFelixNulls {
  input {
    File pheno_cov
    File selected_phenotypes
    File covariate_columns
    File analysis_samples

    # From grm/wdl/BuildSparseGRM.wdl, or any run that built the same GRM.
    File plink_bed
    File plink_bim
    File plink_fam
    File sparse_grm_mtx
    File sparse_grm_sample_ids

    File fit_felix_null_script
    File? export_null_script
    String trait_type = "binary"

    String docker = "us-central1-docker.pkg.dev/broad-dsp-lrma/aou-lr/felix-pilot:0.2.0"
    Int null_cpu = 8
    Int null_memory_gb = 32
    Int null_disk_gb = 100
    Int null_preemptible = 0
  }

  Array[String] phenotypes = read_lines(selected_phenotypes)

  scatter (pheno in phenotypes) {
    call FitFelixNull as Null {
      input:
        pheno_cov = pheno_cov,
        phenotype = pheno,
        covariate_columns = covariate_columns,
        analysis_samples = analysis_samples,
        plink_bed = plink_bed,
        plink_bim = plink_bim,
        plink_fam = plink_fam,
        sparse_grm_mtx = sparse_grm_mtx,
        sparse_grm_sample_ids = sparse_grm_sample_ids,
        fit_null_script = fit_felix_null_script,
        export_null_script = export_null_script,
        trait_type = trait_type,
        docker = docker,
        cpu = null_cpu,
        memory_gb = null_memory_gb,
        disk_gb = null_disk_gb,
        preemptible = null_preemptible
    }
  }

  call CollectNullManifest as Manifest {
    input:
      phenotypes = phenotypes,
      null_rdas = Null.null_rda,
      variance_ratios = Null.variance_ratio,
      samples_used = Null.samples_used,
      null_metas = Null.null_meta,
      docker = docker
  }

  output {
    Array[File] null_rdas = Null.null_rda
    Array[File] variance_ratios = Null.variance_ratio
    Array[File] samples_used = Null.samples_used
    Array[File] null_metas = Null.null_meta
    Array[File] null_export_tars = Null.null_export_tar
    Array[File] pheno_used = Null.pheno_used
    File null_manifest = Manifest.manifest_tsv
    File null_summary = Manifest.summary_tsv
  }

  meta {
    description: "Standalone FELIX step-1 null fits, one per phenotype, for reuse across association runs."
    allowNestedInputs: true
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

# One row per phenotype in input order, plus the case/control counts and the
# covariates each fit actually used, so a later run can be checked against the
# list it was built from rather than trusting array order blindly.
task CollectNullManifest {
  input {
    Array[String] phenotypes
    Array[File] null_rdas
    Array[File] variance_ratios
    Array[File] samples_used
    Array[File] null_metas
    String docker = "us-central1-docker.pkg.dev/broad-dsp-lrma/aou-lr/felix-pilot:0.2.0"
    Int cpu = 1
    Int memory_gb = 4
    Int disk_gb = 20
    Int preemptible = 1
  }

  command <<<
    set -euo pipefail
    cat > phenotypes.txt <<'EOF'
~{sep="\n" phenotypes}
EOF
    cat > null_rdas.txt <<'EOF'
~{sep="\n" null_rdas}
EOF
    cat > variance_ratios.txt <<'EOF'
~{sep="\n" variance_ratios}
EOF
    cat > samples_used.txt <<'EOF'
~{sep="\n" samples_used}
EOF
    cat > null_metas.txt <<'EOF'
~{sep="\n" null_metas}
EOF

    python3 <<'PYEOF'
import csv, os

def lines(path):
    with open(path) as fh:
        return [ln.strip() for ln in fh if ln.strip()]

pheno = lines("phenotypes.txt")
cols = {name: lines(f"{name}.txt")
        for name in ("null_rdas", "variance_ratios", "samples_used", "null_metas")}
for name, vals in cols.items():
    if len(vals) != len(pheno):
        raise SystemExit(f"{name}: {len(vals)} files for {len(pheno)} phenotypes")

with open("null_manifest.tsv", "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t", lineterminator="\n")
    w.writerow(["index", "phenotype", "null_rda", "variance_ratio", "samples_used"])
    for i, p in enumerate(pheno):
        w.writerow([i, p] + [os.path.basename(cols[c][i])
                             for c in ("null_rdas", "variance_ratios", "samples_used")])

rows = []
for i, p in enumerate(pheno):
    with open(cols["null_metas"][i]) as fh:
        meta = next(csv.DictReader(fh, delimiter="\t"), {})
    meta["index"] = i
    meta.setdefault("phenotype", p)
    if meta["phenotype"] != p:
        raise SystemExit(
            f"row {i}: manifest says {p}, null_meta says {meta['phenotype']}"
        )
    rows.append(meta)

fields = ["index"] + [k for k in rows[0] if k != "index"] if rows else ["index"]
with open("null_summary.tsv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t",
                       lineterminator="\n", extrasaction="ignore")
    w.writeheader()
    w.writerows(rows)
print(f"{len(pheno)} null fits")
PYEOF

    head -n 5 null_summary.tsv
  >>>

  output {
    File manifest_tsv = "null_manifest.tsv"
    File summary_tsv = "null_summary.tsv"
  }

  runtime {
    docker: docker
    cpu: cpu
    memory: memory_gb + " GB"
    disks: "local-disk " + disk_gb + " HDD"
    preemptible: preemptible
  }
}
