version 1.0

# Per-sample ExpansionHunter on a Nearline-efficient minicram.
# Extract: broadinstitute/str-analysis make_minicram_for_expansion_hunter
#   (catalog ReferenceRegion + OfftargetRegions + mate pass).
# Binary: bw2/ExpansionHunter (optimized-streaming after the minicram).
# Submit: notebooks/rw/expansion_hunter_01_run.ipynb via `wb workflow`.

workflow ExpansionHunterMinicram {
    input {
        String sample_id

        File cram
        File crai

        File ref_fa
        File ref_fai

        File catalog

        # ExpansionHunter --sex: male | female (m/f/1/2 also accepted).
        String sex = "female"

        # Billing project for requester-pays GCS (AoU CRAMs). Empty → GOOGLE_CLOUD_PROJECT if set.
        String gcloud_project = ""

        Int window_size = 1000
        Int merge_regions_distance = 1000

        Int max_retry = 3
        Int wait_time = 30

        Int minicram_disk_gb = 50
        Int minicram_mem_gb = 16
        Int minicram_n_preemptible = 0

        Int eh_n_cpu = 4
        Int eh_mem_gb = 8
        Int eh_n_preemptible = 2

        String print_reads_docker = "us-central1-docker.pkg.dev/broad-dsp-lrma/aou-lr/aou-locityper-print-reads:0.1.0"
        String eh_docker = "us-central1-docker.pkg.dev/broad-dsp-lrma/aou-lr/aou-expansion-hunter:0.1.0"
    }

    call MakeMinicram {
        input:
            sample_id = sample_id,
            cram = cram,
            crai = crai,
            reference = ref_fa,
            reference_index = ref_fai,
            catalog = catalog,
            gcloud_project = gcloud_project,
            window_size = window_size,
            merge_regions_distance = merge_regions_distance,
            max_retry = max_retry,
            wait_time = wait_time,
            disk_gb = minicram_disk_gb,
            mem_gb = minicram_mem_gb,
            n_preemptible = minicram_n_preemptible,
            docker = print_reads_docker
    }

    call ExpansionHunterGenotype {
        input:
            sample_id = sample_id,
            minicram = MakeMinicram.minicram,
            minicrai = MakeMinicram.minicrai,
            reference = ref_fa,
            reference_index = ref_fai,
            catalog = catalog,
            sex = sex,
            n_cpu = eh_n_cpu,
            mem_gb = eh_mem_gb,
            n_preemptible = eh_n_preemptible,
            docker = eh_docker
    }

    output {
        File eh_json = ExpansionHunterGenotype.eh_json
        File eh_vcf = ExpansionHunterGenotype.eh_vcf
        File minicram = MakeMinicram.minicram
        File minicrai = MakeMinicram.minicrai
        File minicram_transfer_stats = MakeMinicram.transfer_stats
    }
}

task MakeMinicram {
    input {
        String sample_id
        File cram
        File crai
        File reference
        File reference_index
        File catalog
        String gcloud_project
        Int window_size
        Int merge_regions_distance
        Int max_retry
        Int wait_time
        Int disk_gb
        Int mem_gb
        Int n_preemptible
        String docker
    }

    # Do not size disk from the WGS CRAM. localization_optional leaves it in Nearline;
    # Cromwell size() would still see the full object and over-allocate.
    parameter_meta {
        cram: { localization_optional: true }
    }

    String out_cram = sample_id + ".minicram.cram"

    command <<<
        set -euo pipefail
        date

        python3 -c "import json; json.load(open('~{catalog}'))" \
            || { echo "catalog is not valid JSON: ~{catalog}" >&2; exit 1; }

        # Cromwell hashes each File into its own directory; htslib wants FASTA + fai together.
        ln -sf ~{reference} reference.fa
        ln -sf ~{reference_index} reference.fa.fai

        PROJ="~{gcloud_project}"
        if [ -z "${PROJ}" ]; then
            PROJ="${GOOGLE_CLOUD_PROJECT:-}"
        fi
        if [ -n "${PROJ}" ]; then
            echo "requester-pays project: ${PROJ}"
        fi

        run_minicram() {
            local -a cmd=(
                python3 -u -m str_analysis.make_minicram_for_expansion_hunter
                -R reference.fa
                -c ~{catalog}
                -i ~{crai}
                -o ~{out_cram}
                -w ~{window_size}
                -d ~{merge_regions_distance}
                --verbose
                --output-data-transfer-stats
            )
            if [ -n "${PROJ}" ]; then
                cmd+=(--gcloud-project "${PROJ}")
            fi
            cmd+=(~{cram})
            "${cmd[@]}"
        }

        MAX_RETRIES=~{max_retry}
        CURR_RETRIES=0
        WAIT_TIME=~{wait_time}
        until (( CURR_RETRIES == MAX_RETRIES )); do
            if run_minicram; then
                break
            fi
            echo "make_minicram_for_expansion_hunter failed; retry $(( CURR_RETRIES + 1 ))/${MAX_RETRIES}" >&2
            rm -f ~{out_cram} ~{out_cram}.crai *.data_transfer_stats.tsv
            CURR_RETRIES=$(( CURR_RETRIES + 1 ))
            if (( CURR_RETRIES == MAX_RETRIES )); then
                echo "make_minicram_for_expansion_hunter failed after ${MAX_RETRIES} attempts" >&2
                exit 1
            fi
            sleep "${WAIT_TIME}"
        done

        if [ ! -s ~{out_cram} ]; then
            echo "make_minicram wrote no CRAM" >&2
            exit 1
        fi
        if [ ! -s ~{out_cram}.crai ]; then
            echo "make_minicram did not write ~{out_cram}.crai" >&2
            exit 1
        fi

        ls -lh ~{out_cram} ~{out_cram}.crai *.data_transfer_stats.tsv
        date
    >>>

    output {
        File minicram = out_cram
        File minicrai = out_cram + ".crai"
        File transfer_stats = glob("*.data_transfer_stats.tsv")[0]
    }

    runtime {
        memory: mem_gb + " GB"
        cpu: 4
        disks: "local-disk " + disk_gb + " HDD"
        preemptible: n_preemptible
        docker: docker
    }
}

task ExpansionHunterGenotype {
    input {
        String sample_id
        File minicram
        File minicrai
        File reference
        File reference_index
        File catalog
        String sex
        Int n_cpu
        Int mem_gb
        Int n_preemptible
        String docker
    }

    Int disk_gb = 10 + 3 * ceil(size([minicram, minicrai, reference, reference_index, catalog], "GiB"))
    String prefix = sample_id + ".EH"

    command <<<
        set -euo pipefail
        date

        command -v ExpansionHunter
        ExpansionHunter --version 2>&1 || true
        ExpansionHunter --help > eh.help.txt 2>&1 || true

        sex_lc=$(printf '%s' "~{sex}" | tr '[:upper:]' '[:lower:]' | tr -d ' ')
        case "${sex_lc}" in
            male|m|1)   sex_flag="male"   ;;
            female|f|2) sex_flag="female" ;;
            *)
                echo "WARN: unrecognized sex '~{sex}'; defaulting to female" >&2
                sex_flag="female"
                ;;
        esac
        echo "sex=${sex_flag}"

        ln -sf ~{reference} reference.fa
        ln -sf ~{reference_index} reference.fa.fai
        ln -sf ~{minicram} reads.cram
        ln -sf ~{minicrai} reads.cram.crai

        if grep -q -- '--variant-catalog' eh.help.txt; then
            CAT_FLAG="--variant-catalog"
        else
            CAT_FLAG="--catalog"
        fi

        if grep -q -- 'optimized-streaming' eh.help.txt; then
            MODE="optimized-streaming"
        else
            MODE="seeking"
        fi
        echo "catalog flag=${CAT_FLAG} analysis-mode=${MODE}"

        # Minicram is local; do not stream the WGS CRAM. optimized-streaming
        # (bw2 fork) walks the extracted containers once; seeking is the fallback.
        extra=()
        grep -q -- '--dont-output-consensus-sequences' eh.help.txt \
            && extra+=(--dont-output-consensus-sequences)
        grep -q -- '--dont-output-quality-metrics' eh.help.txt \
            && extra+=(--dont-output-quality-metrics)
        grep -q -- '--disable-all-plots' eh.help.txt \
            && extra+=(--disable-all-plots)

        ExpansionHunter \
            --reads reads.cram \
            --reads-index reads.cram.crai \
            --reference reference.fa \
            "${CAT_FLAG}" ~{catalog} \
            --output-prefix ~{prefix} \
            --sex "${sex_flag}" \
            --analysis-mode "${MODE}" \
            "${extra[@]}"

        if [ ! -s ~{prefix}.json ]; then
            echo "ExpansionHunter wrote no ~{prefix}.json" >&2
            ls -la
            exit 1
        fi
        if [ ! -s ~{prefix}.vcf ]; then
            echo "ExpansionHunter wrote no ~{prefix}.vcf" >&2
            ls -la
            exit 1
        fi

        python3 -c "import json; json.load(open('~{prefix}.json'))" \
            || { echo "~{prefix}.json is not valid JSON" >&2; exit 1; }

        ls -lh ~{prefix}.json ~{prefix}.vcf
        date
    >>>

    output {
        File eh_json = prefix + ".json"
        File eh_vcf = prefix + ".vcf"
    }

    runtime {
        memory: mem_gb + " GB"
        cpu: n_cpu
        disks: "local-disk " + disk_gb + " HDD"
        preemptible: n_preemptible
        docker: docker
    }
}
