version 1.0

workflow SubsetAndSplitVcf {

    input {
        File joint_vcf
        File joint_vcf_idx
        String region
        String gcs_output_dir
        String output_tag           # per-sample VCFs land at gcs_output_dir/{sample}.{output_tag}.vcf.gz

        # optional hierarchical split
        File? sample_batches_tsv      # see bcftools +split --help; e.g., contains rows: {sample_name_1},{sample_name_2},...\t-\tbatch-0
    }

    call SubsetAndSplit { input:
        vcf = joint_vcf,
        vcf_idx = joint_vcf_idx,
        region = region,
        gcs_output_dir = if defined(sample_batches_tsv) then gcs_output_dir + "/batches" else gcs_output_dir,
        sample_batches_tsv = sample_batches_tsv,
        output_tag = output_tag
    }
    
    if (defined(sample_batches_tsv)) {
        Int num_batches = length(read_lines(select_first([sample_batches_tsv])))
        scatter (i in range(num_batches)) {
            call SubsetAndSplit as SubsetAndSplitBatch { input:
                vcf = gcs_output_dir + "/batches/batch-" + i + "." + output_tag + ".vcf.gz",
                gcs_output_dir = gcs_output_dir,
                wait_for_me = SubsetAndSplit.number_of_sites,
                output_tag = output_tag
            }
        }
    }

    output {
        Int number_of_sites = SubsetAndSplit.number_of_sites
        String phased_fraction = SubsetAndSplit.phased_fraction

        # Deliberately a COUNT, not the paths. Terra saves every workflow output
        # as an entity attribute and caps the count at 20,000; returning one
        # string per sample per chromosome exceeded that and failed the workflow
        # after every task had succeeded. The paths are reconstructible anyway --
        # <gcs_output_dir>/<sample>.<output_tag>.vcf.gz is fixed by construction,
        # which is what make_trgt_phase_table.py --split-base relies on, so the
        # list was never load-bearing. The per-task manifests remain as files.
        Int number_of_split_files = length(select_first([flatten(select_first([SubsetAndSplitBatch.split_vcf_paths])),
                                                         SubsetAndSplit.split_vcf_paths]))
        Array[File]? split_vcf_manifests = SubsetAndSplitBatch.split_vcf_paths_file
    }
}

struct RuntimeAttr {
    Float? mem_gb
    Int? cpu_cores
    Int? disk_gb
    Int? boot_disk_gb
    Boolean? use_ssd
    Int? preemptible_tries
    Int? max_retries
    String? docker
}

task SubsetAndSplit {

    input {
        File vcf
        File? vcf_idx       # only needed for initial stream, not for batches
        String? region      # only needed for initial stream, not for batches
        String gcs_output_dir
        String output_tag
        File? sample_batches_tsv
        Int? wait_for_me    # dummy input to gate scatter until initial split is done
        Int view_verbosity = 3
        Int index_jobs = 8
        Int phase_check_sites = 10000
        Float min_phased_fraction = 0.5
        RuntimeAttr? runtime_attr_override
    }

    parameter_meta {
        # NOT localization_optional any more. Streaming through htslib broke on
        # chr14 with "Libcurl reported error 92 (Stream error in the HTTP/2
        # framing layer)" after 649,510 of its records -- the GOAWAY/libcurl
        # flakiness this task's own comment already warns about, which is why it
        # writes an intermediate rather than piping into +split.
        #
        # localization_optional was there to avoid copying a huge joint callset,
        # but these per-chromosome files are about 1.6-2 GiB (chr20 measured at
        # 1,696,914,957 bytes), so Cromwell's own localization -- which retries
        # and verifies -- costs roughly a minute and removes the failure mode
        # entirely. Restore the attribute if a future input is large enough that
        # copying it is the problem again.
        vcf: {
            description: "Localized, not streamed; see above."
        }
    }

    String view_input_arg = if defined(vcf_idx) then "\"~{vcf}##idx##~{vcf_idx}\"" else "~{vcf}"

    command <<<
        set -euxo pipefail

        # see https://github.com/samtools/htslib/issues/803#issuecomment-444514336, https://github.com/broadinstitute/bcftools-patched
        mkfifo /tmp/token_fifo
        ( while true ; do curl -s -H "Metadata-Flavor: Google" http://metadata.google.internal/computeMetadata/v1/instance/service-accounts/default/token > /tmp/token_fifo ; done ) &
        export HTS_AUTH_LOCATION="/tmp/token_fifo"

        # we stream to an intermediate file, since piping directly to bcftools +split still results in GOAWAY/Libcurl issues
        bcftools view --no-version ~{view_input_arg} \
            ~{"--regions " + region} \
            --regions-overlap 0 \
            --verbosity ~{view_verbosity} \
            -Ob | \
        bcftools norm --rm-dup all -Ob -o ~{output_tag}.bcf

        # check number of sites to guard against streaming errors
        bcftools view -H ~{output_tag}.bcf | wc -l > number_of_sites.txt

        # The SHAPEIT4 phase is the entire point of these files: TrgtHaplotypePhase
        # reads haplotype order out of them. A split that dropped the '|'
        # separators would yield files that look valid while silently destroying
        # every orientation downstream. Measured on the intermediate, before
        # 12,000 files are written from it.
        # The subshell disables pipefail for this one pipeline, deliberately.
        # `head` closes the pipe after N lines, bcftools takes SIGPIPE and exits
        # 141, and under `set -o pipefail` that fails the whole task -- which is
        # exactly how the first chr20 submission died. Reading all 3M records
        # just to sample 10k would be the alternative.
        ( set +o pipefail
          bcftools query -f '[%GT]\n' ~{output_tag}.bcf \
            | head -n ~{phase_check_sites} > gt_sample.txt )
        N_GT=$(wc -l < gt_sample.txt | tr -d ' ')
        N_PHASED=$(grep -c '|' gt_sample.txt || true)
        echo "phase check: ${N_PHASED} of ${N_GT} sampled genotypes carry '|'" >&2
        if [[ "${N_GT}" -eq 0 ]]; then
            echo "0" > phased_fraction.txt
            echo "error: no genotypes read from the intermediate; refusing to split." >&2
            exit 1
        fi
        FRAC=$(awk -v a="${N_PHASED}" -v b="${N_GT}" 'BEGIN{printf "%.4f", a/b}')
        echo "${FRAC}" > phased_fraction.txt
        if awk -v f="${FRAC}" -v m=~{min_phased_fraction} 'BEGIN{exit !(f < m)}'; then
            echo "error: only ${FRAC} of sampled genotypes are phased, below" >&2
            echo "       min_phased_fraction ~{min_phased_fraction}. These files feed" >&2
            echo "       haplotype-resolved analysis; unphased input would be used as" >&2
            echo "       though it were phased." >&2
            exit 1
        fi

        # bgzipped VCF, not BCF: whatshap reads these directly and VCF is the
        # format it documents, and .tbi pairs with the ##idx## syntax downstream.
        bcftools +split ~{output_tag}.bcf ~{"--samples-file " + sample_batches_tsv} -Oz -o output

        for f in output/*.vcf.gz; do
            base=$(basename "$f" .vcf.gz)
            mv "$f" "output/${base}.~{output_tag}.vcf.gz"
        done

        # Index every output: a region read downstream and whatshap --regions
        # both need random access, and a missing index fails inside whatshap
        # rather than here. Parallel -- at cohort scale this is 12,261 calls.
        find output -name '*.vcf.gz' -print0 \
            | xargs -0 -P ~{index_jobs} -I@@ bcftools index -t -f @@

        N_OUT=$(find output -name '*.vcf.gz' | wc -l | tr -d ' ')
        N_IDX=$(find output -name '*.vcf.gz.tbi' | wc -l | tr -d ' ')
        echo "wrote ${N_OUT} per-sample VCFs and ${N_IDX} indexes" >&2
        if [[ "${N_OUT}" -ne "${N_IDX}" ]]; then
            echo "error: ${N_OUT} VCFs but ${N_IDX} indexes." >&2
            exit 1
        fi

        gcloud storage cp output/*.vcf.gz output/*.vcf.gz.tbi ~{gcs_output_dir}/

        # List THIS task's own outputs, from the local directory. `gsutil ls` on
        # the destination returns everything any batch has written so far, so
        # each of the 25 batch tasks reported a growing share of all 12,261
        # files -- which is how the workflow came to claim 225,092 paths.
        ( cd output && ls -1 ./*.vcf.gz ) \
            | sed "s|^\./|~{gcs_output_dir}/|" > output_vcf_paths.txt
        wc -l < output_vcf_paths.txt | tr -d ' ' > n_split_files.txt
        echo "this task contributed $(cat n_split_files.txt) per-sample VCFs" >&2
    >>>

    output {
        Array[String] split_vcf_paths = read_lines("output_vcf_paths.txt")
        Int n_split_files = read_int("n_split_files.txt")
        File split_vcf_paths_file = "output_vcf_paths.txt"
        Int number_of_sites = read_int("number_of_sites.txt")
        String phased_fraction = read_string("phased_fraction.txt")
    }
    ###################
    RuntimeAttr default_attr = object {
        cpu_cores:          8,
        mem_gb:             8,
        disk_gb:            20 + 4*ceil(size(vcf, "GB")),
        boot_disk_gb:       10,
        use_ssd:            true,
        preemptible_tries:  5,
        max_retries:        2,
        docker:             "us.gcr.io/broad-dsde-methods/slee/lr-bcftools-patched-gcloud/lr-bcftools-patched-gcloud:1.23"      # see https://github.com/broadinstitute/bcftools-patched
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + if select_first([runtime_attr.use_ssd, default_attr.use_ssd]) then " SSD" else " HDD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}