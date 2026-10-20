# TRGT haplotype-phasing image: whatshap + TRGT + samtools/bcftools.
#
# Used by wdl/TrgtHaplotypePhase.wdl. whatshap haplotags reads against the
# phased callset so HP is in *our* haplotype frame, TRGT re-genotypes the
# repeats from those tagged reads (it has consumed HP during genotyping since
# 0.7.0 and written read-to-allele assignments into the spanning BAM since
# 0.5.0), and scripts/phase_trgt_by_haplotag.py joins the two.
#
# Build from repo root (linux/amd64):
#   docker build --platform linux/amd64 -f repeat_loci/docker/Dockerfile \
#     -t aou-trgt-phase:0.1.0 .
# Bump the WDL runtime.docker tag when rebuilding.

FROM mambaorg/micromamba:1.5.10

USER root

ARG TRGT_VERSION=1.5.1
ARG WHATSHAP_VERSION=2.3
ARG SAMTOOLS_VERSION=1.21
ARG BCFTOOLS_VERSION=1.21

ENV DEBIAN_FRONTEND=noninteractive \
    MAMBA_ROOT_PREFIX=/opt/conda \
    PATH=/opt/conda/bin:${PATH}

RUN micromamba install -y -n base -c conda-forge -c bioconda \
        python=3.11 \
        trgt=${TRGT_VERSION} \
        whatshap=${WHATSHAP_VERSION} \
        samtools=${SAMTOOLS_VERSION} \
        bcftools=${BCFTOOLS_VERSION} \
        htslib \
        tabix \
        pip \
    && micromamba clean --all --yes \
    && trgt --version \
    && whatshap --version \
    && samtools --version | head -1 \
    && bcftools --version | head -1 \
    && python --version \
    # The task strips HiPhase's HP/PS before re-tagging, so -x is load-bearing.
    # Fail the build rather than a 12,000-sample run if the flag ever moves.
    && samtools view --help 2>&1 | grep -q -- "--remove-tag" \
    && whatshap haplotag --help >/dev/null

WORKDIR /work
CMD ["/bin/bash"]
