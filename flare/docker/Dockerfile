# bcftools/htslib for FlareByPopulation GCS streams.
#
# aou-sv-annotation:0.1.6 pins bcftools 1.20, which has no libcurl retry.
# HTSlib 1.23+ reconnects mid-transfer on curl 92 (HTTP/2 GOAWAY) — see
# https://github.com/samtools/htslib/pull/1987
# Do not bump the SV image for this; 1.22+ changed VCF GT phasing storage.

FROM mambaorg/micromamba:1.5.10

USER root

ARG BCFTOOLS_VERSION=1.24

ENV DEBIAN_FRONTEND=noninteractive \
    MAMBA_ROOT_PREFIX=/opt/conda \
    PATH=/opt/conda/bin:${PATH} \
    HTS_RETRY_MAX=8 \
    HTS_RETRY_DELAY=500 \
    HTS_RETRY_MAX_DELAY=60000

RUN micromamba install -y -n base -c conda-forge -c bioconda \
        python=3.12 \
        bcftools=${BCFTOOLS_VERSION} \
        htslib=${BCFTOOLS_VERSION} \
        tabix \
    && micromamba clean --all --yes \
    && bcftools --version \
    && python3 --version

WORKDIR /work
CMD ["/bin/bash"]
