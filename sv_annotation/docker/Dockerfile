FROM mambaorg/micromamba:1.5.10

USER root

ARG CADD_SV_VERSION=2.0.2
ARG BCFTOOLS_VERSION=1.20

ENV DEBIAN_FRONTEND=noninteractive \
    MAMBA_ROOT_PREFIX=/opt/conda \
    PATH=/opt/conda/bin:/opt/aou_sv/scripts:${PATH} \
    CADD_SV_CONDA_PREFIX=/caddsv-conda

RUN micromamba install -y -n base -c conda-forge -c bioconda \
        python=3.12 \
        bcftools=${BCFTOOLS_VERSION} \
        htslib \
        tabix \
        bedtools \
        samtools \
        pip \
        conda \
        pandas \
        matplotlib \
        pytest \
    && micromamba clean --all --yes \
    && pip install --no-cache-dir "caddsv==${CADD_SV_VERSION}" \
    && caddsv --help >/dev/null \
    && bcftools --version >/dev/null \
    && python --version \
    && mkdir -p /caddsv-conda /opt/aou_sv/scripts \
    && chmod 777 /caddsv-conda

COPY .docker_scripts/ /opt/aou_sv/scripts/

RUN chmod 755 /opt/aou_sv/scripts/* \
    && python -m py_compile /opt/aou_sv/scripts/*.py

WORKDIR /work
CMD ["/bin/bash"]
