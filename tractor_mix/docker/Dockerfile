# Tractor-Mix + SAIGE Terra pilot image.
# Base on the official SAIGE image (precompiled) to avoid brittle SAIGE source
# builds against newer oneTBB / BH on rocker. Then add PLINK2, bcftools,
# Tractor / Tractor-Mix, GMMAT, and the Rust FLARE tract extractor.
#
# Bump WDL runtime.docker tags when rebuilding.

FROM rust:1.85-alpine AS rust-builder
WORKDIR /src
RUN apk add --no-cache musl-dev

COPY extract_tracts_flare/Cargo.toml extract_tracts_flare/Cargo.lock /src/extract/
COPY extract_tracts_flare/src /src/extract/src
RUN cd extract && cargo build --release --locked --bin extract-tracts-flare \
    && test -x target/release/extract-tracts-flare

COPY tractor_mix_score/Cargo.toml tractor_mix_score/Cargo.lock /src/score/
COPY tractor_mix_score/src /src/score/src
RUN cd score && cargo build --release --locked --bin tractor-mix-score \
    && test -x target/release/tractor-mix-score

FROM wzhou88/saige:1.3.3

ARG TRACTOR_SHA=2860728d22dbf37a816afcf9d41ef42f15facce8
ARG TRACTOR_MIX_SHA=4adb8f1814d9315ecd7868eb729d52ec0c723719
ARG PLINK2_URL=https://s3.amazonaws.com/plink2-assets/alpha7/plink2_linux_x86_64_20260808.zip
ARG BCFTOOLS_VERSION=1.21

ENV DEBIAN_FRONTEND=noninteractive
ENV PATH="/opt/plink2:/opt/bcftools/bin:/opt/Tractor/scripts:/opt/Tractor-Mix:/opt/SAIGE/extdata:/opt/tractor_mix_scripts:${PATH}"

# SAIGE scripts live in /usr/local/bin in the official image; keep the
# /opt/SAIGE/extdata paths our helper scripts expect.
RUN mkdir -p /opt/SAIGE/extdata \
    && ln -sf /usr/local/bin/step1_fitNULLGLMM.R /opt/SAIGE/extdata/step1_fitNULLGLMM.R \
    && ln -sf /usr/local/bin/step2_SPAtests.R /opt/SAIGE/extdata/step2_SPAtests.R \
    && ln -sf /usr/local/bin/createSparseGRM.R /opt/SAIGE/extdata/createSparseGRM.R \
    && Rscript -e 'library(SAIGE); cat("SAIGE OK from official image\n")'

RUN apt-get update && apt-get install -y --no-install-recommends \
        build-essential \
        ca-certificates \
        cmake \
        curl \
        gfortran \
        git \
        gzip \
        libbz2-dev \
        libcurl4-openssl-dev \
        liblzma-dev \
        libssl-dev \
        libxml2-dev \
        python3 \
        python3-pip \
        unzip \
        wget \
        zlib1g-dev \
    && rm -rf /var/lib/apt/lists/*

# PLINK 2
RUN mkdir -p /opt/plink2 \
    && cd /tmp \
    && wget -O plink2.zip "${PLINK2_URL}" \
    && unzip -o plink2.zip -d /opt/plink2 \
    && chmod +x /opt/plink2/plink2 \
    && /opt/plink2/plink2 --version \
    && rm -f plink2.zip

# bcftools
RUN cd /tmp \
    && wget -q "https://github.com/samtools/bcftools/releases/download/${BCFTOOLS_VERSION}/bcftools-${BCFTOOLS_VERSION}.tar.bz2" \
    && tar -xjf "bcftools-${BCFTOOLS_VERSION}.tar.bz2" \
    && cd "bcftools-${BCFTOOLS_VERSION}" \
    && ./configure --prefix=/opt/bcftools \
    && make -j"$(nproc)" \
    && make install \
    && cd / \
    && rm -rf "/tmp/bcftools-${BCFTOOLS_VERSION}" "/tmp/bcftools-${BCFTOOLS_VERSION}.tar.bz2"

# Tractor + Tractor-Mix
RUN git clone --filter=blob:none https://github.com/Atkinson-Lab/Tractor.git /opt/Tractor \
    && cd /opt/Tractor && git checkout "${TRACTOR_SHA}" \
    && git clone --filter=blob:none https://github.com/Atkinson-Lab/Tractor-Mix.git /opt/Tractor-Mix \
    && cd /opt/Tractor-Mix && git checkout "${TRACTOR_MIX_SHA}" \
    && chmod +x /opt/Tractor/scripts/*.py

# GMMAT + helpers for Tractor-Mix (R 3.6 on wzhou88/saige:1.3.3).
# Bioconductor deps required by current GMMAT; pin Bioc 3.10 for R 3.6.
RUN Rscript -e 'install.packages(c("BiocManager","remotes","data.table","foreach","doParallel","abind","Rcpp","RcppArmadillo","CompQuadForm","optparse","jsonlite"), repos="https://cloud.r-project.org")' \
    && Rscript -e 'BiocManager::install(version="3.10", ask=FALSE, update=FALSE); BiocManager::install(c("SeqArray","SeqVarTools"), ask=FALSE, update=FALSE)' \
    && Rscript -e 'install.packages("GMMAT", repos="https://cloud.r-project.org"); library(GMMAT); cat("GMMAT OK\n")'

# Python deps compatible with the SAIGE image Python (~3.8)
RUN pip3 install --no-cache-dir \
        'pandas==1.5.3' \
        'numpy==1.24.4' \
        'matplotlib==3.7.3' \
        'seaborn==0.12.2'

# Build from tractor_mix/ after build_docker.sh stages repo-root scripts/ into
# .docker_scripts/:
#   docker build -f docker/Dockerfile -t tractor-mix-pilot:latest .
COPY .docker_scripts/ /opt/tractor_mix_scripts/
COPY --from=rust-builder /src/extract/target/release/extract-tracts-flare /usr/local/bin/extract-tracts-flare
COPY --from=rust-builder /src/score/target/release/tractor-mix-score /usr/local/bin/tractor-mix-score
RUN chmod 0755 /usr/local/bin/extract-tracts-flare /usr/local/bin/tractor-mix-score \
    && /usr/local/bin/extract-tracts-flare --help >/dev/null \
    && /usr/local/bin/tractor-mix-score --help >/dev/null
ENV PATH="/opt/tractor_mix_scripts:${PATH}"

WORKDIR /data
CMD ["/bin/bash"]
