# FELIX Terra image.
# Consume the published FELIX image (do not compile or fork FELIX). Add PLINK2,
# bcftools, extract-tracts-flare, tractor-mix-score, and repo helper scripts.
#
# Build from repo root (linux/amd64):
#   docker build --platform linux/amd64 -f felix/docker/Dockerfile -t felix-pilot:latest .
# Bump WDL runtime.docker tags when rebuilding.

FROM rust:1.85-alpine AS rust-builder
WORKDIR /src
RUN apk add --no-cache musl-dev

COPY tractor_mix/extract_tracts_flare/Cargo.toml tractor_mix/extract_tracts_flare/Cargo.lock /src/extract/
COPY tractor_mix/extract_tracts_flare/src /src/extract/src
RUN cd extract && cargo build --release --locked --bin extract-tracts-flare \
    && test -x target/release/extract-tracts-flare

COPY tractor_mix/tractor_mix_score/Cargo.toml tractor_mix/tractor_mix_score/Cargo.lock /src/score/
COPY tractor_mix/tractor_mix_score/src /src/score/src
RUN cd score && cargo build --release --locked --bin tractor-mix-score \
    && test -x target/release/tractor-mix-score

FROM lhu1/felix:latest

ARG PLINK2_URL=https://s3.amazonaws.com/plink2-assets/alpha7/plink2_linux_x86_64_20260808.zip
ARG BCFTOOLS_VERSION=1.21

ENV DEBIAN_FRONTEND=noninteractive
ENV PATH="/opt/plink2:/opt/bcftools/bin:/opt/felix_scripts:${PATH}"

RUN mkdir -p /opt/SAIGE/extdata \
    && ln -sf /usr/local/bin/step1_fitNULLGLMM.R /opt/SAIGE/extdata/step1_fitNULLGLMM.R \
    && ln -sf /usr/local/bin/step2_SPAtests.R /opt/SAIGE/extdata/step2_SPAtests.R \
    && ln -sf /usr/local/bin/createSparseGRM.R /opt/SAIGE/extdata/createSparseGRM.R \
    && test -x /usr/local/bin/step1_fitNULLGLMM.R \
    && test -x /usr/local/bin/step2_SPAtests.R \
    && test -x /usr/local/bin/createSparseGRM.R \
    && command -v felixla >/dev/null

RUN apt-get update && apt-get install -y --no-install-recommends \
        build-essential \
        ca-certificates \
        curl \
        gzip \
        libbz2-dev \
        libcurl4-openssl-dev \
        liblzma-dev \
        python3 \
        unzip \
        wget \
        zlib1g-dev \
    && rm -rf /var/lib/apt/lists/*

RUN mkdir -p /opt/plink2 \
    && cd /tmp \
    && wget -O plink2.zip "${PLINK2_URL}" \
    && unzip -o plink2.zip -d /opt/plink2 \
    && chmod +x /opt/plink2/plink2 \
    && /opt/plink2/plink2 --version \
    && rm -f plink2.zip

RUN cd /tmp \
    && wget -q "https://github.com/samtools/bcftools/releases/download/${BCFTOOLS_VERSION}/bcftools-${BCFTOOLS_VERSION}.tar.bz2" \
    && tar -xjf "bcftools-${BCFTOOLS_VERSION}.tar.bz2" \
    && cd "bcftools-${BCFTOOLS_VERSION}" \
    && ./configure --prefix=/opt/bcftools \
    && make -j"$(nproc)" \
    && make install \
    && cd / \
    && rm -rf "/tmp/bcftools-${BCFTOOLS_VERSION}" "/tmp/bcftools-${BCFTOOLS_VERSION}.tar.bz2"

# lhu1/felix:latest puts pixi Python first on PATH (/app/.pixi/envs/default/bin/python3).
# apt python3-pip targets system 3.8; workflow tasks invoke `python3`, so install into pixi.
RUN /app/.pixi/envs/default/bin/python3 -m ensurepip \
    && /app/.pixi/envs/default/bin/python3 -m pip install --no-cache-dir \
        pandas numpy matplotlib seaborn

# Staged by felix/build_docker.sh into felix/.docker_scripts/
COPY felix/.docker_scripts/ /opt/felix_scripts/
COPY --from=rust-builder /src/extract/target/release/extract-tracts-flare /usr/local/bin/extract-tracts-flare
COPY --from=rust-builder /src/score/target/release/tractor-mix-score /usr/local/bin/tractor-mix-score
RUN chmod 0755 /usr/local/bin/extract-tracts-flare /usr/local/bin/tractor-mix-score \
    && /usr/local/bin/extract-tracts-flare --help >/dev/null \
    && /usr/local/bin/tractor-mix-score --help >/dev/null
ENV PATH="/opt/felix_scripts:${PATH}"

RUN felixla --help >/dev/null \
    && step1_fitNULLGLMM.R --help >/dev/null \
    && step2_SPAtests.R --help >/dev/null \
    && createSparseGRM.R --help >/dev/null \
    && plink2 --version \
    && bcftools --version \
    && python3 -c "import pandas, numpy, matplotlib"

WORKDIR /data
CMD ["/bin/bash"]
