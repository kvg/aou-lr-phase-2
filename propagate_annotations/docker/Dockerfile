# Propagate SNV/indel, SV, and FLARE annotations onto a base VCF.

FROM rust:1.85-alpine AS rust-builder
WORKDIR /src
RUN apk add --no-cache musl-dev
COPY Cargo.toml Cargo.lock /src/
COPY src /src/src
RUN cargo build --release --locked --bin propagate-annotations \
    && test -x target/release/propagate-annotations

FROM mambaorg/micromamba:1.5.10

USER root

ARG BCFTOOLS_VERSION=1.21

ENV DEBIAN_FRONTEND=noninteractive \
    MAMBA_ROOT_PREFIX=/opt/conda \
    PATH=/opt/conda/bin:${PATH}

RUN micromamba install -y -n base -c conda-forge -c bioconda \
        bcftools=${BCFTOOLS_VERSION} \
        htslib \
        tabix \
        python=3.12 \
    && micromamba clean --all --yes \
    && bcftools --version >/dev/null \
    && python3 --version >/dev/null

COPY --from=rust-builder /src/target/release/propagate-annotations /usr/local/bin/propagate-annotations
RUN chmod 0755 /usr/local/bin/propagate-annotations \
    && propagate-annotations --help >/dev/null

WORKDIR /work
CMD ["/bin/bash"]
