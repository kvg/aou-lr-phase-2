# bw2/ExpansionHunter for local minicram genotyping.
# FASTA is a WDL input (do not bake Homo_sapiens_assembly38 into the image).
# Build: expansion_hunter/build_docker.sh

FROM debian:bookworm-slim

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        ca-certificates \
        wget \
        git \
        build-essential \
        cmake \
        pkg-config \
        zlib1g-dev \
        libbz2-dev \
        liblzma-dev \
        autoconf \
        libcurl4-openssl-dev \
        libssl-dev \
        libdeflate-dev \
        libboost-all-dev \
        python3 \
    && rm -rf /var/lib/apt/lists/*

RUN wget -q "https://github.com/samtools/htslib/releases/download/1.22/htslib-1.22.tar.bz2" \
    && tar -xjf htslib-1.22.tar.bz2 \
    && cd htslib-1.22 \
    && ./configure \
    && make -j"$(nproc)" \
    && make install \
    && ldconfig \
    && cd .. \
    && rm -rf htslib-1.22 htslib-1.22.tar.bz2

RUN git clone --depth 1 https://github.com/bw2/ExpansionHunter.git \
    && cd ExpansionHunter \
    && mkdir build \
    && cd build \
    && cmake .. \
    && make -j"$(nproc)" \
    && cp install/bin/ExpansionHunter /usr/bin/ \
    && cd ../.. \
    && rm -rf ExpansionHunter

RUN ExpansionHunter --version || ExpansionHunter --help | head -n 20

WORKDIR /work
CMD ["/bin/bash"]
