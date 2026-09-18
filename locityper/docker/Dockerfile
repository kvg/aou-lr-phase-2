# print_reads (broadinstitute/str-analysis) for Nearline-efficient CRAM subsets.
# Build: locityper/build_docker.sh

FROM python:3.11-slim-bookworm

ENV DEBIAN_FRONTEND=noninteractive \
    PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        ca-certificates \
        libcurl4 \
        libssl3 \
    && rm -rf /var/lib/apt/lists/* \
    && pip install --no-cache-dir --upgrade pip \
    && pip install --no-cache-dir \
        'str-analysis' \
        'google-cloud-storage'

RUN python3 -c "import str_analysis.print_reads; import google.cloud.storage; import pysam; import intervaltree"

WORKDIR /work
CMD ["/bin/bash"]
