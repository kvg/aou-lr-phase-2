# print_reads + make_minicram (broadinstitute/str-analysis) for Nearline CRAMs.
# Also installs gcloud so dsub can use this image as DSUB_CLOUD_SDK_IMAGE
# (VPC-SC cannot pull gcr.io/.../cloud-sdk:slim).
# Build: locityper/build_docker.sh

FROM python:3.11-slim-bookworm

ENV DEBIAN_FRONTEND=noninteractive \
    PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    CLOUDSDK_CORE_DISABLE_PROMPTS=1

# Pin a commit that includes make_minicram_for_expansion_hunter (not on PyPI).
ARG STR_ANALYSIS_GIT=https://github.com/broadinstitute/str-analysis.git
ARG STR_ANALYSIS_REF=41198f2fd8840f182531ba699a6817648fbf4d0b

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        ca-certificates \
        curl \
        gnupg \
        apt-transport-https \
        git \
        libcurl4 \
        libssl3 \
    && curl -fsSL https://packages.cloud.google.com/apt/doc/apt-key.gpg \
        | gpg --dearmor -o /usr/share/keyrings/cloud.google.gpg \
    && echo "deb [signed-by=/usr/share/keyrings/cloud.google.gpg] https://packages.cloud.google.com/apt cloud-sdk main" \
        > /etc/apt/sources.list.d/google-cloud-sdk.list \
    && apt-get update \
    && apt-get install -y --no-install-recommends google-cloud-cli \
    && rm -rf /var/lib/apt/lists/* \
    && pip install --no-cache-dir --upgrade pip \
    && pip install --no-cache-dir \
        "google-cloud-storage" \
        "git+${STR_ANALYSIS_GIT}@${STR_ANALYSIS_REF}"

RUN python3 -c "import str_analysis.print_reads; import str_analysis.make_minicram_for_expansion_hunter; import google.cloud.storage; import pysam" \
    && gcloud version \
    && command -v gcloud

WORKDIR /work
CMD ["/bin/bash"]
