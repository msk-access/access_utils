################## BASE IMAGE ######################

FROM python:3.6-slim

################## ARGUMENTS/Environments ##########

ARG BUILD_DATE
ARG BUILD_VERSION
ARG LICENSE="Apache-2.0"
ARG ACCESS_UTILS_VERSION
ARG VCS_REF

################## METADATA ########################
LABEL org.opencontainers.image.vendor="MSKCC"
LABEL org.opencontainers.image.authors="Ian Johnson (johnsoni@mskcc.org)"

LABEL org.opencontainers.image.created=${BUILD_DATE} \
    org.opencontainers.image.version=${BUILD_VERSION} \
    org.opencontainers.image.licenses=${LICENSE} \
    org.opencontainers.image.version.sequence_qc=${ACCESS_UTILS_VERSION} \
#    org.opencontainers.image.source.sequence_qc="https://pypi.org/project/access_utils/" \
    org.opencontainers.image.vcs-url="https://github.com/msk-access/access_utils.git" \
    org.opencontainers.image.vcs-ref=${VCS_REF}

LABEL org.opencontainers.image.description="This container uses python3.6 as the base image to build \
    access_utils version ${ACCESS_UTILS_VERSION}"

################## INSTALL ##########################

RUN apt-get -y update \
    && apt-get -y install wget unzip \
    && cd /tmp \
    && wget https://github.com/msk-access/access_utils/archive/refs/heads/multiqc-utils.zip \
    && unzip multiqc-utils.zip \
    && cd access_utils-multiqc-utils \
    && pip install -r requirements.txt \
    && mv ./multiqc_prep/general_stats_parse.py /usr/bin \
    && chmod 755 /usr/bin/general_stats_parse.py \
    && rm -rf /tmp/*
