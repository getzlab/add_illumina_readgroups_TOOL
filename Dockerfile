FROM google/cloud-sdk:slim

RUN set -ex \
    && curl https://packages.cloud.google.com/apt/doc/apt-key.gpg | apt-key add - \
    && apt-get update -y --allow-releaseinfo-change && apt-get install -y --no-install-recommends \
        autoconf \
        automake \
        gcc \
        perl \
        bzip2 \
        zlib1g-dev \
        libbz2-dev \
        liblzma-dev \
        libcurl4-gnutls-dev \
        libssl-dev \
        libncurses5-dev \
    && rm -rf /var/lib/apt/lists/* \
    && pip3 --no-cache-dir install pysam

COPY add_illumina_readgroups.py /usr/local/bin/
