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
    && HTSLIB_CONFIGURE_OPTIONS="--enable-libcurl" pip3 -v --no-cache-dir install pysam --no-binary :all:

COPY add_illumina_readgroups.py /usr/local/bin/
