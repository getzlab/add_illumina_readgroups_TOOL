FROM google/cloud-sdk:slim

RUN set -ex \
    && apt-get update -y --allow-releaseinfo-change && apt-get install -y --no-install-recommends \
        autoconf \
        automake \
        gcc \
        make \
        perl \
        bzip2 \
        zlib1g-dev \
        libbz2-dev \
        liblzma-dev \
        libcurl4-gnutls-dev \
        libssl-dev \
        libncurses-dev \
    && rm -rf /var/lib/apt/lists/* \
    && HTSLIB_CONFIGURE_OPTIONS="--enable-libcurl" pip3 -v --no-cache-dir install --break-system-packages pysam --no-binary :all:

COPY add_illumina_readgroups.py /usr/local/bin/
