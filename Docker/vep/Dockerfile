FROM ensemblorg/ensembl-vep:release_113.2

ENV PERL5LIB=${PERL5LIB_TMP}:/vep/loftee

#1. Set build arguments
ARG SAMTOOLS_VERSION=1.21
ARG LOFTEE_VERSION=1.0.3

#2. Install samtools for loftee
USER root
RUN set -ex \
    # Install dependencies
    && apt-get update > /dev/null \
    && apt-get install --no-install-recommends --yes \
        libncurses5-dev \
        > /dev/null \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* \
    # Install samtools
    && mkdir -p /opt/tmp \    
    && curl -sLk --output /opt/tmp/samtools-${SAMTOOLS_VERSION}.tar.bz2 \
        https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2 \
    && tar -jxvf /opt/tmp/samtools-${SAMTOOLS_VERSION}.tar.bz2 -C /opt/tmp \
    && cd /opt/tmp/samtools-${SAMTOOLS_VERSION} \
    && ./configure --prefix=${OPT_SRC}/samtools-${SAMTOOLS_VERSION}  \
    && make \
    && make install \
    && cp ${OPT_SRC}/samtools-${SAMTOOLS_VERSION}/bin/* /usr/local/bin \
    # Install loftee
    && curl -sLk --output /opt/tmp/master.zip \
        https://github.com/konradjk/loftee/archive/refs/tags/v${LOFTEE_VERSION}.zip \
    && unzip /opt/tmp/master.zip -d /opt/tmp \
    && rm -rf /opt/tmp/loftee-${LOFTEE_VERSION}/LICENSE /opt/tmp/loftee-${LOFTEE_VERSION}/README.md \
    # Rename loftee folder
    && mv /opt/tmp/loftee-${LOFTEE_VERSION}/ /opt/tmp/loftee \
    # Move loftee to vep folder
    && mkdir -p /vep/loftee \
    && cp -r /opt/tmp/loftee/* /vep/loftee \
    # Clean up
    && rm -rf /opt/tmp

#3. Finalize
# RUN ln -s $OPT/.vep /data
USER vep
WORKDIR /data