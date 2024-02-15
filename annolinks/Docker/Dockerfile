FROM python:3.12.1-slim-bullseye

ENV PYTHONPATH=${PYTHONPATH}:/app
ARG UCSC_BASE_URL="ftp://hgdownload.soe.ucsc.edu/goldenPath"

WORKDIR /app

RUN set -ex \
    && apt-get update --allow-releaseinfo-change > /dev/null \
    && apt-get install --no-install-recommends --yes \
        build-essential \
        curl \
        procps \
        tini \
        zlib1g-dev \
        > /dev/null \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* \
    && mkdir -p /app/resources \
    && curl -sSL -o /app/resources/hg19ToHg38.over.chain.gz \
        "${UCSC_BASE_URL}"/hg19/liftOver/hg19ToHg38.over.chain.gz \
    && pip install --upgrade pip \
    && pip install --no-cache-dir \
        numpy openpyxl pandas pandarallel tqdm pathlib2 pyyaml liftover

COPY . /app/
ENTRYPOINT ["tini", "--"]
CMD [ "/bin/bash" ]
