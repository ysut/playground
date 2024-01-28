FROM python:3.12.1-slim-bullseye

ENV PYTHONPATH=${PYTHONPATH}:/app

WORKDIR /app

RUN set -ex \
    && apt-get update --allow-releaseinfo-change > /dev/null \
    && apt-get install --no-install-recommends --yes \
        build-essential \
        procps \
        tini \
        zlib1g-dev \
        > /dev/null \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* \
    && pip install --upgrade pip \
    && pip install --no-cache-dir \
        numpy openpyxl pandas pandarallel tqdm pathlib2 pyyaml liftover

COPY . /app/
ENTRYPOINT ["tini", "--"]
CMD [ "/bin/bash" ]
