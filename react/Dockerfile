FROM node:23.9.0-bullseye-slim

RUN set -ex \
    && apt-get update > /dev/null \
    && apt-get install --no-install-recommends --yes \
    git \
    vim \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

USER node
WORKDIR /home/node/app

EXPOSE 5173