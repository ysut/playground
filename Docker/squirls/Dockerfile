FROM eclipse-temurin:23.0.1_11-jre

ARG SQUIRLS_VERSION=2.0.1
ARG BASE_URL="https://github.com/monarch-initiative/Squirls/releases/download"

WORKDIR /opt/tmp

RUN set -ex \
    && apt-get update --allow-releaseinfo-change > /dev/null \
    && apt-get install --no-install-recommends --yes \
        unzip \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* 

RUN set -ex \
    && mkdir -p /opt/squirls-cli-${SQUIRLS_VERSION} \
    && wget -c ${BASE_URL}/v${SQUIRLS_VERSION}/squirls-cli-${SQUIRLS_VERSION}-distribution.zip \
    && unzip squirls-cli-${SQUIRLS_VERSION}-distribution.zip \
    && mv squirls-cli-${SQUIRLS_VERSION}/squirls-cli-${SQUIRLS_VERSION}.jar \
        /opt/squirls-cli-${SQUIRLS_VERSION}/ \
    && rm -rf /opt/tmp

WORKDIR /

# docker run -it --rm -v $(pwd):/squirls -e SQUIRLS_DATA=/squirls/2203_hg19 \
# c016e5eec61b \
# java -jar /squirls/squirls-cli-2.0.1/squirls-cli-2.0.1.jar A \
# --data-directory=/squirls/2203_hg19 \
# --output-format=vcf \
# --out-dir=/squirls \
# /squirls/squirls-cli-2.0.1/examples/example.vcf docker_test_output