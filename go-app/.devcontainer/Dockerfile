# Use official Golang image
FROM golang:1.24.1

RUN set -ex \
	&& apt-get update > /dev/null \
	&& apt-get install --no-install-recommends --yes \
		git \
	&& apt-get clean \
	&& rm -rf /var/lib/apt/lists/*

WORKDIR /workspace