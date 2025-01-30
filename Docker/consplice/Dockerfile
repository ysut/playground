FROM mambaorg/micromamba:2.0.4-debian12-slim

RUN set -ex \
	&& micromamba config append channels defaults \
	&& micromamba config append channels ggd-genomics \
	&& micromamba config append channels bioconda \
	&& micromamba config append channels conda-forge \
	&& micromamba create --name ConSplice python=3.8 --yes

SHELL ["micromamba", "run", "-n", "ConSplice", "/bin/bash", "-c"]

RUN set -ex \
    && micromamba install -c bioconda consplice ggd --yes \
	# Data installation for ConSplice using ggd
	&& ggd install -c genomics grch38-reference-genome-gencode-v1 \
	&& ggd install -c genomics grch38-canonical-transcript-features-gencode-v1 \
	&& ggd install -c genomics grch38-segmental-dups-ucsc-v1 \
	&& ggd install -c genomics grch38-self-chain-ucsc-v1 \
	&& micromamba clean --all --yes

WORKDIR /


# bash -c ". /usr/local/bin/_activate_current_env.sh && micromamba activate ConSplice && ......"