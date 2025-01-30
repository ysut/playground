FROM continuumio/miniconda3:24.9.2-0

WORKDIR /opt

# RUN
RUN set -ex \
    && apt-get update --allow-releaseinfo-change > /dev/null \
    && apt-get install --no-install-recommends --yes \
        ca-certificates \
        git \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* \
    && conda update --yes --name base --channel defaults conda \
    && conda create --yes --name pangolin python=3.9 

SHELL ["conda", "run", "-n", "pangolin", "/bin/bash", "-c"]

RUN set -ex \
    && conda install --yes --channel conda-forge pyvcf \
    && pip3 install -U pip \
    && pip3 install biopython==1.84 gffutils==0.13 pyfastx==2.1.0 \
    # && conda install --yes --channel bioconda --channel conda-forge \
        # biopython gffutils pyfastx \
    && pip3 install pandas==2.2.3 numpy==2.0.2 \
    && pip3 install torch torchvision torchaudio --index-url \
        https://download.pytorch.org/whl/cu121 \
    && git clone https://github.com/tkzeng/Pangolin.git \
    && cd Pangolin \
    && pip3 install . \
    && pip3 cache purge \
    && conda clean --yes --tarballs --index-cache --packages \
    && conda clean --yes --force-pkgs-dirs --all \
    && rm -rf ~/.cache/pip \
    && rm -rf /opt/Pangolin

WORKDIR /