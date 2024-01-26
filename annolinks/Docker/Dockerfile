FROM condaforge/miniforge3:23.3.1-1

ENV LANG=C.UTF-8 LC_ALL=C.UTF-8 
ENV CONDA_DIR=/opt/conda 
ENV PATH=${CONDA_DIR}/bin:${PATH}
ENV PYTHONPATH=${PYTHONPATH}:/app

WORKDIR /app

RUN set -ex \
    && apt-get update --allow-releaseinfo-change > /dev/null \
    && apt-get install --no-install-recommends --yes \
        procps \
        > /dev/null \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* \
    && conda update -n base conda --yes \
    && conda install --yes \
         numpy openpyxl pandas pandarallel pathlib2 pyyaml tqdm \
    && conda clean --tarballs --index-cache --packages --yes

COPY . /app/
ENTRYPOINT ["tini", "--"]
CMD [ "/bin/bash" ]
