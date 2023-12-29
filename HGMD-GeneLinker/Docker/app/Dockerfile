FROM debian:bookworm-slim


# #0. Set sql.cnf path
ENV SQL_CNF=/home/sql.cnf

#1. Set build arguments
#1-1. miniforge
ARG MINIFORGE_NAME=Miniforge3
ARG MINIFORGE_VERSION=23.1.0-4
ARG MINIFORGE_GITHUB=https://github.com/conda-forge/miniforge/releases/download/

#1-2. task
ARG TASK_VERSION=v3.31.0
ARG TASK_GITHUB=https://github.com/go-task/task/releases/download

#2. Set environment variables
ENV LANG=C.UTF-8 LC_ALL=C.UTF-8
ENV CONDA_DIR=/opt/conda 
ENV PATH=${CONDA_DIR}/bin:${PATH}

WORKDIR /usr/tools/tmp

RUN set -ex \
    && apt-get update > /dev/null \
    && apt-get install --no-install-recommends --yes \
        build-essential \
        bzip2 \
        ca-certificates \
        curl \
        default-mysql-client \
        git \
        sudo \
        tini \
        > /dev/null \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* \
    && curl -sSL \
         ${MINIFORGE_GITHUB}/${MINIFORGE_VERSION}/${MINIFORGE_NAME}-${MINIFORGE_VERSION}-Linux-$(uname -m).sh \
         -o /tmp/miniforge.sh \
    && /bin/bash /tmp/miniforge.sh -b -p ${CONDA_DIR} \
    && rm /tmp/miniforge.sh \
    && conda clean --tarballs --index-cache --packages --yes \
    && find ${CONDA_DIR} -follow -type f -name '*.a' -delete \
    && find ${CONDA_DIR} -follow -type f -name '*.pyc' -delete \
    && conda clean --force-pkgs-dirs --all --yes \
    && echo ". ${CONDA_DIR}/etc/profile.d/conda.sh && conda activate base" >> /etc/skel/.bashrc \
    && echo ". ${CONDA_DIR}/etc/profile.d/conda.sh && conda activate base" >> ~/.bashrc

# Install conda packages using environment.yml
COPY ./environment.yaml ./
RUN set -ex \
    && conda update --name base --channel conda-forge conda \
    && conda env create --file environment.yaml \
    && conda clean --force-pkgs-dirs --all --yes

# Install Task
RUN set -ex \
    && UNAME_M=$(uname -m) \
    && curl -ssL ${TASK_GITHUB}/${TASK_VERSION}/task_checksums.txt \
        -o /usr/tools/tmp/task_checksums.txt \
    && if [ "${UNAME_M}" = "x86_64" ]; then \
         url="${TASK_GITHUB}/${TASK_VERSION}/task_linux_amd64.deb"; \
         SHA256SUM=$(grep task_linux_amd64.deb task_checksums.txt \
           | awk '{print $1}'); \
       elif [ "${UNAME_M}" = "aarch64" ]; then \
         url="${TASK_GITHUB}/${TASK_VERSION}/task_linux_arm64.deb"; \
         SHA256SUM=$(grep task_linux_arm64.deb task_checksums.txt \
           | awk '{print $1}'); \
       fi \
    && curl -sSL ${url} -o /usr/tools/tmp/task.deb \
    && echo "${SHA256SUM} /usr/tools/tmp/task.deb" >shasum \
    && sha256sum --check --status shasum \
    && dpkg -i /usr/tools/tmp/task.deb

RUN rm -rf /usr/tools/tmp

USER nobody
COPY --chown=nobody:nogroup ./* /home/

WORKDIR /home
ENTRYPOINT ["tini", "--"]
CMD [ "/bin/bash" ]