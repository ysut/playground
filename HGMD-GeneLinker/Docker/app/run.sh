#!/bin/bash

source ${CONDA_DIR}/etc/profile.d/conda.sh
conda activate maria
python /home/summarize.py ${HGMD_VERSION} ${BIND_DIR}