#!/bin/sh

set -eu

nextflow -log log/nflog run nf_script_1.nf --ped_file myfiles/cohort.ped --fastq_dir myfiles/fastq
