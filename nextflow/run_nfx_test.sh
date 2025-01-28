#!/bin/sh

set -eu

nextflow run nf_script_1.nf --ped_file myfiles/cohort.ped --fastq_dir myfiles/fastq
