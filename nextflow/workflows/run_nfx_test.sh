#!/bin/sh

set -eu

nextflow -log log/nflog run nf_script_1.nf --ped_file myfiles/cohort.ped --fastq_dir myfiles/fastq

nextflow -log log/nfx.log run -resume workflows/main.nf --ped_file workflows/myfiles/server_test.ped --fastq_dir files/fastqs/

nextflow -log log/nfx.log run workflows/main.nf --ped_file workflows/myfiles/server_test.ped --fastq_dir files/fastqs/