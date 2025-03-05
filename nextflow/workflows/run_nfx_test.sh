#!/bin/sh

set -eu

source /betelgeuse07/analysis/utsu/envs/nextflow.sh
scp -r utsu@192.168.2.97:/home/utsu/Desktop/GitHub/playground/nextflow/workflows ./
nextflow -log log/nfx.log run -resume workflows/main.nf --ped_file workflows/myfiles/server_test.ped --fastq_dir files/fastqs/ --bam_dir files/bams/


## Trio -bg
nextflow -log log/nfx.log run -bg -resume workflows/main.nf --ped_file workflows/myfiles/server_test_trio.ped --fastq_dir files/fastqs/ --bam_dir files/bams/ > bg.log

nextflow -log log/nfx.log run -resume workflows/main.nf --ped_file workflows/myfiles/server_test.ped --fastq_dir files/fastqs/

