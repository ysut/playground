#!/usr/bin/env nextflow

process ANALYZE {

    input:
     fastq from fastq_ch

    output:
    file "output.txt" into output_ch

    script:
    """
    echo "Hello world!"
    """
}