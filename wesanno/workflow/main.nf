#!/usr/bin/env nextflow
nextflow.enable.dsl=2


process ANALYZE {
    """
    python -m wesanno \
    --input ${params.input} \
    --output ${params.output} \
    --resources ${params.resources} \
    """
}

workflow {
    ANALYZE()
}