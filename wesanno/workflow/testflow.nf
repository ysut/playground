#!/usr/bin/env nextflow

process ANALYZE {
    script:
    """
    python -m wesanno \
    --input ${input} \
    --output ${output} \
    --resources ${resources}    
    """
}

workflow {
    ANALYZE()
}