#!/usr/bin/env nextflow

process ANALYZE {
    script:
    """
    cd ../ &&
    python -m wesanno \
    --input ${input} \
    --output ${output} \
    --resources ${resources}    
    """
}

workflow {
    ANALYZE()
}