#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.workflow_WES = '/Users/utsu/work/Github/playground/wesanno/scripts'

process SPLIT {
    input:
    path vcf

    output:
    path '*.part*'

    script:
    """
    ${params.workflow_WES}/split.sh $vcf
    """
}

process PRINT {
    input:
    path split_vcf

    output:
    stdout

    script:
    """
    echo "Printing $split_vcf"
    ${params.workflow_WES}/printhead.sh $split_vcf
    """
}

workflow {
    input_vcf = Channel.fromPath(params.inputvcf)

    splitFiles = SPLIT(input_vcf)
    
    PRINT(splitFiles)
        .view { "Result: ${it.text}" }
}