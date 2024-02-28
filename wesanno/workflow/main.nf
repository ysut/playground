#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.workflow_WES = '/betelgeuse07/analysis/utsu/MyTools/workflow_WES'

process SPLIT {
    input:
    path vcf

    output:
    path 'split_*'

    script:
    """
    ${params.workflow_WES}/scripts/extract_chr_split_vcf.sh $vcf
    """
}

process SPLICEAI {
    input:
    path vcf

    output:
    path 'spliceai_*.vcf'

    container 'betelegeuse:5000/library/utsu:spliceai'

    script:
    """
    conda activate spliceai
    spliceai \\
      -I \$vcf \\
      -O spliceai_\${vcf}.vcf \\
      -R /resouces/ReferenceGenome/exome_pipeline/human_g1k_v37_fix.fasta \\
      -A /resources/SAI10kcalc/SAI10kcalc_spliceai_tx.hg19.tsv \\
      -D 4999 \\
      -M 0
    """
}


workflow {
    input_ch = Channel.fromPath(params.input)

    SPLIT(input_ch)
        .set{ splitFiles_ch }
        
    SPLICEAI(splitFiles_ch)
        .set{ spliceaiFiles_ch }

    spliceaiFiles_ch.view()

}