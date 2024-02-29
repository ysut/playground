#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.input = ''
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
    path '*_splai.vcf'

    script:
    """
    spliceai \\
      -I $vcf \\
      -O ${vcf.baseName}_splai.vcf \\
      -R /resources/ReferenceGenome/exome_pipeline/human_g1k_v37_fix.fasta \\
      -A /resources/SAI10kcalc/SAI10kcalc_spliceai_tx.hg19.tsv \\
      -D 4999 \\
      -M 0
    """
}

process MERGE {
    input:
    path vcf

    output:
    path 'merged_pre_annovar.vcf'

    script:
    """
    cat $vcf > merged_pre_annovar.vcf
    """
}


process ANNOVAR {
    input
    path vcf

    output:
    path 'exome_test_refGene.hg19_multianno.*'



}


workflow {
    input_ch = Channel.fromPath(params.input)

    SPLIT(input_ch)
        .flatMap{ file -> file }
        .set{ splitFiles_ch }

    SPLICEAI(splitFiles_ch)
        .set{ spliceaiFiles_ch }

    MERGE(spliceaiFiles_ch)
        .set{ mergeFile_ch }
    
    mergeFile_ch.view()

}
