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
    path '*.splai.vcf'

    script:
    """
    spliceai \\
      -I $vcf \\
      -O ${vcf.baseName}.splai.vcf \\
      -R /resources/ReferenceGenome/exome_pipeline/human_g1k_v37_fix.fasta \\
      -A /resources/SAI10kcalc/SAI10kcalc_spliceai_tx.hg19.tsv \\
      -D 4999 \\
      -M 0
    """
}

process CONCATNATE {
    input:
    path vcfs

    output:
    path 'splai.concat.pre_annovar.vcf'


    script:
    """
    ${params.bcftools} concat \\
      --output splai.merge.pre_annovar.vcf \\
      --output-type v \\
      --threads 2 \\
      ${vcfs}
    """
}


// process ANNOVAR {
//     input
//     path vcf

//     output:
//     path 'exome_test_refGene.hg19_multianno.*'
// }


workflow {
    input_ch = Channel.fromPath(params.input)

    SPLIT(input_ch)
        .flatMap{ file -> file }
        .set{ splitFiles_ch }

    SPLICEAI(splitFiles_ch)
        .collect()
        .set{ spliceaiFiles_ch }

    spliceaiFiles_ch
        .collect()
        .set{ collectedSpliceaiFiles_ch }    

    CONCATNATE(collectedSpliceaiFiles_ch)
        .set{ mergeFile_ch }
    
    mergeFile_ch.view()

}
