#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// params.input = ''
// params.phenotypes = ''

// The path below is container specific 
params.fasta = '/resources/ReferenceGenome/exome_pipeline/human_g1k_v37_fix.fasta'



process GADO {
    input:
    path phenotypes

    output:
    path 'GADO_results/*.txt'

    script:
    """
    ${WORKFLOW_WES}/scripts/run_gado.sh $phenotypes ./GADO_results/
    mkdir ./GADO_results/log
    mv ./GADO_results/hpoProcessed.txt* ./GADO_results/log/
    mv ./GADO_results/samples.txt ./GADO_results/log/
    cp ./GADO_results/*.txt ${workflow.launchDir}/
    """
}

workflow {
    fasta = "${file(params.fasta).toAbsolutePath()}"

    input_vcf_ch = Channel.fromPath(params.input)
    
    MAVERICK(input_ch)
        .set{ maverickFiles_ch }

    maverickFiles_ch.view()

    if (params.phenotypes) {
        phenotypes_ch = Channel.fromPath(params.phenotypes)
        GADO(phenotypes_ch)
            .set{ gadoFile_ch }
    }


    // ANNOVAR(concatFile_ch)
    //     .set{ annovarFiles_ch }

    // FORMATANNOVAR(annovarFiles_ch)
    //     .set{ renamedFile_ch }

    // HGMDANNOTATOR(renamedFile_ch)
    //     .set{ hgmdAnnotatedFile_ch }

    // hgmdAnnotatedFile_ch.view()    

}
