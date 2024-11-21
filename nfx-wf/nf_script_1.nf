#!/usr/bin/env nextflow

params.aln_dir = ''
params.ped_file = ''

log.info """\
    This is
    a multiline
    message
"""

process aln {
    input:
    tuple val(sample_id), file(fastq_ch1), file(fastq_ch2)

    output:
    stdout

    script:
    """ 
    echo "Sample ID: $sample_id"
    echo "Fastq file 1: $fastq_ch1"
    echo "Fastq file 2: $fastq_ch2"
    """
}

workflow {
    Channel.fromPath( "${params.ped_file}" )
        | splitText
        | map { it.split('\t') } 
        | map { row -> row[1] }
        | set { sample_ids } 

    Channel.fromPath( "${params.aln_dir}" )
        | listFiles( pattern: "*_R{1,2}*.fastq" )
        | map { it.baseName.replaceAll('_R[1,2]', '') }
        | set { fastq_ids }


    // Channel.fromFilePairs( "${params.aln_dir}/${sample_ids}*_R{1,2}*.fastq" )
    //     | map { k, fastqs -> tuple(k, fastqs[0], fastqs[1]) }
    //     | view
    // //     | aln
    // //     | view
    
}