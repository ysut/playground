#!/usr/bin/env nextflow

params.aln_dir = ''
params.ped_file = ''

// params.output = "${launchDir}"
params.output = params.output ?: "${workflow.launchDir}"
params.out_root = "${params.output}/WES-WF_Results_" + new Date().format('yyyyMMdd-HHmmss')


log.info """\
    This is
    a multiline
    message
"""

process MAKE_OUTPUT_ROOT_DIR {
    script:
    """
    mkdir -p ${params.out_root}
    """
}

process MAKE_FAMILY_SAMPLE_DIR {
    input:
    tuple val(family), val(sample_id)

    script:
    """
    mkdir -p ${params.out_root}/${family}/${sample_id}
    """
}

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

    // 1=male, 2=female, 0/-9=unknown/missing
    // 1=unaffected, 2=affected, 0/-9=missing

workflow {
    // Output directory
    MAKE_OUTPUT_ROOT_DIR()

    Channel.fromPath("${params.ped_file}")
        | splitText
        | map { it.split('\t') }
        | filter { it.size() >= 6 }
        | map { columns -> tuple(columns[0], columns[1]) }
        | MAKE_FAMILY_SAMPLE_DIR

    Channel.fromPath("${params.ped_file}")
        | splitText
        // | map { line -> line.trim().split('\t') }
        | map { it.split('\t') }
        | filter { it.size() >= 6 }
        | map { columns -> tuple(
                columns[0], [
                    individualId: columns[1],
                    fatherId    : columns[2],
                    motherId    : columns[3],
                    sex         : columns[4],
                    phenotype   : columns[5]
                    ]
                    )
            }
        | groupTuple()  // group by family ID (= colums[0])
        // | set { families }
        | view

    // families
    //     | map { familyTuple ->
    //     def familyId = familyTuple[0]
    //     def members = familyTuple[1..-1].flatten()
    //     def memberCount = members.size()
    //     def analysisType = ''

    //     if (memberCount == 3) {
    //         analysisType = 'trio'
    //     } else if (memberCount == 2) {
    //         analysisType = 'duo'
    //     } else if (memberCount == 4) {
    //         analysisType = 'quad'
    //     } else {
    //         analysisType = 'singleton'
    //     }

    //     return [familyId: familyId, members: members, analysisType: analysisType]
    //     }
    //     | set { familyAnalysis }

    // Find fastq files from the alignment directory for each sample
    // Channel.frompath
    
     // Channel.fromPath( "${params.aln_dir}" )
    //     | listFiles( pattern: "*_R{1,2}*.fastq" )
    //     | map { it.baseName.replaceAll('_R[1,2]', '') }
    //     | set { fastq_ids }


    // Channel.fromFilePairs( "${params.aln_dir}/${sample_ids}*_R{1,2}*.fastq" )
    //     | map { k, fastqs -> tuple(k, fastqs[0], fastqs[1]) }
    //     | view
    // //     | aln
    // //     | view
    
}