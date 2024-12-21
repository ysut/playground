#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

params.fastq_dir = '/Users/utsu/work/Github/playground/nextflow/myfiles/fastq'
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
    mkdir -p ${params.out_root}/fastq
    """
}

process MAKE_FAMILY_SAMPLE_DIR {
    input:
    tuple val(familyID), val(individual_id)

    output:
    tuple val(familyID), val(individual_id)

    script:
    """
    mkdir -p ${params.out_root}/${familyID}/${individual_id}
    mkdir -p ${params.out_root}/${familyID}/xams
    """
}

process FIND_FASTQ {
    publishDir "${params.out_root}/fastq", mode: 'symlink'

    input:
    tuple val(familyId), val(individualId)

    output:
    tuple val(familyId), val(individualId), path("*_R1*.fastq.gz"), path("*_R2*.fastq.gz")

    script:
    """
    echo "Searching fastq for Family: ${familyId}, Individual: ${individualId}"

    set -e

    find "${params.fastq_dir}" -type f \\
      -name "*${individualId}*_R1*.fastq.gz" -exec ln -s {} . \\; || true
    find "${params.fastq_dir}" -type f \\
      -name "*${individualId}*_R2*.fastq.gz" -exec ln -s {} . \\; || true

    ls -lh ${params.out_root}/fastq/*.fastq.gz 2>/dev/null \\
      || echo "No FASTQ files found for Individual: ${individualId}"
    """
}

process STROBEALIGN {
    publishDir "${params.out_root}/${familyID}/xams", mode: 'symlink'

    input:
    tuple val(familyID), val(individualID), path(fastq_R1), path(fastq_R2)

    output:
    tuple val(familyID), val(individualID), path("*.{b,cr}am"), path("*.*am.{b,cr}ai")

    script:
    """ 
    touch ${individualID}.bam
    touch ${individualID}.bam.bai
    """
}


    // 1=male, 2=female, 0/-9=unknown/missing
    // 1=unaffected, 2=affected, 0/-9=missing

workflow {
    // Output directory
    MAKE_OUTPUT_ROOT_DIR()

    // Read pedigree file
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
        | set { pedigree_info }
    
    // Make output directories for each Family and Sample
    pedigree_info
        | flatMap { familyId, individuals ->
            individuals.collect { indi -> tuple(familyId, indi.individualId)}
            }
        | MAKE_FAMILY_SAMPLE_DIR
        | FIND_FASTQ
        | STROBEALIGN
        | view

    // pedigree_info
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


    // Channel.fromFilePairs( "${params.aln_dir}/${individual_ids}*_R{1,2}*.fastq" )
    //     | map { k, fastqs -> tuple(k, fastqs[0], fastqs[1]) }
    //     | view
    // //     | aln
    // //     | view
    
}