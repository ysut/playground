process EXPANSIONHUNTER {
    publishDir "${params.out_root}/expansionhunter", mode: 'symlink'

    input:
    tuple val(sample_id), path(bam), path(bai)

    output:
    tuple val(sample_id), path("*.expansionhunter.bam"), path("*.expansionhunter.bam.bai")

    script:
    """
    touch ${sample_id}.expansionhunter.bam
    touch ${sample_id}.expansionhunter.bam.bai
    """
}