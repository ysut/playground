process GATK_GCNV {
    publishDir "${params.out_root}/gcnv", mode: 'symlink'

    input:
    tuple val(sample_id), path(bam), path(bai)

    output:
    tuple val(sample_id), path("*.gcnv.bam"), path("*.gcnv.bam.bai")

    script:
    """
    touch ${sample_id}.gcnv.bam
    touch ${sample_id}.gcnv.bam.bai
    """
}

