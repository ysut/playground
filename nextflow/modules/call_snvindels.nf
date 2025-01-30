process DEEPVARIANT {
    publishDir "${params.out_root}/${family_id}/raw_vcfs", mode: 'symlink'

    input:
    tuple val(family_id), val(sample_id), val(anlysis_type)

    output:
    // tuple val(family_id), val(sample_id), path("*.vcf.gz"), path("*.vcf.gz.tbi")
    stdout

    script:
    """
    echo ${params.out_root}/xams/${sample_id}.bam

    touch ${sample_id}.vcf.gz
    touch ${sample_id}.vcf.gz.tbi
    """
}


process DEEPTRIO {
    publishDir "${params.out_root}/${family_id}/raw_vcfs", mode: 'symlink'

    input:
    tuple val(family_id), val(proband_id), val(mother_id), val(father_id), val(analysis_type)

    output:
    tuple val(family_id), path("*.g.vcf.gz"), path("*.g.vcf.gz.tbi")

    script:
    """
    echo "Processing Family ID: ${family_id}"
    echo "Analysis Type: ${analysis_type}"

    # 定義済みの BAM ファイルパス
    PROBAND_BAM="${params.out_root}/xams/${proband_id}.*am"
    MOTHER_BAM="${params.out_root}/xams/${mother_id}.*am"
    FATHER_BAM="${params.out_root}/xams/${father_id}.*am"

    # 分析タイプに基づく処理の分岐
    if [ "${analysis_type}" == "duo" ]; then
        echo "Running DEEPTRIO for Duo"
        echo "Proband BAM: \${PROBAND_BAM}"
        echo "Mother BAM: \${MOTHER_BAM}"
        
        # Docker command for Duo
  
    elif [ "${analysis_type}" == "trio" ]; then
        echo "Running DEEPTRIO for Trio"
        echo "Proband BAM: \${PROBAND_BAM}"
        echo "Mother BAM: \${MOTHER_BAM}"
        echo "Father BAM: \${FATHER_BAM}"
        
        # Docker command for Trio

    else
        echo "Unknown analysis type: ${analysis_type}"
        exit 1
    fi

    # VCFファイルのダミー作成（実際の処理では不要）
    touch ${proband_id}_proband.vcf.gz
    touch ${proband_id}_proband.vcf.gz.tbi
    touch ${proband_id}_proband.g.vcf.gz
    touch ${proband_id}_proband.g.vcf.gz.tbi
    touch ${mother_id}_mother.vcf.gz
    touch ${mother_id}_mother.vcf.gz.tbi
    touch ${mother_id}_mother.g.vcf.gz
    touch ${mother_id}_mother.g.vcf.gz.tbi
    touch ${father_id}_father.vcf.gz
    touch ${father_id}_father.vcf.gz.tbi
    touch ${father_id}_father.g.vcf.gz
    touch ${father_id}_father.g.vcf.gz.tbi
    """
}

process DEEPTRIO_FOR_QUAD {
    publishDir "${params.out_root}/${family_id}/raw_vcfs", mode: 'symlink'

    input:
    tuple val(family_id), val(proband_id), val(mother_id), val(father_id), val(sibling_id), val(analysis_type)

    output:
    tuple val(family_id), path("*.g.vcf.gz"), path("*.g.vcf.gz.tbi"), val(sibling_id)

    script:
    """
    echo "Processing Family ID: ${family_id}"
    echo "Analysis Type: ${analysis_type}"

    # 定義済みの BAM ファイルパス
    PROBAND_BAM="${params.out_root}/xams/${proband_id}.*am"
    MOTHER_BAM="${params.out_root}/xams/${mother_id}.*am"
    FATHER_BAM="${params.out_root}/xams/${father_id}.*am"
    SIBLING_BAM="${sibling_id ? "${params.out_root}/xams/${sibling_id}.*am" : "null"}"

    # 分析タイプに基づく処理の分岐
    # docker run for deeptrio

    # VCFファイルのダミー作成（実際の処理では不要）
    touch ${proband_id}_proband.vcf.gz
    touch ${proband_id}_proband.vcf.gz.tbi
    touch ${proband_id}_proband.g.vcf.gz
    touch ${proband_id}_proband.g.vcf.gz.tbi
    touch ${mother_id}_mother.vcf.gz
    touch ${mother_id}_mother.vcf.gz.tbi
    touch ${mother_id}_mother.g.vcf.gz
    touch ${mother_id}_mother.g.vcf.gz.tbi
    touch ${father_id}_father.vcf.gz
    touch ${father_id}_father.vcf.gz.tbi
    touch ${father_id}_father.g.vcf.gz
    touch ${father_id}_father.g.vcf.gz.tbi
    """
}


process DEEPVARIANT_FOR_SIBLING {
    publishDir "${params.out_root}/${family_id}/raw_vcfs", mode: 'symlink'

    input:
    tuple val(family_id), path(gvcfs_wo_sibling), path(gvcfs_index_wo_sibling), val(sibling_id)

    output:
    tuple val(family_id), path(gvcfs_wo_sibling), path(gvcfs_index_wo_sibling), path("*_sibling.g.vcf.gz"), path("*_sibling.g.vcf.gz.tbi")
    // stdout

    script:
    """
    # docker run for deepvariant

    touch ${sibling_id}_sibling.vcf.gz
    touch ${sibling_id}_sibling.vcf.gz.tbi
    touch ${sibling_id}_sibling.g.vcf.gz
    touch ${sibling_id}_sibling.g.vcf.gz.tbi
    """
}

process GLNEXUS {
    publishDir "${params.out_root}/${family_id}", mode: 'symlink'

    input:
    tuple val(family_id), path(input_gvcfs), path(input_gvcfs_index)

    output:
    tuple val(family_id), path("*.joint.bcf.gz"), path("*.joint.bcf.gz.tbi")
    // stdout

    script:
    """
    echo "Processing Family ID: ${family_id}"
    echo "Processing gVCFs: ${input_gvcfs}"

    touch ${family_id}.joint.bcf.gz
    touch ${family_id}.joint.bcf.gz.tbi
    """
}
