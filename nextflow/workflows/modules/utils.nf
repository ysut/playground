process MAKE_OUTPUT_ROOT_DIR {
    script:
    """
    mkdir -p ${params.out_root}
    """
}

process MAKE_FAMILY_SAMPLE_DIR {
    input:
    tuple val(familyID), val(individual_id)

    output:
    tuple val(individual_id), val(familyID)

    script:
    """
    mkdir -p ${params.out_root}/${familyID}/${individual_id}
    """
}
