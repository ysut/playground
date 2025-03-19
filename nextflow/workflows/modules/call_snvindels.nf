process DEEPVARIANT {
    publishDir "${params.out_root}/${family_id}/raw/gvcfs", 
               mode: 'symlink', pattern: '*.g.vcf.gz*'
    publishDir "${params.out_root}/${family_id}/raw/bcfs", 
               mode: 'symlink', pattern: '*.bcf.gz*'
    publishDir "${params.out_root}/${family_id}/misc/reports/DeepVariant", 
               mode: 'symlink', pattern: '*.html'
    publishDir "${params.out_root}/${family_id}/misc/logs/DeepVariant", 
               mode: 'symlink', pattern: '*.log'

    input:
    tuple val(family_id), path(reference), path(reference_index),
          val(sample_id), path(xam), path(xai),
          val(sex), val(analysis_type),
          path(par_bed), path(par_bed_index)

    output:
    tuple val(family_id), path(reference), path(reference_index),
          val(sample_id),
          path("${sample_id}.bcf.gz"), path("${sample_id}.bcf.gz.csi"), 
          path("${sample_id}.g.vcf.gz"), path("${sample_id}.g.vcf.gz.tbi"), 
          val(analysis_type)

    // 1=male, 2=female, 0/-9=unknown/missing
    script:
    """
    if [ "${sex}" == "1" ]; then
      /opt/deepvariant/bin/run_deepvariant \\
        --model_type=${params.deepvariant_model} \\
        --ref=${reference} \\
        --reads=${xam} \\
        --output_vcf=${sample_id}.vcf.gz \\
        --output_gvcf=${sample_id}.g.vcf.gz \\
        --vcf_stats_report=true \\
        --num_shards=${params.deepvariant_shards} \\
        --haploid_contigs=\\"chrX,chrY\\" \\
        --par_regions_bed=${par_bed} \\
        --dry_run=false && \\
    else
      /opt/deepvariant/bin/run_deepvariant \\
        --model_type=${params.deepvariant_model} \\
        --ref=${reference} \\
        --reads=${xam} \\
        --output_vcf=${sample_id}.vcf.gz \\
        --output_gvcf=${sample_id}.g.vcf.gz \\
        --vcf_stats_report=true \\
        --num_shards=${params.deepvariant_shards} \\
        --dry_run=false && \\
    fi

    /opt/conda/envs/bio/bin/bcftools view \\
      --output-type b \\
      --threads ${params.bcftools_threads} \\
      ${sample_id}.vcf.gz && \\
      > ${sample_id}.bcf 
    /opt/conda/envs/bio/bin/bcftools index ${sample_id}.bcf.gz
    """
}

process GLNEXUS {
    publishDir "${params.out_root}/${family_id}/raw/merged_bcf", 
               mode: 'symlink', pattern: '*.bcf'

    input:
    tuple val(family_id), path(input_vcfs), path(input_vcfs_index), val(analysis_type)

    output:
    tuple val(family_id), path("*.joint.bcf")

    script:
    """
    glnexus_cli \\
      --config DeepVariant_unfiltered \\
      --threads ${params.glnexus_threads} \\
      ${input_vcfs} > ${family_id}.joint.bcf
    """
}


process DEEPTRIO_MALE_CHILD {
    publishDir "${params.out_root}/${family_id}/raw_vcfs/gvcfs", 
               mode: 'symlink', pattern: '*.g.vcf.gz*'
    publishDir "${params.out_root}/${family_id}/raw_vcfs/vcfs", 
               mode: 'symlink', pattern: '*AP.vcf.gz*'
    publishDir "${params.out_root}/${family_id}/raw_vcfs/vcfs", 
               mode: 'symlink', pattern: '*PAR.vcf.gz*'
    publishDir "${params.out_root}/${family_id}/misc/reports/DeepTrio", 
               mode: 'symlink', pattern: '*.html'
    publishDir "${params.out_root}/${family_id}/misc/logs/DeepTrio", 
               mode: 'symlink', pattern: '*.log'

    input:
    tuple val(family_id), path(reference), path(reference_index),
          val(child_id), val(dad_id), val(mom_id), 
          path(child_xam), path(dad_xam), path(mom_xam), 
          path(child_xai), path(dad_xai), path(mom_xai),
          val(sex), val(analysis_type),
          path(par_bed), path(par_bed_index)

    output:
    tuple val(family_id), path(reference), path(reference_index),
          val(child_id), val(dad_id), val(mom_id),
          path("*.AP.g.vcf.gz"), path("*.AP.g.vcf.gz.tbi"), 
          path("*.X_incl_PAR.g.vcf.gz"), path("*.X_incl_PAR.g.vcf.gz.tbi"), 
          path("*.Y_incl_PAR.g.vcf.gz"), path("*.Y_incl_PAR.g.vcf.gz.tbi"), 
          path("dad_*.X_nonPAR.bcf.gz"), path("dad_*.X_nonPAR.bcf.gz.csi"),
          path(par_bed), path(par_bed_index),
          path("*AP.vcf.gz*"), path("*PAR.vcf.gz*"), path("*.html")

    script:
    """
    /opt/deepvariant/bin/deeptrio/run_deeptrio \\
      --regions=${params.autosome_and_par} \\
      --model_type=${params.deepvariant_model} \\
      --ref=${reference} \\
      --reads_child=${child_xam} \\
      --reads_parent1=${dad_xam} \\
      --reads_parent2=${mom_xam} \\
      --output_vcf_child=son_${child_id}.AP.vcf.gz \\
      --output_vcf_parent1=dad_${dad_id}.AP.vcf.gz \\
      --output_vcf_parent2=mom_${mom_id}.AP.vcf.gz \\
      --output_gvcf_child=son_${child_id}.AP.g.vcf.gz \\
      --output_gvcf_parent1=dad_${dad_id}.AP.g.vcf.gz \\
      --output_gvcf_parent2=mom_${mom_id}.AP.g.vcf.gz \\
      --sample_name_child="son_${child_id}" \\
      --sample_name_parent1="dad_${dad_id}" \\
      --sample_name_parent2="mom_${mom_id}" \\
      --postprocess_variants_child_extra_args="--haploid_contigs=\\"chrX,chrY\\",--par_regions_bed=\\"${par_bed}\\"" \\
      --postprocess_variants_parent1_extra_args="--haploid_contigs=\\"chrX,chrY\\",--par_regions_bed=\\"${par_bed}\\"" \\
      --vcf_stats_report=true \\
      --num_shards=${params.deepvariant_shards} \\
      --logging_dir=./ \\
      --dry_run=false && \\

    /opt/deepvariant/bin/deeptrio/run_deeptrio \\
      --regions="chrX" \\
      --model_type=${params.deepvariant_model} \\
      --ref=${reference} \\
      --reads_child=${child_xam} \\
      --reads_parent1=${mom_xam} \\
      --output_vcf_child=son_${child_id}.X_incl_PAR.vcf.gz \\
      --output_vcf_parent1=mom_${mom_id}.X_incl_PAR.vcf.gz \\
      --output_gvcf_child=son_${child_id}.X_incl_PAR.g.vcf.gz \\
      --output_gvcf_parent1=mom_${mom_id}.X_incl_PAR.g.vcf.gz \\
      --sample_name_child="son_${child_id}" \\
      --sample_name_parent1="mom_${mom_id}" \\
      --postprocess_variants_child_extra_args="--haploid_contigs=\\"chrX\\",--par_regions_bed=\\"${par_bed}\\"" \\
      --vcf_stats_report=true \\
      --num_shards=${params.deepvariant_shards} \\
      --logging_dir=./ \\
      --dry_run=false && \\
        
    /opt/deepvariant/bin/deeptrio/run_deeptrio \\
      --regions="chrY" \\
      --model_type=${params.deepvariant_model} \\
      --ref=${reference} \\
      --reads_child=${child_xam} \\
      --reads_parent1=${dad_xam} \\
      --output_vcf_child=son_${child_id}.Y_incl_PAR.vcf.gz \\
      --output_vcf_parent1=dad_${dad_id}.Y_incl_PAR.vcf.gz \\
      --output_gvcf_child=son_${child_id}.Y_incl_PAR.g.vcf.gz \\
      --output_gvcf_parent1=dad_${dad_id}.Y_incl_PAR.g.vcf.gz \\
      --sample_name_child="son_${child_id}" \\
      --sample_name_parent1="dad_${dad_id}" \\
      --postprocess_variants_child_extra_args="--haploid_contigs=\\"chrY\\",--par_regions_bed=\\"${par_bed}\\"" \\
      --postprocess_variants_parent1_extra_args="--haploid_contigs=\\"chrY\\",--par_regions_bed=\\"${par_bed}\\"" \\
      --vcf_stats_report=true \\
      --num_shards=${params.deepvariant_shards} \\
      --logging_dir=./ \\
      --dry_run=false && \\
    
    # For dad's variants in chrX
    /opt/deepvariant/bin/deeptrio/run_deeptrio \\
      --regions="chrX" \\
      --model_type=${params.deepvariant_model} \\
      --ref=${reference} \\
      --reads_child=${child_xam} \\
      --reads_parent1=${dad_xam} \\
      --reads_parent2=${mom_xam} \\
      --output_vcf_child=son_${child_id}.X_incl_PAR.vcf.disuse.gz \\
      --output_vcf_parent1=dad_${dad_id}.X_incl_PAR.vcf.gz \\
      --output_vcf_parent2=mom_${mom_id}.X_incl_PAR.vcf.disuse.gz \\
      --output_gvcf_child=son_${child_id}.X_incl_PAR.g.vcf.disuse.gz \\
      --output_gvcf_parent1=dad_${dad_id}.X_incl_PAR.called_with_trio.g.vcf.gz \\
      --output_gvcf_parent2=mom_${mom_id}.X_incl_PAR.g.vcf.disuse.gz \\
      --sample_name_child="son_${child_id}" \\
      --sample_name_parent1="dad_${dad_id}" \\
      --sample_name_parent2="mom_${mom_id}" \\
      --postprocess_variants_child_extra_args="--haploid_contigs=\\"chrX,chrY\\",--par_regions_bed=\\"${par_bed}\\"" \\
      --postprocess_variants_parent1_extra_args="--haploid_contigs=\\"chrX,chrY\\",--par_regions_bed=\\"${par_bed}\\"" \\
      --vcf_stats_report=true \\
      --num_shards=${params.deepvariant_shards} \\
      --logging_dir=./ \\
      --dry_run=false && \\
    
    /opt/conda/envs/bio/bin/bcftools view \\
      --targets-file ^${par_bed} \\
      --output-type b \\
      --threads ${params.bcftools_threads} \\
      dad_${dad_id}.X_incl_PAR.vcf.gz > dad_${dad_id}.X_nonPAR.bcf.gz && \\

    /opt/conda/envs/bio/bin/bcftools index dad_${dad_id}.X_nonPAR.bcf.gz
    """
}

process GLNEXUS_FOR_MALE_TRIO {
    publishDir "${params.out_root}/${family_id}/misc/intermediate_files/GLnexus", 
               mode: 'symlink', pattern: '*_merged.bcf*'

    input:
    tuple val(family_id), path(reference), path(reference_index),
          val(child_id), val(dad_id), val(mom_id),
          path(input_AP_gvcfs), path(input_AP_gvcfs_tbi), 
          path(input_X_incl_PAR_gvcfs), path(input_X_incl_PAR_gvcfs_tbis), 
          path(input_Y_incl_PAR_gvcfs), path(input_Y_incl_PAR_gvcfs_tbis), 
          path(dad_X_nonPAR_bcf), path(dad_X_nonPAR_bcf_csi),
          path(par_bed), path(par_bed_index),
          val(_), val(_), val(_)
    
    output:
    tuple val(family_id), path(reference), path(reference_index),
          val(child_id), val(dad_id), val(mom_id),
          path("*.AP.trio_merged.bcf.gz"), path("*.AP.trio_merged.bcf.gz.csi"),
          path("*.X_incl_PAR.nonDad_merged.bcf"),
          path("*.Y_incl_PAR.nonMom_merged.bcf"),
          path(dad_X_nonPAR_bcf), path(dad_X_nonPAR_bcf_csi),
          path(par_bed), path(par_bed_index)

    script:
    """
    # Autosome and PAR
    glnexus_cli \\
      --dir ./GLnexus.DB.AP \\
      --config DeepVariant_unfiltered \\
      --threads ${params.glnexus_threads} \\
      ${input_AP_gvcfs} > ${family_id}.AP.trio_merged.bcf && \\
    
    bcftools view \\
      --samples son_${child_id},mom_${mom_id},dad_${dad_id} \\
      --output-type b \\
      --threads ${params.bcftools_threads} \\
      ${family_id}.AP.trio_merged.bcf > ${family_id}.AP.trio_merged.bcf.gz && \\
    
    bcftools index ${family_id}.AP.trio_merged.bcf.gz && \\
    
    # X chromosome
    glnexus_cli \\
      --dir ./GLnexus.DB.X \\
      --config DeepVariant_unfiltered \\
      --threads ${params.glnexus_threads} \\
      ${input_X_incl_PAR_gvcfs} > ${family_id}.X_incl_PAR.nonDad_merged.bcf && \\
    
    # Y chromosome
    glnexus_cli \\
      --dir ./GLnexus.DB.Y \\
      --config DeepVariant_unfiltered \\
      --threads ${params.glnexus_threads} \\
      ${input_Y_incl_PAR_gvcfs} > ${family_id}.Y_incl_PAR.nonMom_merged.bcf

    """
}

process MERGE_BCFS_FOR_MALE_CHILD_TRIO {
    publishDir "${params.out_root}/${family_id}/merged_bcf", 
               mode: 'symlink', pattern: '*.trio_merged.bcf*'

    input:
    tuple val(family_id), path(reference), path(reference_index),
          val(child_id), val(dad_id), val(mom_id),
          path(AP_trio_merged_bcf), path(AP_trio_merged_bcf_csi),
          path(X_incl_PAR_nonDad_merged_bcf),
          path(Y_incl_PAR_nonMom_merged_bcf),
          path(dad_X_nonPAR_bcf), path(dad_X_nonPAR_bcf_csi),
          path(par_bed), path(par_bed_index)

    output:
    tuple val(family_id), 
          path("${family_id}.trio_merged.bcf.gz")

    script:
    """
    set -eux pipefail

    #######################################
    ##      Preprocess X chromosome      ##
    #######################################

    # Exclude PAR regions in X chromosome
    bcftools view \\
      --targets-file ^${par_bed} \\
      --output-type b \\
      --threads ${params.bcftools_threads} \\
      ${X_incl_PAR_nonDad_merged_bcf} \\
      > ${family_id}.X_nonPAR.nonDad_merged.bcf.gz && \\
    
    bcftools index ${family_id}.X_nonPAR.nonDad_merged.bcf.gz && \\
   
    bcftools merge \\
      --output-type b \\
      --merge none \\
      ${family_id}.X_nonPAR.nonDad_merged.bcf.gz ${dad_X_nonPAR_bcf} \\
      > ${family_id}.X_nonPAR.trio_merged_unorder.bcf.gz && \\

    bcftools index ${family_id}.X_nonPAR.trio_merged_unorder.bcf.gz && \\

    bcftools view \\
        --samples son_${child_id},mom_${mom_id},dad_${dad_id} \\
        --output-type b \\
        --threads ${params.bcftools_threads} \\
        ${family_id}.X_nonPAR.trio_merged_unorder.bcf.gz \\
        > ${family_id}.X_nonPAR.trio_merged.bcf.gz && \\

    bcftools index ${family_id}.X_nonPAR.trio_merged.bcf.gz && \\
    
    #######################################
    ##      Preprocess Y chromosome      ##
    #######################################
  
    # Exclude PAR regions in Y chromosome
    bcftools view \\
      --targets-file ^${par_bed} \\
      --output-type b \\
      --threads ${params.bcftools_threads} \\
      ${Y_incl_PAR_nonMom_merged_bcf} \\
      > ${family_id}.Y_nonPAR.nonMom_merged.bcf.gz && \\

    bcftools index ${family_id}.Y_nonPAR.nonMom_merged.bcf.gz && \\

    # Create empty bcf for concatenation
    bcftools view --header-only ${family_id}.Y_nonPAR.nonMom_merged.bcf.gz \\
      | bcftools view --samples son_${child_id} \\
      | bcftools reheader \\
          --samples <(echo "son_${child_id} mom_${mom_id}") \\
          > mom_header_only.bcf && \\
    
    bgzip mom_header_only.bcf && \\
    bcftools index mom_header_only.bcf.gz && \\

    bcftools merge \\
      --merge none \\
      --output-type b \\
      --threads ${params.bcftools_threads} \\
      ${family_id}.Y_nonPAR.nonMom_merged.bcf.gz mom_header_only.bcf.gz \\
      > ${family_id}.Y_nonPAR.trio_merged_unorder.bcf.gz && \\
    
    bcftools index ${family_id}.Y_nonPAR.trio_merged_unorder.bcf.gz && \\

    bcftools view \\
      --samples son_${child_id},mom_${mom_id},dad_${dad_id} \\
      --output-type b \\
      --threads ${params.bcftools_threads} \\
      ${family_id}.Y_nonPAR.trio_merged_unorder.bcf.gz \\
      > ${family_id}.Y_nonPAR.trio_merged.bcf.gz && \\
    
    bcftools index ${family_id}.Y_nonPAR.trio_merged.bcf.gz && \\
    
    #############################################
    ##  Concatenate three bcfs (AP, X, and Y)  ##
    #############################################

    bcftools concat \\
      --output-type u \\
      ${AP_trio_merged_bcf} \\
      ${family_id}.X_nonPAR.trio_merged.bcf.gz \\
      ${family_id}.Y_nonPAR.trio_merged.bcf.gz \\
      | bcftools sort \\
          --output-type b \\
          > ${family_id}.trio_merged.bcf.gz && \\
    
    bcftools index ${family_id}.trio_merged.bcf.gz

    """
}

process DEEPTRIO_FEMALE_CHILD {
    publishDir "${params.out_root}/${family_id}/raw_vcfs/vcfs", 
                mode: 'symlink', pattern: '*.vcf.gz*'
    publishDir "${params.out_root}/${family_id}/misc/reports/DeepTrioreports", 
                mode: 'move', pattern: '*.html'

    input:
    tuple val(family_id), path(reference), path(reference_index),
          val(child_id), val(dad_id), val(mom_id), 
          path(child_xam), path(dad_xam), path(mom_xam), 
          path(child_xai), path(dad_xai), path(mom_xai),
          val(sex), val(analysis_type),
          path(par_bed), path(par_bed_index)

    output:
    tuple val(family_id), path(reference), path(reference_index),
          val(child_id), val(dad_id), val(mom_id),
          path("*.g.vcf.gz"), path("*.g.vcf.gz.tbi"), 
          path(par_bed), path(par_bed_index),
          path("daughter_${child_id}.vcf.gz*"), 
          path("dad_${dad_id}.vcf.gz*"),
          path("mom_${mom_id}.vcf.gz*"),
          path("*.html")

    script:
    """
    /opt/deepvariant/bin/deeptrio/run_deeptrio \\
      --model_type=${params.deepvariant_model} \\
      --ref=${reference} \\
      --reads_child=${child_xam} \\
      --reads_parent1=${dad_xam} \\
      --reads_parent2=${mom_xam} \\
      --output_vcf_child=daughter_${child_id}.vcf.gz \\
      --output_vcf_parent1=dad_${dad_id}.vcf.gz \\
      --output_vcf_parent2=mom_${mom_id}.vcf.gz \\
      --output_gvcf_child=daughter_${child_id}.g.vcf.gz \\
      --output_gvcf_parent1=dad_${dad_id}.g.vcf.gz \\
      --output_gvcf_parent2=mom_${mom_id}.g.vcf.gz \\
      --sample_name_child="daughter_${child_id}" \\
      --sample_name_parent1="father_${dad_id}" \\
      --sample_name_parent2="mother_${mom_id}" \\
      --postprocess_variants_parent1_extra_args="--haploid_contigs=\\"chrX,chrY\\",--par_regions_bed=\\"${par_bed}\\"" \\
      --vcf_stats_report=true \\
      --num_shards=${params.deepvariant_shards} \\
      --logging_dir=./ \\
      --dry_run=false
    """
}

process GLNEXUS_FOR_FEMALE_TRIO {
    publishDir "${params.out_root}/${family_id}/merged_bcf", 
               mode: 'symlink', pattern: '*.trio_merged.bcf*'

    input:
    tuple val(family_id), path(reference), path(reference_index),
          val(child_id), val(dad_id), val(mom_id),
          path(input_gvcfs), path(input_gvcfs_tbi), 
          path(par_bed), path(par_bed_index),
          val(_), val(_), val(_), val(_)
    
    output:
    tuple val(family_id), 
          path("${family_id}.trio_merged.bcf.gz")

    script:
    """
    glnexus_cli \\
      --config DeepVariant_unfiltered \\
      --threads ${params.glnexus_threads} \\
      ${input_gvcfs} > ${family_id}.trio_merged.bcf && \\

    bcftools view \\
      --samples daughter_${child_id},mom_${mom_id},dad_${dad_id} \\
      --output-type b \\
      --threads ${params.bcftools_threads} \\
      ${family_id}.AP.trio_merged.bcf > ${family_id}.AP.trio_merged.bcf.gz && \\
    
    bcftools index ${family_id}.trio_merged.bcf.gz
    """
}

process DEEPTRIO_FOR_DUO {
    publishDir "${params.out_root}/${family_id}/raw_vcfs", mode: 'symlink', pattern: '*.vcf.gz'
    publishDir "${params.out_root}/${family_id}/raw_vcfs", mode: 'symlink', pattern: '*.vcf.gz.tbi'
    publishDir "${params.out_root}/${family_id}/misc/reports/html", mode: 'move', pattern: '*.html'
    
    input:
    tuple val(family_id), 
          val(child_id), 
          val(parent_id),
          path(child_xam), 
          path(child_xai),
          path(pa_xam), 
          path(pa_xai),
          val(sex), 
          val(analysis_type),
          path(reference),
          path(reference_index),
          path(par_bed),
          val(parent_sex)

    output:
    tuple val(family_id), 
          path("*.g.vcf.gz"), 
          path("*.g.vcf.gz.tbi"), 
          val(analysis_type),
          path("rename_map.tsv")

    script:
    """
    echo "${child_id} son_${child_id}" > rename_map.tsv
    if [ "${parent_sex}" == "xx" ]; then
      echo "${parent_id} mother_${parent_id}" >> rename_map.tsv
    else
      echo "${parent_id} father_${parent_id}" >> rename_map.tsv
    fi

    bash -c "
      /opt/deepvariant/bin/deeptrio/run_deeptrio \\
        --model_type=${params.deepvariant_model} \\
        --ref=${reference} \\
        --reads_child=${child_xam} \\
        --reads_parent1=${pa_xam} \\
        --output_vcf_child=${child_id}.vcf.gz \\
        --output_vcf_parent1=${parent_id}.vcf.gz \\
        --output_gvcf_child=${child_id}.g.vcf.gz \\
        --output_gvcf_parent1=${parent_id}.g.vcf.gz \\
        --sample_name_child=${child_id} \\
        --sample_name_parent1=${parent_id} \\
        --vcf_stats_report=true \\
        --num_shards=${params.deepvariant_shards} \\
        --logging_dir=logs \\
        --dry_run=false    
    "
    """
}

process GLNEXUS_FOR_DUO {
    publishDir "${params.out_root}/${family_id}/raw_vcfs", mode: 'symlink', pettern: '*.bcf'

    input:
    tuple val(family_id), 
          path(input_gvcfs), 
          path(input_gvcfs_index),
          val(analysis_type),
          path(rename_map)

    output:
    tuple val(family_id), 
          path("*.duo_joint.bcf")

    script:
    """
    glnexus_cli \\
      --config DeepVariant_unfiltered \\
      --threads ${params.glnexus_threads} \\
      ${input_gvcfs} > ${family_id}.joint.bcf && \\
    
    bcftools reheader \\
      --sample ${rename_map} \\
      ${family_id}.joint.bcf \\
      | bcftools view -O u -o ${family_id}.duo_joint.bcf
    """
}

process DEEPTRIO_FOR_QUAD {
    publishDir "${params.out_root}/${family_id}/raw_vcfs", mode: 'symlink', pattern: '*.vcf.gz'
    publishDir "${params.out_root}/${family_id}/raw_vcfs", mode: 'symlink', pattern: '*.vcf.gz.tbi'
    publishDir "${params.out_root}/${family_id}/misc/reports/html", mode: 'move', pattern: '*.html'

    input:
    tuple val(family_id), val(analysis_type), path(reference), path(reference_index),
          val(proband_id), val(mom_id), val(dad_id), val(sibling_id),
          path(pro_xam), path(mom_xam), path(dad_xam), path(sib_xam), 
          path(pro_xai), path(mom_xai), path(dad_xai), path(sib_xai),
          val(sex), val(sibling_sex), 
          path(par_bed)

    output:
    tuple val(family_id), val(analysis_type), path(reference), path(reference_index),
          path("*.g.vcf.gz"), path("*.g.vcf.gz.tbi"), path("rename_map.tsv"),
          val(proband_id), val(mom_id), val(dad_id), val(sibling_id),
          path(pro_xam), path(mom_xam), path(dad_xam), path(sib_xam),
          path(pro_xai), path(mom_xai), path(dad_xai), path(sib_xai),
          val(sex), val(sibling_sex),
          path(par_bed)

    script:
    """
    echo "${proband_id} proband_${proband_id}" > rename_map.tsv
    echo "${mom_id} mother_${mom_id}" >> rename_map.tsv
    echo "${dad_id} father_${dad_id}" >> rename_map.tsv
    echo "${sibling_id} sibling_${sibling_id}" >> rename_map.tsv

    bash -c "
      /opt/deepvariant/bin/deeptrio/run_deeptrio \\
        --model_type=${params.deepvariant_model} \\
        --ref=${reference} \\
        --reads_child=${pro_xam} \\
        --reads_parent1=${mom_xam} \\
        --reads_parent2=${dad_xam} \\
        --output_vcf_child=${proband_id}.vcf.gz \\
        --output_vcf_parent1=${mom_id}.vcf.gz \\
        --output_vcf_parent2=${dad_id}.vcf.gz \\
        --output_gvcf_child=${proband_id}.g.vcf.gz \\
        --output_gvcf_parent1=${mom_id}.g.vcf.gz \\
        --output_gvcf_parent2=${dad_id}.g.vcf.gz \\
        --sample_name_child=${proband_id} \\
        --sample_name_parent1=${mom_id} \\
        --sample_name_parent2=${dad_id} \\
        --vcf_stats_report=true \\
        --num_shards=${params.deepvariant_shards} \\
        --dry_run=false
    "
    """
}

process DEEPVARIANT_FOR_SIBLING {
    container 'betelgeuse:5000/library/utsu/deepvariant:1.6.0'
    containerOptions "--security-opt seccomp=unconfined"

    publishDir "${params.out_root}/${family_id}/raw_vcfs", mode: 'symlink', pettern: '*.bcf'
    publishDir "${params.out_root}/${family_id}/raw_vcfs", mode: 'symlink', pettern: '*.vcf.gz'
    publishDir "${params.out_root}/${family_id}/raw_vcfs", mode: 'symlink', pettern: '*.vcf.gz.tbi'
    publishDir "${params.out_root}/${family_id}/misc/reports/html", mode: 'move', pattern: '*.html'

    input:
    tuple val(family_id), val(analysis_type), path(reference), path(reference_index),
          path(gvcfs), path(gvcfs_index), path(rename_map),
          val(child_id), val(mom_id), val(dad_id), val(sibling_id),
          path(pro_xam), path(mom_xam), path(dad_xam), path(sib_xam),
          path(pro_xai), path(mom_xai), path(dad_xai), path(sib_xai),
          val(sex), val(sibling_sex),
          path(par_bed)

    output:
    tuple val(family_id), val(analysis_type),
          path(gvcfs), path(gvcfs_index), path(rename_map),
          path("${sibling_id}.g.vcf.gz"), path("${sibling_id}.g.vcf.gz.tbi"),
          val(child_id), val(mom_id), val(dad_id), val(sibling_id),
          path(pro_xam), path(mom_xam), path(dad_xam), path(sib_xam),
          path(pro_xai), path(mom_xai), path(dad_xai), path(sib_xai),
          val(sex), val(sibling_sex),
          path(reference), path(reference_index),
          path(par_bed)

    script:
    """
    if [ "${sibling_sex}" == "1" ]; then
      bash -c " \\
        /opt/deepvariant/bin/run_deepvariant \\
          --model_type=${params.deepvariant_model} \\
          --ref=${reference} \\
          --reads=${sib_xam} \\
          --output_vcf=${sibling_id}.vcf.gz \\
          --output_gvcf=${sibling_id}.g.vcf.gz \\
          --vcf_stats_report=true \\
          --num_shards=${params.deepvariant_shards} \\
          --sample_name=sibling_${sibling_id} \\
          --haploid_contigs=\"chrX,chrY" \\
          --par_regions_bed=${par_bed} \\
          --dry_run=false && \\
        /opt/conda/envs/bio/bin/bcftools view \\
          -O u -o ${sibling_id}.bcf ${sibling_id}.vcf.gz
      "
    else
      bash -c " \\
        /opt/deepvariant/bin/run_deepvariant \\
          --model_type=${params.deepvariant_model} \\
          --ref=${reference} \\
          --reads=${sib_xam} \\
          --output_vcf=${sibling_id}.vcf.gz \\
          --output_gvcf=${sibling_id}.g.vcf.gz \\
          --sample_name=sibling_${sibling_id} \\
          --vcf_stats_report=true \\
          --num_shards=${params.deepvariant_shards} \\
          --dry_run=false && \\
        /opt/conda/envs/bio/bin/bcftools view \\
          -O u -o ${sibling_id}.bcf ${sibling_id}.vcf.gz
      "
    fi
    """
}


process GLNEXUS_FOR_QUAD {
    containerOptions "--security-opt seccomp=unconfined"

    publishDir "${params.out_root}/${family_id}/raw_vcfs", mode: 'symlink', pettern: '*.bcf'

    input:
    tuple val(family_id), 
          path(input_gvcfs), 
          path(input_gvcfs_index),
          val(analysis_type),
          path(rename_map)

    output:
    tuple val(family_id), 
          path("*.quad_joint.bcf")

    script:
    """
    glnexus_cli \\
      --config DeepVariant_unfiltered \\
      --threads ${params.glnexus_threads} \\
      ${input_gvcfs} > ${family_id}.joint.bcf && \\
    
    bcftools reheader \\
      --sample ${rename_map} \\
      ${family_id}.joint.bcf \\
      | bcftools view -O u -o ${family_id}.quad_joint.bcf
    """
}

process DEEPVARIANT_FOR_OTHER {
    container 'betelgeuse:5000/library/utsu/deepvariant:1.6.0'
    containerOptions "--security-opt seccomp=unconfined"

    publishDir "${params.out_root}/${family_id}/raw_vcfs", mode: 'symlink', pettern: '*.bcf'
    publishDir "${params.out_root}/${family_id}/raw_vcfs", mode: 'symlink', pettern: '*.vcf.gz'
    publishDir "${params.out_root}/${family_id}/raw_vcfs", mode: 'symlink', pettern: '*.vcf.gz.tbi'
    publishDir "${params.out_root}/${family_id}/misc/reports/html", mode: 'move', pattern: '*.html'

    input:
    tuple val(family_id), 
          val(sample_id), 
          path(xam), 
          path(xai),
          val(sex), 
          val(analysis_type),
          path(reference),
          path(reference_index),
          path(par_bed)

    output:
    tuple val(family_id), 
          path("${sample_id}.vcf.gz"), 
          path("${sample_id}.vcf.gz.tbi"), 
          val(analysis_type)

    script:
    """
    echo "TEST"
    """
}

process GLNEXUS_FOR_OTHER {
    container 'betelgeuse:5000/library/utsu/glnexus:1.4.3'
    containerOptions "--security-opt seccomp=unconfined"

    publishDir "${params.out_root}/${family_id}", mode: 'symlink'

    input:
    tuple val(family_id), 
          path(input_gvcfs), 
          path(input_gvcfs_index),
          val(analysis_type),
          path(rename_map)

    output:
    tuple val(family_id), 
          path("*.other_joint.bcf")
      //     path("*.joint.vcf.gz.tbi")

    script:
    """
    glnexus_cli \\
      --config DeepVariant_unfiltered \\
      --threads ${params.glnexus_threads} \\
      ${input_gvcfs} > ${family_id}.other_joint.bcf
    """
}