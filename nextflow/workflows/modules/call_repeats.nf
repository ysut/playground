process EXPANSIONHUNTER {
    container "betelgeuse:5000/library/utsu/expansionhunter:5.0.0"
    containerOptions "--security-opt seccomp=unconfined"

    publishDir "${params.out_root}/${family_id}/${sample_id}/repeat", mode: 'symlink'

    input:
    tuple val(family_id), val(sample_id), path(reference), path(reference_index),
          path(xam), path(xai),
          val(xxxy), 
          path(par_bed),
          val(anlysis_type)

    output:
    tuple val(family_id), 
          val(sample_id), 
          path("*.eh_realigned.bam"), 
      //     path("*.eh_realigned.bam.bai"),
          path("*.vcf"),
          path("*.json")
        
    script:
    """
    ExpansionHunter \\
      --reads ${xam} \\
      --reference ${reference} \\
      --variant-catalog /RepeatCatalogs/hg38/variant_catalog.json \\
      --output-prefix ${sample_id}.eh \\
      --sex ${(xxxy == "1") ? "male" : "female"} \\
      --analysis-mode ${params.eh_mode} \\
      --region-extension-length ${params.eh_region_extension_length} \\
      --threads ${params.eh_threads}
    """
}

process INDEX_EH_BAM {
    publishDir "${params.out_root}/${family_id}/${sample_id}/repeat", mode: 'symlink', pattern: "*.bam.bai"

    input:
    tuple val(family_id), val(sample_id), 
          path(eh_bam),
          path("*.vcf"),
          path("*.json")

    output:
    tuple val(family_id),
          val(sample_id), 
          path("*.eh_realigned.bam"), 
          path("*.eh_realigned.bam.bai")
    
    script:
    """
    samtools index ${eh_bam}
    """
}

process REVIEWER {
      publishDir "${params.out_root}/${family_id}/${sample_id}/repeat", mode: 'symlink'
      
      input:
      tuple val(family_id), val(sample_id), 
              path(xam), 
              path(xai),
              path(fasta)
      
      output:
      tuple val(family_id), val(sample_id), 
              path("*.reviewer.bam"), 
              path("*.reviewer.bam.bai")
      
      script:
      """
      reviewer \\
        --reads ${xai} \\
        --reference ${fasta} \\
        --output-prefix ${sample_id}.reviewer
      """
}