// The rule of FASTQ file naming is as follows:
// 1. The FASTQ file name should contain the individual ID
// 2. The FASTQ file name should contain lane information (e.g. L001, L002, etc.)
// 3. The FASTQ file name should contain read information (e.g. R1, R2, etc.)
// 4. The FASTQ file name should contain the file extension (.fastq.gz)
// 5. The FASTQ file name should NOT necessarily contain the family ID
// 6. The FASTQ file name should NOT necessarily contain the sample number information (e.g. S1, S2, etc.)
// Example: Rare_disease_cohort_99999_L003_R2_001.fastq.gz

/*
Using docker container with conda environment for the process, 
"bash -c" is used to run the commands.
e.g. 
    script:
    """
    bash -c "source /opt/conda/etc/profile.d/conda.sh ...... "
    """
*/

process STROBEALIGN_INDEX {
    input:
    tuple val(reference)

    output:
    tuple val(reference), path("*.fai")

    script:
    """
    bash -c " \\
    source /opt/conda/etc/profile.d/conda.sh && \\
    conda activate strobealign && \\

    strobealign index ${reference} \\
    
    samtools faidx ${reference}
    """
}



process STROBEALIGN {
    input:
    tuple val(fileID), val(key), path(reference), path(reference_index),
          path(fastq_R1), path(fastq_R2)

    output:
    tuple val(fileID), val(key), path("*.sorted.bam"), path("*.sorted.bam.bai")

    script:
    """
    bash -c " \\
    source /opt/conda/etc/profile.d/conda.sh && \\
    conda activate strobealign && \\
    
    strobealign \\
      ${reference} ${fastq_R1} ${fastq_R2} \\
      --threads=${params.strobealign_threads} \\
      --rg-id=tmp_id --rg=SM:tmp_id \\
      --rg=LB:${params.rg_library} --rg=PL:${params.rg_platform} \\
      | samtools sort -@ ${params.samtools_threads} -o ${key}.sorted.bam && \\

    samtools index ${key}.sorted.bam \\
    "
    """
}

process MERGE_MULTIPLE_LANE_XAMS {
    input:
    tuple val(fileID), path(xams), path(xais)

    output:
    tuple val(fileID), path("*_merged.bam"), path("*_merged.bam.bai")

    script:
    """
    input_xams="${xams.join(' ')}"
    input_xais="${xais.join(' ')}"
    samtools merge \\
      --threads ${params.samtools_threads} \\
      -o ${fileID}_merged.bam \\
      \${input_xams} && \\

    samtools index ${fileID}_merged.bam
    """
}

process MARKDUP {
    publishDir "${params.out_root}/reports/MarkdupMetrics", 
               mode: 'move', pattern: '*.marked_dup_metrics.txt'
    publishDir "${params.out_root}/intermediate_files/MarkDuplicates", 
               mode: 'symlink', pattern: '*.bam'
    publishDir "${params.out_root}/intermediate_files/MarkDuplicates", 
               mode: 'symlink', pattern: '*.bai'

    input:
    tuple val(fileID), path(xam), path(xai)

    output:
    tuple val(fileID), 
          path("*.sorted.marked.bam"), path("*.sorted.marked.bai"), 
          path("*.marked_dup_metrics.txt")

    script:
    """
    /opt/java/openjdk/bin/java -jar /usr/picard/picard.jar \\
      MarkDuplicates \\
      -I ${xam} \\
      -O ${fileID}.sorted.marked.bam \\
      -M ${fileID}.marked_dup_metrics.txt \\
      --CREATE_INDEX true
    """
}

process RENAME_XAM {
    publishDir "${params.out_root}/intermediate_files", mode: 'symlink'

    input:
    val(sample_id)

    output:
    tuple val(sample_id),
          path("${sample_id}.sorted.marked.bam"),
          path("${sample_id}.sorted.marked.bam.bai")

    script:
    """
    tmp_dir="${params.out_root}/intermediate_files/MarkDuplicates"

    for ext in bam bai; do
      files=( "\${tmp_dir}"/*${sample_id}*.\${ext} )
      if [ -e "\${files[0]}" ]; then
        case "\${ext}" in
          bam)
            new_name="${sample_id}.sorted.marked.bam"
            ;;
          bai)
            new_name="${sample_id}.sorted.marked.bam.bai"
            ;;
        esac
        mv "\${files[0]}" "\${new_name}"
      fi
    done
    """
}

process EDIT_RG {
    input:
    tuple val(sample_id), path(xam), path(xai)

    output:
    tuple val(sample_id), path("${sample_id}.processed.bam")

    script:
    """ 
    /opt/java/openjdk/bin/java -jar /usr/picard/picard.jar \\
      AddOrReplaceReadGroups \\
      -I ${xam} \\
      -O ${sample_id}.processed.bam \\
      -ID ${sample_id} \\
      -LB ${params.rg_library} \\
      -PL ${params.rg_platform} \\
      -PU unit1 \\
      -SM 20 \\
      --CREATE_INDEX false
    """
}

process INDEX_XAM {
    publishDir "${params.out_root}/xams", mode: 'symlink'
    
    input:
    tuple val(sample_id), path(processed_xam)

    output:
    tuple val(sample_id), path(processed_xam), path("*.*ai")

    script:
    """
    samtools index -@ ${params.samtools_threads} ${processed_xam}
    """
}