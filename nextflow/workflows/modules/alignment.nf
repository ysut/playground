// The rule of FASTQ file naming is as follows:
// 1. The FASTQ file name should contain the individual ID
// 2. The FASTQ file name should contain lane information (e.g. L001, L002, etc.)
// 3. The FASTQ file name should contain read information (e.g. R1, R2, etc.)
// 4. The FASTQ file name should contain the file extension (.fastq.gz)
// 5. The FASTQ file name should NOT necessarily contain the family ID
// 6. The FASTQ file name should NOT necessarily contain the sample number information (e.g. S1, S2, etc.)
// Example: Rare_disease_cohort_99999_L003_R2_001.fastq.gz (99999 is the individual ID)

/*
Using docker container with conda environment for the process, 
"bash -c" is used to run the commands.
In the script section, the commands are written in bash.
e.g. 
    script:
    """
    bash -c "source /opt/conda/etc/profile.d/conda.sh ...... "
    """
*/


process STROBEALIGN {
    // container 'betelgeuse:5000/library/utsu/strobealign:0.14.0'
    
    input:
    tuple val(fileName), val(laneID), path(fastq_R1), path(fastq_R2), path(reference)
    
    output:
    tuple val(fileName), val(laneID), path("*.sorted.*am"), path("*.sorted.*ai")

    script:
    """
    bash -c "source /opt/conda/etc/profile.d/conda.sh && \\
    conda activate strobealign && \\
    
    strobealign \\
      --threads=8 ${reference} ${fastq_R1} ${fastq_R2} \\
      --rg-id=tmp_id --rg=SM:tmp_id --rg=LB:mylibrary --rg=PL:Illumina \\
        | samtools sort -o ${fileName}_${laneID}.sorted.bam && \\

    samtools index ${fileName}_${laneID}.sorted.bam"
    """
}

process MERGE_MULTIPLE_LANE_XAMS {
    input:
    tuple val(fileID), val(laneID), path(xams), path(xais)

    output:
    tuple val(fileID), path("*_lane_merged.*am"), path("*_lane_merged.*ai")

    script:
    """
    input_bams="${xams.join(' ')}"
    input_xais="${xais.join(' ')}"
    samtools merge --threads 4 \\
      -o ${fileID}_lane_merged.bam \\
      \${input_bams} && \\

    samtools index ${fileID}_lane_merged.bam
    """
}

process MARKDUP {
    publishDir "${params.out_root}/tmp", mode: 'symlink'
    // publishDir "${params.output}", mode: 'copy', pattern: '*.txt'

    input:
    tuple val(fileID), path(xam), path(xai)

    output:
    tuple val(fileID), path("*.markdup.*am")
    // tuple val(fileID), path("*.marked.*am"), path("*.marked.*ai")

    script:
    """
    /opt/java/openjdk/bin/java -jar /usr/picard/picard.jar MarkDuplicates \\
      -I ${xam} \\
      -O ${fileID}_sorted.markdup.bam \\
      -M ${fileID}_marked_dup_metrics.txt \\
      --CREATE_INDEX true

    # samtools index ${fileID}_sorted.markdup.bam
    """
}

process RENAME_XAM {
    publishDir "${params.out_root}/tmp", mode: 'symlink'

    input:
    val(sample_id)

    output:
    tuple val(sample_id), path("${sample_id}.*am"), path("${sample_id}.*ai")

    script:
    """
    # Find the .cram and .bam files in the temporary directory

    # Rename the .cram files
    cram_files=( "${params.out_root}/tmp/"*${sample_id}*.cram )
    if [ -e "\${cram_files[0]}" ]; then
        mv "\${cram_files[0]}" "${sample_id}.sorted.marked.cram"
    fi

    # Rename the .bam files
    bam_files=( "${params.out_root}/tmp/"*${sample_id}*.bam )
    echo "bam_files: \${bam_files}"
    if [ -e "\${bam_files[0]}" ]; then
        mv "\${bam_files[0]}" "${sample_id}.sorted.marked.bam"
    fi

    # Rename the .crai files
    crai_files=( "${params.out_root}/tmp/"*${sample_id}*.crai )
    if [ -e "\${crai_files[0]}" ]; then
        mv "\${crai_files[0]}" "${sample_id}.sorted.marked.cram.crai"
    fi

    # Rename the .bai files
    bai_files=( "${params.out_root}/tmp/"*${sample_id}*.bai )
    if [ -e "\${bai_files[0]}" ]; then
        mv "\${bai_files[0]}" "${sample_id}.sorted.marked.bam.bai"
    fi
    """
}


process EDIT_RG {
    publishDir "${params.out_root}/xams", mode: 'symlink'

    input:
    tuple val(sample_id), path(xam), path(xai)

    output:
    tuple val(sample_id), path("${sample_id}.*am"), path("${sample_id}.*ai")

    script:
    """
    touch ${sample_id}.bam
    touch ${sample_id}.bam.bai
    """
}



process STROBEALIGN2 {
    container 'betelgeuse:5000/library/utsu/strobealign:0.14.0'
    
    input:
    tuple path(read1), path(read2), path(reference), val(sample_id)
    
    output:
    tuple path('sorted.bam'), path('sorted.bam.bai'), val(sample_id)

    script:
    """
    conda activate strobealign && \\
    strobealign \\
      --threads=8 ${reference} ${read1} ${read2} 
      --rg-id=tmp_id --rg=SM:tmp_id --rg=LB:mylibrary --rg=PL:Illumina
        | samtools sort -o sorted.bam && \\
    samtools index sorted.bam
    """
}