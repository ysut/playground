#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

params.fastq_dir = '/Users/utsu/work/Github/playground/nextflow/myfiles/fastq'
params.ped_file = ''

if (!params.ped_file) {
    error "Enter the path to the pedigree file using the --ped_file parameter."
}

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
    mkdir -p ${params.out_root}/xams
    mkdir -p ${params.out_root}/tmp
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
    mkdir -p ${params.out_root}/${familyID}/raw_VCFs
    """
}

//     // The rule of FASTQ file naming is as follows:
//     // 1. The FASTQ file name should contain the individual ID
//     // 2. The FASTQ file name should contain lane information (e.g. L001, L002, etc.)
//     // 3. The FASTQ file name should contain read information (e.g. R1, R2, etc.)
//     // 4. The FASTQ file name should contain the file extension (.fastq.gz)
//     // 5. The FASTQ file name should NOT necessarily contain the family ID
//     // 6. The FASTQ file name should NOT necessarily contain the sample number information (e.g. S1, S2, etc.)
//     // Example: Rare_disease_cohort_99999_L003_R2_001.fastq.gz (99999 is the individual ID)


process STROBEALIGN {
    input:
    tuple val(fileName), val(laneID), path(fastq_R1), path(fastq_R2)

    output:
    tuple val(fileName), val(laneID), path("*.{b,cr}am"), path("*.*am.{b,cr}ai")

    script:
    """ 
    touch ${fileName}_${laneID}.bam
    touch ${fileName}_${laneID}.bam.bai
    """
}

process MERGE_MULTIPLE_LANE_BAMS {
    input:
    tuple val(fileID), val(laneID), path(bams), path(bais)

    output:
    tuple val(fileID), path("*_lane_merged.bam"), path("*_lane_merged.bam.bai")

    script:
    """
    # echo "${bams.join(' ')}"
    touch ${fileID}_lane_merged.bam
    touch ${fileID}_lane_merged.bam.bai
    """
}


process MARKDUP {
    publishDir "${params.out_root}/tmp", mode: 'symlink'

    input:
    tuple val(fileID), path(bam), path(bai)

    output:
    tuple val(fileID), path("*.marked.bam"), path("*.marked.bam.bai")

    script:
    """
    touch ${fileID}.marked.bam
    touch ${fileID}.marked.bam.bai
    touch ${fileID}.marked_dup_metrics.txt
    """
}

process FIND_AND_RENAME_BAM {
    input:
    val(sample_id)

    output:
    tuple val(sample_id), path("*.rg.bam"), path("*.rg.bam.bai")

    script:
    """
    find "${params.out_root}/tmp" -name "*${sample_id}*.bam" \\
        -exec mv {} "${sample_id}.rg.bam" \\; &&
    find "${params.out_root}/tmp" -name "*${sample_id}*.bai" \\
        -exec mv {} "${sample_id}.rg.bam.bai" \\;
    """
}


process EDIT_READGROUP {
    publishDir "${params.out_root}/xams", mode: 'symlink'

    input:
    tuple val(sample_id), path(bam), path(bai)

    output:
    tuple val(sample_id), path("*.rg.bam"), path("*.rg.bam.bai")

    script:
    """
    touch ${sample_id}.rg.bam
    touch ${sample_id}.rg.bam.bai
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

// 1=male, 2=female, 0/-9=unknown/missing
// 1=unaffected, 2=affected, 0/-9=missing

workflow {
    // STEP 1: Create the root output directory
    MAKE_OUTPUT_ROOT_DIR()

    // STEP 2: Read the input pedigree file and create directories for each family and sample
    /* 
    Pedigree fileの読み込み．
    pedigree_infoとしてまとめておく 
    */
    Channel.fromPath("${params.ped_file}")
        | splitText
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

    pedigree_info
        | map { familyTuple ->
                def familyId = familyTuple[0]
                def members = familyTuple[1..-1].flatten()
                def memberCount = members.size()
                def analysisType = ''

                if (memberCount == 3) {
                    analysisType = 'trio'
                } else if (memberCount == 2) {
                    analysisType = 'duo'
                } else if (memberCount == 4) {
                    analysisType = 'quad'
                } else {
                    analysisType = 'singleton'
                }
                return [familyId: familyId, members: members, analysisType: analysisType]
                }
        | flatMap { family -> 
                    family.members.collect {
                        indi -> tuple(
                            family.familyId,
                            indi.individualId,
                            indi.fatherId,
                            indi.motherId,
                            indi.sex,
                            indi.phenotype,
                            family.analysisType
                            )
                        }
                    }
        | map { indiTuple -> indiTuple[1] }
        | set { sampleIDs }

    /*
    pedigree_infoを使って，familyIdとindividualId毎のディレクトリを作成する．
    この結果のディレクトリ構造は以下の通り．
    out_root/
        xams/
        familyId/
            individualId/
    */
    pedigree_info
        | flatMap { familyId, individuals ->
            individuals.collect { indi -> tuple(familyId, indi.individualId)}
            }
        | MAKE_FAMILY_SAMPLE_DIR


    // STEP 3: Align the FASTQ files
    /*
    fastq_dir内の*.fastq.gzを探してきて，アライメントする．
    
    1. fastq_dir内のファイル名からBAMを探索してくる．
    ファイル名からサンプルIDとレーンIDを抽出する．サンプルによっては，
    複数ペアのfastqがある場合もある．
    
    2. 1で集めたfastqからBAMを生成する．
    複数ペアのfastqがある場合は，BAMを作った後にmergeが必要．
    生成したBAMは，各Family/Sampleディレクトリに保存する．

    3. 2で生成したBAMに対して，MarkDuplicatesを実行する．
    MarkDuplicatesで生成したBAMは，各Family/Sampleディレクトリに保存する．
    */
    Channel.fromFilePairs("${params.fastq_dir}/*_R{1,2}*.fastq.gz", 
                            flat: true, size: 2)
        | map { pair ->
                def file1 = pair[1]
                def file2 = pair[2]

                def filename1 = file1.getName()

                def sampleMatch = filename1 =~ /^(.*?)_L\d+_R[12]/
                def fileID = sampleMatch ? sampleMatch[0][1] : pair[0]

                def laneMatch = filename1 =~ /_(L\d+)_R\d+/
                def laneID = laneMatch ? laneMatch[0][1] : 'unknown'

                tuple(fileID, laneID, file1, file2)
            }
        | STROBEALIGN
        | groupTuple()  // Group by fileID. Then, merge BAM files if there are multiple lanes.
        | set { groupedBAMs }

    groupedBAMs
        | filter { fileID, laneID, bam_paths, bai_paths -> laneID.size() == 1 }
        | map { fileID, laneID, bam_paths, bai_paths -> tuple(fileID, bam_paths[0], bai_paths[0]) }
        | set { singleBAMs }

    groupedBAMs
        | filter { fileID, laneID, bam_paths, bai_paths -> laneID.size() > 1 }
        | MERGE_MULTIPLE_LANE_BAMS
        | concat(singleBAMs)
        | MARKDUP
        | collect
        | set { markdup_done }
    
    markdup_done
        | map { true }
        | combine(sampleIDs)
        | map { _, sampleID -> sampleID }
        | FIND_AND_RENAME_BAM
        | view

    // Generate VCF 
        // | set { flat_pedigree_info }
}
