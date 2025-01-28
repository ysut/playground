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

// process FIND_FASTQ {
//     // The rule of FASTQ file naming is as follows:
//     // 1. The FASTQ file name should contain the individual ID
//     // 2. The FASTQ file name should contain lane information (e.g. L001, L002, etc.)
//     // 3. The FASTQ file name should contain read information (e.g. R1, R2, etc.)
//     // 4. The FASTQ file name should contain the file extension (.fastq.gz)
//     // 5. The FASTQ file name should NOT necessarily contain the family ID
//     // 6. The FASTQ file name should NOT necessarily contain the sample number information (e.g. S1, S2, etc.)
//     // Example: Rare_disease_cohort_99999_L003_R2_001.fastq.gz (99999 is the individual ID)

//     publishDir "${params.out_root}/fastq", mode: 'symlink'

//     input:
//     tuple val(familyId), val(individualId)

//     output:
//     tuple val(familyId), val(individualId), 
//           val(laneID), path("*_R1*.fastq.gz"), path("*_R2*.fastq.gz")

//     script:
//     """
//     echo "Searching FASTQ for Family: ${familyId}, Individual: ${individualId}"

//     set -e

//     cd ${params.out_root}/fastq

//     # Find all R1 files for the individual
//     R1_FILES=(\$(find "${params.fastq_dir}" -type f -name "*${individualId}*_R1*.fastq.gz"))

//     for R1 in "\${R1_FILES[@]}"; do
//         # Extract laneID from the filename
//         LANE=\$(echo \$(basename "\$R1") | grep -oP '(?<=_L)\d{3}(?=_)')

//         # Find the corresponding R2 file
//         R2=\$(find "${params.fastq_dir}" -type f -name "*${individualId}*_L\${LANE}_R2*.fastq.gz")

//         if [[ -f "\$R2" ]]; then
//             # Create symlinks with lane information
//             ln -sf "\$R1" "${individualId}_L\${LANE}_R1.fastq.gz"
//             ln -sf "\$R2" "${individualId}_L\${LANE}_R2.fastq.gz"

//             # Output the tuple
//             echo "\${familyId} \${individualId} \${LANE} \${individualId}_L\${LANE}_R1.fastq.gz \${individualId}_L\${LANE}_R2.fastq.gz"
//         else
//             echo "R2 file not found for \${R1}" >&2
//         fi
//     done
//     """
// }

// process STROBEALIGN {
//     input:
//     tuple val(familyID), val(individualID), path(fastq_R1), path(fastq_R2)

//     output:
//     tuple val(familyID), val(individualID), path("*.{b,cr}am"), path("*.*am.{b,cr}ai")

//     script:
//     """ 
//     touch ${individualID}.bam
//     touch ${individualID}.bam.bai
//     """
// }

process MARKDUP {
    publishDir "${params.out_root}/${familyID}/xams", mode: 'symlink'

    input:
    tuple val(familyID), val(individualID), path(bam), path(bai)

    output:
    tuple val(familyID), val(individualID), path("*.marked.bam"), path("*.marked.bam.bai")

    script:
    """
    touch ${individualID}.marked.bam
    touch ${individualID}.marked.bam.bai
    touch ${individualID}.marked_dup_metrics.txt
    """
}

// process DEEPVARIANT {
//     input:
//     path sequences

//     input:
//     tuple val(familyID), val(individualID), path(bam), path(bai)

//     script:
//     if( mode == 'tcoffee' )        
//         t_coffee -in $sequences > out_file

//     else if( mode == 'mafft' )
//         mafft --anysymbol --parttree --quiet $sequences > out_file

//     else if( mode == 'clustalo' )
//         clustalo -i $sequences -o out_file

//     else
//         error "Invalid alignment mode: ${mode}"

// }




// process DEEPTRIO {
//     input:
//     tuple val(familyID), val(individualID), path(bam), path(bai)

//     output:
//     tuple val(familyID), val(individualID), path("*.deeptrio.vcf.gz"), path("*.deeptrio.vcf.gz.tbi")

//     script:
//     """
// sudo docker run \
//   -v "${INPUT_DIR}":"/input" \
//   -v "${OUTPUT_DIR}":"/output" \
//   google/deepvariant:deeptrio-"${BIN_VERSION}" \
//   /opt/deepvariant/bin/deeptrio/run_deeptrio \
//   --model_type=WGS \
//   --ref=/input/GRCh38_no_alt_analysis_set.fasta \
//   --reads_child=/input/HG002.chr20.10_10p1mb.bam \
//   --reads_parent1=/input/HG003.chr20.10_10p1mb.bam \
//   --reads_parent2=/input/HG004.chr20.10_10p1mb.bam \
//   --output_vcf_child /output/HG002.output.vcf.gz \
//   --output_vcf_parent1 /output/HG003.output.vcf.gz \
//   --output_vcf_parent2 /output/HG004.output.vcf.gz \
//   --sample_name_child 'HG002' \
//   --sample_name_parent1 'HG003' \
//   --sample_name_parent2 'HG004' \
//   --num_shards $(nproc)  \
//   --regions "chr20:10,000,000-10,010,000" \
//   --intermediate_results_dir /output/intermediate_results_dir \
//   --output_gvcf_child /output/HG002.g.vcf.gz \
//   --output_gvcf_parent1 /output/HG003.g.vcf.gz \
//   --output_gvcf_parent2 /output/HG004.g.vcf.gz
//     """
// }


// 1=male, 2=female, 0/-9=unknown/missing
// 1=unaffected, 2=affected, 0/-9=missing

// process GLNEXUS {

//     input:
//     tuple val(familyID), val(individualID), path(bam), path(bai)

//     output:
//     tuple val(familyID), val(individualID), path("*.glnexus.vcf.gz"), path("*.glnexus.vcf.gz.tbi")

//     script:
//     """
//     touch ${individualID}.glnexus.vcf.gz
//     touch ${individualID}.glnexus.vcf.gz.tbi
//     """
// }


// process GCNV {
//     input:
//     tuple val(familyID), val(individualID), path(bam), path(bai)

//     output:
//     tuple val(familyID), val(individualID), path("*.gcnv.vcf.gz"), path("*.gcnv.vcf.gz.tbi")

//     script:
//     """
//     touch ${individualID}.gcnv.vcf.gz
//     touch ${individualID}.gcnv.vcf.gz.tbi
//     """
// }

// process AUTOMAP {
//     input:
//     tuple val(familyID), val(individualID), path(bam), path(bai)

//     output:
//     tuple val(familyID), val(individualID), path("*.automap.vcf.gz"), path("*.automap.vcf.gz.tbi")

//     script:
//     """
//     touch ${individualID}.automap.vcf.gz
//     touch ${individualID}.automap.vcf.gz.tbi
//     """
// }



workflow {
    /// Output directory    ========================================
    MAKE_OUTPUT_ROOT_DIR()

    // Read pedigree file   ========================================
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
    

    // =========================================================================
    // Create BAMs and output directories for each Family and Sample
    pedigree_info
        | flatMap { familyId, individuals ->
            individuals.collect { indi -> tuple(familyId, indi.individualId)}
            }
        | MAKE_FAMILY_SAMPLE_DIR
    //     | FIND_FASTQ
    //     | STROBEALIGN
    //     | MARKDUP
    

    Channel.fromFilePairs(
        "${params.fastq_dir}/**/*_L{lane}_R{read}*.fastq.gz", flat: true, size: 2)
        | map { pair ->
                def file1 = pair[0]
                def file2 = pair[1]

                // ファイル名からサンプルIDとレーンIDを抽出
                def filename1 = file1.getName()
                // def filename2 = file2.getName()

                def sampleMatch = filename1 =~ /(S\d+)/
                def laneMatch = filename1 =~ /_L(\d{3})_R\d+/

                def sampleID = sampleMatch ? sampleMatch[0][1] : 'unknown_sample'
                def laneID = laneMatch ? laneMatch[0][1] : 'unknown_lane'

                tuple(sampleID, laneID, file1, file2)
            }
        | view
        // | set { fastq_pairs }

    // =========================================================================
    // Generate VCF 
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
            family.members.collect { indi ->
                tuple(
                    family.familyId,
                    indi.individualId,
                    indi.fatherId,
                    indi.motherId,
                    indi.sex,
                    indi.phenotype,
                    family.analysisType  // ここでanalysisTypeも含める
                )
            }
        }
        | view
    
}




// [
//     familyId:FAM_trio, 
//     members:[
//         [individualId:32975, fatherId:32976, motherId:32977, sex:1, phenotype:2], 
//         [individualId:32976, fatherId:0, motherId:0, sex:1, phenotype:1], 
//         [individualId:32977, fatherId:0, motherId:0, sex:2, phenotype:1]
//         ],
//     analysisType:trio
// ]


// [familyId:FAM_quad, members:[[individualId:40001, fatherId:40002, motherId:40003, sex:2, phenotype:2
// ], [individualId:40002, fatherId:40002, motherId:40003, sex:1, phenotype:1
// ], [individualId:40003, fatherId:0, motherId:0, sex:1, phenotype:1
// ], [individualId:40004, fatherId:0, motherId:0, sex:2, phenotype:1
// ]], analysisType:quad]


// [familyId:FAM_pro, members:[[individualId:99999, fatherId:0, motherId:0, sex:0, phenotype:2
// ]], analysisType:singleton]


// [
//     familyId:FAM_duo, 
//     members:[
//         [individualId:32800, fatherId:0, motherId:32801, sex:1, phenotype:2], 
//         [individualId:32801, fatherId:0, motherId:0, sex:2, phenotype:1]
//         ], 
//     analysisType:duo
// ]

