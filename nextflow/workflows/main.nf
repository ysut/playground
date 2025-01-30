#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

params.fastq_dir = '/Users/utsu/work/Github/playground/nextflow/myfiles/fastq'
params.ped_file = ''

params.output = params.output ?: "${workflow.launchDir}"
params.out_root = "${params.output}/WES-WF_Results_" + new Date().format('yyyyMMdd-HHmmss')


// === Include modules ===

include { 
    STROBEALIGN; MERGE_MULTIPLE_LANE_XAMS; MARKDUP; RENAME_XAM; EDIT_RG 
    } from './modules/alignment.nf'
include { 
    DEEPVARIANT; DEEPTRIO; DEEPTRIO_FOR_QUAD; DEEPVARIANT_FOR_SIBLING; GLNEXUS
    } from './modules/call_snvindels.nf'
// include { 
//     EXPANSIONHUNTER 
//     } from './modules/call_repeat.nf'

// Debugging includes (suffix 2)
// include { 
//     STROBEALIGN; MERGE_MULTIPLE_LANE_XAMS; MARKDUP; RENAME_XAM; EDIT_RG 
//     } from './modules/alignment2.nf'
// include { 
//     DEEPVARIANT; DEEPTRIO; DEEPTRIO_FOR_QUAD; DEEPVARIANT_FOR_SIBLING; GLNEXUS
//     } from './modules/call_snvindels2.nf'
// include { 
//     EXPANSIONHUNTER 
//     } from './modules/call_repeat2.nf'



// === Processes ===
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

// 1=male, 2=female, 0/-9=unknown/missing
// 1=unaffected, 2=affected, 0/-9=missing

// === Main workflow ===
workflow {
    // STEP 0. Check the input parameters
    log.info """\
    This is
    a multiline
    message
    """
    
    if (!params.ped_file) {
    error "Enter the path to the pedigree file using the --ped_file parameter."
    }

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
        | flatMap { familyTuple ->
                // familyTupleの構造は [familyId, [member1, member2, ...]]
                def members = familyTuple[1]
                // 各メンバーのindividualIdを収集
                members.collect { member ->
                    member.individualId
                }
            }
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
        | filter { it -> it[1].size() > 1 }
        | MERGE_MULTIPLE_LANE_XAMS
        | concat(singleBAMs)
        | MARKDUP
        | collect
        | set { markdup_done }
    
    markdup_done
        | map { true }
        | combine(sampleIDs)
        | map { _, sampleID -> sampleID }
        | RENAME_XAM
        | EDIT_RG
        | collect
        | set { align_process_done }


    // STEP 4: Call SNVs and indels
    pedigree_info
        | map { familyTuple ->
                def familyId = familyTuple[0]
                def members = familyTuple[1].collect { 
                    [
                        individualId: it.individualId,
                        fatherId    : it.fatherId,
                        motherId    : it.motherId,
                        sex         : it.sex,
                        phenotype   : it.phenotype
                    ]
                }
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
        | set { family_ch }

    /*
    一人だけdeepvariantに入力する場合は，singleton_chに入れる．
    deepvariantでcallして，vcf作る．
    */
    // Generate the singleton channel
    align_process_done
        | map { true }
        | combine(family_ch)
        | map { _, family -> family }
        | filter { it.analysisType == 'singleton' }
        | map { family -> 
                return [
                    family_id    : family.familyId,
                    sample_id    : family.members[0].individualId,
                    analysisType : family.analysisType
                ]
            }
        | DEEPVARIANT
    
    /*
    複数人（2-4人）の場合は，family_chに入れる．
    == family_chのサンプル抽出基準 ==
    [duo]  - proband は，fatherId または motherId が 0 でない．
           - parent は，proband の fatherId または motherId から抽出．
    
    [trio]  - proband は，fatherId と motherId が 0 でない
            - mother は，proband の motherId から抽出．
            - father は，proband の fatherId から抽出．
    
    [quad]  - proband は，fatherId と motherId が 0 でない
    　　　    かつ，phenotype が 2 で，複数人いたら，IDが若い方 を抽出．
            - mother は，proband の motherId から抽出．
            - father は，proband の fatherId から抽出．
            - sibling は，fatherId と motherId が 0 でない
              かつ proband 以外の候補を抽出．sibling_idは，DeepVariantへ．
              あとで，GLNEXUSでジョイント．
    */

    // Generate the family channel
    align_process_done
        | map { true }
        | combine(family_ch)
        | map { _, family -> family }
        | filter { it.analysisType != 'singleton' && it.analysisType != 'quad' }
        | map { family -> 
                def proband = null
                switch (family.analysisType) {
                    case 'trio':
                        proband = family.members.find {
                            (it.motherId as Integer) != 0 
                            && (it.phenotype as Integer) == 2
                        }
                        break
                    case 'duo':
                        proband = family.members.find { 
                            (it.fatherId as Integer) != 0 
                            || (it.motherId as Integer) != 0 
                        }
                        break
                    default:
                        println "Unknown analysisType: ${family.analysisType}"
                }

                // proband が見つからなかった場合のハンドリング
                if (!proband) {
                    println "Proband not found for Family ID: ${family.familyId}"
                    return null
                }

                // proband が見つかった場合、親を探す
                def mother = family.members.find { 
                    (it.individualId as Integer) == (proband.motherId as Integer) 
                }
                def father = family.members.find { 
                    (it.individualId as Integer) == (proband.fatherId as Integer) 
                }

                return [
                    family_id    : family.familyId,
                    proband_id   : proband.individualId,
                    mother_id    : mother ? mother.individualId : null,
                    father_id    : father ? father.individualId : null,
                    analysisType : family.analysisType
                ]
            }
        | DEEPTRIO
        | set { variant_call_results }

    align_process_done
        | map { true }
        | combine(family_ch)
        | map { _, family -> family }
        | filter { it.analysisType == 'quad' }
        | map { family ->
                def proband = null
                proband_candidate = family.members.findAll {
                    (it.fatherId as Integer) != 0 
                    && (it.motherId as Integer) != 0 
                    && (it.phenotype as Integer) == 2
                    }
                
                if (proband_candidate.size() == 2) {
                    proband = proband_candidate.min { it.individualId as Integer }
                    sibling = family.proband_candidate.find {
                        (it.individualId as Integer) != (proband.individualId as Integer)
                    }
                } else {
                    proband = proband_candidate[0]
                    sibling = family.members.find {
                        (it.fatherId as Integer) != 0 
                        && (it.motherId as Integer) != 0 
                        &&  (it.individualId as Integer) != (proband.individualId as Integer)
                    }
                }
                
                // proband が見つからなかった場合のハンドリング
                if (!proband) {
                    println "Proband not found for Family ID: ${family.familyId}"
                    return null
                }

                // proband が見つかった場合、親を探す
                def mother = family.members.find { 
                    (it.individualId as Integer) == (proband.motherId as Integer) 
                }
                def father = family.members.find { 
                    (it.individualId as Integer) == (proband.fatherId as Integer) 
                }

                return [
                    family_id    : family.familyId,
                    proband_id   : proband.individualId,
                    mother_id    : mother ? mother.individualId : null,
                    father_id    : father ? father.individualId : null,
                    sibling_id   : sibling ? sibling.individualId : null,
                    analysisType : family.analysisType
                ]
            }
        | DEEPTRIO_FOR_QUAD
        | DEEPVARIANT_FOR_SIBLING
        | map { quad_dv_results -> 
                def familyId = quad_dv_results[0]
                def gvcf_paths = quad_dv_results[1] + [quad_dv_results[3]]
                def tbi_paths = quad_dv_results[2] + [quad_dv_results[4]]
                tuple (familyId, gvcf_paths, tbi_paths)
            }
        | set { variant_call_for_quad_results }

    variant_call_results.concat(variant_call_for_quad_results)
        | GLNEXUS
        | view
    

}
