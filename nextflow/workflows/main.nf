#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

params.fastq_dir = ''
params.output = params.output ?: "${workflow.launchDir}"
params.out_root = "${params.output}/WES_Results_" + new Date().format('yyyyMMdd-HHmmss')
params.bam_dir = params.bam_dir ?: "${params.out_root}/xams"
params.ped_file = params.ped_file ?: "${workflow.launchDir}/*.ped"
params.use_deepvariant = params.use_deepvariant ?: false
params.relationship = params.relationship ?: true

// === Include modules ===
include { 
    STROBEALIGN; MERGE_MULTIPLE_LANE_XAMS; 
    MARKDUP; RENAME_XAM; EDIT_RG; INDEX_XAM
    } from './modules/alignment.nf'

include {
    DEEPVARIANT as DEEPVARIANT_SINGLETON;
    DEEPVARIANT as DEEPVARIANT_OTHERS; 
    GLNEXUS as GLNEXUS_OTHERS;
    DEEPTRIO_MALE_CHILD;
    GLNEXUS_FOR_MALE_TRIO; 
    MERGE_BCFS_FOR_MALE_CHILD_TRIO;
    DEEPTRIO_FEMALE_CHILD; 
    GLNEXUS_FOR_FEMALE_TRIO; 
    DEEPTRIO_FOR_DUO; GLNEXUS_FOR_DUO;
    DEEPTRIO_FOR_QUAD; DEEPVARIANT_FOR_SIBLING; GLNEXUS_FOR_QUAD;
    } from './modules/call_snvindels.nf'

include { 
    EXPANSIONHUNTER; INDEX_EH_BAM
    } from './modules/call_repeats.nf'

// 1=male, 2=female, 0/-9=unknown/missing
// 1=unaffected, 2=affected, 0/-9=missing

// === Main workflow ===

def find_children(String family_id, String sample_id, List family_ch) {
    def family = family_ch.find { it ->
        it.familyId == family_id && it.members.find { member ->
            member.fatherId == sample_id || member.motherId == sample_id
        }
    }
    return family ? family.members.findAll { member ->
        member.fatherId == sample_id || member.motherId == sample_id
    } : []
}

workflow {
    // STEP 0. Check the input parameters

    log.info """\

    WES analysis pipeline started with the following parameters:
    
    - FASTQ directory       : ${params.fastq_dir}
    - BAM directory         : ${params.bam_dir}
    - Pedigree file         : ${params.ped_file}
    - Output directory      : ${params.out_root}
    - Use DeepVariant       : ${params.use_deepvariant}                                                                                          
    """
    
    if (!params.ped_file) {
    error "Enter the path to the pedigree file using the --ped_file parameter."
    }
    if (!params.fastq_dir) {
    error "Enter the path to the FASTQ files using the --fastq_dir parameter."
    }

    pedigree_info = Channel.fromPath("${params.ped_file}")
                        .splitText()
                        .map { it.trim() } 
                        .filter { it } 
                        .map { it.split('\t') }
                        .filter { it.size() >= 6 }
                        .map { columns -> 
                            tuple(
                                columns[0], 
                                [
                                    individualId: columns[1],
                                    fatherId    : columns[2],
                                    motherId    : columns[3],
                                    sex         : columns[4],
                                    phenotype   : columns[5]
                                ]
                            )
                        }
                        .groupTuple()

    // The structure of familyTuple is [familyId, [member1, member2, ...]]
    sampleIDs_ch = pedigree_info
                    .flatMap { familyTuple ->
                        def members = familyTuple[1]
                        members.collect { member -> member.individualId
                        }
                    }
    sampleIDs = sampleIDs_ch.toList()
    
    family_ch = pedigree_info
                    .map { familyTuple ->
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
                            if (memberCount == 1) {
                                analysisType = 'singleton'
                            } else if (memberCount == 3) {
                                analysisType = 'trio'
                            } else if (memberCount == 2) {
                                analysisType = 'duo'
                            } else {
                                analysisType = 'others'
                            }
                            return [
                                familyId: familyId, 
                                members: members, 
                                analysisType: analysisType
                            ]
                        }

    // STEP 3: Align the FASTQ files
    Channel.fromPath("${params.fastq_dir}/*_R*.fastq.gz")
        | filter { fastq ->
                   sampleIDs.get().any { sampleID -> fastq.name.contains(sampleID) }
                   }
        | map { file ->
                /* Match the file name with the expected pattern
                 * General pattern: FileID_LaneID_R[12]_PairID.fastq.gz ("R" is optional).
                 * LaneID and PairID are also optional. */
                def m = file.name =~ /^(.*?)(?:_(L\d+))?_R?([12])(?:_(\d+))?\.fastq\.gz$/
                if( !m ) {
                    error "Not compatible with the expected format. ${file.name}"
                }

                // Extract the sample ID, lane ID, read number, and pair ID
                def fileID  = m[0][1]
                def laneID  = m[0][2] ?: "" 
                def readID  = m[0][3]
                def pairID  = m[0][4] ?: "" 
            
                /* Generate a key based on the sample ID, lane ID, and pair ID, 
                 * not including the read number, because pair reads will be 
                 * grouped by the read number (e.g. R1 and R2). */
                def key = fileID + (laneID ? "_" + laneID : "") + (pairID ? "_" + pairID : "")
            
                return tuple(key, fileID, readID, file)
            }
        | groupTuple
        | map { it ->
                def fileID = it[1][0]
                def key = it[0]
                def read1 = it[3][0]
                def read2 = it[3][1]
                return tuple(
                    fileID, key, params.fasta, params.fasta_index, read1, read2
                    )
            }
        | STROBEALIGN
        | groupTuple(by: 0)  // Group by fileID. 
        | set { groupedBAMs }

    groupedBAMs
        // If there is only one key, it means that there is only one pair of FASTQ files.
        | filter { it -> it[1].size() == 1 } 
        | map { it -> 
                def fileID = it[0]
                def xam_path = it[2][0]
                def xai_path = it[3][0]
                return tuple(fileID, xam_path, xai_path) }
        | set { singleBAMs }

    groupedBAMs
        | filter { it -> it[1].size() > 1 }
        | map { it -> 
                def fileID = it[0]
                def xam_paths = it[2]
                def xai_paths = it[3]
                return tuple(fileID, xam_paths, xai_paths) }
        | MERGE_MULTIPLE_LANE_XAMS
        | concat(singleBAMs)
        // | MARKDUP
        | collect
        | set { markdup_done }

    
    // STEP 4: Modify the read group information with individual IDs
    markdup_done
        // | view
        | map { true }
        | combine(sampleIDs_ch)
        | view
        | map { _, sampleID -> sampleID }
        // | RENAME_XAM
        // | EDIT_RG
        // | INDEX_XAM
        | collect
        | map { true }
        | combine(family_ch)
        | map { _, family -> family }
        | set { align_process_done }

    // STEP 5: Call SNVs and indels
    /*  This step was divided into five parts:
     *  1. Singleton cases
     *  2. Trio cases
     *  3. Duo cases
     *  4. Others
     */

    // 1. Singleton cases
    align_process_done
        | filter { it.analysisType == 'singleton' }
        | map { family -> 
                def sample_id = family.members[0].individualId

                def xams = file("${params.bam_dir}/${sample_id}.processed.bam")
                def xais = file("${params.bam_dir}/${sample_id}.processed.bam.bai")

                return [
                    family_id    : family.familyId,
                    fasta        : params.fasta,
                    fai          : params.fasta_index,
                    sample_id    : sample_id,
                    xam          : xams,
                    xai          : xais,
                    sex          : family.members[0].sex,
                    analysisType : family.analysisType,
                    par_bed      : params.deepvariant_par,
                    par_bed_index: params.deepvariant_par_index
                ]
            }
        | set { singleton_ch }
    
    // singleton_ch
    //     | view
    //     | EXPANSIONHUNTER
    //     | INDEX_EH_BAM
    //     | view

    singleton_ch
        // | DEEPVARIANT_SINGLETON
        | view
        | set { singleton_results }

    
    // /*
    // 複数人（2-4人）の場合は，family_chに入れる．
    // == family_chのサンプル抽出基準 ==
    // [duo]  - proband は，fatherId または motherId が 0 でない．
    //        - parent は，proband の fatherId または motherId から抽出．
    
    // [trio]  - proband は，fatherId と motherId が 0 でない
    //         - mother は，proband の motherId から抽出．
    //         - father は，proband の fatherId から抽出．
    
    // [quad]  - proband は，fatherId と motherId が 0 でない
    // 　　　    かつ，phenotype が 2 で，複数人いたら，IDが若い方 を抽出．
    //         - mother は，proband の motherId から抽出．
    //         - father は，proband の fatherId から抽出．
    //         - sibling は，fatherId と motherId が 0 でない
    //           かつ proband 以外の候補を抽出．sibling_idは，DeepVariantへ．
    //           あとで，GLNEXUSでジョイント．
    // */
    // 2. Trio cases
    align_process_done
        | filter { it.analysisType == 'trio' }
        | map { family -> 
                def proband = family.members.find {
                    (it.motherId as Integer) != 0 
                    && (it.fatherId as Integer) != 0
                    }
                // Error handling when proband is not found
                if (!proband) {
                    println "Proband not found for Family ID: ${family.familyId}"
                    return null
                    }
                // When proband is found, find the parents
                def mother = family.members.find { 
                    (it.individualId as Integer) == (proband.motherId as Integer) 
                    }
                def father = family.members.find { 
                    (it.individualId as Integer) == (proband.fatherId as Integer) 
                    }
                // Error handling when mother or father is not found
                if (!mother || !father) {
                    println "Mother or father not found for Family ID: ${family.familyId}"
                    return null
                    }

                def child_id = proband.individualId
                def child_xam = file("${params.bam_dir}/${child_id}.processed.*am")
                def child_xai = file("${params.bam_dir}/${child_id}.processed.*ai")
                def mother_id = mother.individualId
                def mother_xam = file("${params.bam_dir}/${mother_id}.processed.*am")
                def mother_xai = file("${params.bam_dir}/${mother_id}.processed.*ai")
                def father_id = father.individualId
                def father_xam = file("${params.bam_dir}/${father_id}.processed.*am")
                def father_xai = file("${params.bam_dir}/${father_id}.processed.*ai")

                return [
                    family_id    : family.familyId,
                    fasta        : params.fasta,
                    fai          : params.fasta_index,
                    child_id     : child_id,
                    father_id    : father_id,
                    mother_id    : mother_id,
                    child_xam    : child_xam,
                    father_xam   : father_xam,
                    mother_xam   : mother_xam,
                    child_xai    : child_xai,
                    father_xai   : father_xai,
                    mother_xai   : mother_xai,
                    child_xxxy   : proband.sex,
                    analysisType : family.analysisType,
                    par_bed      : params.deepvariant_par,
                    par_bed_index: params.deepvariant_par_index
                ]
            }
        | set { trio_ch }

    // Separate the trio by xx and xy
    trio_ch
        | filter {it.child_xxxy == '1' }
        | DEEPTRIO_MALE_CHILD
        | GLNEXUS_FOR_MALE_TRIO
        | MERGE_BCFS_FOR_MALE_CHILD_TRIO
        | set { trio_male_results }
    
    trio_ch
        | filter {it.child_xxxy != '1' }
        | DEEPTRIO_FEMALE_CHILD
        | GLNEXUS_FOR_FEMALE_TRIO
        | set { trio_female_results }


    // 3. Duo cases
    align_process_done
        | filter { it.analysisType == 'duo' }
        | map { family -> 
                def proband = family.members.find { 
                    (it.fatherId as Integer) != 0 
                    || (it.motherId as Integer) != 0 
                    }
                // Error handling when proband is not found
                if (!proband) {
                    println "Proband not found for Family ID: ${family.familyId}"
                    return null
                    }
                // When proband is found, find the parents
                def mother = family.members.find { 
                    (it.individualId as Integer) == (proband.motherId as Integer) 
                    }
                def father = family.members.find { 
                    (it.individualId as Integer) == (proband.fatherId as Integer) 
                    }
                def parent_xxxy = "xx"
                if (father) {parent_xxxy = "xy"}

                def child_id = proband.individualId
                def child_xam = file("${params.bam_dir}/${child_id}.processed.*am")
                def child_xai = file("${params.bam_dir}/${child_id}.processed.*ai")
                def parent_id = [ mother?.individualId, father?.individualId ].findAll { it != null }[0]
                def parent_xam = file("${params.bam_dir}/${parent_id}.processed.*am")
                def parent_xai = file("${params.bam_dir}/${parent_id}.processed.*ai")

                return [
                    family_id    : family.familyId,
                    child_id     : proband.individualId,
                    parent_id    : parent_id,
                    child_xam    : child_xam,
                    child_xai    : child_xai,
                    parent_xam   : parent_xam,
                    parent_xai   : parent_xai,
                    child_xxxy   : proband.sex,
                    analysisType : family.analysisType,
                    fasta        : params.fasta,
                    fai          : params.fasta_index,
                    par_bed      : params.deepvariant_par,
                    parent_xxxy  : parent_xxxy
                ]
            }
        | DEEPTRIO_FOR_DUO
        | GLNEXUS_FOR_DUO
        | set { duo_results }

    // More than 4 individuals
    align_process_done
        | filter { it.analysisType == 'others' }
        | view
        | flatMap { family ->
                    def family_id = family.familyId
                    family.members.collect { member ->
                        return member << [ family_id: family_id ]
                    }
                }
        | view
        | map { individual ->
            def family_id = individual.family_id
            def sample_id = individual.individualId
            // def sex = individual.sex
            // def xams = file("${params.bam_dir}/${sample_id}.processed.bam")
            // def xais = file("${params.bam_dir}/${sample_id}.processed.bam.bai")
            // def father_id = individual.fatherId as Integer
            // def mother_id = individual.motherId as Integer

            // def relationship = null
            if (individual.fatherId == "0" && individual.motherId == "0") { 
                def familyList = family_ch.toList().get()
                def children_list = find_children(family_id, sample_id, familyList)
                println "${children_list}"
            } else {
                println "Parent"
            }
            return [
                family_id    : family_id,
                fasta        : params.fasta,
                fai          : params.fasta_index,
                sample_id    : sample_id
                ]
        }   
            //         if (has_children) {
            //             def has_grandchildren = family_ch.find { 
            //                 it.familyId == family_id 
            //                 && it.members.find { 
            //                     it.fatherId == sample_id 
            //                     || it.motherId == sample_id 
            //                 }
            //             }

            //         relationship = 'proband'
                
                
            //     } else {
            //         relationship = 'sibling'
            //     }

            //     return [
            //         family_id    : family_id,
            //         fasta        : params.fasta,
            //         fai          : params.fast,
            //         sample_id    : sample_id,
            //         xam          : xams,
            //         xai          : xais,
            //         sex          : family.members[0].sex,
            //         analysisType : family.analysisType,
            //         par_bed      : params.deepvariant_par,
            //         par_bed_index: params.deepvariant_par_index
            //         relationship : individual.phenotype
            //     ]
            // }
        
        | view
        | set { dummy }



        // | map { family ->
        //         def sample_id = null
        //         if family.members

        //         def proband = null
        //         def proband_candidate = family.members.findAll {
        //             (it.fatherId as Integer) != 0 
        //             && (it.motherId as Integer) != 0 
        //             && (it.phenotype as Integer) == 2
        //             }
                
        //         def sibling = null
        //         if (proband_candidate.size() == 2) {
        //             proband = proband_candidate.min { it.individualId as Integer }
        //             sibling = family.proband_candidate.find {
        //                 (it.individualId as Integer) != (proband.individualId as Integer)
        //             }
        //         } else if (proband_candidate.size() == 1) {
        //             proband = proband_candidate[0]
        //             sibling = family.members.find {
        //                 (it.fatherId as Integer) != 0 
        //                 && (it.motherId as Integer) != 0 
        //                 &&  (it.individualId as Integer) != (proband.individualId as Integer)
        //             }
        //         } else {
        //             // Error handling when proband is not found
        //             println "Proband not found for Family ID: ${family.familyId}"
        //             return null
        //         }
        //         println "Proband: ${proband.individualId}, Sibling: ${sibling.individualId}"
        //         // proband が見つかった場合、親を探す
        //         def mother = family.members.find { 
        //             (it.individualId as Integer) == (proband.motherId as Integer) 
        //         }
        //         def father = family.members.find { 
        //             (it.individualId as Integer) == (proband.fatherId as Integer) 
        //         }
        //         // Error handling when mother or father is not found
        //         if (!mother || !father) {
        //             println "Mother or father not found for Family ID: ${family.familyId}"
        //             return null
        //         }
        //         println "Mother: ${mother.individualId}, Father: ${father.individualId}"
        //         println "Proband: ${proband.individualId}, Sibling: ${sibling.individualId}"

        //         def proband_id = proband.individualId
        //         def proband_xam = file("${params.bam_dir}/${proband_id}.processed.*am")
        //         def proband_xai = file("${params.bam_dir}/${proband_id}.processed.*ai")
        //         println "Child: ${proband_id}, ${proband_xam}, ${proband_xai}"
        //         def mother_id = mother.individualId
        //         def mother_xam = file("${params.bam_dir}/${mother_id}.processed.*am")
        //         def mother_xai = file("${params.bam_dir}/${mother_id}.processed.*ai")
        //         def father_id = father.individualId
        //         def father_xam = file("${params.bam_dir}/${father_id}.processed.*am")
        //         def father_xai = file("${params.bam_dir}/${father_id}.processed.*ai")
        //         def sibling_id = sibling.individualId
        //         def sibling_xam = file("${params.bam_dir}/${sibling_id}.processed.*am")
        //         def sibling_xai = file("${params.bam_dir}/${sibling_id}.processed.*ai")

        //         return [
        //             family_id    : family.familyId,
        //             analysisType : family.analysisType,
        //             fasta        : params.fasta,
        //             fai          : params.fasta_index,
        //             proband_id   : proband.individualId,
        //             mother_id    : mother ? mother.individualId : null,
        //             father_id    : father ? father.individualId : null,
        //             sibling_id   : sibling ? sibling.individualId : null,
        //             child_xam    : proband_xam,
        //             mother_xam   : mother_xam,
        //             father_xam   : father_xam,
        //             sibling_xam  : sibling_xam,
        //             child_xai    : proband_xai,
        //             mother_xai   : mother_xai,
        //             father_xai   : father_xai,
        //             sibling_xai  : sibling_xai,
        //             child_xxxy   : proband.sex,
        //             sibling_xxxy : sibling.sex,
        //             par_bed      : params.deepvariant_par
        //         ]
        //     }
        // | view
        // | map { quad_dv_results -> 
        //         def familyId = quad_dv_results[0]
        //         def gvcf_paths = quad_dv_results[1] + [quad_dv_results[3]]
        //         def tbi_paths = quad_dv_results[2] + [quad_dv_results[4]]
        //         tuple (familyId, gvcf_paths, tbi_paths)
        //     }
        // | set { quad_deepvariant_results }

    // trio_deepvariant_results
    //     | GLNEXUS
    //     | view
    //     | combine(duo_deeptrio_results)
    //     | view

    // duo_deeptrio_results
        // | view
        // | combine(quad_deepvariant_results)
    // variant_call_results.concat(variant_call_for_quad_results)
        // | GLNEXUS
    //     | view
}

workflow.onComplete {
    println ""
    println "~*~*~*~*~~*~*~*~*~~*~*~*~*~~*~*~*~*~*~*~*~~*~*~*~*~*~*~*~*~"
    println "Pipeline completed at: $workflow.complete"
    println "Execution time       : $workflow.duration"
    println "Execution status     : ${ workflow.success ? 'OK' : 'failed' }"
    println "~*~*~*~*~~*~*~*~*~~*~*~*~*~~*~*~*~*~*~*~*~~*~*~*~*~*~*~*~*~"
}




