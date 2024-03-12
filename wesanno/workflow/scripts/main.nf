#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.input = ''

process VCFANNO {
    input:
    path vcf

    output:
    path 'vcfanno.vcf'

    script:
    """
    vcfanno \\
      -p 8 \\
      ${WORKFLOW_WES}/config/vcfanno_${ASSEMBLY}.toml \\
      $vcf \\
      > vcfanno.vcf
    """
}

process SPLIT {
    input:
    path vcf

    output:
    path 'chunk_*'

    script:
    """
    ${WORKFLOW_WES}/scripts/chrsplit.sh $vcf
    """
}

process SPLICEAI {
    input:
    path vcf

    output:
    path '*.splai.vcf'

    script:
    """
    spliceai \\
      -I $vcf \\
      -O ${vcf.baseName}.splai.vcf \\
      -R /resources/ReferenceGenome/exome_pipeline/human_g1k_v37_fix.fasta \\
      -A /resources/SAI10kcalc/SAI10kcalc_spliceai_tx.hg19.tsv \\
      -D 4999 \\
      -M 0
    """
}

process CONCATENATESORT {
    input:
    path collected_vcfs

    output:
    path 'splai.concat.sort.pre_annovar.vcf'

    script:
    """
    bcftools concat \\
      --output splai.concat.pre_annovar.vcf \\
      --output-type v \\
      --threads 4 \\
      ${collected_vcfs}

    bcftools sort \\
      --output-file splai.concat.sort.pre_annovar.vcf \\
      --output-type v \\
      --max-mem 8G \\
      splai.concat.pre_annovar.vcf
    """
}


process ANNOVAR {
    input:
    path pre_annovar_vcf

    output:
    tuple path( '*.txt' ), path( '*.vcf' ), emit: annovar_files

    script:
    """
    ${WORKFLOW_WES}/scripts/run_annovar_wes_${ASSEMBLY}.sh \\
      ${pre_annovar_vcf} \\
      ${ANNOVAR_DB} \\
      ${SPLICING_THRESHOLD}
    """
}


process FORMATANNOVAR {
    input:
    tuple path(txt), path(vcf)

    output:
    path "${txt.baseName}.renamed.txt"

    script:
    """
    ${WORKFLOW_WES}/bin/anvrformatter \\
      ${txt} \\
      ${WORKFLOW_WES}/config/rename.toml
    """
}


process HGMDANNOTATOR {
    input:
    path txt

    output:
    path 'exome_summary_*.txt'

    script:
    """
    run_hgmd_annotator ${txt} exome_summary.txt
    DATE=\$(date +'%Y%m%d_%H%M%S')
    mv exome_summary.txt exome_summary_\${DATE}.txt
    cp exome_summary_\${DATE}.txt ${workflow.launchDir}/
    """
}


workflow {
    input_ch = Channel.fromPath(params.input)

    VCFANNO(input_ch)
        .set{ vcfannoFile_ch }

    SPLIT(vcfannoFile_ch)
        .flatMap{ file -> file }
        .set{ splitFiles_ch }

    SPLICEAI(splitFiles_ch)
        .collect()
        .set{ spliceaiFiles_ch }

    spliceaiFiles_ch
        .collect()
        .set{ collectedSpliceaiFiles_ch }

    CONCATENATESORT(collectedSpliceaiFiles_ch)
        .set{ concatFile_ch }

    ANNOVAR(concatFile_ch)
        .set{ annovarFiles_ch }

    FORMATANNOVAR(annovarFiles_ch)
        .set{ renamedFile_ch }

    HGMDANNOTATOR(renamedFile_ch)
        .set{ hgmdAnnotatedFile_ch }

    hgmdAnnotatedFile_ch.view()    

}
