#!/usr/bin/env nextflow

params.reads = "$baseDir/reads/*_{1,2}.fq.gz"
params.assembly = "$baseDir/assembly/assembly.fasta"
params.reference = "$baseDir/reference/reference.fasta"
params.gff = "$baseDir/reference/reference.gff"
params.genes = "$baseDir/reference/genes.fasta"
params.genes_gff = "$baseDir/reference/genes.gff"
params.genes_gtf = "$baseDir/reference/genes.gtf"
params.genes_bed = "$baseDir/reference/genes.bed"
params.genes_fasta = "$baseDir/reference/genes.fasta"
params.genes_fasta_index = "$baseDir/reference/genes.fasta.fai"
