#!/usr/bin/env nextflow

// Module INCLUDE statements
include { FASTQC } from './modules/fastqc.nf'
include { TRIM_GALORE } from './modules/trim_galore.nf'
include { HISAT2_ALIGN } from './modules/hisat2_align.nf'

/*
 * Pipeline parameters
 */

// Primary input
params.reads = "data/reads/ENCSR000COQ1_1.fastq.gz"
params.hisat2_index_zip = file('data/genome_index.tar.gz')

workflow {

    // Create input channel
    read_ch = channel.fromPath(params.reads)

    // Initial Quality Control
    FASTQC(read_ch)

    // Adapter Trimming and post-trimming QC
    TRIM_GALORE(read_ch)

    // Alignment to a refence genome
    HISAT2_ALIGN(TRIM_GALORE.output.trimmed_reads, params.hisat2_index_zip)

}