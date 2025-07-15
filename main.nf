#!/usr/bin/env nextflow

// Module INCLUDE statements
include { FASTQC } from './modules/fastqc.nf'
include { TRIM_GALORE } from './modules/trim_galore.nf'
include { HISAT2_ALIGN } from './modules/hisat2_align.nf'

/*
 * Pipeline parameters
 */

// Primary input
params.reads = "$baseDir/data/reads/ENCSR000COQ1_1.fastq.gz"
params.hisat2_index_zip = "$baseDir/data/genome_index.tar.gz"

workflow {
    log.info """\
        R N A S E Q - N F   P I P E L I N E
        ===================================
        transcriptome: ${params.hisat2_index_zip}
        reads        : ${params.reads}
    """
    // Create input channel
    read_ch = channel.fromPath(params.reads)

    // Initial Quality Control
    FASTQC(read_ch)

    // Adapter Trimming and post-trimming QC
    TRIM_GALORE(read_ch)

    // Alignment to a refence genome
    HISAT2_ALIGN(TRIM_GALORE.output.trimmed_reads, file(params.hisat2_index_zip))

}