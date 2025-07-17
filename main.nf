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
params.transcriptome = "$baseDir/data/aligned/genome_index.tar.gz"
params.outdir = "results"

workflow {
    log.info """\
        R N A S E Q - N F   P I P E L I N E
        ===================================
        transcriptome: ${params.transcriptome}
        reads        : ${params.reads}
        outdir       : ${params.outdir}
    """
    // Create input channel
    read_ch = channel.fromPath(params.reads)

    // Initial Quality Control
    FASTQC(read_ch)

    // Adapter Trimming and post-trimming QC
    TRIM_GALORE(read_ch)

    // Alignment to a refence genome
    HISAT2_ALIGN(TRIM_GALORE.output.trimmed_reads, file(params.transcriptome))

}