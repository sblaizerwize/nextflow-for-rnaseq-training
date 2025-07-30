#!/usr/bin/env nextflow

process KALLISTO_QUAN {

    container "community.wave.seqera.io/library/kallisto:d61e08588cd20f0f"
    publishDir params.outdir, mode: 'copy'

    input:
    path reads
    path index_zip

    output:
    path "${reads.simpleName}-kallisto", emit: abundance 
    path "${reads.simpleName}-kallisto.log", emit: log

    script:
    """
    tar -xzvf $index_zip
    kallisto quant -i ${index_zip.simpleName} \
    -o ${reads.simpleName}-kallisto \
    -t ${task.cpus} \
    --single \
    -l 275 \
    -s 30 \
    $reads \
    &> ${reads.simpleName}-kallisto.log
    """
}
