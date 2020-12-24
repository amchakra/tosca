#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl=2

process FILTER_SPLICED_READS {

    tag "${sample_id}"
    publishDir "${params.outdir}/filtered", mode: 'copy', overwrite: true

    time '12h'

    input:
        tuple val(sample_id), path(bam), path(bai)

    output:
        tuple val(sample_id), path("${sample_id}.unspliced.fastq.gz"), emit: fastq
        path("*.filter_spliced_reads.log"), emit: log

    script:

    cmd = "filter_spliced_reads.py $bam ${sample_id} > ${sample_id}.filter_spliced_reads.log"

    if(params.verbose) { println ("[MODULE] FILTERSPLICEDREADS: " + cmd) }

    """
    $cmd
    """
    
}