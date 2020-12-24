#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl=2

process filtersplicedreads {

    tag "${meta.sample_id}"
    publishDir "${params.outdir}/filtered", mode: 'copy', overwrite: true

    time '12h'

    input:
        tuple val(meta), path(bam), path(bai)

    output:
        tuple val(meta), path("${meta.sample_id}.unspliced.fastq.gz")

    script:
    """
    remove_spliced_reads.py $bam ${meta.sample_id}
    """
    
}