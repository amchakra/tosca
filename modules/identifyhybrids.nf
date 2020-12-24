#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl=2

process IDENTIFY_HYBRIDS {

    tag "${sample_id}"
    // publishDir "${params.outdir}/hybrids/split/", mode: 'copy', overwrite: true

    time '24h'
    cpus = 8
    memory '64 G'

    input:
        tuple val(sample_id), path(blast8), path(reads)

    output:
        tuple val(sample_id), path("${sample_id}.hybrids.tsv.gz"), emit: hybrids

    script:

    """
    identify_hybrids.R -t ${task.cpus} -b $blast8 -f $reads -o ${sample_id}.hybrids.tsv.gz
    """

}