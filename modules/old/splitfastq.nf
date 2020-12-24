#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl=2

process SPLIT_FASTQ {
    tag "${sample_id}"

    cpus 1
    time '24h'

    input:
        tuple val(sample_id), path(reads)

    output:
        path("${sample_id}_*.fastq.gz"), emit: fastq

    script:

    split_size = params.split_size * 4

    cmd = "gunzip -c $reads | split -l $split_size --additional-suffix .fastq - ${sample_id}_ && pigz *.fastq"

    if(params.verbose) { println ("[MODULE] CUTADAPT: " + cmd) }

    """
    $cmd
    """

}