#!/usr/bin/env nextflow

// Specify DSL2
nextflow.preview.dsl = 2

process mapblat {
    tag "${sample_id}"
    publishDir "${params.outdir}/split", mode: 'copy', overwrite: true

    cpus 1
    time '24h'

    input:
        tuple val(sample_id), path(reads)

    output:
        tuple val(sample_id), path("${sample_id}_*.fastq.gz")

    shell:
    """
    zcat $i | split -l 4000000 -d --additional-suffix .fastq - ${sample_id}_
    pigz *.fastq
    """

}