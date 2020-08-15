#!/usr/bin/env nextflow

// Specify DSL2
nextflow.preview.dsl = 2

process convert_fastq_to_fasta {
    tag "${sample_id}"
    cache false

    cpus 8
    time '24h'

    input:
        tuple val(sample_id), path(reads)

    output:
        tuple val(sample_id), path("${sample_id}.fasta")

    shell:
    """
    ml FASTX-Toolkit

    zcat $reads | fastq_to_fasta -n -o "${sample_id}.fasta"
    """
}
