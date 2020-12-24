#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl=2

process convert_fastq_to_fasta {
    tag "${sample_id}"
    cache true

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

process FASTQ_TO_FASTA {
    tag "${sample_id}"
    cache true

    cpus 8
    time '24h'

    input:
        tuple val(sample_id), path(reads)

    output:
        tuple val(sample_id), path("${sample_id}.fasta"), emit: fasta

    script:
    """
    reformat.sh in1=$reads out1=${sample_id}.fasta
    """
}