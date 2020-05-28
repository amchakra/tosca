#!/usr/bin/env nextflow

// Specify DSL2
nextflow.preview.dsl = 2

process trim {
    tag "${sample_id}"
    publishDir 'results/trimmed', mode: 'copy', overwrite: false

    cpus 8
    time '24h'

    input:
        tuple val(sample_id), path(reads)

    output:
        tuple val(sample_id), path("${sample_id}.trimmed.fastq.gz")

    shell:
    """
    cutadapt -j ${task.cpus} --minimum-length 16 -q 10 -a AGATCGGAAGAGC -o "${sample_id}.trimmed.fastq.gz" $reads > "${sample_id}_cutadapt.log"
    """
}