#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl=2

process CUTADAPT {
    tag "${sample_id}"
    publishDir "${params.outdir}/trimmed", mode: 'copy', overwrite: true

    cpus 8
    time '24h'

    input:
        tuple val(sample_id), path(reads)

    output:
        tuple val(sample_id), path("${sample_id}.trimmed.fastq.gz"), emit: fastq

    args = " -j ${task.cpus}"
    args += " -a " + params.adapter 
    args += " -q " + params.min_quality
    args += " --minimum-length " + params.min_readlength
    args += " -o ${sample_id}.trimmed.fastq.gz"

    cmd = "cutadapt $args $reads > ${sample_id}_cutadapt.log"

    println ("[MODULE] cutadapt command: " + cmd)

    shell:
    """
    $cmd
    """
}
