#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl=2

process CUTADAPT {

    tag "${sample_id}"
    // publishDir "${params.outdir}/trimmed", mode: 'copy', overwrite: true

    input:
        tuple val(sample_id), path(reads)

    output:
        tuple val(sample_id), path("${sample_id}.trimmed.fastq.gz"), emit: fastq
        path("*.cutadapt.log"), emit: log

    script:
    args = " -j ${task.cpus}"
    args += " -a " + params.adapter 
    args += " -q " + params.min_quality
    args += " --minimum-length " + params.min_readlength
    args += " -o ${sample_id}.trimmed.fastq.gz"

    """
    cutadapt $args $reads > ${sample_id}.cutadapt.log
    """
    
}
