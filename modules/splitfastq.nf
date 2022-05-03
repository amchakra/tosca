#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl=2

process SPLIT_FASTQ {
    tag "${sample_id}"
    if(!params.keep_cache) cache false
    
    // cpus 8
    // time '1h'

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

process FASTQ_TO_FASTA {
    tag "${sample_id}"
    if(!params.keep_cache) cache false

    cpus 8
    time '1h'

    input:
        tuple val(sample_id), path(reads)

    output:
        tuple val(sample_id), path("${sample_id}.fasta"), emit: fasta

    script:
    """
    reformat.sh in1=$reads out1=${sample_id}.fasta
    """
}

