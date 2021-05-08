#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl=2

process BLAT {
    tag "${sample_id}"
    // if(params.keep_intermediates) cache true

    cpus 8
    memory '16 G'
    time '4h'

    input:
        tuple val(sample_id), path(reads)
        path(transcript_fa)

    output:
        tuple val(sample_id), path("${sample_id}.blast8.gz"), emit: blast8

    script:

    args = " -threads=${task.cpus} "
    args += " -stepSize=5 -tileSize=11 -minScore=15 -out=blast8 "
    args += " -dots=1000000 "

    """
    pblat $args $transcript_fa $reads ${sample_id}.blast8
    pigz ${sample_id}.blast8
    """

}

process FILTER_BLAT {

    tag "${sample_id}"
    if(params.keep_intermediates) publishDir "${params.outdir}/mapped", mode: 'copy', overwrite: true
    
    memory '16 G'
    time '4h'

    input:
        tuple val(sample_id), path(blast8)

    output:
        tuple val(sample_id), path("${sample_id}.filtered.blast8.gz"), emit: blast8

    script:

    evalue = params.evalue
    maxhits = params.maxhits

    """
    filter_blat.py $blast8 ${sample_id}.filtered.blast8.gz $evalue $maxhits
    """

}