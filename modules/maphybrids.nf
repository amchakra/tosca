#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl=2

process BLAT {
    tag "${sample_id}"
    if(!params.keep_cache) cache false

    input:
        tuple val(sample_id), path(reads)
        path(transcript_fa)

    output:
        tuple val(sample_id), path("${sample_id}.blast8.gz"), emit: blast8

    script:

    args = " -threads=${task.cpus} "
    args += " -stepSize=" + params.step_size 
    args += " -tileSize=" + params.tile_size
    args += " -minScore=" + params.min_score 
    args += " -out=blast8 "
    args += " -dots=1000000 "

    """
    pblat $args $transcript_fa $reads ${sample_id}.blast8
    pigz ${sample_id}.blast8
    """

}

process FILTER_BLAT {

    tag "${sample_id}"
    if(params.keep_intermediates) publishDir "${params.outdir}/mapped", mode: 'copy', overwrite: true

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
