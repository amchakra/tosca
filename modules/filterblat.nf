#!/usr/bin/env nextflow

// Specify DSL2
nextflow.preview.dsl = 2

params.evalue = params.evalue
params.maxhits = params.maxhits

process filterblat {

    tag "${sample_id}"
    publishDir "${params.outdir}/mapped", mode: 'copy', overwrite: true

    time '24h'

    input:
        tuple val(sample_id), path(blast8)

    output:
        tuple val(sample_id), path("${sample_id}.filtered.blast8.gz")

    script:

    evalue = "$params.evalue"
    maxhits = "$params.maxhits"

    """
    filter_blat.py $blast8 ${sample_id}.filtered.blast8.gz $evalue $maxhits
    """

}