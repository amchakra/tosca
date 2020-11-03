#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl=2

params.evalue = params.evalue
params.maxhits = params.maxhits

process identifyhybrids {

    tag "${sample_id}"
    // publishDir "${params.outdir}/hybrids", mode: 'copy', overwrite: true

    time '24h'
    cpus = 8
    memory '64 G'

    input:
        tuple val(sample_id), path(blast8), path(reads)

    output:
        tuple val(sample_id), path("${sample_id}.hybrids.tsv.gz")

    script:

    identify_hybrids = "$baseDir/bin/identify_hybrids.R"

    """
    Rscript --vanilla $identify_hybrids -t ${task.cpus} -b $blast8 -f $reads -o ${sample_id}.hybrids.tsv.gz
    """

}