#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl=2

process ANALYSE_STRUCTURE {

    tag "${sample_id}"
    publishDir "${params.outdir}/hybrids", mode: 'copy', overwrite: true

    // memory 16G
    // time '12h'

    input:
        tuple val(sample_id), path(hybrids)
        path(transcript_fa)

    output:
        tuple val(sample_id), path("${sample_id}.hybrids.mfe.tsv.gz"), emit: hybrids

    script:

    args = ''
    if ( params.shuffled_mfe ) args += ' --shuffled_mfe '
    if ( params.clusters_only) args += ' --clusters_only'

    """
    analyse_structure.R --hybrids $hybrids --fasta $transcript_fa --output ${sample_id}.hybrids.mfe.tsv.gz $args
    """

}