#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl=2

process GET_BINDING_ENERGY {

    tag "${sample_id}"
    publishDir "${params.outdir}/hybrids", mode: 'copy', overwrite: true

    time '24h'

    input:
        tuple val(sample_id), path(hybrids)
        path(transcript_fa)

    output:
        tuple val(sample_id), path("${sample_id}.hybrids.mfe.tsv.gz")

    script:

    args = ''
    if ( params.shuffled_mfe ) args += '--shuffled_mfe'

    """
    get_binding_energy.R --hybrids $hybrids --fasta $transcript_fa --output ${sample_id}.hybrids.mfe.tsv.gz $args
    """

}