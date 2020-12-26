#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl=2

process GET_CONTACT_MAPS {

    tag "${sample_id}"
    publishDir "${params.outdir}/maps", mode: 'copy', overwrite: true

    time '24h'
    memory '32 G'

    input:
        tuple val(sample_id), path(hybrids)
        path(fai)
        path(genes)

    output:
        tuple val(sample_id), path("${sample_id}.*.mat.rds"), emit: map
        tuple val(sample_id), path("${sample_id}.*.binned_map.tsv.gz"), emit: binned_map

    script:
    """
    get_contact_map.R --hybrids $hybrids --genes $genes --fai $fai --output ${sample_id}
    """

}