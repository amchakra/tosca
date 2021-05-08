#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl=2

process GET_CONTACT_MAPS {

    tag "${sample_id}"
    publishDir "${params.outdir}/maps", mode: 'copy', overwrite: true

    time '6h'
    memory '64 G'

    input:
        tuple val(sample_id), path(hybrids)
        path(fai)
        path(genes)

    output:
        tuple val(sample_id), path("${sample_id}.*.mat.rds"), emit: map
        tuple val(sample_id), path("${sample_id}.*_binned_map.tsv.gz"), emit: binned_map

    script:

    bin_size = params.bin_size

    """
    get_contact_map.R --hybrids $hybrids --genes $genes --fai $fai --bin_size $bin_size --output ${sample_id}
    """

}