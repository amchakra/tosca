#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl=2

process DEDUPLICATE {

    tag "${sample_id}"
    publishDir "${params.outdir}/hybrids", mode: 'copy', overwrite: false

    memory '64G'
    time '6h'

    input:
        tuple val(sample_id), path(hybrids)

    output:
        tuple val(sample_id), path("${sample_id}.hybrids.dedup.tsv.gz"), emit: hybrids
        path("*.dedup.log"), emit: log

    script:

    umi_separator = params.umi_separator
    dedup_method = params.dedup_method

    """
    deduplicate_hybrids.py $hybrids ${sample_id}.hybrids.dedup.tsv.gz $umi_separator $dedup_method > ${sample_id}.dedup.log
    """

}