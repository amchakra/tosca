#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl=2

process AGGREGATE_LOGS {

    tag "${sample_id}"
    label 'process_medium'

    publishDir "${params.outdir}/hybrids/logs", mode: 'copy', overwrite: true, pattern: "${sample_id}.*.log"

    input:
        tuple val(sample_id), path(filter_blat_logs), path(identify_hybrids)

    output:
        tuple val(sample_id), path("${sample_id}.filter_blat.log"), path("${sample_id}.identify_hybrids.log"), emit: logs

    script:

    """
    aggregate_logs.R -l . -t ${task.cpus} -o ${sample_id}
    """

}

process TRACK_READ_FATE {

    tag "${sample_id}"
    label 'process_low'

    container 'iraiosub/nf-riboseq-qc:latest'

    publishDir "${params.outdir}/hybrids/logs", mode: 'copy', overwrite: true

    input:
        tuple val(sample_id), path(logs)

    output:
        tuple val(sample_id), path("${sample_id}_sankey.html"), path("*_sankey_files"), emit: sankey

    script:

    """
    plot_sankey.R -l . -o ${sample_id}_sankey.html
    """

}