#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl=2

process AGGREGATE_LOGS {

    tag "${sample_id}"
    label 'process_medium'

    publishDir "${params.outdir}/logs", mode: 'copy', overwrite: true, pattern: "${sample_id}.*.log"

    input:
        tuple val(sample_id), path(filter_blat_logs), path(identify_hybrids_logs)

    output:
        tuple val(sample_id), path("${sample_id}.filter_blat.log"), path("${sample_id}.identify_hybrids.log"), emit: logs

    script:

    """
    aggregate_logs.R -l . -t ${task.cpus} -o ${sample_id}
    """

}