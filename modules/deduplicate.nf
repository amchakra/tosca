#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl=2

process DEDUPLICATE {

    tag "${sample_id}"
    label 'process_high'

    publishDir "${params.outdir}/hybrids", mode: 'copy', overwrite: true, pattern: '*.hybrids.dedup.tsv.gz'
    publishDir "${params.outdir}/hybrids/logs", mode: 'copy', overwrite: true, pattern: '*.dedup.log'


    input:
        tuple val(sample_id), path(hybrids)

    output:
        tuple val(sample_id), path("${sample_id}.hybrids.dedup.tsv.gz"), emit: hybrids
        tuple val(sample_id), path("*.dedup.log"), emit: log

    script:

    umi_separator = params.umi_separator
    dedup_method = params.dedup_method

    """
    deduplicate_hybrids.py $hybrids ${sample_id}.hybrids.dedup.tsv.gz $umi_separator $dedup_method > ${sample_id}.dedup.log
    """

}

process FLAG_REPEAT_SHIFTS {

    // Adds repeat_shift to the deduplicated hybrids (bin/flag_repeat_shifts.py): genomic-order
    // hybrids whose short transcript gap is an exact copy of its flank, i.e. a linear read with
    // fewer tandem-repeat copies than the reference. Flag only, unless drop_repeat_shifts.

    tag "${sample_id}"
    label 'process_medium'

    publishDir "${params.outdir}/hybrids", mode: 'copy', overwrite: true, pattern: '*.hybrids.dedup.flagged.tsv.gz'
    publishDir "${params.outdir}/hybrids/logs", mode: 'copy', overwrite: true, pattern: '*.repeat_shift.log'

    input:
        tuple val(sample_id), path(hybrids)
        path(transcript_fa)
        path(transcript_fai)

    output:
        tuple val(sample_id), path("${sample_id}.hybrids.dedup.flagged.tsv.gz"), emit: hybrids
        tuple val(sample_id), path("*.repeat_shift.log"), emit: log

    script:

    drop = params.drop_repeat_shifts ? '--drop' : ''

    """
    flag_repeat_shifts.py $hybrids $transcript_fa ${sample_id}.hybrids.dedup.flagged.tsv.gz $drop > ${sample_id}.repeat_shift.log
    """

}
