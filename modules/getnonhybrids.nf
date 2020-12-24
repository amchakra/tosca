#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl=2

process GET_NON_HYBRIDS {

    tag "${sample_id}"
    publishDir "${params.outdir}/nonhybrids", mode: 'copy', overwrite: true

    time '24h'

    input:
        tuple val(sample_id), path(hybrids), path(reads)

    output:
        tuple val(sample_id), path("${sample_id}.nonhybrid.fastq.gz"), emit: nonhybrids

    script:
    
    """
    awk '{ if(\$2 ~ "single|multi_overlap") { print \$1 } }' > ${sample_id}.hits.txt

    filterbyname.sh in=$reads out=${sample_id}.nonhybrid.fastq.gz names=${sample_id}.hits.txt include=f
    """

}