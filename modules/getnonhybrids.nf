#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl=2

process GET_NON_HYBRIDS {

    tag "${sample_id}"
    label 'process_low'
    publishDir "${params.outdir}/nonhybrids", mode: 'copy', overwrite: true

    input:
        tuple val(sample_id), path(hybrids), path(reads)

    output:
        tuple val(sample_id), path("${sample_id}.nonhybrid.fastq.gz"), emit: nonhybrids

    script:
    
    """
    gunzip -c $hybrids | \
    awk -v col=name 'NR==1{for(i=1;i<=NF;i++){if(\$i==col){colnum=i;break}}} {print \$colnum}'  \
    > ${sample_id}.hits.txt 

    filterbyname.sh in=$reads out=${sample_id}.nonhybrid.fastq.gz names=${sample_id}.hits.txt include=f
    """

}