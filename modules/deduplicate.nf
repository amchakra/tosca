#!/usr/bin/env nextflow

// Specify DSL2
nextflow.preview.dsl = 2

process deduplicate {

    tag "${sample_id}"
    publishDir 'results/dedup', mode: 'copy', overwrite: false

    memory '64G'
    time '24h'

    input:
        tuple val(sample_id), path(reads), path(bai)

    output:
        tuple val(sample_id), path("${sample_id}.dedup.bam"), path("${sample_id}.dedup.bam.bai")

    shell:
    """
    if `samtools view $reads | cut -f 1 | head -1 | grep -q ':rbc:'`; then
        echo iCount
        umi_tools dedup --method directional --umi-separator=":" -I $reads -S ${sample_id}.dedup.bam
    else
        echo UMItools
        umi_tools dedup --method directional -I $reads -S ${sample_id}.dedup.bam
    fi
    
    sambamba index ${sample_id}.dedup.bam
    """

}

process deduplicate_unique {

    tag "${sample_id}"
    publishDir 'results/dedup', mode: 'copy', overwrite: false

    memory '64G'
    time '24h'

    input:
        tuple val(sample_id), path(reads), path(bai)

    output:
        tuple val(sample_id), path("${sample_id}.dedup.bam"), path("${sample_id}.dedup.bam.bai")

    shell:
    """
    if `samtools view $reads | cut -f 1 | head -1 | grep -q ':rbc:'`; then
        echo iCount
        umi_tools dedup --method unique --umi-separator=":" -I $reads -S ${sample_id}.dedup.bam
    else
        echo UMItools
        umi_tools dedup --method unique -I $reads -S ${sample_id}.dedup.bam
    fi
    
    sambamba index ${sample_id}.dedup.bam
    """

}