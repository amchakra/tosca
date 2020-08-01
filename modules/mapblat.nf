#!/usr/bin/env nextflow

// Specify DSL2
nextflow.preview.dsl = 2

process mapblat {
    tag "${sample_id}"
    cache false

    cpus 8
    memory '32 G'
    time '24h'

    input:
        tuple val(sample_id), path(reads), path(transcript_fa)

    output:
        tuple val(sample_id), path("${sample_id}.blast8.gz")

    shell:
    """
    pblat -threads=${task.cpus} -dots=1000000 -stepSize=5 -tileSize=11 -minScore=15 -out=blast8 $transcript_fa $reads ${sample_id}.blast8
    pigz ${sample_id}.blast8
    """

}