#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl=2

process CLUSTER_HYBRIDS {

    tag "${sample_id}"
    publishDir "${params.outdir}/hybrids", mode: 'copy', overwrite: true

    cpus 8
    memory 32G
    time '24h'

    input:
        tuple val(sample_id), path(hybrids)

    output:
        tuple val(sample_id), path("${sample_id}.mfe.clusters.tsv.gz"), emit: hybrids

    script:

    percent_overlap = params.percent_overlap
    sample_size = params.sample_size

    """
    #!/usr/bin/env Rscript

    suppressPackageStartupMessages(library(data.table))
    suppressPackageStartupMessages(library(primavera))
    suppressPackageStartupMessages(library(parallel))

    setDTthreads(${task.cpus})

    # Load hybrids
    hybrids.dt <- fread("$hybrids")
    clusters.dt <- cluster_hybrids(hybrids.dt, percent_overlap = $percent_overlap, sample_size = $sample_size, cores = ${task.cpus})
    fwrite(clusters.dt, "${sample_id}.mfe.clusters.tsv.gz", sep = "\t")
    """

}

process COLLAPSE_CLUSTERS {

    tag "${sample_id}"
    publishDir "${params.outdir}/clusters", mode: 'copy', overwrite: true

    cpus 4
    memory 16G
    time '12h'

    input:
        tuple val(sample_id), path(hybrids)

    output:
        tuple val(sample_id), path("${sample_id}.clusters.tsv.gz"), emit: clusters

    script:

    """
    #!/usr/bin/env Rscript

    suppressPackageStartupMessages(library(data.table))
    suppressPackageStartupMessages(library(primavera))

    hybrids.dt <- fread("$hybrids")

    # Collapse clusters
    clusters.dt <- collapse_clusters(hybrids.dt)
    fwrite(clusters.dt, "${sample_id}.clusters.tsv.gz", sep = "\t")

    message("Completed!")
    """
}
