#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl=2

process CLUSTER_HYBRIDS {

    tag "${sample_id}"
    publishDir "${params.outdir}/${type}", mode: 'copy', overwrite: true

    // cpus 8
    // memory 32G
    // time '12h'

    input:
        val(type)
        tuple val(sample_id), path(hybrids)

    output:
        tuple val(sample_id), path("${sample_id}.${type}.clustered.tsv.gz"), emit: hybrids

    script:

    percent_overlap = params.percent_overlap
    sample_size = params.sample_size

    args = ""
    // if(workflow.profile.contains("crick")) { args += " --slurm" }
    if(params.slurm) { args += " --slurm" }

    """
    cluster_hybrids.R \
        --sample ${sample_id} \
        --hybrids ${hybrids} \
        --type ${type} \
        --percent_overlap ${percent_overlap} \
        --sample_size ${sample_size} \
        ${args}
    """

}

process COLLAPSE_CLUSTERS {

    tag "${sample_id}"
    publishDir "${params.outdir}/${type}", mode: 'copy', overwrite: true

    // cpus 4
    // memory 16G
    // time '1h'

    input:
        val(type)
        tuple val(sample_id), path(hybrids)

    output:
        tuple val(sample_id), path("${sample_id}.clusters.tsv.gz"), emit: clusters

    script:

    """
    #!/usr/bin/env Rscript

    suppressPackageStartupMessages(library(data.table))
    suppressPackageStartupMessages(library(toscatools))

    hybrids.dt <- fread("$hybrids")

    # Collapse clusters
    clusters.dt <- collapse_clusters(hybrids.dt)
    fwrite(clusters.dt, "${sample_id}.clusters.tsv.gz", sep = "\t")

    message("Completed!")
    """
}