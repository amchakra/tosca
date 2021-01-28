#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl=2

process ANNOTATE_HYBRIDS {

    tag "${sample_id}"
    publishDir "${params.outdir}/${type}", mode: 'copy', overwrite: true

    cpus 8
    memory 32G
    time '24h'

    input:
        val(type)
        tuple val(sample_id), path(hybrids)
        path(regions_gtf)

    output:
        tuple val(sample_id), path("${sample_id}.annotated.tsv.gz"), emit: hybrids

    script:

    """
    #!/usr/bin/env Rscript

    suppressPackageStartupMessages(library(data.table))
    suppressPackageStartupMessages(library(primavera))

    setDTthreads(${task.cpus})

    hybrids.dt <- fread("$hybrids")
    regions.gr <- rtracklayer::import.gff2("$regions_gtf")
    hybrids.dt <- annotate_hybrids(hybrids.dt, regions.gr)

    fwrite(hybrids.dt, "${sample_id}.annotated.tsv.gz", sep = "\t")
    """

}