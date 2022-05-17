#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl=2

process CONVERT_COORDINATES {

    tag "${sample_id}"
    publishDir "${params.outdir}/${type}", mode: 'copy', overwrite: true

    input:
        val(type)
        tuple val(sample_id), path(hybrids)
        path(transcript_gtf)

    output:
        tuple val(sample_id), path("${sample_id}.${type}.gc.tsv.gz"), emit: hybrids

    script:
    """
    #!/usr/bin/env Rscript

    suppressPackageStartupMessages(library(data.table))
    suppressPackageStartupMessages(library(toscatools))

    # Load data
    genes.gr <- rtracklayer::import.gff2("$transcript_gtf")
    hybrids.dt <- fread("$hybrids")
    
    # hybrids.dt <- hybrids.dt[grep("^rRNA", L_seqnames, invert = TRUE)] # Remove rRNA
    # hybrids.dt <- hybrids.dt[grep("Mt", L_seqnames, invert = TRUE)] # Remove MT for now

    coord.dt <- convert_coordinates(hybrids.dt = hybrids.dt, genes.gr = genes.gr)
    fwrite(coord.dt, "${sample_id}.${type}.gc.tsv.gz", sep = "\t")
    """
    
}