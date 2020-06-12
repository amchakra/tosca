#!/usr/bin/env nextflow

// Specify DSL2
nextflow.preview.dsl = 2

process convertcoordinates {

    tag "${sample_id}"
    publishDir "${params.outdir}/bed", mode: 'copy', overwrite: true

    cpus 8
    time '24h'

    input:
        tuple val(sample_id), path(hybrids), path(transcript_gtf)

    output:
        tuple val(sample_id), path("${sample_id}.intragenic_hybrids.mfe.clusters.bed.gz")

    script:
    """
    #!/usr/bin/env Rscript

    suppressPackageStartupMessages(library(optparse))
    suppressPackageStartupMessages(library(data.table))
    suppressPackageStartupMessages(library(primavera))
    suppressPackageStartupMessages(library(tictoc))

    setDTthreads(8)

    genes.gr <- rtracklayer::import.gff2("$transcript_gtf")

    hybrids.dt <- fread("$hybrids")
    seq.dt <- hybrids.dt[overlapping_hybrid %in% c(NA, FALSE)] # Remove overlapping hybrids, NA is genomic orientation
    seq.dt <- seq.dt[grep("Mt", L_seqnames, invert = TRUE)] # Remove MT for now

    # New method
    message("Converting coordinates...")
    tic()
    f_out <- gsub("tsv.gz", "bed", "$hybrids")
    ExportGenomicBED(seq.dt = seq.dt, genes.gr = genes.gr, sam_tag = TRUE, filename = f_out)
    system(paste("pigz", f_out))
    toc()

    message("Completed!")
    """
}