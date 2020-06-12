#!/usr/bin/env nextflow

// Specify DSL2
nextflow.preview.dsl = 2

process collapseclusters {

    tag "${sample_id}"
    publishDir "${params.outdir}/clusters", mode: 'copy', overwrite: true

    cpus 8
    time '24h'

    input:
        tuple val(sample_id), path(hybrids), path(transcript_gtf)

    output:
        tuple val(sample_id), path("${sample_id}.clusters.tsv.gz"), path("${sample_id}.clusters.bed.gz")

    script:
    """
    #!/usr/bin/env Rscript

    suppressPackageStartupMessages(library(optparse))
    suppressPackageStartupMessages(library(data.table))
    suppressPackageStartupMessages(library(primavera))
    suppressPackageStartupMessages(library(tictoc))

    setDTthreads(8)

    hybrids.dt <- fread("$hybrids")

    clusters.dt <- CollapseClusters(hybrids.dt)
    f_out <- gsub(".intragenic_hybrids.mfe.clusters.tsv.gz", ".clusters.tsv.gz", "$hybrids")
    fwrite(clusters.dt, f_out, sep = "\t")

    message(nrow(clusters.dt[L_end >= R_start]), " clusters removed")
    clusters.dt <- clusters.dt[L_end < R_start]

    # ========== SPLIT OUT?

    # New method
    message("Converting coordinates...")
    genes.gr <- rtracklayer::import.gff2("$transcript_gtf")
    tic()
    f_out <- gsub("tsv.gz", "bed", f_out)
    ExportGenomicBED(seq.dt = clusters.dt, genes.gr = genes.gr, sam_tag = TRUE, filename = gsub("tsv.gz", "bed", f_out))
    system(paste("pigz", f_out))
    toc()

    message("Completed!")
    """
}
