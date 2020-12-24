#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl=2

process MERGE_HYBRIDS {

    tag "${sample_id}"
    publishDir "${params.outdir}/hybrids", mode: 'copy', overwrite: true

    time '24h'

    input:
        tuple val(sample_id), path(hybrids)

    output:
        tuple val(sample_id), path("${sample_id}.hybrids.tsv.gz"), emit: hybrids

    script:
    
    // zcat $hybrids | pigz > ${sample_id}.hybrids.tsv.gz
    // This doesn't account for empty table files 

    """
    #!/usr/bin/env Rscript

    suppressPackageStartupMessages(library(data.table))

    # print("$hybrids")
    hybrids.files <- strsplit("$hybrids", " ")[[1]]
    hybrids.list <- lapply(hybrids.files, fread)
    hybrids.dt <- rbindlist(hybrids.list, use.names = TRUE)

    fwrite(hybrids.dt, file = paste0("$sample_id", ".hybrids.tsv.gz"), sep = "\t")
    """

}