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
    suppressPackageStartupMessages(library(primavera))

    # print("$hybrids")
    hybrids.files <- strsplit("$hybrids", " ")[[1]]
    hybrids.list <- lapply(hybrids.files, fread)
    hybrids.dt <- rbindlist(hybrids.list, use.names = TRUE)

    # Get orientations where appropriate
    hybrids.dt <- reorient_hybrids(unique.hybrid.dt)
    hybrids.dt[L_seqnames == R_seqnames, orientation := ifelse(L_read_start <= R_read_start, "genomic", "reverse")]

    # Refashion hybrid.dt
    hybrids.dt[, `:=` (L_width = L_end - L_start + 1,
                        L_strand = "+",
                        R_width = R_end - R_start + 1,
                        R_strand = "+")]

    fwrite(hybrids.dt, file = paste0("$sample_id", ".hybrids.tsv.gz"), sep = "\t")
    """

}