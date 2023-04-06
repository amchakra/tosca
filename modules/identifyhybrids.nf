#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl=2

process IDENTIFY_HYBRIDS {

    tag "${sample_id}"
    label 'process_high'

    input:
        tuple val(sample_id), path(blast8), path(reads)

    output:
        tuple val(sample_id), path("${sample_id}.hybrids.tsv.gz"), emit: hybrids

    script:

    """
    identify_hybrids.R -t ${task.cpus} -b $blast8 -f $reads -o ${sample_id}.hybrids.tsv.gz
    """

}

process MERGE_HYBRIDS {

    tag "${sample_id}"
    label 'process_medium'
    
    publishDir "${params.outdir}/${type}", mode: 'copy', overwrite: true

    input:
        val(type)
        tuple val(sample_id), path(hybrids)

    output:
        tuple val(sample_id), path("${sample_id}.hybrids.tsv.gz"), emit: hybrids

    script:
    
    // zcat $hybrids | pigz > ${sample_id}.hybrids.tsv.gz
    // This doesn't account for empty table files 

    """
    #!/usr/bin/env Rscript

    suppressPackageStartupMessages(library(data.table))
    suppressPackageStartupMessages(library(toscatools))

    # print("$hybrids")
    hybrids.files <- strsplit("$hybrids", " ")[[1]]
    hybrids.list <- lapply(hybrids.files, fread)
    hybrids.dt <- rbindlist(hybrids.list, use.names = TRUE)

    if("$type" == "atlas") hybrids.dt\$sample <- rep(basename(hybrids.files), S4Vectors::elementNROWS(hybrids.list))

    fwrite(hybrids.dt, file = "${sample_id}.hybrids.tsv.gz", sep = "\t")
    """

}