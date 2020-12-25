#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl=2

process IDENTIFY_HYBRIDS {

    tag "${sample_id}"
    // publishDir "${params.outdir}/hybrids/split/", mode: 'copy', overwrite: true

    time '24h'
    cpus = 8
    memory '64 G'

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

    fwrite(hybrids.dt, file = paste0("$sample_id", ".hybrids.tsv.gz"), sep = "\t")
    """

}

process IDENTIFY_HYBRIDS_2 {

    tag "${sample_id}"
    // publishDir "${params.outdir}/hybrids/split/", mode: 'copy', overwrite: true

    time '24h'
    cpus = 8
    memory '64 G'

    input:
        tuple val(sample_id), path(blast8), path(reads)

    output:
        tuple val(sample_id), path("${sample_id}.hybrids.tsv.gz"), emit: hybrids

    script:

    """
    #!/usr/bin/env Rscript

    suppressPackageStartupMessages(library(optparse))
    suppressPackageStartupMessages(library(data.table))
    suppressPackageStartupMessages(library(primavera))
    suppressPackageStartupMessages(library(parallel))

    opt <- list(blast = "$blast8",
            fasta = "$reads",
            output = "${sample_id}.hybrids.tsv.gz",
            threads = 8)

    # Load blast and get read lengths
    blast.dt <- load_blast8(opt\$blast)
    blast.dt <- add_read_lengths(blast.dt = blast.dt, fasta = opt\$fasta)
    blast.dt <- calculate_blast8_metrics(blast.dt = blast.dt)

    blast.list <- split(blast.dt, blast.dt\$query)
    cl <- makeForkCluster(opt\$threads) # otherwise really slow...
    hybrids.list <- parLapply(cl = cl, blast.list, function(x) get_valid_hybrids(blast.query.dt = x))
    stopCluster(cl)

    message(sum(S4Vectors::elementNROWS(hybrids.list) == 0), " out of ", length(hybrids.list), " reads did not have hybrids")
    message(round(sum(S4Vectors::elementNROWS(hybrids.list) != 0)/length(hybrids.list), 4) * 100, "% of reads had hybrids")

    # Filter multi hits
    hybrids.dt <- rbindlist(hybrids.list)
    hybrids.dt <- reorient_hybrids(hybrids.dt)

    hybrids.dt[, multi := .N, by = query]
    single.hybrids.dt <- hybrids.dt[multi == 1]
    multi.hybrids.dt <- hybrids.dt[multi > 1]

    ol <- find_hybrid_overlaps(multi.hybrids.dt, single.hybrids.dt)
    multi.hybrids.dt <- multi.hybrids.dt[unique(ol\$queryHits)]
    multi.hybrids.dt <- multi.hybrids.dt[, N := .N, by = query][N == 1] # Keep only those with a unique mathc
    multi.hybrids.dt[, N := NULL]

    single.hybrids.dt[, hybrid_selection := "single"]
    multi.hybrids.dt[, hybrid_selection := "multi_overlap"]
    valid.hybrids.dt <- rbindlist(list(single.hybrids.dt, multi.hybrids.dt), use.names = TRUE)

    other.hybrids.dt <- hybrids.dt[!query %in% valid.hybrids.dt\$query]
    other.hybrids.dt[, hybrid_selection := "multi"]

    fwrite(valid.hybrids.dt, file = opt\$output, sep = "\t", col.names = TRUE)

    """

}

process MERGE_HYBRIDS_2 {

    tag "${sample_id}"
    publishDir "${params.outdir}/hybrids_2", mode: 'copy', overwrite: true

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

    fwrite(hybrids.dt, file = paste0("$sample_id", ".hybrids.tsv.gz"), sep = "\t")
    """

}