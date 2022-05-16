#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl=2

process ANALYSE_STRUCTURE {

    tag "${sample_id}"
    publishDir "${params.outdir}/hybrids", mode: 'copy', overwrite: true

    // memory 16G
    // time '12h'

    input:
        tuple val(sample_id), path(hybrids)
        path(transcript_fa)

    output:
        tuple val(sample_id), path("${sample_id}.hybrids.mfe.tsv.gz"), emit: hybrids

    script:

    args = ''
    if ( params.shuffled_mfe ) args += ' --shuffled_mfe '
    if ( params.clusters_only) args += ' --clusters_only'

    """
    analyse_structure.R --hybrids $hybrids --fasta $transcript_fa --output ${sample_id}.hybrids.mfe.tsv.gz $args
    """

}

process CHUNK_SEQUENCES {

    tag "${sample_id}"
    // publishDir "${params.outdir}/${type}", mode: 'copy', overwrite: true

    // memory 16G
    // time '12h'

    input:
        val(type)
        tuple val(sample_id), path(hybrids)
        path(transcript_fa)

    output:
       path("${sample_id}_*.structure.rds"), emit: rds

    script:

    clusters_only = params.clusters_only
    chunk_number = params.chunk_number

    """
    #!/usr/bin/env Rscript

    suppressPackageStartupMessages(library(data.table))
    suppressPackageStartupMessages(library(toscatools))
    suppressPackageStartupMessages(library(parallel))

    setDTthreads(${task.cpus})

    # Load genome
    message("Loading genome...")
    genome.fa <- Biostrings::readDNAStringSet(opt$fasta)
    genome.dt <- data.table(gene_id = names(genome.fa),
                            sequence = as.character(genome.fa))

    hybrids.dt <- fread(opt\$hybrids)
    setkey(hybrids.dt, name)
    stopifnot(!any(duplicated(hybrids.dt\$name))) # Check no duplicates

    # Get sequences
    message("Getting sequences...")
    tic()
    hybrids.dt <- get_sequence(hybrids.dt = hybrids.dt, genome.dt = genome.dt)
    toc()

    stopifnot(!any(is.na(c(hybrids.dt\$L_sequence, hybrids.dt\$R_sequence))))

    # Remove rRNA
    sel.hybrids.dt <- hybrids.dt[!(L_seqnames == "rRNA_45S" & R_seqnames == "rRNA_45S")]
    sel.hybrids.dt <- sel.hybrids.dt[!(L_seqnames == "rDNA" & R_seqnames == "rDNA")]
    sel.hybrids.dt <- sel.hybrids.dt[!(L_seqnames == "rRNA_5S" & R_seqnames == "rRNA_5S")]

    # if(clusters_only == "true") sel.hybrids.dt <- sel.hybrids.dt[!is.na(cluster)][cluster != "."] # Some are "" if unclustered as too many
    if(clusters_only == "true") sel.hybrids.dt <- sel.hybrids.dt[!is.na(cluster)][grep("^C", cluster)]

    structure.dt <- sel.hybrids.dt[, .(name, L_sequence, R_sequence)]

    if($chunk_number > nrow(structure.dt)) {
        structure.list.chunks <- split(structure.dt, cut(seq_len(nrow(structure.dt)), nrow(structure.dt), label = FALSE))
        lapply(seq_len(nrow(structure.dt)), function(i) { saveRDS(atlas.hybrids.list.chunks[[i]], paste0("${sample_id}", "_", i, ".rds")) })
    }
    if($chunk_number > 1) {
        structure.list.chunks <- split(structure.dt, cut(seq_len(nrow(structure.dt)), $chunk_number, label = FALSE))
        lapply(seq_len($chunk_number), function(i) { saveRDS(structure.list.chunks[[i]], paste0("${sample_id}", "_", i, "structure.rds")) })
    } else {
        structure.list.chunks <- structure.dt
        saveRDS(structure.list.chunks, paste0("${sample_id}", "_", 1, ".structure.rds"))
    }
    """

}

process CALCULATE_STRUCTURES {

    tag "${sample_id}"
    // publishDir "${params.outdir}/${type}", mode: 'copy', overwrite: true

    // memory 16G
    // time '12h'

    input:
        val(type)
        tuple val(sample_id), path(rds)

    output:
       tuple val(sample_id), path("${sample_id}.structures.tsv.gz"), emit: tsv

    script:

    """
    #!/usr/bin/env Rscript

    suppressPackageStartupMessages(library(data.table))
    suppressPackageStartupMessages(library(toscatools))
    suppressPackageStartupMessages(library(parallel))

    setDTthreads(${task.cpus})

    sel.hybrids.dt <- readRDS("$rds")

    structure.list <- parallel::mclapply(seq_len(nrow(sel.hybrids.dt)), function(i) {

        analyse_structure(name = sel.hybrids.dt\$name[i], L_sequence = sel.hybrids.dt\$L_sequence[i], R_sequence = sel.hybrids.dt\$L_sequence[i])
    
    }, mc.cores = ${task.cpus})

    structure.dt <- rbindlist(structure.list, use.names = TRUE)
    fwrite(atlas.clusters.dt, paste0("${sample_id}", ".structures.tsv.gz"), sep = "\t")

    """

}

process CALCULATE_SHUFFLED_ENERGIES {

    tag "${sample_id}"
    // publishDir "${params.outdir}/${type}", mode: 'copy', overwrite: true

    // memory 16G
    // time '12h'

    input:
        val(type)
        tuple val(sample_id), path(rds)

    output:
       tuple val(sample_id), path("${sample_id}.shuffled.tsv.gz"), emit: tsv

    script:

    """
    #!/usr/bin/env Rscript

    suppressPackageStartupMessages(library(data.table))
    suppressPackageStartupMessages(library(toscatools))
    suppressPackageStartupMessages(library(parallel))

    setDTthreads(${task.cpus})

    sel.hybrids.dt <- readRDS("$rds")

    shuffled.list <- parallel::mclapply(seq_len(nrow(sel.hybrids.dt)), function(i) {

        get_shuffled_mfe(name = sel.hybrids.dt\$name[i], L_sequence = sel.hybrids.dt\$L_sequence[i], R_sequence = sel.hybrids.dt\$L_sequence[i])
    
    }, mc.cores = ${task.cpus})

    shuffled.dt <- rbindlist(shuffled.list, use.names = TRUE)
    fwrite(atlas.clusters.dt, paste0("${sample_id}", ".shuffled.tsv.gz"), sep = "\t")

    """

}

process MERGE_STRUCTURES {

    tag "${sample_id}"
    publishDir "${params.outdir}/${type}", mode: 'copy', overwrite: true

    // memory 16G
    // time '12h'

    input:
        val(type)
        tuple val(sample_id), path(hybrids), path(structures)

    output:
       tuple val(sample_id), path("${sample_id}.${type}.gc.annotated.mfe.tsv.gz"), emit: tsv

    script:

    """
    #!/usr/bin/env Rscript

    suppressPackageStartupMessages(library(data.table))
    suppressPackageStartupMessages(library(toscatools))

    setDTthreads(${task.cpus})

    hybrids.dt <- fread("$hybrids")
    structures.files <-  strsplit("$structures", " ")[[1]]
    structures.list <- lapply(structures.files, fread)
    structures.dt <- rbindlist(structures.list, use.names = TRUE, fill = TRUE)

    hybrids.dt <- merge(hybrids.dt, structures.dt, by = "name", all.x = TRUE)
    fwrite(atlas.clusters.dt, paste0("${sample_id}", ".", "${type}", ".gc.annotated.mfe.tsv.gz"), sep = "\t")

    """

}

process MERGE_SHUFFLED {

    tag "${sample_id}"
    publishDir "${params.outdir}/${type}", mode: 'copy', overwrite: true

    // memory 16G
    // time '12h'

    input:
        val(type)
        tuple val(sample_id), path(hybrids), path(shuffled)

    output:
       tuple val(sample_id), path("${sample_id}.${type}.gc.annotated.mfe.shuffled.tsv.gz"), emit: tsv

    script:

    """
    #!/usr/bin/env Rscript

    suppressPackageStartupMessages(library(data.table))
    suppressPackageStartupMessages(library(toscatools))

    setDTthreads(${task.cpus})

    hybrids.dt <- fread("$hybrids")
    shuffled.files <-  strsplit("$shuffled", " ")[[1]]
    shuffled.list <- lapply(shuffled.files, fread)
    shuffled.dt <- rbindlist(shuffled.list, use.names = TRUE, fill = TRUE)

    hybrids.dt <- merge(hybrids.dt, shuffled.dt, by = "name", all.x = TRUE)
    fwrite(atlas.clusters.dt, paste0("${sample_id}", ".", "${type}", ".gc.annotated.mfe.shuffled.tsv.gz"), sep = "\t")

    """

}