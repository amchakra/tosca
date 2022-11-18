#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl=2

process CLUSTER_HYBRIDS_SLURM {

    tag "${sample_id}"
    publishDir "${params.outdir}/${type}", mode: 'copy', overwrite: true

    input:
        val(type)
        tuple val(sample_id), path(hybrids)

    output:
        tuple val(sample_id), path("${sample_id}.${type}.slurm.clustered.tsv.gz"), emit: hybrids

    script:

    percent_overlap = params.percent_overlap
    sample_size = params.sample_size

    args = ""
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

process CHUNK_HYBRIDS {

    tag "${sample_id}"

    input:
        val(type)
        tuple val(sample_id), path(hybrids)

    output:
        tuple val(sample_id), path("${sample_id}.unclustered.tsv.gz"), emit: tsv
        path("${sample_id}_*.rds"), emit: rds

    script:

    sample_size = params.sample_size
    chunk_number = params.chunk_number

    """
    #!/usr/bin/env Rscript

    suppressPackageStartupMessages(library(data.table))
    suppressPackageStartupMessages(library(toscatools))
    suppressPackageStartupMessages(library(parallel))

    setDTthreads("${task.cpus}")
    set.seed(42)

    # Load hybrids
    hybrids.dt <- fread("$hybrids")
    setkey(hybrids.dt, name)

    # Filter hybrids
    message("Number of hybrids: ", nrow(hybrids.dt))
    message("Removing rRNA-rRNA")
    atlas.hybrids.dt <- hybrids.dt[, total_count := .N, by = .(L_seqnames, R_seqnames)]
    atlas.hybrids.dt <- atlas.hybrids.dt[!(L_seqnames == "rRNA_45S" & R_seqnames == "rRNA_45S")]
    atlas.hybrids.dt <- atlas.hybrids.dt[!(L_seqnames == "rDNA" & R_seqnames == "rDNA")]
    atlas.hybrids.dt <- atlas.hybrids.dt[!(L_seqnames == "rRNA_5S" & R_seqnames == "rRNA_5S")]

    message(nrow(atlas.hybrids.dt[total_count > 1e4]), " high incidence (>10,000) gene pairs not clustered")
    # atlas.hybrids.dt <- atlas.hybrids.dt[total_count < 1e4 & total_count > 1]
    atlas.hybrids.dt <- atlas.hybrids.dt[total_count > 1]

    # Subsample as indicated
    if($sample_size != -1) atlas.hybrids.dt <- atlas.hybrids.dt[sample(1:nrow(atlas.hybrids.dt), $sample_size)]
    message("Number of hybrids to cluster: ", nrow(atlas.hybrids.dt))

    # Keep ones not clustered to add back in later
    unclustered.hybrids.dt <- hybrids.dt[!name %in% atlas.hybrids.dt\$name]
    stopifnot(nrow(unclustered.hybrids.dt) + nrow(atlas.hybrids.dt) == nrow(hybrids.dt))

    fwrite(unclustered.hybrids.dt, paste0("${sample_id}", ".unclustered.tsv.gz"), sep = "\t")

    # Split into list to parallelise
    atlas.hybrids.list <- split(atlas.hybrids.dt, by = c("L_seqnames", "R_seqnames"))
    message("Gene pairs to cluster: ", length(atlas.hybrids.list))

    # Split into chunks and write out
    if($chunk_number < length(atlas.hybrids.list)) {
        atlas.hybrids.list.chunks <-  split(atlas.hybrids.list, cut(seq_along(atlas.hybrids.list), length(atlas.hybrids.list), label = FALSE))
        lapply(seq_len(length(atlas.hybrids.list)), function(i) { saveRDS(atlas.hybrids.list.chunks[[i]], paste0("${sample_id}", "_", i, ".rds")) })
    } else if($chunk_number > 1) {
        atlas.hybrids.list.chunks <- split(atlas.hybrids.list, cut(seq_along(atlas.hybrids.list), $chunk_number, label = FALSE))
        lapply(seq_len($chunk_number), function(i) { saveRDS(atlas.hybrids.list.chunks[[i]], paste0("${sample_id}", "_", i, ".rds")) })
    } else {
        atlas.hybrids.list.chunks <- atlas.hybrids.list
        saveRDS(atlas.hybrids.list.chunks, paste0("${sample_id}", "_", 1, ".rds"))
    }
    """
}

process IDENTIFY_CLUSTERS {

    tag "${sample_id}"

    input:
        val(type)
        tuple val(sample_id), path(rds)

    output:
        tuple val(sample_id), path("${sample_id}.clustered.tsv.gz"), emit: tsv
    
    script:

    percent_overlap = params.percent_overlap

    """
    #!/usr/bin/env Rscript

    suppressPackageStartupMessages(library(data.table))
    suppressPackageStartupMessages(library(toscatools))
    suppressPackageStartupMessages(library(parallel))

    setDTthreads(${task.cpus})
    set.seed(42)

    atlas.hybrids.list <- readRDS("$rds")
    atlas.clusters.list <- parallel::mclapply(atlas.hybrids.list, cluster_hybrids, percent_overlap = $percent_overlap, mc.cores = ${task.cpus})
    atlas.clusters.dt <- rbindlist(atlas.clusters.list, use.names = TRUE, fill = TRUE)

    fwrite(atlas.clusters.dt, paste0("${sample_id}", ".clustered.tsv.gz"), sep = "\t")

    """
}

process MERGE_CLUSTERS {

    tag "${sample_id}"
    publishDir "${params.outdir}/${type}", mode: 'copy', overwrite: true

    input:
        val(type)
        tuple val(sample_id), path(hybrids), path(unclustered), path(clustered)

    output:
        tuple val(sample_id), path("${sample_id}.${type}.clustered.tsv.gz"), emit: hybrids
    
    script:

    """
    #!/usr/bin/env Rscript

    suppressPackageStartupMessages(library(data.table))
    suppressPackageStartupMessages(library(toscatools))

    setDTthreads(${task.cpus})
    
    hybrids.dt <- fread("$hybrids")
    unclustered.hybrids.dt <- fread("$unclustered")

    clusters.files <- strsplit("$clustered", " ")[[1]]
    clusters.list <- lapply(clusters.files, fread)
    clusters.dt <- rbindlist(clusters.list, use.names = TRUE, fill = TRUE)

    clusters.dt <- rbindlist(list(clusters.dt, unclustered.hybrids.dt), use.names = TRUE, fill = TRUE)
    setorder(clusters.dt, name)

    stopifnot(nrow(clusters.dt) == nrow(hybrids.dt))
    stopifnot(all(clusters.dt\$name %in% hybrids.dt\$name))

    fwrite(clusters.dt, paste0("${sample_id}", ".", "${type}", ".clustered.tsv.gz"), sep = "\t")

    """
}