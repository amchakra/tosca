#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl=2

process CLUSTER_HYBRIDS_SLURM {

    tag "${sample_id}"
    label 'process_medium'
    
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
    label 'process_low'

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
    label 'process_medium'

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

    # Inter-transcript only
    inter_only = TRUE
    if(inter_only) {
        atlas.hybrids.dt <- atlas.hybrids.dt[L_seqnames != R_seqnames]
        print(unique(atlas.hybrids.dt[, .(L_seqnames, R_seqnames, total_count)]))
    }

    # Keep ones not clustered to add back in later
    unclustered.hybrids.dt <- hybrids.dt[!name %in% atlas.hybrids.dt\$name]
    stopifnot(nrow(unclustered.hybrids.dt) + nrow(atlas.hybrids.dt) == nrow(hybrids.dt))

    fwrite(unclustered.hybrids.dt, paste0("${sample_id}", ".unclustered.tsv.gz"), sep = "\t")

    # Split into list to parallelise
    atlas.hybrids.list <- split(atlas.hybrids.dt, by = c("L_seqnames", "R_seqnames"))
    message("Gene pairs to cluster: ", length(atlas.hybrids.list))

    # Split into chunks and write out
    if($chunk_number > length(atlas.hybrids.list)) {
        atlas.hybrids.list.chunks <-  split(atlas.hybrids.list, seq_along(atlas.hybrids.list))
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
    label 'process_high'

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

    # ==========
    # Functions
    # ==========

    find_hybrid_overlaps_fraction <- function(hybrids.dt, fraction_overlap) {

        hybrids.dt <- toscatools::reorient_hybrids(hybrids.dt)

        # Create BEDPE and get overlaps
        bedpe.colnames <- c("L_seqnames", "L_start", "L_end", "R_seqnames", "R_start", "R_end", "name", "total_count", "L_strand", "R_strand")
        bedpe.dt <- hybrids.dt[, ..bedpe.colnames]
        bedpe.dt[, `:=`(
            L_start = L_start - 1,
            R_start = R_start - 1
        )]

        bedpe <- tempfile(tmpdir = getwd(), fileext = ".bedpe")
        ol <- tempfile(tmpdir = getwd(), fileext = ".bedpe")

        fwrite(bedpe.dt, file = bedpe, sep = "\t", col.names = FALSE)
        cmd <- paste("bedtools pairtopair -rdn -f ", fraction_overlap, " -a", bedpe, "-b", bedpe, ">", ol)
        message(cmd)
        system(cmd)

        # Check if there are no overlaps
        if (file.size(ol) != 0) {
            bedpe.dt <- fread(ol, col.names = c(paste0(bedpe.colnames, ".x"), paste0(bedpe.colnames, ".y")))
            # Delete temporary files
            invisible(file.remove(bedpe))
            invisible(file.remove(ol))
        } else {

            # Delete temporary files
            invisible(file.remove(bedpe))
            invisible(file.remove(ol))
            return(data.table())
        }

        # Get calculations and filter
        bedpe.dt[, `:=`(
            L_ol = min(L_end.x, L_end.y) - max(L_start.x, L_start.y) + 1,
            R_ol = min(R_end.x, R_end.y) - max(R_start.x, R_start.y) + 1
        ),
        by = .(name.x, name.y)
        ]

        bedpe.dt[, `:=`(
            L_sp = max(L_end.x, L_end.y) - min(L_start.x, L_start.y) + 1,
            R_sp = max(R_end.x, R_end.y) - min(R_start.x, R_start.y) + 1
        ),
        by = .(name.x, name.y)
        ]

        bedpe.dt[, `:=`(
            L_p = L_ol / L_sp,
            R_p = R_ol / R_sp
        ),
        by = .(name.x, name.y)
        ]

        bedpe.dt[, mean_p := mean(c(L_p, R_p)), by = .(name.x, name.y)]
        return(bedpe.dt)
    }

    cluster_hybrids_fraction <- function(hybrids.dt, percent_overlap = 0.75, verbose = FALSE) {

        hybrids.bedpe.dt <- find_hybrid_overlaps_fraction(hybrids.dt, fraction_overlap = percent_overlap)

        if (nrow(hybrids.bedpe.dt) == 0) {
            return(hybrids.dt[, cluster := as.character(NA)])
        }

        sel.bedpe.dt <- hybrids.bedpe.dt[L_p > percent_overlap & R_p > percent_overlap]

        # igraph
        g <- igraph::graph_from_edgelist(el = as.matrix(sel.bedpe.dt[, .(name.x, name.y)]), directed = FALSE)
        igraph::E(g)\$weight <- sel.bedpe.dt\$mean_p # weight by percent overlap

        c <- igraph::components(g)
        if (verbose) message(c\$no, " clusters")

        clusters.dt <- data.table(
            name = names(c\$membership),
            cluster = c\$membership
        )
        setorder(clusters.dt, cluster)

        # Merge back
        if (nrow(clusters.dt) == 0) clusters.dt[, name := character()] # In case there are no clusters
        setkey(clusters.dt, name)
        if (nrow(clusters.dt) != 0) stopifnot(any(!duplicated(clusters.dt\$name))) # Make sure no hybrid is in more than one cluster, but only if there are clusters
        clusters.dt[, tempcluster := paste0("C", cluster)][, cluster := NULL]

        # Order cluster names by number of hybrids
        clusters.order.dt <- clusters.dt[, .N, by = tempcluster]
        setorder(clusters.order.dt, -N, tempcluster)[, cluster := paste0("C", stringr::str_pad(1:.N, width = 3, pad = 0))]
        setnames(clusters.order.dt, "N", "cluster_hybrid_count")
        clusters.dt <- merge(clusters.dt, clusters.order.dt, by = "tempcluster")
        clusters.dt[, tempcluster := NULL]

        # Merge back
        setkey(hybrids.dt, name)
        if ("cluster" %in% names(hybrids.dt)) hybrids.dt[, cluster := NULL] # if clusters already assigned, remove them
        hybrids.clustered.dt <- merge(hybrids.dt, clusters.dt, by = "name", all.x = TRUE)
        hybrids.clustered.dt[is.na(cluster), cluster := "."]

        return(hybrids.clustered.dt)
    }

    atlas.hybrids.list <- readRDS("$rds")
    # atlas.clusters.list <- parallel::mclapply(atlas.hybrids.list, cluster_hybrids, percent_overlap = $percent_overlap, mc.cores = ${task.cpus})
    atlas.clusters.list <- parallel::mclapply(atlas.hybrids.list, cluster_hybrids_fraction, percent_overlap = ${percent_overlap}, mc.cores = ${task.cpus})
    atlas.clusters.dt <- rbindlist(atlas.clusters.list, use.names = TRUE, fill = TRUE)

    fwrite(atlas.clusters.dt, paste0("${sample_id}", ".clustered.tsv.gz"), sep = "\t")

    """
}

process MERGE_CLUSTERS {

    tag "${sample_id}"
    label 'process_medium'
    
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