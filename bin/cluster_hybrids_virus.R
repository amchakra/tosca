#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(primavera))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(rslurm))


option_list <- list(make_option(c("", "--hybrids"), action = "store", type = "character", help = "Hybrids data.table"),
            make_option(c("", "--type"), action = "store", type = "character", help = "Type of hybrids data.table"),
            make_option(c("", "--sample"), action = "store", type = "character", help = "Sample id"),
            make_option(c("", "--percent_overlap"), action = "store", type = "numeric", help = "Percent overlap"),
            make_option(c("", "--sample_size"), action = "store", type = "integer", help = "Sample size"),
            make_option(c("", "--threads"), action = "store", type = "integer", default = 8, help = "Number of threads [default: %default]"),
            make_option(c("", "--slurm"), action = "store_true", type = "logical", default = FALSE, help = "Use SLURM [default: %default]"))
opt_parser = OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

setDTthreads(opt$threads)
set.seed(42)

# ==========

find_hybrid_overlaps2 <- function(hybrids.dt, percent_overlap, verbose = TRUE) {
  hybrids.dt <- primavera::reorient_hybrids(hybrids.dt)

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
  cmd <- paste("bedtools pairtopair -rdn -is -f", percent_overlap, "-a", bedpe, "-b", bedpe, ">", ol)
  if (verbose) message(cmd)
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

# ==========

cluster_hybrids2 <- function(hybrids.dt, percent_overlap = 0.75, verbose = TRUE) {
  hybrids.bedpe.dt <- find_hybrid_overlaps2(hybrids.dt, percent_overlap = percent_overlap, verbose = verbose)

  if (nrow(hybrids.bedpe.dt) == 0) {
    return(hybrids.dt[, cluster := as.character(NA)])
  }

  sel.bedpe.dt <- hybrids.bedpe.dt[L_p > percent_overlap & R_p > percent_overlap]

  # igraph
  g <- igraph::graph_from_edgelist(el = as.matrix(sel.bedpe.dt[, .(name.x, name.y)]), directed = FALSE)
  igraph::E(g)$weight <- sel.bedpe.dt$mean_p # weight by percent overlap

  c <- igraph::components(g)
  if (verbose) message(c$no, " clusters")

  clusters.dt <- data.table(
    name = names(c$membership),
    cluster = c$membership
  )
  setorder(clusters.dt, cluster)

  # Merge back
  if (nrow(clusters.dt) == 0) clusters.dt[, name := character()] # In case there are no clusters
  setkey(clusters.dt, name)
  if (nrow(clusters.dt) != 0) stopifnot(any(!duplicated(clusters.dt$name))) # Make sure no hybrid is in more than one cluster, but only if there are clusters
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

# ==========

# Load hybrids
hybrids.dt <- fread(opt$hybrids)
setkey(hybrids.dt, name)

# Filter hybrids
message("Number of hybrids: ", nrow(hybrids.dt))
atlas.hybrids.dt <- hybrids.dt[, total_count := .N, by = .(L_seqnames, R_seqnames)]
atlas.hybrids.dt <- atlas.hybrids.dt[total_count > 1]

# Subsample as indicated
if(opt$sample_size != -1) atlas.hybrids.dt <- atlas.hybrids.dt[sample(1:nrow(atlas.hybrids.dt, opt$sample_size))]
message("Number of hybrids to cluster: ", nrow(atlas.hybrids.dt))

# Keep ones not clustered to add back in later
unclustered.hybrids.dt <- hybrids.dt[!name %in% atlas.hybrids.dt$name]
stopifnot(nrow(unclustered.hybrids.dt) + nrow(atlas.hybrids.dt) == nrow(hybrids.dt))

# Split into list to parallelise
atlas.hybrids.list <- split(atlas.hybrids.dt, by = c("L_seqnames", "R_seqnames"))
message("Gene pairs to cluster: ", length(atlas.hybrids.list))

if(opt$slurm) {

    # Submit to cluster
    sjob <- slurm_map(atlas.hybrids.list, f = cluster_hybrids2, percent_overlap = opt$percent_overlap, global_objects = c("cluster_hybrids2", "find_hybrid_overlaps2"), jobname = sapply(strsplit(basename(opt$hybrids), "\\."), "[[", 1), nodes = 100, cpus_per_node = 1, slurm_options = list(time = "12:00:00", mem = "64G"))
    status <- FALSE
    while(status == FALSE) {
        squeue.out <- system(paste("squeue -n", sjob$jobname), intern = TRUE) # Get contents of squeue for this job
        if(length(squeue.out) == 1) status <- TRUE # i.e. only the header left
        Sys.sleep(60)
    }

    atlas.clusters.list <- get_slurm_out(sjob)
    cleanup_files(sjob) 

} else {

    atlas.clusters.list <- parallel::mclapply(atlas.hybrids.list, cluster_hybrids, percent_overlap = opt$percent_overlap, mc.cores = opt$threads)
    
}

atlas.clusters.dt <- rbindlist(atlas.clusters.list, use.names = TRUE, fill = TRUE)
stopifnot(nrow(atlas.clusters.dt) == nrow(atlas.hybrids.dt))

# Add back in unclustered
atlas.clusters.dt <- rbindlist(list(atlas.clusters.dt, unclustered.hybrids.dt), use.names = TRUE, fill = TRUE)

stopifnot(nrow(atlas.clusters.dt) == nrow(hybrids.dt))
fwrite(atlas.clusters.dt, paste0(opt$sample, ".", opt$type, ".clustered.tsv.gz"), sep = "\t")