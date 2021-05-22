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

# Load hybrids
hybrids.dt <- fread(opt$hybrids)
setkey(hybrids.dt, name)

# Filter hybrids
message("Number of hybrids: ", nrow(hybrids.dt))
message("Removing rRNA-rRNA")
atlas.hybrids.dt <- hybrids.dt[, total_count := .N, by = .(L_seqnames, R_seqnames)]
atlas.hybrids.dt <- atlas.hybrids.dt[!(L_seqnames == "rRNA_45S" & R_seqnames == "rRNA_45S")]
atlas.hybrids.dt <- atlas.hybrids.dt[!(L_seqnames == "rDNA" & R_seqnames == "rDNA")]
atlas.hybrids.dt <- atlas.hybrids.dt[!(L_seqnames == "rRNA_5S" & R_seqnames == "rRNA_5S")]

message(nrow(atlas.hybrids.dt[total_count > 1e4]), " high incidence (>10,000) gene pairs not clustered")
atlas.hybrids.dt <- atlas.hybrids.dt[total_count < 1e4 & total_count > 1]

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
    sjob <- slurm_map(atlas.hybrids.list, f = cluster_hybrids, percent_overlap = opt$percent_overlap, jobname = sapply(strsplit(basename(opt$hybrids), "\\."), "[[", 1), nodes = 100, cpus_per_node = 8, slurm_options = list(time = "12:00:00", mem = "64G"))
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