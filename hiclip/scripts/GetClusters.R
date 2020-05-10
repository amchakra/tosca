# Script to get intragenic clusters
# A. M. Chakrabarti
# 8th May 2020


suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(primavera))
# suppressPackageStartupMessages(library(rslurm))
suppressPackageStartupMessages(library(tictoc))
# suppressPackageStartupMessages(library(parallel))

option_list <- list(make_option(c("-i", "--input"), action = "store", type = "character", help = "Input hybrids file"),
                    make_option(c("-p", "--percent"), action = "store", type = "numeric", help = "Percentage overlap"),
                    make_option(c("-o", "--output"), action = "store", type = "character", help = "Output hybrids file"))

opt_parser = OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

ptm <- proc.time()

# Load hybrids
hybrids.dt <- fread(opt$input)

# Get intragenic hybrids
intragenic.hybrids.dt <- hybrids.dt[L_seqnames == R_seqnames][grep("ENSMUSG", L_seqnames)]
fwrite(intragenic.hybrids.dt, gsub("\\.hybrids\\.", "\\.intragenic_hybrids\\.", opt$input), sep = "\t")

# Get Cluster
message("Clustering...")
tic()

# Split out by gene
intragenic.hybrids.list <- split(intragenic.hybrids.dt, intragenic.hybrids.dt$L_seqnames)
solo.intragenic.hybrids.list <- intragenic.hybrids.list[S4Vectors::elementNROWS(intragenic.hybrids.list) == 1] # Remove solos to add in later
message(length(solo.intragenic.hybrids.list), " genes with one hybrid")
toomany.intragenic.hybrids.list <- intragenic.hybrids.list[S4Vectors::elementNROWS(intragenic.hybrids.list) > 5000] # Remove too many
message(length(toomany.intragenic.hybrids.list), " genes with >5000 hybrids")
intragenic.hybrids.list <- intragenic.hybrids.list[S4Vectors::elementNROWS(intragenic.hybrids.list) > 1]
intragenic.hybrids.list <- intragenic.hybrids.list[S4Vectors::elementNROWS(intragenic.hybrids.list) <= 5000]
message(length(intragenic.hybrids.list), " genes to cluster")

# TODO: add in check for length

library(tictoc)
tic()
intragenic.hybrids.clusters.list <- lapply(1:length(intragenic.hybrids.list), function(i) {

  # message(i)
  ClusterHybrids(intragenic.hybrids.list[[i]], percent_overlap = opt$percent)

})

# Name and id order flipped for genes without clusters, because of merging clusters back in, hence use.names = TRUE
intragenic.hybrids.clusters.dt <- rbindlist(intragenic.hybrids.clusters.list, use.names = TRUE)
solo.intragenic.hybrids.dt <- rbindlist(solo.intragenic.hybrids.list, use.names = TRUE) # Add solos back in
toomany.intragenic.hybrids.dt <- rbindlist(toomany.intragenic.hybrids.list, use.names = TRUE)[, cluster := Inf]
intragenic.hybrids.clusters.dt <- rbind(intragenic.hybrids.clusters.dt, solo.intragenic.hybrids.dt, toomany.intragenic.hybrids.dt, use.names = TRUE, fill = TRUE)
toc()

stopifnot(nrow(intragenic.hybrids.clusters.dt) == nrow(intragenic.hybrids.dt))
fwrite(intragenic.hybrids.clusters.dt, opt$output, sep = "\t")

message("Completed!")

