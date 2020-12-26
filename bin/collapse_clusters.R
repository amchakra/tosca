#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(primavera))
suppressPackageStartupMessages(library(parallel))

option_list <- list(make_option(c("", "--hybrids"), action = "store", type = "character", help = "Hybrids file"),
            make_option(c("-g", "--gtf"), action = "store", type = "character", help = "Transcript gtf file"),
            make_option(c("-o", "--output"), action = "store", type = "character", help = "Output stem"))
opt_parser = OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Load data
genes.gr <- rtracklayer::import.gff2(opt$gtf)
hybrids.dt <- fread(opt$hybrids)

clusters.dt <- collapse_clusters(hybrids.dt)
fwrite(clusters.dt, paste0(opt$output, ".clusters.tsv.gz"), sep = "\t")

# Select intragenic only
# Adjust ones where there is some s_overlap (e.g. clusters)
clusters.dt <- clusters.dt[L_seqnames == R_seqnames]
clusters.dt[R_start < L_end, R_start := L_end + 1]
clusters.dt <- clusters.dt[grep("^rRNA", L_seqnames, invert = TRUE)] # Remove rRNA
clusters.dt <- clusters.dt[grep("Mt", L_seqnames, invert = TRUE)] # Remove MT for now
setnames(clusters.dt, "cluster", "name") # Adjust for BED

convert_coordinates(clusters.dt, genes.gr = genes.gr, filename = paste0(opt$output, ".clusters.bed.gz"), sam_tag = FALSE)

message("Completed!")