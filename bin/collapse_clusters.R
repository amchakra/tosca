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

# Collapse clusters
clusters.dt <- collapse_clusters(hybrids.dt)
fwrite(clusters.dt, paste0(opt$output, ".clusters.tsv.gz"), sep = "\t")

# Convert to genomic coordinates
clusters.dt <- clusters.dt[grep("^rRNA", L_seqnames, invert = TRUE)] # Remove rRNA
clusters.dt <- clusters.dt[grep("Mt", L_seqnames, invert = TRUE)] # Remove MT for now
clusters.dt <- convert_coordinates(hybrids.dt = clusters.dt, genes.gr = genes.gr)
fwrite(clusters.dt, paste0(opt$output, ".clusters.gc.tsv.gz"), sep = "\t")

# Export genomic bed
intragenic.dt <- clusters.dt[L_seqnames == R_seqnames]
export_genomic_bed(hybrids.dt = intragenic.dt, sam_tag = TRUE, filename = paste0(opt$output, ".intragenic_clusters.bed.gz"))

message("Completed!")