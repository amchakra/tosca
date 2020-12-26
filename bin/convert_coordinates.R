#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(primavera))
suppressPackageStartupMessages(library(tictoc))

option_list <- list(make_option(c("", "--hybrids"), action = "store", type = "character", help = "Hybrids file"),
        make_option(c("-g", "--gtf"), action = "store", type = "character", help = "Transcript gtf file"),
        make_option(c("-o", "--output"), action = "store", type = "character", help = "Output bed file"))
opt_parser = OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Load data
genes.gr <- rtracklayer::import.gff2(opt$gtf)
hybrids.dt <- fread(opt$hybrids)

# Select only intragenic
hybrids.dt <- hybrids.dt[type == "intragenic"]

# Adjust ones where there is some s_overlap (e.g. clusters)
# message(nrow(hybrids.dt[R_start < L_end]))
hybrids.dt[R_start < L_end, R_start := L_end + 1]    
hybrids.dt <- hybrids.dt[grep("^rRNA", L_seqnames, invert = TRUE)] # Remove rRNA
hybrids.dt <- hybrids.dt[grep("Mt", L_seqnames, invert = TRUE)] # Remove MT for now

convert_coordinates(hybrids.dt = hybrids.dt, genes.gr = genes.gr, sam_tag = TRUE, filename = opt$output)

message("Completed!")