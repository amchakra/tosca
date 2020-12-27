#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(primavera))
suppressPackageStartupMessages(library(parallel))

option_list <- list(make_option(c("", "--hybrids"), action = "store", type = "character", help = "Hybrids file"),
            make_option(c("-p", "--percent_overlap"), action = "store", type = "character", help = "Percentage overlap"),
            make_option(c("-s", "--sample_size"), action = "store", type = "integer", default = -1, help = "Sample size [default: %default]"),
            make_option(c("-o", "--output"), action = "store", type = "character", help = "Output file"),
            make_option(c("-t", "--threads"), action = "store", type = "integer", default = 8, help = "Number of threads [default: %default]"))
opt_parser = OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

setDTthreads(opt$threads)

# Load hybrids
hybrids.dt <- fread(opt$hybrids)
clusters.dt <- cluster_hybrids(hybrids.dt, percent_overlap = opt$percent_overlap, sample_size = opt$sample_size, cores = opt$threads)
fwrite(clusters.dt, opt$output, sep = "\t")

message("Completed!")
