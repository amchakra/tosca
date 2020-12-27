#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(primavera))

option_list <- list(make_option(c("", "--hybrids"), action = "store", type = "character", help = "Hybrids file"),
        make_option(c("-g", "--gtf"), action = "store", type = "character", help = "Transcript gtf file"),
        make_option(c("-o", "--output"), action = "store", type = "character", help = "Output stem"))
opt_parser = OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Load data
genes.gr <- rtracklayer::import.gff2(opt$gtf)
hybrids.dt <- fread(opt$hybrids)
hybrids.dt <- hybrids.dt[grep("^rRNA", L_seqnames, invert = TRUE)] # Remove rRNA
hybrids.dt <- hybrids.dt[grep("Mt", L_seqnames, invert = TRUE)] # Remove MT for now

coord.dt <- convert_coordinates(hybrids.dt = hybrids.dt, genes.gr = genes.gr)
fwrite(coord.dt, paste0(opt$output, ".gc.tsv.gz"), sep = "\t")

intragenic.dt <- coord.dt[type == "intragenic"]
export_genomic_bed(hybrids.dt = intragenic.dt, sam_tag = TRUE, filename = paste0(opt$output, ".intragenic.bed.gz"))

message("Completed!")