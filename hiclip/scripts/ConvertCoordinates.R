# Script to convert to genomic coordinates
# A. M. Chakrabarti
# 8th May 2020

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(primavera))
# suppressPackageStartupMessages(library(rslurm))
suppressPackageStartupMessages(library(tictoc))
# suppressPackageStartupMessages(library(parallel))

option_list <- list(make_option(c("-i", "--input"), action = "store", type = "character", help = "Input hybrids file"),
                    make_option(c("-a", "--annotation"), action = "store", type = "character", help = "Annotation GTF"),
                    make_option(c("-o", "--output"), action = "store", type = "character", help = "Output BED file"))

opt_parser = OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

genes.gr <- rtracklayer::import.gff2(opt$annotation)

hybrids.dt <- fread(opt$input)
seq.dt <- hybrids.dt[overlapping_hybrid %in% c(NA, FALSE)]
seq.dt <- seq.dt[grep("Mt", L_seqnames, invert = TRUE)] # Remove MT for now

"Converting coordinates..."
tic()
g.grl <- ConvertCoordinates(seq.dt = seq.dt, genes.gr = genes.gr, cores = 8)
toc()

saveRDS(g.grl, gsub("bed$", "ggrl.rds", opt$output))

ExportBED(g.grl = g.grl, hybrids.dt = seq.dt, filename = opt$output, sam_tag = TRUE)