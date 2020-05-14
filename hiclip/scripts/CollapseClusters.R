# Script to get collapse clusters
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
                    make_option(c("-o", "--output"), action = "store", type = "character", help = "Output clusters file"))

opt_parser = OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

hybrids.dt <- fread(opt$input)

clusters.dt <- CollapseClusters(hybrids.dt)
fwrite(clusters.dt, opt$output, sep = "\t")

message(nrow(clusters.dt[L_end >= R_start]), " clusters removed")
clusters.dt <- clusters.dt[L_end < R_start]

# message("Converting coordinates (old)...")
# genes.gr <- rtracklayer::import.gff2(opt$annotation)
# tic()
# clusters.grl <- ConvertCoordinates(clusters.dt, genes.gr = genes.gr, cores = 8)
# ExportBED(clusters.grl, clusters.dt, filename = gsub("tsv$", "old.bed", opt$output), sam_tag = TRUE)
# toc()

# New method
message("Converting coordinates...")
genes.gr <- rtracklayer::import.gff2(opt$annotation)
tic()
ExportGenomicBED(seq.dt = clusters.dt, genes.gr = genes.gr, sam_tag = TRUE, filename = gsub("tsv$", "bed", opt$output))
toc()