# Script to get MFE
# A. M. Chakrabarti
# 6th May 2020


suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(primavera))
suppressPackageStartupMessages(library(tictoc))
library(parallel)

option_list <- list(make_option(c("-i", "--input"), action = "store", type = "character", help = "Input hybrids file"),
                    make_option(c("-r", "--ref"), action = "store", type = "character", help = "Reference fasta"),
                    make_option(c("-o", "--output"), action = "store", type = "character", help = "Output hybrids file"))

opt_parser = OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

ptm <- proc.time()

# Load genome
message("Loading genome...")
tic()
genome.dt <- fread(opt$ref)
toc()

hybrids.dt <- fread(opt$input)

# Extract hybrids
message("Getting sequences...")
tic()
hybrids.dt <- GetSequence(hybrids.dt = hybrids.dt, genome.dt = genome.dt)
toc()

# Getting MFE SJ motifs
message("Calculating MFE...")
tic()
cl <- makeForkCluster(8)
hybrids.dt$mfe <- parSapply(cl = cl, 1:nrow(hybrids.dt), function(i) GetMFE(hybrids.dt$L_sequence[i], hybrids.dt$R_sequence[i]))
toc()

fwrite(hybrids.dt, opt$output, sep = "\t")

message("Completed!")
print(proc.time() - ptm)