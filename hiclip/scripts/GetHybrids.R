# Script to get hybrids
# A. M. Chakrabarti
# 6th May 2020


suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(primavera))
suppressPackageStartupMessages(library(tictoc))

option_list <- list(make_option(c("-a", "--aligned"), action = "store", type = "character", help = "Input aligned file"),
					make_option(c("-c", "--chimeric"), action = "store", type = "character", help = "Input chimeric file"),
                    make_option(c("-r", "--ref"), action = "store", type = "character", help = "Reference fasta"),
                    make_option(c("-o", "--output"), action = "store", type = "character", help = "Output hybrids file"))

opt_parser = OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Load genome
message("Loading genome...")
tic()
genome.dt <- fread(opt$ref)
toc()

# Extract hybrids
message("Extracting hybrids...")
tic()
hybrids.dt <- ExtractHybrids(aligned.bam = opt$aligned, chimeric.junction = opt$chimeric)
hybrids.dt <- ReorientHybrids(hybrids.dt)
toc()

# Get SJ motifs
message("Getting SJ motifs...")
tic()
hybrids.dt <- GetSJMotifs(hybrids.dt = hybrids.dt, genome.dt = genome.dt)
toc()

fwrite(hybrids.dt, opt$output, sep = "\t")

message("Completed!")