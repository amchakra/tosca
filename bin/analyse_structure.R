#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(primavera))
suppressPackageStartupMessages(library(rslurm))
suppressPackageStartupMessages(library(tictoc))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(optparse))

option_list <- list(make_option(c("", "--hybrids"), action = "store", type = "character", help = "Hybrids file"),
                make_option(c("", "--fasta"), action = "store", type = "character", help = "Transcript fasta"),
                make_option(c("", "--output"), action = "store", type = "character", help = "Output file"),
                make_option(c("", "--nodes"), action = "store", type = "integer", default = 100, help = "Number of nodes to allocate [default: %default]"),
                make_option(c("", "--shuffled_mfe"), action = "store_true", type = "logical", help = "Calculate shuffled binding energy (100 iterations)", default = FALSE),
                make_option(c("", "--clusters_only"), action = "store_true", type = "logical", help = "Analyse structure for hybrids in clusters only", default = FALSE))
opt_parser = OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Load genome
message("Loading genome...")
tic()
genome.fa <- Biostrings::readDNAStringSet(opt$fasta)
genome.dt <- data.table(gene_id = names(genome.fa),
                        sequence = as.character(genome.fa))
toc()

hybrids.dt <- fread(opt$hybrids)
setkey(hybrids.dt, name)
stopifnot(!any(duplicated(hybrids.dt$name))) # Check no duplicates

# Get sequences
message("Getting sequences...")
tic()
hybrids.dt <- get_sequence(hybrids.dt = hybrids.dt, genome.dt = genome.dt)
toc()

stopifnot(!any(is.na(c(hybrids.dt$L_sequence, hybrids.dt$L_sequence))))

sel.hybrids.dt <- hybrids.dt[!(L_seqnames == "rRNA_45S" & R_seqnames == "rRNA_45S")]
sel.hybrids.dt <- sel.hybrids.dt[!(L_seqnames == "rDNA" & R_seqnames == "rDNA")]
sel.hybrids.dt <- sel.hybrids.dt[!(L_seqnames == "rRNA_5S" & R_seqnames == "rRNA_5S")]

if(opt$clusters_only) sel.hybrids.dt <- sel.hybrids.dt[!is.na(cluster)][cluster != "."]

# Getting MFE and structure
message("Analysing structure...")
tic()

# Cluster jobs
sjob <- slurm_apply(analyse_structure, sel.hybrids.dt[, .(name, L_sequence, R_sequence)], 
                    jobname = sapply(strsplit(basename(opt$hybrids), "\\."), "[[", 1), 
                    nodes = opt$nodes, 
                    cpus_per_node = 1, 
                    slurm_options = list(time = "24:00:00"), 
                    submit = TRUE)

Sys.sleep(60)
status <- FALSE
while(status == FALSE) {

    squeue.out <- system(paste("squeue -n", sjob$jobname), intern = TRUE) # Get contents of squeue for this job
    if(length(squeue.out) == 1) status <- TRUE # i.e. only the header left
    Sys.sleep(60)

}
structure.list <- get_slurm_out(sjob)
cleanup_files(sjob) # Remove temporary files

structure.dt <- rbindlist(structure.list)
hybrids.dt <- merge(hybrids.dt, structure.dt, by = "name")

# Do shuffled depending on flag
if(opt$shuffled_mfe) {
    sjob <- slurm_apply(get_shuffled_mfe, sel.hybrids.dt[, .(name, L_sequence, R_sequence)], 
                        jobname = sapply(strsplit(basename(opt$hybrids), "\\."), "[[", 1), 
                        nodes = opt$nodes, 
                        cpus_per_node = 1, 
                        slurm_options = list(time = "24:00:00"), 
                        submit = TRUE)

    Sys.sleep(60) # To give it enough time to submit before the first check
    status <- FALSE
    while(status == FALSE) {

    squeue.out <- system(paste("squeue -n", sjob$jobname), intern = TRUE) # Get contents of squeue for this job
    if(length(squeue.out) == 1) status <- TRUE # i.e. only the header left
    Sys.sleep(60)

    }

    mfe.list <- get_slurm_out(sjob)
    cleanup_files(sjob) # Remove temporary files

    mfe.dt <- rbindlist(mfe.list)
    hybrids.dt <- merge(hybrids.dt, mfe.dt, by = "name")

}

fwrite(hybrids.dt, opt$output, sep = "\t")

message("Completed!")