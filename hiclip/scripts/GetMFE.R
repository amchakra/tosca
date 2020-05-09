# Script to get MFE
# A. M. Chakrabarti
# 6th May 2020


suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(primavera))
suppressPackageStartupMessages(library(rslurm))
suppressPackageStartupMessages(library(tictoc))
suppressPackageStartupMessages(library(parallel))

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

# Get sequences
message("Getting sequences...")
tic()
hybrids.dt <- GetSequence(hybrids.dt = hybrids.dt, genome.dt = genome.dt)
toc()

# Add ID
hybrids.dt[, id := paste0("H", 1:.N)]

# Getting MFE SJ motifs
message("Calculating MFE...")
# tic()
# cl <- makeForkCluster(8)
# hybrids.dt$mfe <- parSapply(cl = cl, 1:nrow(hybrids.dt), function(i) GetMFE(hybrids.dt$L_sequence[i], hybrids.dt$R_sequence[i]))
# toc()

tic()
# Split out relevant part of hybrids.dt
seq.df <- data.frame(hybrids.dt[, .(id, L_sequence, R_sequence)])

# Cluster jobs
sjob <- slurm_apply(.slurm_GetMFE, seq.df, jobname = sapply(strsplit(basename(opt$input), "\\."), "[[", 1), nodes = 100, cpus_per_node = 8, slurm_options = list(time = "24:00:00"), submit = TRUE)
Sys.sleep(10) # To give it enough time to submit before the first check
status <- get_job_status(sjob)$completed
while(status == FALSE) status <- get_job_status(sjob)$completed
mfe <- get_slurm_out(sjob, outtype = 'raw')

# Remove temporary files
cleanup_files(sjob) 

# Merge back
names(mfe) <- NULL # Not sure why it is names with L_sequence
mfe.dt <- as.data.table(unlist(mfe), keep.rownames = TRUE)
setnames(mfe.dt, c("id", "mfe"))
setkey(mfe.dt, id)
mfe.dt[mfe > 0, mfe := 0] # Get rid of positives
mfe.dt[is.na(mfe), mfe := 0] # Fix for positives < 10?

setkey(hybrids.dt, id)

hybrids.dt <- merge(hybrids.dt, mfe.dt, by = "id")
toc()



fwrite(hybrids.dt, opt$output, sep = "\t")

message("Completed!")
print(proc.time() - ptm)