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
sjob <- slurm_apply(.slurm_GetMFE, seq.df, jobname = sapply(strsplit(basename(opt$input), "\\."), "[[", 1), nodes = 100, cpus_per_node = 1, slurm_options = list(time = "24:00:00"), submit = TRUE)
Sys.sleep(60) # To give it enough time to submit before the first check

status <- FALSE
while(status == FALSE) {

	squeue.out <- system(paste("squeue -n", sjob$jobname), intern = TRUE) # Get contents of squeue for this job
	if(length(squeue.out) == 1) status <- TRUE # i.e. only the header left
	Sys.sleep(60)

}

# status <- get_job_status(sjob)$completed # This fails if not all the jobs have been submitted yet (i.e. in an array)
# while(status == FALSE) {
# 	status <- try(get_job_status(sjob)$completed)
# 	Sys.sleep(10)
# }

mfe <- get_slurm_out(sjob, outtype = 'raw')

# Remove temporary files
cleanup_files(sjob) 

# Merge back
names(mfe) <- NULL # Not sure why it is names with L_sequence - not anymore...?
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