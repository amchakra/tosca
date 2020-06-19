#!/usr/bin/env nextflow

// Specify DSL2
nextflow.preview.dsl = 2

process getbindingenergy {

    tag "${sample_id}"
    publishDir "${params.outdir}/hybrids", mode: 'copy', overwrite: true

    time '24h'

    input:
        tuple val(sample_id), path(hybrids), path(transcript_fa)

    output:
        tuple val(sample_id), path("${sample_id}.hybrids.mfe.tsv.gz")

    script:
    """
    #!/usr/bin/env Rscript

    suppressPackageStartupMessages(library(data.table))
    suppressPackageStartupMessages(library(primavera))
    suppressPackageStartupMessages(library(rslurm))
    suppressPackageStartupMessages(library(tictoc))
    suppressPackageStartupMessages(library(parallel))

    # Load genome
    message("Loading genome...")
    tic()
    genome.fa <- Biostrings::readDNAStringSet("$transcript_fa")
    genome.dt <- data.table(gene_id = names(genome.fa),
                            sequence = as.character(genome.fa))
    toc()

    hybrids.dt <- fread("$hybrids")

    # Get sequences
    message("Getting sequences...")
    tic()
    hybrids.dt <- GetSequence(hybrids.dt = hybrids.dt, genome.dt = genome.dt)
    toc()

    # Add ID
    hybrids.dt[, id := paste0("H", 1:.N)]

    # Getting MFE
    message("Calculating MFE...")
    tic()
    
    seq.df <- data.frame(hybrids.dt[, .(id, L_sequence, R_sequence)]) # Split out relevant part of hybrids.dt

    # Cluster jobs
    sjob <- slurm_apply(.slurm_GetMFE, seq.df, jobname = sapply(strsplit(basename("$hybrids"), "\\\\."), "[[", 1), nodes = 100, cpus_per_node = 1, slurm_options = list(time = "24:00:00"), submit = TRUE)
    Sys.sleep(60) # To give it enough time to submit before the first check

    status <- FALSE
    while(status == FALSE) {

        squeue.out <- system(paste("squeue -n", sjob\$jobname), intern = TRUE) # Get contents of squeue for this job
        if(length(squeue.out) == 1) status <- TRUE # i.e. only the header left
        Sys.sleep(60)

    }

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

    f_out <- paste0("$sample_id", "hybrids.mfe.tsv.gz")
    fwrite(hybrids.dt, f_out, sep = "\t")

    message("Completed!")
    """

}