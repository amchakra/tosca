#!/usr/bin/env nextflow

// Specify DSL2
nextflow.preview.dsl = 2

process clusterbindingenergy {

    tag "${sample_id}"
    publishDir "${params.outdir}/clusters", mode: 'copy', overwrite: true

    cpus 8
    time '24h'

    input:
        tuple val(sample_id), path(clusters), path(clustersbed), path(transcript_fa)

    output:
        tuple val(sample_id), path("${sample_id}.clusters.mfe.tsv.gz")

    script:
    """
    #!/usr/bin/env Rscript

    suppressPackageStartupMessages(library(data.table))
    suppressPackageStartupMessages(library(primavera))
    suppressPackageStartupMessages(library(rslurm))
    suppressPackageStartupMessages(library(tictoc))
    suppressPackageStartupMessages(library(parallel))

    setDTthreads(8)

    # Load genome
    message("Loading genome...")
    tic()
    genome.fa <- Biostrings::readDNAStringSet("$transcript_fa")
    genome.dt <- data.table(gene_id = names(genome.fa),
                            sequence = as.character(genome.fa))
    toc()

    clusters.dt <- fread("$clusters")

    # Get sequences
    message("Getting sequences...")
    tic()
    clusters.dt <- GetSequence(hybrids.dt = clusters.dt, genome.dt = genome.dt)
    toc()

    # Getting MFE
    message("Calculating MFE...")
    tic()
    
    seq.df <- data.frame(clusters.dt[, .(name, L_sequence, R_sequence)]) # Split out relevant part of hybrids.dt
    colnames(seq.df) <- c("id", "L_sequence", "R_sequence")

    # Cluster jobs
    sjob <- slurm_apply(.slurm_AnalyseMFE, seq.df, jobname = sapply(strsplit(basename("$clusters"), "\\\\."), "[[", 1), nodes = 100, cpus_per_node = 1, slurm_options = list(time = "24:00:00"), submit = TRUE)
    Sys.sleep(60) # To give it enough time to submit before the first check

    status <- FALSE
    while(status == FALSE) {

        squeue.out <- system(paste("squeue -n", sjob\$jobname), intern = TRUE) # Get contents of squeue for this job
        if(length(squeue.out) == 1) status <- TRUE # i.e. only the header left
        Sys.sleep(60)

    }

    mfe <- get_slurm_out(sjob, outtype = 'table')
    mfe.dt <- as.data.table(mfe, keep.rownames = TRUE)
    setnames(mfe.dt, "rn", "name")
    setkey(mfe.dt, name)
    mfe.dt[mfe > 0, mfe := 0] # Get rid of positives
    mfe.dt[is.na(mfe), mfe := 0] # Fix for positives < 10?

    # Remove temporary files
    cleanup_files(sjob)

    # Merge back
    setkey(clusters.dt, name)
    clusters.dt <- merge(clusters.dt, mfe.dt, by = "name")
    toc()

    f_out <- paste0("$sample_id", ".clusters.mfe.tsv.gz")
    fwrite(clusters.dt, f_out, sep = "\t")

    message("Completed!")
    """
}
