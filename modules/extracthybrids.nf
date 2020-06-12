#!/usr/bin/env nextflow

// Specify DSL2
nextflow.preview.dsl = 2

process extracthybrids {
    
    tag "${sample_id}"
    publishDir "${params.outdir}/hybrids", mode: 'copy', overwrite: true

    cpus 4
    memory '32G'
    time '12h'

    input:
        tuple val(sample_id), path(reads), path(bai), path(transcript_fa)

    output:
        tuple val(sample_id), path("${sample_id}.hybrids.tsv.gz")

    script:
    """
    #!/usr/bin/env Rscript

    suppressPackageStartupMessages(library(data.table))
    suppressPackageStartupMessages(library(primavera))
    suppressPackageStartupMessages(library(tictoc))

    # Extract hybrids
    message("Extracting hybrids...")
    tic()
    hybrids.dt <- ExtractHybridsWithinBAM(aligned.bam = "$reads")
    hybrids.dt <- ReorientHybrids(hybrids.dt)
    toc()

    # Get SJ motifs
    message("Loading genome...")
    genome.fa <- Biostrings::readDNAStringSet("$transcript_fa")
    genome.dt <- data.table(gene_id = names(genome.fa),
                            sequence = as.character(genome.fa))

    message("Getting SJ motifs...")
    tic()
    hybrids.dt <- GetSJMotifs(hybrids.dt = hybrids.dt, genome.dt = genome.dt)
    toc()

    f_out <- gsub("dedup.bam", "hybrids.tsv.gz", "$reads")
    fwrite(hybrids.dt, f_out, sep = "\t")

    message("Completed!")
    """

}