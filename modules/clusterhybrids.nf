#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl=2

process CLUSTER_HYBRIDS {

    tag "${sample_id}"
    publishDir "${params.outdir}/${type}", mode: 'copy', overwrite: true

    // cpus 8
    // memory 32G
    // time '12h'

    input:
        val(type)
        tuple val(sample_id), path(hybrids)

    output:
        tuple val(sample_id), path("${sample_id}.${type}.clustered.tsv.gz"), emit: hybrids

    script:

    percent_overlap = params.percent_overlap
    sample_size = params.sample_size

    args = ""
    if(workflow.profile.contains("crick")) { args += " --slurm" }

    """
    cluster_hybrids.R \
        --sample ${sample_id} \
        --hybrids ${hybrids} \
        --type ${type} \
        --percent_overlap ${percent_overlap} \
        --sample_size ${sample_size} \
        ${args}
    """

}

process CLUSTER_HYBRIDS_VIRUS {

    tag "${sample_id}"
    publishDir "${params.outdir}/${type}", mode: 'copy', overwrite: true

    // cpus 8
    // memory 32G
    // time '12h'

    input:
        val(type)
        tuple val(sample_id), path(hybrids)

    output:
        tuple val(sample_id), path("${sample_id}.${type}.clustered.tsv.gz"), emit: hybrids

    script:

    percent_overlap = params.percent_overlap
    sample_size = params.sample_size

    args = ""
    if(workflow.profile.contains("crick")) { args += " --slurm" }

    """
    cluster_hybrids_virus.R \
        --sample ${sample_id} \
        --hybrids ${hybrids} \
        --type ${type} \
        --percent_overlap ${percent_overlap} \
        --sample_size ${sample_size} \
        ${args}
    """

}

// process CLUSTER_HYBRIDS {

//     tag "${sample_id}"
//     publishDir "${params.outdir}/${type}", mode: 'copy', overwrite: true

//     cpus 8
//     memory 32G
//     time '12h'

//     input:
//         val(type)
//         tuple val(sample_id), path(hybrids)

//     output:
//         tuple val(sample_id), path("${sample_id}.${type}.clustered.tsv.gz"), emit: hybrids

//     script:

//     percent_overlap = params.percent_overlap
//     sample_size = params.sample_size

//     """
//     #!/usr/bin/env Rscript

//     suppressPackageStartupMessages(library(data.table))
//     suppressPackageStartupMessages(library(primavera))
//     suppressPackageStartupMessages(library(parallel))
//     suppressPackageStartupMessages(library(rslurm))

//     setDTthreads(${task.cpus})
//     set.seed(42)

//     # Load hybrids
//     hybrids.dt <- fread("$hybrids")
//     setkey(hybrids.dt, name)

//     # Filter hybrids
//     message("Number of hybrids: ", nrow(hybrids.dt))
//     message("Removing rRNA-rRNA")
//     atlas.hybrids.dt <- hybrids.dt[, total_count := .N, by = .(L_seqnames, R_seqnames)]
//     atlas.hybrids.dt <- atlas.hybrids.dt[!(L_seqnames == "rRNA_45S" & R_seqnames == "rRNA_45S")]
//     atlas.hybrids.dt <- atlas.hybrids.dt[!(L_seqnames == "rDNA" & R_seqnames == "rDNA")]
//     atlas.hybrids.dt <- atlas.hybrids.dt[!(L_seqnames == "rRNA_5S" & R_seqnames == "rRNA_5S")]
    
//     message(nrow(atlas.hybrids.dt[total_count > 1e4]), " high incidence (>10,000) gene pairs not clustered")
//     atlas.hybrids.dt <- atlas.hybrids.dt[total_count < 1e4 & total_count > 1]

//     # Subsample as indicated
//     if($sample_size != -1) atlas.hybrids.dt <- atlas.hybrids.dt[sample(1:nrow(atlas.hybrids.dt, $sample_size))]
//     message("Number of hybrids to cluster: ", nrow(atlas.hybrids.dt))

//     # Keep ones not clustered to add back in later
//     unclustered.hybrids.dt <- hybrids.dt[!name %in% atlas.hybrids.dt\$name]
//     stopifnot(nrow(unclustered.hybrids.dt) + nrow(atlas.hybrids.dt) == nrow(hybrids.dt))

//     # Split into list to parallelise
//     atlas.hybrids.list <- split(atlas.hybrids.dt, by = c("L_seqnames", "R_seqnames"))
//     message("Gene pairs to cluster: ", length(atlas.hybrids.list))

//     atlas.clusters.list <- parallel::mclapply(atlas.hybrids.list, cluster_hybrids, percent_overlap = $percent_overlap, mc.cores = ${task.cpus})
        
//     # Submit to cluster
//     #sjob <- slurm_map(atlas.hybrids.list, f = cluster_hybrids, percent_overlap = $percent_overlap, jobname = sapply(strsplit(basename("$hybrids"), "\\\\."), "[[", 1), nodes = 100, cpus_per_node = 8, slurm_options = list(time = "12:00:00", mem = "64G"))
//     #status <- FALSE
//     #while(status == FALSE) {
//     #    squeue.out <- system(paste("squeue -n", sjob\$jobname), intern = TRUE) # Get contents of squeue for this job
//     #    if(length(squeue.out) == 1) status <- TRUE # i.e. only the header left
//     #    Sys.sleep(60)
//     #}

//     #atlas.clusters.list <- get_slurm_out(sjob)
//     #cleanup_files(sjob) 

//     atlas.clusters.dt <- rbindlist(atlas.clusters.list, use.names = TRUE, fill = TRUE)
//     stopifnot(nrow(atlas.clusters.dt) == nrow(atlas.hybrids.dt))

//     # Add back in unclustered
//     atlas.clusters.dt <- rbindlist(list(atlas.clusters.dt, unclustered.hybrids.dt), use.names = TRUE, fill = TRUE)

//     stopifnot(nrow(atlas.clusters.dt) == nrow(hybrids.dt))
//     fwrite(atlas.clusters.dt, "${sample_id}.${type}.clustered.tsv.gz", sep = "\t")
//     """

// }

process COLLAPSE_CLUSTERS {

    tag "${sample_id}"
    publishDir "${params.outdir}/${type}", mode: 'copy', overwrite: true

    // cpus 4
    // memory 16G
    // time '1h'

    input:
        val(type)
        tuple val(sample_id), path(hybrids)

    output:
        tuple val(sample_id), path("${sample_id}.clusters.tsv.gz"), emit: clusters

    script:

    """
    #!/usr/bin/env Rscript

    suppressPackageStartupMessages(library(data.table))
    suppressPackageStartupMessages(library(toscatools))

    hybrids.dt <- fread("$hybrids")

    # Collapse clusters
    clusters.dt <- collapse_clusters(hybrids.dt)
    fwrite(clusters.dt, "${sample_id}.clusters.tsv.gz", sep = "\t")

    message("Completed!")
    """
}

// process CLUSTER_HYBRIDS_ATLAS {

//     tag "${sample_id}"
//     publishDir "${params.outdir}/atlas", mode: 'copy', overwrite: true

//     cpus 8
//     memory 32G
//     time '24h'

//     input:
//         tuple val(sample_id), path(hybrids)

//     output:
//         tuple val(sample_id), path("${sample_id}.mfe.clusters.tsv.gz"), emit: hybrids

//     script:

//     percent_overlap = params.percent_overlap
//     sample_size = params.sample_size

//     """
//     #!/usr/bin/env Rscript

//     suppressPackageStartupMessages(library(data.table))
//     suppressPackageStartupMessages(library(primavera))
//     suppressPackageStartupMessages(library(parallel))
//     suppressPackageStartupMessages(library(rslurm))

//     setDTthreads(${task.cpus})

//     # Load hybrids
//     hybrids.dt <- fread("$hybrids")

//     # hybrids.list <- split(hybrids.dt, by = c("L_seqnames", "R_seqnames"))
//     # message(length(hybrids.list), " gene pairs to cluster")
//     # message("Distribution of hybrids per gene pair:")
//     # table(S4Vectors::elementNROWS(hybrids.list))

//     message("Number of hybrids: ", nrow(hybrids.dt))
//     message("Removing rRNA-rRNA and very high incidence genes (>100,000)...")
//     atlas.hybrids.dt <- hybrids.dt[L_seqnames != "rRNA_45S"][R_seqnames != "rRNA_45S"]
//     atlas.hybrids.dt <- hybrids.dt[L_seqnames != "rDNA"][R_seqnames != "rDNA"]
//     atlas.hybrids.dt <- atlas.hybrids.dt[L_seqnames != "rRNA_5S"][R_seqnames != "rRNA_5S"]
//     atlas.hybrids.dt <- atlas.hybrids.dt[, total_count := .N, by = .(L_seqnames, R_seqnames)][total_count < 1e4]
//     atlas.hybrids.dt[, total_count := NULL]
//     message("Number of hybrids remaining: ", nrow(atlas.hybrids.dt))
  
//     # Split into list for rslurm
//     atlas.hybrids.list <- split(atlas.hybrids.dt, by = c("L_seqnames", "R_seqnames"))

//     sjob <- slurm_map(atlas.hybrids.list, f = cluster_hybrids, percent_overlap = 0.75, jobname = sapply(strsplit(basename("$hybrids"), "\\\\."), "[[", 1), nodes = 100, cpus_per_node = 8, slurm_options = list(time = "12:00:00", mem = "64G"))
//     status <- FALSE
//     while(status == FALSE) {

//         squeue.out <- system(paste("squeue -n", sjob\$jobname), intern = TRUE) # Get contents of squeue for this job
//         if(length(squeue.out) == 1) status <- TRUE # i.e. only the header left
//         Sys.sleep(60)

//     }
    
//     atlas.clusters.list <- get_slurm_out(sjob)
//     cleanup_files(sjob) 

//     atlas.clusters.dt <- rbindlist(atlas.clusters.list, use.names = TRUE, fill = TRUE)
//     stopifnot(nrow(atlas.clusters.dt) == nrow(atlas.hybrids.dt))

//     fwrite(clusters.dt, "${sample_id}.clusters.tsv.gz", sep = "\t")
//     """

// }