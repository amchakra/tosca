#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl=2

process EXPORT_GENOMIC_BED {

    tag "${sample_id}"
    label 'process_low'

    publishDir "${params.outdir}/igv", mode: 'copy', overwrite: true

    input:
        val(type)
        tuple val(sample_id), path(hybrids)

    output:
        tuple val(sample_id), path("${sample_id}.${type}.intragenic.bed.gz"), emit: bed

    script:
    """
    #!/usr/bin/env Rscript

    suppressPackageStartupMessages(library(data.table))
    suppressPackageStartupMessages(library(toscatools))

    hybrids.dt <- fread("$hybrids")
    intragenic.dt <- hybrids.dt[L_seqnames == R_seqnames]
    intragenic.dt <- intragenic.dt[grep("^rRNA|^rDNA", L_seqnames, invert = TRUE)] # Remove rRNA

    export_genomic_bed(hybrids.dt = intragenic.dt, sam_tag = TRUE, filename = "${sample_id}.${type}.intragenic.bed.gz")
    """

}

process EXPORT_GENOMIC_BAM {

    tag "${sample_id}"
    label 'process_low'

    publishDir "${params.outdir}/igv", mode: 'copy', overwrite: true

    input:
        tuple val(sample_id), path(bed)
        path(genome_fai)

    output:
        tuple val(sample_id), path("${sample_id}.*.bam"), path("${sample_id}.*.bam.bai"), emit: bam

    script:
    """
    convert_hybrid_bed_to_bam.py $bed $genome_fai
    """

}

process GET_CONTACT_MAPS {

    tag "${sample_id}"
    label 'process_medium'

    publishDir "${params.outdir}/maps", mode: 'copy', overwrite: true

    input:
        tuple val(sample_id), path(hybrids)
        path(fai)
        path(genes)

    output:
        tuple val(sample_id), path("${sample_id}.*.mat.rds"), emit: map
        tuple val(sample_id), path("${sample_id}.*_binned_map.tsv.gz"), emit: binned_map

    script:

    bin_size = params.bin_size

    """
    get_contact_map.R --hybrids $hybrids --genes $genes --fai $fai --bin_size $bin_size --output ${sample_id}
    """

}

process GET_ARCS {

    tag "${sample_id}"
    label 'process_low'
    
    publishDir "${params.outdir}/igv", mode: 'copy', overwrite: true

    input:
        tuple val(sample_id), path(clusters)
        path(genes)

    output:
        tuple val(sample_id), path("${sample_id}.*.bp"), emit: arcs

    script:

    breaks = params.breaks

    """
    get_arcs.R --clusters $clusters --genes $genes --breaks $breaks --output ${sample_id}
    """

}

process EXPORT_BEDPE {

    tag "${sample_id}"
    label 'process_low'

    publishDir "${params.outdir}/igv", mode: 'copy', overwrite: true

    input:
        val(type)
        tuple val(sample_id), path(hybrids)

    output:
        tuple val(sample_id), path("${sample_id}.${type}.bedpe.gz"), emit: bedpe

    script:
    """
    #!/usr/bin/env Rscript

    suppressPackageStartupMessages(library(data.table))
    suppressPackageStartupMessages(library(toscatools))

    hybrids.dt <- fread("$hybrids")
    hybrids.dt <- toscatools::reorient_hybrids(hybrids.dt)

    if("$type" == "hybrids") {
        bedpe.colnames <- c("L_seqnames", "L_start", "L_end", "R_seqnames", "R_start", "R_end", "name", "total_count", "L_strand", "R_strand")
        bedpe.dt <- hybrids.dt[, ..bedpe.colnames]
    } else if($type == "clusters") {
        bedpe.colnames <- c("L_seqnames", "L_start", "L_end", "R_seqnames", "R_start", "R_end", "name", "cluster_hybrid_count", "L_strand", "R_strand")
        bedpe.dt <- hybrids.dt[, ..bedpe.colnames]
    }

    bedpe.dt[, `:=` (L_start = L_start - 1, 
                     R_start = R_start - 1)]

    fwrite(bedpe.dt, "${sample_id}.${type}.bedpe.gz", sep = "\t", col.names = FALSE, quote = FALSE)
    """

}