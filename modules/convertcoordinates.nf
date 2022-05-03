#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl=2

process CONVERT_COORDINATES {

    tag "${sample_id}"
    publishDir "${params.outdir}/${type}", mode: 'copy', overwrite: true

    // cpus 4
    // memory 16G
    // time '1h'

    input:
        val(type)
        tuple val(sample_id), path(hybrids)
        path(transcript_gtf)

    output:
        tuple val(sample_id), path("${sample_id}.${type}.gc.tsv.gz"), emit: hybrids

    script:
    """
    #!/usr/bin/env Rscript

    suppressPackageStartupMessages(library(data.table))
    suppressPackageStartupMessages(library(primavera))

    # Load data
    genes.gr <- rtracklayer::import.gff2("$transcript_gtf")
    hybrids.dt <- fread("$hybrids")
    
    # hybrids.dt <- hybrids.dt[grep("^rRNA", L_seqnames, invert = TRUE)] # Remove rRNA
    # hybrids.dt <- hybrids.dt[grep("Mt", L_seqnames, invert = TRUE)] # Remove MT for now

    coord.dt <- convert_coordinates(hybrids.dt = hybrids.dt, genes.gr = genes.gr)
    fwrite(coord.dt, "${sample_id}.${type}.gc.tsv.gz", sep = "\t")
    """
}

process EXPORT_GENOMIC_BED {

    tag "${sample_id}"
    publishDir "${params.outdir}/igv", mode: 'copy', overwrite: true

    // cpus 4
    // memory 16G
    // time '1h'

    input:
        val(type)
        tuple val(sample_id), path(hybrids)

    output:
        tuple val(sample_id), path("${sample_id}.${type}.intragenic.bed.gz"), emit: bed

    script:
    """
    #!/usr/bin/env Rscript

    suppressPackageStartupMessages(library(data.table))
    suppressPackageStartupMessages(library(primavera))

    hybrids.dt <- fread("$hybrids")
    intragenic.dt <- hybrids.dt[L_seqnames == R_seqnames]
    intragenic.dt <- intragenic.dt[grep("^rRNA|^rDNA", L_seqnames, invert = TRUE)] # Remove rRNA

    export_genomic_bed(hybrids.dt = intragenic.dt, sam_tag = TRUE, filename = "${sample_id}.${type}.intragenic.bed.gz")
    """

}

process CONVERT_BED_TO_BAM {

    tag "${sample_id}"
    publishDir "${params.outdir}/igv", mode: 'copy', overwrite: true

    // cpus 4
    // memory 16G
    // time '1h'

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