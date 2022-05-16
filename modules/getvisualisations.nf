#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl=2

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
    suppressPackageStartupMessages(library(toscatools))

    hybrids.dt <- fread("$hybrids")
    intragenic.dt <- hybrids.dt[L_seqnames == R_seqnames]
    intragenic.dt <- intragenic.dt[grep("^rRNA|^rDNA", L_seqnames, invert = TRUE)] # Remove rRNA

    export_genomic_bed(hybrids.dt = intragenic.dt, sam_tag = TRUE, filename = "${sample_id}.${type}.intragenic.bed.gz")
    """

}

process EXPORT_GENOMIC_BAM {

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

process GET_CONTACT_MAPS {

    tag "${sample_id}"
    publishDir "${params.outdir}/maps", mode: 'copy', overwrite: true

    // time '6h'
    // memory '64 G'

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