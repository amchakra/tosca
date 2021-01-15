#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl=2

include { EXPORT_GENOMIC_BED as EXPORT_HYBRID_GENOMIC_BED; EXPORT_GENOMIC_BED as EXPORT_CLUSTERS_GENOMIC_BED } from '../modules/convertcoordinates.nf'
include { CONVERT_BED_TO_BAM } from '../modules/convertcoordinates.nf'

workflow EXPORT_INTRAGENIC {

    take:
    hybrids         // channel: hybrids
    clusters        // channel: clusters
    genome_fai      // channel: transcript_fa

    main:
    
    EXPORT_HYBRID_GENOMIC_BED("hybrids", hybrids)
    CONVERT_BED_TO_BAM(EXPORT_HYBRID_GENOMIC_BED.out.bed, genome_fai.collect())
    EXPORT_CLUSTERS_GENOMIC_BED("clusters", clusters)

    emit:
    hybrid_bed = EXPORT_HYBRID_GENOMIC_BED.out.bed
    hybrid_bam = CONVERT_BED_TO_BAM.out.bam
    cluster_bed = EXPORT_CLUSTERS_GENOMIC_BED.out.bed

}