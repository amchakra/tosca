#!/usr/bin/env nextflow

// Define DSL2
nextflow.enable.dsl=2

include { MERGE_HYBRIDS as MERGE_ATLAS_HYBRIDS } from '../modules/identifyhybrids.nf'
include { ANNOTATE_HYBRIDS } from '../modules/annotatehybrids.nf'
include { CLUSTER_HYBRIDS as CLUSTER_ATLAS_HYBRIDS; COLLAPSE_CLUSTERS as COLLAPSE_ATLAS_CLUSTERS } from '../modules/clusterhybrids.nf'
include { CONVERT_COORDINATES as CONVERT_CLUSTER_COORDINATES; EXPORT_GENOMIC_BED as EXPORT_ATLAS_BED; EXPORT_GENOMIC_BED as EXPORT_CLUSTER_BED; CONVERT_BED_TO_BAM as CONVERT_ATLAS_BED_TO_BAM } from '../modules/convertcoordinates.nf'

workflow GET_ATLAS {

    take:
    hybrids         // channel: hybrids (merged)
    transcript_gtf  // channel: transcript_gtf
    regions_gtf     // channel: regions_gtf
    genome_fai      // channel: genome_fai

    main:

    ch_all_hybrids = hybrids
        .map { [ 'all', it[1] ] }
        .groupTuple(by: 0)
        // .view()

    MERGE_ATLAS_HYBRIDS("atlas", ch_all_hybrids)
    CLUSTER_ATLAS_HYBRIDS("atlas", MERGE_HYBRIDS.out.hybrids) // Get clusters
    EXPORT_ATLAS_BED("atlas",  CLUSTER_ATLAS_HYBRIDS.out.hybrids)
    CONVERT_ATLAS_BED_TO_BAM(EXPORT_ATLAS_BED.out.bed, genome_fai.collect())

    COLLAPSE_ATLAS_CLUSTERS("atlas_clusters", CLUSTER_ATLAS_HYBRIDS.out.hybrids) // Collapse clusters
    CONVERT_CLUSTER_COORDINATES("atlas_clusters", COLLAPSE_ATLAS_CLUSTERS.out.clusters, transcript_gtf.collect())
    ANNOTATE_HYBRIDS("atlas_clusters", CONVERT_CLUSTER_COORDINATES.out.hybrids, regions_gtf.collect())
    EXPORT_CLUSTER_BED("atlas_clusters",  CLUSTER_HYBRIDS.out.hybrids)

}