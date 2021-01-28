#!/usr/bin/env nextflow

// Define DSL2
nextflow.enable.dsl=2

include { MERGE_HYBRIDS } from '../modules/identifyhybrids.nf'
include { ANNOTATE_HYBRIDS } from '../modules/annotatehybrids.nf'
include { CLUSTER_HYBRIDS_ATLAS; COLLAPSE_CLUSTERS } from '../modules/clusterhybrids.nf'
include { CONVERT_COORDINATES; EXPORT_GENOMIC_BED } from '../modules/convertcoordinates.nf'

workflow GET_ATLAS {

    take:
    hybrids         // channel: hybrids (merged)
    transcript_gtf  // channel: transcript_gtf
    regions_gtf     // channel: regions_gtf

    main:

    ch_all_hybrids = hybrids
        .map { [ 'all', it[1] ] }
        .groupTuple(by: 0)
        // .view()

    MERGE_HYBRIDS("atlas", ch_all_hybrids)
    CLUSTER_HYBRIDS_ATLAS(MERGE_HYBRIDS.out.hybrids) // Get clusters
    COLLAPSE_CLUSTERS("atlas", CLUSTER_HYBRIDS_ATLAS.out.hybrids) // Collapse clusters
    CONVERT_COORDINATES("atlas", COLLAPSE_CLUSTERS.out.clusters, transcript_gtf.collect()) // Get genomic coordinates for hybrids
    ANNOTATE_HYBRIDS("atlas", CONVERT_COORDINATES.out.hybrids, regions_gtf.collect()) // Annotate    
    EXPORT_GENOMIC_BED("atlas",  CONVERT_COORDINATES.out.hybrids)

}