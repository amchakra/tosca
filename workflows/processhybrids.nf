#!/usr/bin/env nextflow

// Define DSL2
nextflow.enable.dsl=2

include { ANNOTATE_HYBRIDS as ANNOTATE_HYBRIDS; ANNOTATE_HYBRIDS as ANNOTATE_CLUSTERS } from '../modules/annotatehybrids.nf'
include { CLUSTER_HYBRIDS_SLURM; COLLAPSE_CLUSTERS } from '../modules/clusterhybrids.nf'
include { CLUSTER_HYBRIDS } from '../subworkflows/clusterhybrids.nf'
include { CONVERT_COORDINATES as CONVERT_HYBRID_COORDINATES; CONVERT_COORDINATES as CONVERT_CLUSTER_COORDINATES } from '../modules/convertcoordinates.nf'
include { ANALYSE_STRUCTURES } from '../subworkflows/analysestructures.nf'

workflow PROCESS_HYBRIDS {

    take:
        hybrids         // channel: hybrids
        transcript_fa   // channel: transcript_fa
        transcript_gtf  // channel: transcript_gtf
        regions_gtf     // channel: regions_gtf

    main:
        if(params.cluster_old) {
            CLUSTER_HYBRIDS_SLURM("hybrids", hybrids)
        }

        CLUSTER_HYBRIDS("hybrids", hybrids)
        CONVERT_HYBRID_COORDINATES("hybrids", CLUSTER_HYBRIDS.out.hybrids, transcript_gtf.collect()) // Get genomic coordinates for hybrids
        ANNOTATE_HYBRIDS("hybrids", CONVERT_HYBRID_COORDINATES.out.hybrids, regions_gtf.collect()) // Annotate

        COLLAPSE_CLUSTERS("clusters", CLUSTER_HYBRIDS.out.hybrids) // Collapse clusters
        CONVERT_CLUSTER_COORDINATES("clusters", COLLAPSE_CLUSTERS.out.clusters, transcript_gtf.collect()) // Get genomic coordinates for clusters
        ANNOTATE_CLUSTERS("clusters", CONVERT_CLUSTER_COORDINATES.out.hybrids, regions_gtf.collect()) // Annotate      

        if(params.analyse_structures) {

            ANALYSE_STRUCTURES("hybrids", ANNOTATE_HYBRIDS.out.hybrids, transcript_fa.collect())
            
            output_hybrids = ANALYSE_STRUCTURES.out.hybrids
            output_clusters = ANNOTATE_CLUSTERS.out.hybrids
        
        } else {

            output_hybrids = ANNOTATE_HYBRIDS.out.hybrids
            output_clusters = ANNOTATE_CLUSTERS.out.hybrids
        }

    emit:
        hybrids = output_hybrids
        clusters = output_clusters

}