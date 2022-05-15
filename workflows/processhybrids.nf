#!/usr/bin/env nextflow

// Define DSL2
nextflow.enable.dsl=2

// include { GET_BINDING_ENERGY } from '../modules/getbindingenergy.nf'
include { ANNOTATE_HYBRIDS as ANNOTATE_HYBRIDS; ANNOTATE_HYBRIDS as ANNOTATE_CLUSTERS } from '../modules/annotatehybrids.nf'
include { CLUSTER_HYBRIDS_SLURM; COLLAPSE_CLUSTERS } from '../modules/clusterhybrids.nf'
include { CLUSTER_HYBRIDS } from '../subworkflows/clusterhybrids.nf'
include { CONVERT_COORDINATES as CONVERT_HYBRID_COORDINATES; CONVERT_COORDINATES as CONVERT_CLUSTER_COORDINATES } from '../modules/convertcoordinates.nf'

workflow PROCESS_HYBRIDS {

    take:
    hybrids         // channel: hybrids
    transcript_fa   // channel: transcript_fa
    transcript_gtf  // channel: transcript_gtf
    regions_gtf     // channel: regions_gtf

    main:
    
    // GET_BINDING_ENERGY(hybrids, transcript_fa.collect()) // Get binding energies

    // CLUSTER_HYBRIDS(GET_BINDING_ENERGY.out.hybrids) // Get clusters
    if(params.cluster_old) {
        CLUSTER_HYBRIDS_SLURM("hybrids", hybrids)
    }

    CLUSTER_HYBRIDS("hybrids", hybrids)
    CONVERT_HYBRID_COORDINATES("hybrids", CLUSTER_HYBRIDS.out.hybrids, transcript_gtf.collect()) // Get genomic coordinates for hybrids
    ANNOTATE_HYBRIDS("hybrids", CONVERT_HYBRID_COORDINATES.out.hybrids, regions_gtf.collect()) // Annotate

    COLLAPSE_CLUSTERS("clusters", CLUSTER_HYBRIDS.out.hybrids) // Collapse clusters
    CONVERT_CLUSTER_COORDINATES("clusters", COLLAPSE_CLUSTERS.out.clusters, transcript_gtf.collect()) // Get genomic coordinates for clusters
    ANNOTATE_CLUSTERS("clusters", CONVERT_CLUSTER_COORDINATES.out.hybrids, regions_gtf.collect()) // Annotate      

    emit:
    // mfe = GET_BINDING_ENERGY.out.hybrids
    // clustered = CLUSTER_HYBRIDS.out.hybrids
    hybrids = ANNOTATE_HYBRIDS.out.hybrids
    clusters = ANNOTATE_CLUSTERS.out.hybrids

}