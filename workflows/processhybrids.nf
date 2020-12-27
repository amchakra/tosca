#!/usr/bin/env nextflow

// Define DSL2
nextflow.enable.dsl=2

include { DEDUPLICATE } from '../modules/deduplicate.nf'
include { GET_BINDING_ENERGY } from '../modules/getbindingenergy.nf'
include { CLUSTER_HYBRIDS; COLLAPSE_CLUSTERS } from '../modules/clusterhybrids.nf'
include { CONVERT_COORDINATES as CONVERT_HYBRID_COORDINATES; CONVERT_COORDINATES as CONVERT_CLUSTER_COORDINATES } from '../modules/convertcoordinates.nf'

workflow PROCESS_HYBRIDS {

    take:
    hybrids         // channel: hybrids
    transcript_fa   // channel: transcript_fa
    transcript_gtf   // channel: transcript_gtf

    main:
    
    DEDUPLICATE(hybrids) // Remove PCR duplicates
    GET_BINDING_ENERGY(DEDUPLICATE.out.hybrids, transcript_fa.collect()) // Get binding energies
    CLUSTER_HYBRIDS(GET_BINDING_ENERGY.out.hybrids) // Get clusters
    COLLAPSE_CLUSTERS(CLUSTER_HYBRIDS.out.hybrids) // Collapse clusters

    CONVERT_HYBRID_COORDINATES("hybrids", CLUSTER_HYBRIDS.out.hybrids, transcript_gtf.collect()) // Get genomic coordinates for hybrids
    CONVERT_CLUSTER_COORDINATES("clusters", COLLAPSE_CLUSTERS.out.clusters, transcript_gtf.collect()) // Get genomic coordinates for clusters


    emit:
    dedup = DEDUPLICATE.out.hybrids
    mfe = GET_BINDING_ENERGY.out.hybrids
    clustered = CLUSTER_HYBRIDS.out.hybrids

    hybrids = CONVERT_HYBRID_COORDINATES.out.hybrids
    clusters = CONVERT_CLUSTER_COORDINATES.out.hybrids

}