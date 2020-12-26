#!/usr/bin/env nextflow

// Define DSL2
nextflow.enable.dsl=2

include { DEDUPLICATE } from '../modules/deduplicate.nf'
include { GET_BINDING_ENERGY } from '../modules/getbindingenergy.nf'
include { CLUSTER_HYBRIDS } from '../modules/clusterhybrids.nf'
include { CONVERT_COORDINATES } from '../modules/convertcoordinates.nf'

workflow PROCESS_HYBRIDS {

    take:
    hybrids         // channel: hybrids
    transcript_fa   // channel: transcript_fa

    main:
    DEDUPLICATE(hybrids) // Remove PCR duplicates
    GET_BINDING_ENERGY(DEDUPLICATE.out.hybrids, transcript_fa.collect()) // Get binding energies
    CLUSTER_HYBRIDS(GET_BINDING_ENERGY.out.hybrids) // Get clusters

    emit:
    mfe = GET_BINDING_ENERGY.out.hybrids
    hybrids = CLUSTER_HYBRIDS.out.hybrids

}