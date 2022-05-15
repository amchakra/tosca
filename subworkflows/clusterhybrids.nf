#!/usr/bin/env nextflow

// Define DSL2
nextflow.enable.dsl=2

include { CHUNK_HYBRIDS; IDENTIFY_CLUSTERS; MERGE_CLUSTERS } from '../modules/clusterhybrids.nf'


workflow CLUSTER_HYBRIDS {

    take:
    
    type
    hybrids

    main:

    CHUNK_HYBRIDS(type, hybrids)
    IDENTIFY_CLUSTERS(type, CHUNK_HYBRIDS.out.rds
                                           .flatten()
                                           .map { file -> tuple(file.simpleName, file) })
    MERGE_CLUSTERS(type, hybrids.join(CHUNK_HYBRIDS.out.tsv)
                                .join(IDENTIFY_CLUSTERS.out.tsv
                                                           .map { [ it[0].split('_')[0..-2].join('_'), it[1] ] }
                                                           .groupTuple(by: 0)))
    
    emit:
    hybrids = MERGE_CLUSTERS.out.hybrids

}