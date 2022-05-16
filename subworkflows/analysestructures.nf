#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl=2

include { CHUNK_SEQUENCES; CALCULATE_STRUCTURES; CALCULATE_SHUFFLED_ENERGIES; MERGE_STRUCTURES; MERGE_SHUFFLED } from '../modules/analysestructures.nf'

workflow ANALYSE_STRUCTURES {

    take:
    
    type
    hybrids
    transcript_fa   // channel: transcript_fa

    main:

    CHUNK_SEQUENCES(type, hybrids)
    CALCULATE_STRUCTURES(type, CHUNK_SEQUENCES.out.rds
                                           .flatten()
                                           .map { file -> tuple(file.simpleName, file) })
    MERGE_STRUCTURES(type, hybrids.join(CALCULATE_STRUCTURES.out.tsv
                                                           .map { [ it[0].split('_')[0..-2].join('_'), it[1] ] }
                                                           .groupTuple(by: 0)))

    if(params.shuffled_mfe) {

        CALCULATE_SHUFFLED_ENERGIES(type, CHUNK_SEQUENCES.out.rds
                                           .flatten()
                                           .map { file -> tuple(file.simpleName, file) })
        MERGE_SHUFFLED(type, MERGE_STRUCTURES.out.tsv.join(CALCULATE_STRUCTURES.out.tsv
                                                           .map { [ it[0].split('_')[0..-2].join('_'), it[1] ] }
                                                           .groupTuple(by: 0)))

        emit:
        hybrids = MERGE_SHUFFLED.out.hybrids

    } else {
    
        emit:
        hybrids = MERGE_STRUCTURES.out.hybrids

    }

}