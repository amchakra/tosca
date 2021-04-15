#!/usr/bin/env nextflow

// Define DSL2
nextflow.enable.dsl=2

include { SPLIT_FASTQ; FASTQ_TO_FASTA } from '../modules/splitfastq.nf'
include { BLAT; FILTER_BLAT } from '../modules/maphybrids.nf'
include { IDENTIFY_HYBRIDS; MERGE_HYBRIDS } from '../modules/identifyhybrids.nf'
include { DEDUPLICATE } from '../modules/deduplicate.nf'

workflow GET_HYBRIDS {

    take:
    reads       // channel: fastq
    fasta       // channel: /path/to/transcript/fasta/

    main:

    // Split
    SPLIT_FASTQ(reads)

    ch_split_fastq = SPLIT_FASTQ.out.fastq
        .flatten()
        .map { file -> tuple(file.simpleName, file) }

    // Convert to fasta
    FASTQ_TO_FASTA(ch_split_fastq)

    // Map hybrids
    BLAT(FASTQ_TO_FASTA.out.fasta, fasta.collect())
    FILTER_BLAT(BLAT.out.blast8)

    // Identify hybrids
    IDENTIFY_HYBRIDS(FILTER_BLAT.out.blast8.join(FASTQ_TO_FASTA.out.fasta))

    // Merge hybrids
    ch_merge_hybrids = IDENTIFY_HYBRIDS.out.hybrids
        .map { [ it[0].split('_')[0..-2].join('_'), it[1] ] }
        .groupTuple(by: 0)

    MERGE_HYBRIDS("hybrids", ch_merge_hybrids)
    DEDUPLICATE(MERGE_HYBRIDS.out.hybrids) // Remove PCR duplicates

    emit:
    hybrids = DEDUPLICATE.out.hybrids

}


