#!/usr/bin/env nextflow

// Define DSL2
nextflow.enable.dsl=2

include { STAR; FILTER_SPLICED_READS } from '../modules/premap.nf'

workflow PREMAP {

    take:
    reads       // channel: fastq
    index       // channel: /path/to/star/index/

    main:
    STAR(reads, index.collect())
    FILTER_SPLICED_READS(STAR.out.bam)

    emit:
    fastq = FILTER_SPLICED_READS.out.fastq

}