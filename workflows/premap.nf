#!/usr/bin/env nextflow

// Define DSL2
nextflow.enable.dsl=2

include { STAR; STAR_TX_INDEX; STAR_TX } from '../modules/premap.nf'
include { FILTER_PREMAP_READS; FILTER_PREMAP_READS as FILTER_PREMAP_READS_TX } from '../modules/premap.nf'

// Which alignment the linear-read filter reads: the genome premap if it runs, otherwise an
// unspliced STAR alignment to the transcript FASTA (linear_reference = 'auto'), or as forced.
def linearReference() {
    if (params.skip_linear_filter) return 'none'
    if (params.linear_reference == 'auto') return params.skip_premap ? 'transcript' : 'genome'
    return params.linear_reference
}

workflow PREMAP {

    // Called when the genome premap and/or the linear-read filter is enabled.
    //   premap + linear (genome)      STAR -> filter (spliced + linear)
    //   premap, no linear filter      STAR -> filter (spliced only; previous behaviour)
    //   premap + linear (transcript)  STAR -> filter (spliced) -> STAR_TX -> filter (linear)
    //   no premap + linear            STAR_TX -> filter (linear only)

    take:
        reads           // channel: [ sample_id, fastq ]
        star_genome     // channel: /path/to/star/genome/index/ (empty if skip_premap)
        transcript_fa   // channel: /path/to/transcript/fasta/

    main:
        linear_ref = linearReference()

        if (linear_ref == 'transcript') {
            ch_tx_index = params.linear_star_index ?
                Channel.fromPath(params.linear_star_index, checkIfExists: true) :
                STAR_TX_INDEX(transcript_fa.collect()).index
        }

        if (!params.skip_premap) {

            STAR(reads, star_genome.collect())
            ch_genome_bam = STAR.out.bam.map { id, bam, bai -> [ id, bam, bai, [] ] }

            if (linear_ref == 'transcript') {

                FILTER_PREMAP_READS(ch_genome_bam, '--no-linear-filter', '.stage1')
                STAR_TX(FILTER_PREMAP_READS.out.fastq, ch_tx_index.collect())
                ch_tx_bam = STAR_TX.out.bam.join(FILTER_PREMAP_READS.out.log) // [ id, bam, bai, stage-1 log ]
                FILTER_PREMAP_READS_TX(ch_tx_bam, '--no-spliced-filter', '')
                ch_fastq = FILTER_PREMAP_READS_TX.out.fastq
                ch_logs = FILTER_PREMAP_READS_TX.out.log

            } else {

                filter_args = linear_ref == 'none' ? '--no-linear-filter' : ''
                FILTER_PREMAP_READS(ch_genome_bam, filter_args, '')
                ch_fastq = FILTER_PREMAP_READS.out.fastq
                ch_logs = FILTER_PREMAP_READS.out.log

            }

        } else {

            if (linear_ref != 'transcript') {
                error "linear_reference = '${params.linear_reference}' needs the genome premap; use 'auto' or 'transcript' with --skip_premap"
            }
            STAR_TX(reads, ch_tx_index.collect())
            ch_tx_bam = STAR_TX.out.bam.map { id, bam, bai -> [ id, bam, bai, [] ] }
            FILTER_PREMAP_READS_TX(ch_tx_bam, '--no-spliced-filter', '')
            ch_fastq = FILTER_PREMAP_READS_TX.out.fastq
            ch_logs = FILTER_PREMAP_READS_TX.out.log

        }

    emit:
        fastq = ch_fastq
        logs = ch_logs

}
