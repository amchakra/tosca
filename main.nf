#!/usr/bin/env nextflow

/*
========================================================================================
                                    amchakra/tosca
========================================================================================
hiCLIP analysis pipeline.
 #### Homepage / Documentation
 https://github.com/amchakra/tosca
----------------------------------------------------------------------------------------
*/

// Define DSL2
nextflow.enable.dsl=2

// Need to set these before module is loaded else not propagated
params.intronmin = 10
params.evalue = 0.001
params.maxhits = 100

params.adapter = 'AGATCGGAAGAGC'
params.min_quality = 10
params.min_readlength = 16

params.split_size = 100000

// Processes
include { hiclipheader } from './modules/utils.nf'
include { METADATA } from './modules/metadata.nf'
include { CUTADAPT } from './modules/cutadapt.nf'
include { PREMAP } from './modules/premap.nf'
include { FILTER_SPLICED_READS } from './modules/premap.nf'
include { SPLIT_FASTQ } from './modules/splitfastq.nf'
include { FASTQ_TO_FASTA } from './modules/convert_fastq_to_fasta.nf'
include { BLAT } from './modules/mapblat.nf'
include { FILTER_BLAT } from './modules/mapblat.nf'

// include { mapchimeras } from './modules/mapchimeras.nf'
// include { deduplicate } from './modules/deduplicate.nf'
// include { deduplicate_unique } from './modules/deduplicate.nf'
// include { extracthybrids } from './modules/extracthybrids.nf'
// include { getbindingenergy } from './modules/getbindingenergy.nf'
// include { clusterhybrids } from './modules/clusterhybrids.nf'
// include { collapseclusters } from './modules/collapseclusters.nf'
// include { clusterbindingenergy } from './modules/clusterbindingenergy.nf'
// include { convertcoordinates } from './modules/convertcoordinates.nf'
// include { hybridbedtohybridbam } from './modules/hybridbedtohybridbam.nf'

// include { splitfastq } from './modules/splitfastq.nf'
// include { convert_fastq_to_fasta } from './modules/convert_fastq_to_fasta.nf'
// include { mapblat } from './modules/mapblat.nf'
// include { filterblat } from './modules/filterblat.nf'
// include { identifyhybrids } from './modules/identifyhybrids.nf'
// include { mergehybrids } from './modules/mergehybrids.nf'
// include { deduplicate_blat } from './modules/deduplicate.nf'
// include { getnonhybrids } from './modules/getnonhybrids.nf'

// Main workflow

// General variables
params.quickdedup = true
params.premap = true


// Input variables
// params.input='metadata.csv'
// params.org='mouse'

// Genome variables
params.genome_fa = params.genomes[ params.org ].genome_fa
params.genome_fai = params.genomes[ params.org ].genome_fai
params.genome_gtf = params.genomes[ params.org ].genome_gtf 
params.transcript_fa = params.genomes[ params.org ].transcript_fa
params.transcript_gtf = params.genomes[ params.org ].transcript_gtf
params.star_genome = params.genomes[ params.org ].star_genome
params.star_transcript = params.genomes[ params.org ].star_transcript

// Create channels for static files
ch_star_genome = Channel.fromPath(params.star_genome, checkIfExists: true)
ch_star_transcript = Channel.fromPath(params.star_transcript, checkIfExists: true)
ch_transcript_fa = Channel.fromPath(params.transcript_fa, checkIfExists: true)
ch_genome_fai = Channel.fromPath(params.genome_fai, checkIfExists: true)
ch_transcript_gtf = Channel.fromPath(params.transcript_gtf, checkIfExists: true)

// Show header
log.info hiclipheader()
def summary = [:]
summary['Output directory'] = params.outdir
summary['Trace directory'] = params.tracedir
// summary['Genome fasta'] = params.genome_fa
// summary['Genome fasta index'] = params.genome_fai
summary['Genome annotation'] = params.genome_gtf
summary['Transcriptome fasta'] = params.transcript_fa
summary['Transcriptome annotation'] = params.transcript_gtf
summary['STAR genome'] = params.star_genome
// summary['STAR transcriptome'] = params.star_transcript
// summary['Deduplicate quickly'] = params.quickdedup
// summary['Minimum intron length'] = params.intronmin

log.info summary.collect { k,v -> "${k.padRight(25)}: $v" }.join("\n")
log.info "-\033[2m---------------------------------------------------------------\033[0m-"

def settings = [:]
settings['Adapter sequence'] = params.adapter
settings['Minimum read quality'] = params.min_quality
settings['Minimum read length'] = params.min_readlength
settings['FASTQ split size'] = params.split_size
settings['Minimum e-value'] = params.evalue
settings['Maximum hits/read'] = params.maxhits

log.info settings.collect { k,v -> "${k.padRight(25)}: $v" }.join("\n")
log.info "-\033[2m---------------------------------------------------------------\033[0m-"

// Pipeline
workflow {

    // Get fastq paths 
    METADATA(params.input)

    // Trim adapters
    CUTADAPT(METADATA.out)

    // Filter spliced reads
    PREMAP(CUTADAPT.out.fastq, ch_star_genome.collect())
    FILTER_SPLICED_READS(PREMAP.out.bam)

    // Split
    SPLIT_FASTQ(FILTER_SPLICED_READS.out.fastq)

    ch_split_fastq = SPLIT_FASTQ.out.fastq
        .flatten()
        .map { file -> tuple(file.simpleName, file) }

    // Convert to fasta
    FASTQ_TO_FASTA(ch_split_fastq)

    // Map hybrids
    BLAT(FASTQ_TO_FASTA.out.fasta, ch_transcript_fa.collect())
    FILTER_BLAT(BLAT.out.blast8)

    // // Merge back test
    // // ch_comb = convert_fastq_to_fasta.out
    // //     .map { [ it[0].split('_')[0..-2].join('_'), it[1] ] }
    // //     .groupTuple(by: 0)
    // //     .view()

    // // Map chimeras
    // mapblat(convert_fastq_to_fasta.out.combine(ch_transcript_fa))
    // filterblat(mapblat.out)

    // // Identify hybrids
    // identifyhybrids(filterblat.out.join(convert_fastq_to_fasta.out))

    // // Merge hybrids
    // ch_comb = identifyhybrids.out
    //     .map { [ it[0].split('_')[0..-2].join('_'), it[1] ] }
    //     .groupTuple(by: 0)
    //     // .view()

    // mergehybrids(ch_comb)

    // // Get non-hybrid reads for later
    // getnonhybrids(mergehybrids.out.join(filtersplicedreads.out))
    // Map chimerias
    // mapchimeras(filtersplicedreads.out.combine(ch_star_transcript))

    // Remove PCR duplicates
    // if ( params.quickdedup ) {
    //     deduplicate_unique(mapchimeras.out)
    // } else {
    //     deduplicate(mapchimeras.out)
    // }

    // deduplicate_blat(mergehybrids.out)

    // // // Extract hybrids
    // // if ( params.quickdedup ) {
    // //     extracthybrids(deduplicate_unique.out.combine(ch_transcript_fa))
    // // } else {
    // //     extracthybrids(deduplicate.out.combine(ch_transcript_fa))
    // // }
    // // // Get binding energies
    // // getbindingenergy(extracthybrids.out.combine(ch_transcript_fa))
    // getbindingenergy(deduplicate_blat.out.combine(ch_transcript_fa))

    // // // Get clusters
    // clusterhybrids(getbindingenergy.out)

    // // // Convert coordinates
    // // // Write hybrid BAM
    // convertcoordinates(clusterhybrids.out.combine(ch_transcript_gtf))
    // hybridbedtohybridbam(convertcoordinates.out.combine(ch_genome_fai))

    // // // Collapse clusters
    // collapseclusters(clusterhybrids.out.combine(ch_transcript_gtf))
    // clusterbindingenergy(collapseclusters.out.combine(ch_transcript_fa))

}

workflow.onComplete {

    log.info "\nPipeline complete!\n"

}