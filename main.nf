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
params.adapter = 'AGATCGGAAGAGC'
params.min_quality = 10
params.min_readlength = 16

params.split_size = 100000

params.evalue = 0.001
params.maxhits = 100

params.shuffled_mfe = false

// Processes
include { hiclipheader } from './modules/utils.nf'
include { METADATA } from './modules/metadata.nf'
include { CUTADAPT } from './modules/cutadapt.nf'
include { PREMAP } from './workflows/premap.nf'
include { GET_HYBRIDS } from './workflows/gethybrids.nf'

include { GET_NON_HYBRIDS } from './modules/getnonhybrids.nf'
include { deduplicate_blat } from './modules/deduplicate.nf'
include { GET_BINDING_ENERGY } from './modules/getbindingenergy.nf'

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
settings['Premapping'] = params.premap
settings['FASTQ split size'] = params.split_size
settings['Minimum e-value'] = params.evalue
settings['Maximum hits/read'] = params.maxhits
settings['Shuffled binding energy'] = params.shuffled_mfe

log.info settings.collect { k,v -> "${k.padRight(25)}: $v" }.join("\n")
log.info "-\033[2m---------------------------------------------------------------\033[0m-"
log.info "-----------------------------------------------------------------"

// Pipeline
workflow {


    /* 
    PREPARE INPUTS
    */
    METADATA(params.input) // Get fastq paths 
    CUTADAPT(METADATA.out) // Trim adapters

    /* 
    IDENTIFY HYBRIDS
    */
    if(params.premap) {

        PREMAP(CUTADAPT.out.fastq, ch_star_genome) // Filter spliced reads
        GET_HYBRIDS(PREMAP.out.fastq, ch_transcript_fa) // Identify hybrids

    } else {

        GET_HYBRIDS(CUTADAPT.out.fastq, ch_transcript_fa) // Identify hybrids

    }

    /* 
    IDENTIFY NON-HYBRIDS
    */
    GET_NON_HYBRIDS(GET_HYBRIDS.out.hybrids.join(METADATA.out))


    /* 
    PROCESS HYBRIDS
    */
    // Remove PCR duplicates
    if ( params.quickdedup ) {
        deduplicate_blat(GET_HYBRIDS.out.hybrids)
    } else {
        deduplicate_blat(GET_HYBRIDS.out.hybrids)
    }

    // Get binding energies
    GET_BINDING_ENERGY(deduplicate_blat.out.hybrids, ch_transcript_fa.collect())

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