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

// Parameters
// Need to set these before module is loaded else not propagated
params.keep_intermediates = true

params.adapter = 'AGATCGGAAGAGC'
params.min_quality = 10
params.min_readlength = 16

params.split_size = 100000

params.evalue = 0.001
params.maxhits = 100

params.umi_separator = '_'
params.dedup_method = 'directional'

params.shuffled_mfe = false

params.percent_overlap = 0.75
params.sample_size = -1

params.goi = false

// Processes
include { hiclipheader } from './modules/utils.nf'
include { METADATA } from './modules/metadata.nf'
include { CUTADAPT } from './modules/cutadapt.nf'
include { PREMAP } from './workflows/premap.nf'
include { GET_HYBRIDS } from './workflows/gethybrids.nf'
include { GET_NON_HYBRIDS } from './modules/getnonhybrids.nf'
include { PROCESS_HYBRIDS } from './workflows/processhybrids.nf'
include { EXPORT_INTRAGENIC } from './workflows/exportbedbam.nf'
include { GET_CONTACT_MAPS } from './modules/getcontactmaps.nf'

// General variables
params.premap = true
if(params.org == 'rSARS-CoV-2' | params.org == 'SARS-CoV-2-England-2-2020') { params.virus = true }

// Genome variables
params.genome_fai = params.genomes[ params.org ].genome_fai
params.transcript_fa = params.genomes[ params.org ].transcript_fa
params.transcript_fai = params.genomes[ params.org ].transcript_fai
params.transcript_gtf = params.genomes[ params.org ].transcript_gtf
params.star_genome = params.genomes[ params.org ].star_genome

// Create channels for static files
ch_star_genome = Channel.fromPath(params.star_genome, checkIfExists: true)
ch_transcript_fa = Channel.fromPath(params.transcript_fa, checkIfExists: true)
ch_transcript_fai = Channel.fromPath(params.transcript_fai, checkIfExists: true)
ch_genome_fai = Channel.fromPath(params.genome_fai, checkIfExists: true)
ch_transcript_gtf = Channel.fromPath(params.transcript_gtf, checkIfExists: true)

// Channels for optional inputs
if(params.goi) ch_goi = Channel.fromPath(params.goi, checkIfExists: true)

// Show header
log.info hiclipheader()

def summary = [:]
summary['Output directory'] = params.outdir
summary['Trace directory'] = params.tracedir
if(workflow.repository) summary['Pipeline repository'] = workflow.repository
if(workflow.revision) summary['Pipeline revision'] = workflow.revision
summary['Pipeline directory'] = workflow.projectDir
summary['Run name'] = workflow.runName
summary['Profile'] = workflow.profile
if(workflow.container) summary['Container'] = workflow.container
if(params.keep_intermediates) summary['Keep intermediates'] = params.keep_intermediates
log.info summary.collect { k,v -> "${k.padRight(25)}: $v" }.join("\n")
log.info "-\033[2m---------------------------------------------------------------\033[0m-"

def settings = [:]
settings['Organism'] = params.org
settings['Adapter sequence'] = params.adapter
settings['Minimum read quality'] = params.min_quality
settings['Minimum read length'] = params.min_readlength
settings['Premapping'] = params.premap
settings['FASTQ split size'] = params.split_size
settings['Minimum e-value'] = params.evalue
settings['Maximum hits/read'] = params.maxhits
settings['Deduplication method'] = params.dedup_method
if(params.dedup_method != 'none') settings['UMI separator'] = params.umi_separator
settings['Shuffled binding energy'] = params.shuffled_mfe
settings['Clustering sample size'] = params.sample_size
settings['Clustering overlap'] = params.percent_overlap
if(params.goi) { settings['Genes for contact maps'] = params.goi } else { settings['Genes for contact maps'] = "none" }
log.info settings.collect { k,v -> "${k.padRight(25)}: $v" }.join("\n")
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
    PROCESS_HYBRIDS(GET_HYBRIDS.out.hybrids, ch_transcript_fa, ch_transcript_gtf)

    /* 
    EXPORT INTRAGENIC
    */
    EXPORT_INTRAGENIC(PROCESS_HYBRIDS.out.hybrids, PROCESS_HYBRIDS.out.clusters, ch_genome_fai)

    /* 
    GET CONTACT MAPS
    */
    if(params.goi) GET_CONTACT_MAPS(DEDUPLICATE.out.hybrids, ch_transcript_fai.collect(), ch_goi.collect())

}

workflow.onComplete {

    log.info "\nPipeline complete!\n"

}