#!/usr/bin/env nextflow

/*
========================================================================================
                                    amchakra/tosca
========================================================================================
hiCLIP/proximity ligation analysis pipeline.
 #### Homepage / Documentation
 https://github.com/amchakra/tosca
----------------------------------------------------------------------------------------
*/

// Define DSL2
nextflow.enable.dsl=2

// Processes
include { hiclipheader } from './modules/utils.nf'
include { METADATA } from './modules/metadata.nf'
include { CUTADAPT } from './modules/cutadapt.nf'
include { PREMAP } from './workflows/premap.nf'
include { GET_HYBRIDS } from './workflows/gethybrids.nf'
include { GET_NON_HYBRIDS } from './modules/getnonhybrids.nf'
include { PROCESS_HYBRIDS } from './workflows/processhybrids.nf'
include { GET_VISUALISATIONS } from './workflows/getvisualisations.nf'
include { GET_ATLAS } from './workflows/getatlas.nf'
include { MAKE_REPORT } from './workflows/makereport.nf'

// Genome variables
if(params.org && params.genomesdir) {

    params.genome_fai = params.genomes[ params.org ].genome_fai
    params.transcript_fa = params.genomes[ params.org ].transcript_fa
    params.transcript_fai = params.genomes[ params.org ].transcript_fai
    params.transcript_gtf = params.genomes[ params.org ].transcript_gtf
    params.star_genome = params.genomes[ params.org ].star_genome
    params.regions_gtf = params.genomes[ params.org ].regions_gtf

} else {

    if(!params.genome_fai) { exit 1, "--genome_fai is not specified." } 
    if(!params.transcript_fa) { exit 1, "--transcript_fa is not specified." } 
    if(!params.transcript_fai) { exit 1, "--transcript_fai is not specified." } 
    if(!params.transcript_gtf) { exit 1, "--transcript_gtf is not specified." } 
    if(!params.regions_gtf) { exit 1, "--regions_gtf is not specified." } 

}


// Create channels for static files
ch_transcript_fa = Channel.fromPath(params.transcript_fa, checkIfExists: true)
ch_transcript_fai = Channel.fromPath(params.transcript_fai, checkIfExists: true)
ch_genome_fai = Channel.fromPath(params.genome_fai, checkIfExists: true)
ch_transcript_gtf = Channel.fromPath(params.transcript_gtf, checkIfExists: true)
ch_regions_gtf = Channel.fromPath(params.regions_gtf, checkIfExists: true)

// Channels for optional inputs
if(!params.skip_premap) {
    ch_star_genome = Channel.fromPath(params.star_genome, checkIfExists: true)
} else {
    ch_star_genome = Channel.empty()
}

if(params.goi) {
    ch_goi = Channel.fromPath(params.goi, checkIfExists: true) 
} else {
    ch_goi = Channel.empty()
}

// Channel for MultiQC config
ch_multiqc_config = Channel.fromPath(params.multiqc_config, checkIfExists: true)

// Show header
log.info hiclipheader()

def summary = [:]
summary['Output directory'] = params.outdir
summary['Trace directory'] = params.tracedir
if(workflow.repository) summary['Pipeline repository'] = workflow.repository
if(workflow.revision) summary['Pipeline revision'] = workflow.revision
summary['Pipeline directory'] = workflow.projectDir
summary['Working dir'] = workflow.workDir
summary['Run name'] = workflow.runName
summary['Profile'] = workflow.profile
if(workflow.container) summary['Container'] = workflow.container
if(params.keep_intermediates) summary['Keep intermediates'] = params.keep_intermediates
if(!params.keep_cache) summary['Keep cache'] = params.keep_cache
log.info summary.collect { k,v -> "${k.padRight(25)}: $v" }.join("\n")
log.info "-\033[2m---------------------------------------------------------------\033[0m-"

def settings = [:]
settings['Organism'] = params.org
if(params.skip_qc) { settings['Skip QC'] = params.skip_qc } 
// if(params.skip_atlas) { settings['Skip atlas generation'] = params.skip_atlas }
if(params.skip_premap) { settings['Skip premapping'] = params.skip_premap } 
settings['Adapter sequence'] = params.adapter
settings['Minimum read quality'] = params.min_quality
settings['Minimum read length'] = params.min_readlength
settings['FASTQ split size'] = params.split_size
settings['Minimum e-value'] = params.evalue
settings['Maximum hits/read'] = params.maxhits
settings['Deduplication method'] = params.dedup_method
if(params.dedup_method != 'none') settings['UMI separator'] = params.umi_separator
if(params.slurm) settings['Use SLURM'] = params.slurm
settings['Clustering chunk number'] = params.chunk_number
settings['Clustering sample size'] = params.sample_size
settings['Clustering overlap'] = params.percent_overlap
settings['Analyse structures'] = params.analyse_structures
if(params.analyse_structures) settings['Analyse clusters only'] = params.clusters_only
if(params.analyse_structures) settings['Analyse shuffled energies'] = params.shuffled_energies

if(params.goi) { settings['Genes of interest'] = params.goi }
if(params.goi) { settings['Bin size for contact maps'] = params.bin_size } 
if(params.goi) { settings['Breaks for arcs'] = params.breaks } 
log.info settings.collect { k,v -> "${k.padRight(25)}: $v" }.join("\n")
log.info "-----------------------------------------------------------------"

if(atlas) {

    ch_all_hybrids = Channel.fromPath(params.input)
    GET_ATLAS(ch_all_hybrids, ch_transcript_gtf, ch_regions_gtf, ch_genome_fai)

} else {

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
        if(!params.skip_premap) {
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
        PROCESS_HYBRIDS(GET_HYBRIDS.out.hybrids, ch_transcript_fa, ch_transcript_gtf, ch_regions_gtf)

        // /* 
        // GET ATLAS
        // */
        // if(!params.skip_atlas) {
        //     GET_ATLAS(PROCESS_HYBRIDS.out.hybrids, ch_transcript_gtf, ch_regions_gtf, ch_genome_fai)
        // }

        /* 
        GET VISUALISATIONS
        */
        GET_VISUALISATIONS(PROCESS_HYBRIDS.out.hybrids, PROCESS_HYBRIDS.out.clusters, ch_genome_fai, ch_transcript_fai, ch_goi)

        /* 
        MAKE REPORT
        */
        if(!params.skip_qc) {
            MAKE_REPORT(PREMAP.out.logs.collect(), GET_HYBRIDS.out.logs.collect(), GET_HYBRIDS.out.raw_hybrids.collect{it[1]}, PROCESS_HYBRIDS.out.hybrids.collect{it[1]}, PROCESS_HYBRIDS.out.clusters.collect{it[1]}, ch_multiqc_config)
        }

    }

}

workflow.onComplete {

    if (workflow.success) {
        log.info "-\033[0;34m[Tosca]\033[0;32m Pipeline completed successfully\033[0m-\n"
    } else {
        log.info "-\033[0;34m[Tosca]\033[1;91m Pipeline completed with errors\033[0m\n"
    }

}