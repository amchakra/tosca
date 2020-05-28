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
nextflow.preview.dsl=2

// Processes

include hiclipheader from './modules/utils.nf'
include metadata from './modules/metadata.nf'
include trim from './modules/trim.nf'
include premap from './modules/premap.nf'
include filtersplicedreads from './modules/filtersplicedreads.nf'
include mapchimeras from './modules/mapchimeras.nf'
include deduplicate from './modules/deduplicate.nf'
include extracthybrids from './modules/extracthybrids.nf'
include getbindingenergy from './modules/getbindingenergy.nf'
include clusterhybrids from './modules/clusterhybrids.nf'
include collapseclusters from './modules/collapseclusters.nf'
include convertcoordinates from './modules/convertcoordinates.nf'
include hybridbedtohybridbam from './modules/hybridbedtohybridbam.nf'


// Main workflow

params.input='metadata.csv'
params.star_genome_index = '/camp/lab/luscomben/home/users/chakraa2/projects/flora/mouse/ref/STAR_GRCm38_GencodeM24'
params.star_transcript_index = '/camp/lab/luscomben/home/users/chakraa2/projects/flora/mouse/ref/Mm_GencodeM24_rRNA_MT_genes'
params.transcript_fasta_csv = '/camp/lab/luscomben/home/users/chakraa2/projects/flora/mouse/ref/Mm_GencodeM24_rRNA_MT_genes.csv.gz'
params.genome_fai = '/camp/lab/luscomben/home/users/chakraa2/projects/flora/mouse/ref/GRCm38.primary_assembly.genome.fa.fai'
params.transcript_gtf = '/camp/lab/luscomben/home/users/chakraa2/projects/flora/mouse/ref/Mm_GencodeM24_rRNA_MT_genes.gtf.gz'
params.genome_csv = '/camp/lab/luscomben/home/users/chakraa2/projects/flora/mouse/ref/Mm_GencodeM24_rRNA_MT_genes.csv.gz'

// Create channels for static files
ch_star_genome_index = Channel.fromPath(params.star_genome_index, checkIfExists: true)
ch_star_transcript_index = Channel.fromPath(params.star_transcript_index, checkIfExists: true)
ch_transcript_fasta_csv = Channel.fromPath(params.transcript_fasta_csv, checkIfExists: true)
ch_genome_fai = Channel.fromPath(params.genome_fai, checkIfExists: true)
ch_transcript_gtf = Channel.fromPath(params.transcript_gtf, checkIfExists: true)
ch_genome_csv = Channel.fromPath(params.genome_csv, checkIfExists: true)

// Show banner
log.info hiclipheader()

workflow {

    // Get fastq paths 
    metadata(params.input)

    // Trim
    trim(metadata.out)

    // Filter spliced reads
    premap(trim.out.combine(ch_star_genome_index))
    filtersplicedreads(premap.out)

    // Map chimerias
    mapchimeras(filtersplicedreads.out.combine(ch_star_transcript_index))

    // Remove PCR duplicates
    deduplicate(mapchimeras.out)

    // Extract hybrids
    // Get binding energies
    extracthybrids(deduplicate.out.combine(ch_genome_csv))
    getbindingenergy(extracthybrids.out.combine(ch_transcript_fasta_csv))

    // Get clusters
    clusterhybrids(getbindingenergy.out)

    // Convert coordinates
    // Write hybrid BAM
    convertcoordinates(clusterhybrids.out.combine(ch_transcript_gtf))
    hybridbedtohybridbam(convertcoordinates.out.combine(ch_genome_fai))

    // Collapse clusters
    collapseclusters(clusterhybrids.out.combine(ch_transcript_gtf))

}

workflow.onComplete {

    // def msg = """\
    //         Pipeline execution summary
    //         ---------------------------
    //         Completed at: ${workflow.complete}
    //         Duration    : ${workflow.duration}
    //         Success     : ${workflow.success}
    //         workDir     : ${workflow.workDir}
    //         exit status : ${workflow.exitStatus}
    //         """
    //         .stripIndent()

    //     sendMail(to: 'anob.chakrabarti@crick.ac.uk', subject: 'hiCLIP pipeline execution', body: msg)

    log.info "\nPipeline complete!\n"

}