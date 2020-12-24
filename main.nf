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

// Processes
include { hiclipheader } from './modules/utils.nf'
include { metadata } from './modules/metadata.nf'
include { trim } from './modules/trim.nf'
include { premap } from './modules/premap.nf'
include { filtersplicedreads } from './modules/filtersplicedreads.nf'
include { mapchimeras } from './modules/mapchimeras.nf'
include { deduplicate } from './modules/deduplicate.nf'
include { deduplicate_unique } from './modules/deduplicate.nf'
include { extracthybrids } from './modules/extracthybrids.nf'
include { getbindingenergy } from './modules/getbindingenergy.nf'
include { clusterhybrids } from './modules/clusterhybrids.nf'
include { collapseclusters } from './modules/collapseclusters.nf'
include { clusterbindingenergy } from './modules/clusterbindingenergy.nf'
include { convertcoordinates } from './modules/convertcoordinates.nf'
include { hybridbedtohybridbam } from './modules/hybridbedtohybridbam.nf'

include { splitfastq } from './modules/splitfastq.nf'
include { convert_fastq_to_fasta } from './modules/convert_fastq_to_fasta.nf'
include { mapblat } from './modules/mapblat.nf'
include { filterblat } from './modules/filterblat.nf'
include { identifyhybrids } from './modules/identifyhybrids.nf'
include { mergehybrids } from './modules/mergehybrids.nf'
include { deduplicate_blat } from './modules/deduplicate.nf'
include { getnonhybrids } from './modules/getnonhybrids.nf'

// Main workflow

// General variables
params.quickdedup = true

// Input variables
params.input='metadata.csv'
params.org='mouse'

// Genome variables
params.genome_fa = params.genomes[ params.org ].genome_fa
params.genome_fai = params.genomes[ params.org ].genome_fai
params.genome_gtf = params.genomes[ params.org ].genome_gtf 
params.transcript_fa = params.genomes[ params.org ].transcript_fa
params.transcript_gtf = params.genomes[ params.org ].transcript_gtf
params.star_genome = params.genomes[ params.org ].star_genome
params.star_transcript = params.genomes[ params.org ].star_transcript

// params.star_genome_index = '/camp/lab/luscomben/home/users/chakraa2/projects/flora/mouse/ref/STAR_GRCm38_GencodeM24'
// params.star_transcript_index = '/camp/lab/luscomben/home/users/chakraa2/projects/flora/mouse/ref/Mm_GencodeM24_rRNA_MT_genes'
// params.transcript_fasta_csv = '/camp/lab/luscomben/home/users/chakraa2/projects/flora/mouse/ref/Mm_GencodeM24_rRNA_MT_genes.csv.gz'
// params.genome_fai = '/camp/lab/luscomben/home/users/chakraa2/projects/flora/mouse/ref/GRCm38.primary_assembly.genome.fa.fai'
// params.transcript_gtf = '/camp/lab/luscomben/home/users/chakraa2/projects/flora/mouse/ref/Mm_GencodeM24_rRNA_MT_genes.gtf.gz'

// Create channels for static files
ch_star_genome = Channel.fromPath(params.star_genome, checkIfExists: true)
ch_star_transcript = Channel.fromPath(params.star_transcript, checkIfExists: true)
ch_transcript_fa = Channel.fromPath(params.transcript_fa, checkIfExists: true)
ch_genome_fai = Channel.fromPath(params.genome_fai, checkIfExists: true)
ch_transcript_gtf = Channel.fromPath(params.transcript_gtf, checkIfExists: true)

// Show header
log.info hiclipheader()
// def summary = [:]
// summary['Output directory'] = params.outdir
// summary['Trace directory'] = params.tracedir
// summary['Genome fasta'] = params.genome_fa
// summary['Genome fasta index'] = params.genome_fai
// summary['Genome annotation'] = params.genome_gtf
// summary['Transcriptome fasta'] = params.transcript_fa
// summary['Transcriptome annotation'] = params.transcript_gtf
// summary['STAR genome'] = params.star_genome
// summary['STAR transcriptome'] = params.star_transcript
// summary['Deduplicate quickly'] = params.quickdedup
// summary['Minimum intron length'] = params.intronmin

// log.info summary.collect { k,v -> "${k.padRight(25)}: $v" }.join("\n")
// log.info "-\033[2m---------------------------------------------------------------\033[0m-"

// Pipeline
workflow {

    // Get fastq paths 
    metadata(params.input)

    // metadata.out.view()

    // // Split
    // splitfastq(metadata.out)

    // ch_spl = splitfastq.out
    //     .flatten()
    //     .map { file -> tuple(file.simpleName, file) }

    // Trim
    // trim(ch_spl)
    trim(metadata.out)

    // Filter spliced reads
    premap(trim.out.combine(ch_star_genome))
    filtersplicedreads(premap.out)

    // Split
    splitfastq(filtersplicedreads.out)

    ch_spl = splitfastq.out
        .flatten()
        .map { file -> tuple(file.simpleName, file) }

    // Convert to fasta
    // convert_fastq_to_fasta(filtersplicedreads.out)
    convert_fastq_to_fasta(ch_spl)

    // Merge back test
    // ch_comb = convert_fastq_to_fasta.out
    //     .map { [ it[0].split('_')[0..-2].join('_'), it[1] ] }
    //     .groupTuple(by: 0)
    //     .view()

    // Map chimeras
    mapblat(convert_fastq_to_fasta.out.combine(ch_transcript_fa))
    filterblat(mapblat.out)

    // Identify hybrids
    identifyhybrids(filterblat.out.join(convert_fastq_to_fasta.out))

    // Merge hybrids
    ch_comb = identifyhybrids.out
        .map { [ it[0].split('_')[0..-2].join('_'), it[1] ] }
        .groupTuple(by: 0)
        // .view()

    mergehybrids(ch_comb)

    // Get non-hybrid reads for later
    getnonhybrids(mergehybrids.out.join(filtersplicedreads.out))
    // Map chimerias
    // mapchimeras(filtersplicedreads.out.combine(ch_star_transcript))

    // Remove PCR duplicates
    // if ( params.quickdedup ) {
    //     deduplicate_unique(mapchimeras.out)
    // } else {
    //     deduplicate(mapchimeras.out)
    // }

    deduplicate_blat(mergehybrids.out)

    // // Extract hybrids
    // if ( params.quickdedup ) {
    //     extracthybrids(deduplicate_unique.out.combine(ch_transcript_fa))
    // } else {
    //     extracthybrids(deduplicate.out.combine(ch_transcript_fa))
    // }
    // // Get binding energies
    // getbindingenergy(extracthybrids.out.combine(ch_transcript_fa))
    getbindingenergy(deduplicate_blat.out.combine(ch_transcript_fa))

    // // Get clusters
    clusterhybrids(getbindingenergy.out)

    // // Convert coordinates
    // // Write hybrid BAM
    convertcoordinates(clusterhybrids.out.combine(ch_transcript_gtf))
    hybridbedtohybridbam(convertcoordinates.out.combine(ch_genome_fai))

    // // Collapse clusters
    collapseclusters(clusterhybrids.out.combine(ch_transcript_gtf))
    clusterbindingenergy(collapseclusters.out.combine(ch_transcript_fa))

}

workflow.onComplete {

    log.info "\nPipeline complete!\n"

}