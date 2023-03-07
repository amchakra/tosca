#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl=2

process BLAT {
    tag "${sample_id}"
    if(!params.keep_cache) cache false

    input:
        tuple val(sample_id), path(reads)
        path(transcript_fa)

    output:
        tuple val(sample_id), path("${sample_id}.blast8.gz"), emit: blast8

    script:

    args = " -threads=${task.cpus} "
    args += " -stepSize=" + params.step_size 
    args += " -tileSize=" + params.tile_size
    args += " -minScore=" + params.min_score 
    args += " -out=blast8 "
    args += " -dots=1000000 "

    """
    pblat $args $transcript_fa $reads ${sample_id}.blast8
    pigz ${sample_id}.blast8
    """

}

process FILTER_BLAT {

    tag "${sample_id}"
    if(params.keep_intermediates) publishDir "${params.outdir}/mapped", mode: 'copy', overwrite: true

    input:
        tuple val(sample_id), path(blast8)

    output:
        tuple val(sample_id), path("${sample_id}.filtered.blast8.gz"), emit: blast8

    script:

    evalue = params.evalue
    maxhits = params.maxhits

    """
    filter_blat.py $blast8 ${sample_id}.filtered.blast8.gz $evalue $maxhits
    """

}

process BLAT_ALL_IN_ONE {

    tag "${sample_id}"
    if(!params.keep_cache) cache false

    cpus '8'
    memory '32G'
    time '8h'

    input:
        tuple val(sample_id), path(reads)
        path(transcript_fa)

    output:
        tuple val(sample_id), path("${sample_id}.hybrids.tsv.gz"), emit: hybrids

    script:

    args = " -threads=${task.cpus} "
    args += " -stepSize=" + params.step_size 
    args += " -tileSize=" + params.tile_size
    args += " -minScore=" + params.min_score 
    args += " -out=blast8 "
    args += " -dots=1000000 "

    evalue = params.evalue
    maxhits = params.maxhits

    """
    reformat.sh in1=$reads out1=${sample_id}.fasta

    pblat $args $transcript_fa ${sample_id}.fasta ${sample_id}.blast8
    pigz ${sample_id}.blast8
    
    filter_blat.py ${sample_id}.blast8.gz ${sample_id}.filtered.blast8.gz $evalue $maxhits
    rm ${sample_id}.blast8.gz

    identify_hybrids.R -t ${task.cpus} -b ${sample_id}.filtered.blast8.gz -f ${sample_id}.fasta -o ${sample_id}.hybrids.tsv.gz
    rm ${sample_id}.fasta
    """    

}
