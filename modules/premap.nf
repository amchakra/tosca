#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl=2

process STAR {

    tag "${sample_id}"
    if(params.keep_intermediates) publishDir "${params.outdir}/premap", mode: 'copy', overwrite: true

    // cpus '8'
    // memory '48G'
    // time '2h'

    input:
        tuple val(sample_id), path(reads)
        file(star_genome_index)
    
    output:
        tuple val(sample_id), path("${sample_id}.Aligned.sortedByCoord.out.bam"), path("${sample_id}.Aligned.sortedByCoord.out.bam.bai"), emit: bam
        path("*.Log.final.out"), emit: log

    script:

    args = " --runThreadN ${task.cpus} "
    args += " --genomeDir $star_genome_index --genomeLoad NoSharedMemory " 
    args += " --readFilesIn $reads --readFilesCommand gunzip -c "
    args += " --outFileNamePrefix ${sample_id}. "
    args += " --outSAMattributes All --outSAMtype BAM SortedByCoordinate --outSAMunmapped Within " // need to keep unmapped for later filtering
    args += " --outFilterMultimapNmax 20 --outSAMmultNmax 1 "
    args += " --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outFilterType BySJout "
    args += " --alignIntronMin 20 --alignIntronMax 100000 "
    args += " --limitBAMsortRAM 60000000000"
    args += " " + params.star_args

    cmd = "STAR $args && samtools index -@ ${task.cpus} ${sample_id}.Aligned.sortedByCoord.out.bam"

    if(params.verbose) { println ("[MODULE] PREMAP: " + cmd) }
    
    """
    $cmd
    """
}

process FILTER_SPLICED_READS {

    tag "${sample_id}"
    if(params.keep_intermediates) publishDir "${params.outdir}/filtered", mode: 'copy', overwrite: true

    time '2h'

    input:
        tuple val(sample_id), path(bam), path(bai)

    output:
        tuple val(sample_id), path("${sample_id}.unspliced.fastq.gz"), emit: fastq
        path("*.filter_spliced_reads.log"), emit: log

    script:

    cmd = "filter_spliced_reads.py $bam ${sample_id} > ${sample_id}.filter_spliced_reads.log"

    if(params.verbose) { println ("[MODULE] FILTERSPLICEDREADS: " + cmd) }

    """
    $cmd
    """
    
}