#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl=2

process PREMAP {

    tag "${sample_id}"
    publishDir "${params.outdir}/premap", mode: 'copy', overwrite: true

    cpus 8
    memory 48G
    time '12h'

    input:
        tuple val(sample_id), path(reads)
        file(star_genome_index)
    
    output:
        tuple val(sample_id), path("${sample_id}.Aligned.sortedByCoord.out.bam"), path("${sample_id}.Aligned.sortedByCoord.out.bam.bai"), emit: bam

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

    cmd = "STAR $args && samtools index -@ ${task.cpus} ${sample_id}.Aligned.sortedByCoord.out.bam"

    if(params.verbose) { println ("[MODULE] PREMAP: " + cmd) }
    
    """
    $cmd
    """
}