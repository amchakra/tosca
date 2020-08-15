#!/usr/bin/env nextflow

// Specify DSL2
nextflow.preview.dsl = 2

process premap {

    tag "${sample_id}"
    publishDir "${params.outdir}/premap", mode: 'copy', overwrite: true

    cpus 8
    time '12h'

    input:
        tuple val(sample_id), path(reads), path(star_genome_index)
    
    output:
        tuple val(sample_id), path("${sample_id}.Aligned.sortedByCoord.out.bam")

    shell:
    """
    STAR --runThreadN ${task.cpus} \
        --genomeDir $star_genome_index --genomeLoad NoSharedMemory \
        --readFilesIn $reads --readFilesCommand zcat \
        --outFileNamePrefix ${sample_id}. \
        --outFilterMultimapNmax 20 --outSAMmultNmax 1 \
        --outSAMunmapped Within \
        --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outFilterType BySJout \
        --alignIntronMin 20 --alignIntronMax 100000 \
        --outSAMattributes All --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 60000000000

    sambamba index -t ${task.cpus} ${sample_id}.Aligned.sortedByCoord.out.bam
    """
}