#!/usr/bin/env nextflow

// Specify DSL2
nextflow.preview.dsl = 2

process mapchimeras {

    tag "${sample_id}"
    publishDir "${params.outdir}/mapped", mode: 'copy', overwrite: true

    cpus 8
    time '12h'

    input:
        tuple val(sample_id), path(reads), path(star_transcript_index)

    output:
        tuple val(sample_id), path("${sample_id}.Aligned.sortedByCoord.out.bam"), path("${sample_id}.Aligned.sortedByCoord.out.bam.bai")

    shell:

    intronmin = params.intronmin

    """
    STAR --runThreadN ${task.cpus} \
    --genomeDir $star_transcript_index --genomeLoad NoSharedMemory \
    --readFilesIn $reads --readFilesCommand zcat \
    --outFileNamePrefix ${sample_id}. \
    --outFilterMultimapNmax 1 \
    --alignIntronMin ${intronmin} --scoreGapNoncan 0 --scoreGapATAC 0 --scoreGapGCAG 0 --scoreGap 0 \
    --chimSegmentMin 12 --chimJunctionOverhangMin 12  --chimScoreJunctionNonGTAG 0 \
    --chimNonchimScoreDropMin 10 --chimOutType WithinBAM \
    --alignSJoverhangMin 12 --alignSJDBoverhangMin 12 --alignSJstitchMismatchNmax 5 -1 5 5 \
    --outSAMattributes All --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 9839723217

    sambamba index -t ${task.cpus} ${sample_id}.Aligned.sortedByCoord.out.bam
    """

}