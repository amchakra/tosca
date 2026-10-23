#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl=2

process STAR {

    tag "${sample_id}"
    label 'process_high'
    if(params.keep_intermediates) publishDir "${params.outdir}/premap", mode: 'copy', overwrite: true

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

    """
    STAR $args && samtools index -@ ${task.cpus} ${sample_id}.Aligned.sortedByCoord.out.bam
    """
}

process STAR_TX_INDEX {

    // STAR index of the transcript FASTA that BLAT maps to, for the linear-read filter when the
    // genome premap is skipped (e.g. viral runs) or linear_reference = 'transcript'.

    tag "${transcript_fa}"
    label 'process_high'

    input:
        path(transcript_fa)

    output:
        path("star_tx_index"), emit: index

    script:

    // Small references (viral genomes, test data) need genomeSAindexNbases and genomeChrBinNbits
    // scaled down, as recommended in the STAR manual.
    """
    stats=\$(awk '/^>/ {n++; next} {len += length(\$0)} END {print len, n}' $transcript_fa)
    genome_length=\$(echo \$stats | cut -d' ' -f1)
    n_refs=\$(echo \$stats | cut -d' ' -f2)
    sa_nbases=\$(awk -v L=\$genome_length 'BEGIN {x = int(log(L)/log(2)/2 - 1); print (x < 14 ? x : 14)}')
    chr_nbits=\$(awk -v L=\$genome_length -v n=\$n_refs 'BEGIN {x = int(log((L/n > 100 ? L/n : 100))/log(2)); print (x < 18 ? x : 18)}')

    mkdir star_tx_index
    STAR --runMode genomeGenerate --runThreadN ${task.cpus} \\
        --genomeDir star_tx_index --genomeFastaFiles $transcript_fa \\
        --genomeSAindexNbases \$sa_nbases --genomeChrBinNbits \$chr_nbits
    """
}

process STAR_TX {

    // Unspliced STAR alignment of the reads to the transcript FASTA, for the linear-read filter.

    tag "${sample_id}"
    label 'process_high'
    if(params.keep_intermediates) publishDir "${params.outdir}/premap", mode: 'copy', overwrite: true

    input:
        tuple val(sample_id), path(reads)
        path(star_tx_index)

    output:
        tuple val(sample_id), path("${sample_id}.tx.Aligned.sortedByCoord.out.bam"), path("${sample_id}.tx.Aligned.sortedByCoord.out.bam.bai"), emit: bam
        path("*.Log.final.out"), emit: log

    script:

    args = " --runThreadN ${task.cpus} "
    args += " --genomeDir $star_tx_index --genomeLoad NoSharedMemory "
    args += " --readFilesIn $reads --readFilesCommand gunzip -c "
    args += " --outFileNamePrefix ${sample_id}.tx. "
    args += " --outSAMattributes All --outSAMtype BAM SortedByCoordinate --outSAMunmapped Within " // need to keep unmapped for later filtering
    args += " --outFilterMultimapNmax 20 --outSAMmultNmax 1 "
    args += " --alignIntronMax 1 --alignSJoverhangMin 999 " // no splicing on a transcript/viral reference
    args += " --limitBAMsortRAM 60000000000"
    args += " " + params.linear_star_args

    """
    STAR $args && samtools index -@ ${task.cpus} ${sample_id}.tx.Aligned.sortedByCoord.out.bam
    """
}

process FILTER_PREMAP_READS {

    // Removes annotated-spliced and/or linear reads (bin/filter_premap_reads.py) before BLAT.
    // suffix distinguishes the first stage of a two-stage run (genome spliced filter, then
    // transcriptome linear filter); the final stage always writes <sample>.filter_premap_reads.log.

    tag "${sample_id}"
    label 'process_low'
    if(params.keep_intermediates) publishDir "${params.outdir}/filtered", mode: 'copy', overwrite: true

    input:
        tuple val(sample_id), path(bam), path(bai), path(upstream_log)
        val(filter_args)
        val(suffix)

    output:
        tuple val(sample_id), path("${sample_id}${suffix}.filtered.fastq.gz"), emit: fastq
        tuple val(sample_id), path("${sample_id}.filter_premap_reads${suffix}.log"), emit: log
        path("*.{spliced,linear}.bam"), optional: true, emit: discarded

    script:

    linear_args = "--linear-min-frac ${params.linear_min_frac} --linear-max-del ${params.linear_max_del}"
    upstream_arg = upstream_log ? "--upstream-log ${upstream_log}" : ""
    discarded_arg = params.keep_intermediates ? "--write-discarded" : ""

    """
    filter_premap_reads.py $bam ${sample_id}${suffix} $filter_args $linear_args $upstream_arg $discarded_arg > ${sample_id}.filter_premap_reads${suffix}.log
    """

}
