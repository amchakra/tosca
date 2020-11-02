#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl=2

process splitfastq {
    tag "${sample_id}"
    // publishDir "${params.outdir}/split", mode: 'copy', overwrite: true

    cpus 8
    time '24h'

    input:
        tuple val(sample_id), path(reads)

    output:
        path("${sample_id}_*.fastq.gz")

    shell:
    """
    zcat $reads | split -l 4000000 --additional-suffix .fastq - ${sample_id}_
    pigz *.fastq
    """

}

// workflow splitfastq {

//     main:
//         split_fastq()
//         Channel
//             .fromPath(split_fastq.out)
//             .map { file -> tuple(file.baseName, file)}
//             .set { data }
//     emit: data

// }