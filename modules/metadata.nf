#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl=2

workflow METADATA {
    take: csv
    main:
        Channel
            .fromPath( csv )
            .splitCsv(header:true)
            .map { row -> [ row.sample, file(row.fastq, checkIfExists: true) ]  }
            .set { data }
    emit:
        data
}