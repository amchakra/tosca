#!/usr/bin/env nextflow

// Define DSL2
nextflow.enable.dsl=2

include { TOSCA_QC; MULTIQC } from '../modules/multiqc.nf'

workflow MAKE_REPORT {

    take:
        filter_spliced_reads_logs
        dedup_logs
        raw_hybrids
        hybrids
        clusters
        multiqc_config

    main:
        TOSCA_QC(filter_spliced_reads_logs, dedup_logs, raw_hybrids, hybrids, clusters)
        MULTIQC(multiqc_config, TOSCA_QC.out)

}
