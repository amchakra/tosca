#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl=2

ch_multiqc_config = file("$projectDir/assets/multiqc_config.yaml", checkIfExists: true)


process TOSCA_QC {
    tag "${workflow.runName}"

    cpus 4
    memory '32 G'
    time '1h'

    publishDir "${params.outdir}/multiqc", mode: 'copy', overwrite: true

    input:
    path(filter_spliced_reads_logs)
    path(dedup_logs)
    path(raw_hybrids)
    path(hybrids)
    path(clusters)

    output:
    path("*.tsv")

    script:
    """
    tosca_qc.R
    """

}

process MULTIQC {
    tag "${workflow.runName}"

    cpus 1
    memory '16 G'
    time '1h'

    publishDir "${params.outdir}/multiqc", mode: 'copy', overwrite: true

    input:
    path(multiqc_config)
    path(tsv)

    output:
    path("*multiqc_report.html")
    path("*_data")
    path("multiqc_plots")

    script:

    rtitle = ''
    rfilename = ''
    if (!(workflow.runName ==~ /[a-z]+_[a-z]+/)) {
        rtitle = "--title \"${workflow.runName}\""
        rfilename = "--filename " + workflow.runName.replaceAll('\\W','_').replaceAll('_+','_') + "_multiqc_report"
    }

    custom_config_file = "--config ${params.multiqc_config}"

    """
    multiqc -f -p $rtitle $rfilename $custom_config_file .
    """

}