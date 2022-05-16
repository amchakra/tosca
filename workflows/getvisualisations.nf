#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl=2

include { EXPORT_GENOMIC_BED as EXPORT_HYBRID_GENOMIC_BED; 
          EXPORT_GENOMIC_BED as EXPORT_CLUSTERS_GENOMIC_BED;
          EXPORT_GENOMIC_BAM;
          GET_CONTACT_MAPS;
          GET_ARCS } from '../modules/getvisualisations.nf'


workflow GET_VISUALISATIONS {

    take:
    hybrids         // channel: hybrids
    clusters        // channel: clusters
    genome_fai      // channel: genome_fai
    transcript_fai      // channel: transcript_fai
    goi

    main:
    
    EXPORT_HYBRID_GENOMIC_BED("hybrids", hybrids)
    EXPORT_GENOMIC_BAM(EXPORT_HYBRID_GENOMIC_BED.out.bed, genome_fai.collect())
    EXPORT_CLUSTERS_GENOMIC_BED("clusters", clusters)

    if(params.goi) {

        GET_CONTACT_MAPS(ch_hybrids, ch_transcript_fai.collect(), ch_goi.collect())
        GET_ARCS(ch_clusters, ch_goi.collect())

    }

    emit:
    hybrid_bed = EXPORT_HYBRID_GENOMIC_BED.out.bed
    hybrid_bam = CONVERT_BED_TO_BAM.out.bam
    cluster_bed = EXPORT_CLUSTERS_GENOMIC_BED.out.bed

}