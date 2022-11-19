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

            GET_CONTACT_MAPS(hybrids, transcript_fai.collect(), goi.collect())

            if(!params.skip_arcs) {
                GET_ARCS(clusters, goi.collect())
            }

        }

}