#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl=2

process deduplicate {

    tag "${sample_id}"
    publishDir "${params.outdir}/dedup", mode: 'copy', overwrite: true

    memory '64G'
    time '24h'

    input:
        tuple val(sample_id), path(reads), path(bai)

    output:
        tuple val(sample_id), path("${sample_id}.dedup.bam"), path("${sample_id}.dedup.bam.bai")

    shell:
    """
    if `samtools view $reads | cut -f 1 | head -1 | grep -q ':rbc:'`; then
        echo iCount
        umi_tools dedup --method directional --umi-separator=":" -I $reads -S ${sample_id}.dedup.bam
    else
        echo UMItools
        umi_tools dedup --method directional -I $reads -S ${sample_id}.dedup.bam
    fi
    
    sambamba index ${sample_id}.dedup.bam
    """

}

process deduplicate_unique {

    tag "${sample_id}"
    publishDir 'results/dedup', mode: 'copy', overwrite: false

    memory '64G'
    time '24h'

    input:
        tuple val(sample_id), path(reads), path(bai)

    output:
        tuple val(sample_id), path("${sample_id}.dedup.bam"), path("${sample_id}.dedup.bam.bai")

    shell:
    """
    if `samtools view $reads | cut -f 1 | head -1 | grep -q ':rbc:'`; then
        echo iCount
        umi_tools dedup --method unique --umi-separator=":" -I $reads -S ${sample_id}.dedup.bam
    else
        echo UMItools
        umi_tools dedup --method unique -I $reads -S ${sample_id}.dedup.bam
    fi
    
    sambamba index ${sample_id}.dedup.bam
    """

}

process deduplicate_blat {

    tag "${sample_id}"
    publishDir "${params.outdir}/dedup", mode: 'copy', overwrite: false

    memory '32G'
    time '24h'

    input:
        tuple val(sample_id), path(hybrids)

    output:
        tuple val(sample_id), path("${sample_id}.uniquehybrids.tsv.gz")

    script:

    """
    #!/usr/bin/env Rscript

    suppressPackageStartupMessages(library(data.table))
    suppressPackageStartupMessages(library(primavera))

    hybrid.dt <- fread("$hybrids")
    hybrid.dt\$hybrid_selection <- factor(hybrid.dt\$hybrid_selection, levels = c("single", "multi_overlap", "qual", "intra"))

    # Remove PCR duplicates
    hybrid.dt[, rbc := sub(".*\\\\:", "", read)]
    setorder(hybrid.dt, hybrid_selection) # Order so keep single if present
    unique.hybrid.dt <- unique(hybrid.dt, by = c("L_seqnames", "L_start", "L_end", "R_seqnames", "R_start", "R_end", "rbc"))

    message(paste0("PCR duplication ratio: ", round(nrow(hybrid.dt)/nrow(unique.hybrid.dt), 2)))

    # Get orientations where appropriate
    unique.hybrid.dt <- ReorientHybrids(unique.hybrid.dt)
    unique.hybrid.dt[L_seqnames == R_seqnames, orientation := ifelse(L_read_start <= R_read_start, "genomic", "reverse")]

    # Refashion hybrid.dt
    unique.hybrid.dt[, `:=` (L_width = L_end - L_start + 1,
                            L_strand = "+",
                            R_width = R_end - R_start + 1,
                            R_strand = "+")]

    unique.hybrid.dt <- unique.hybrid.dt[, .(read, L_seqnames, L_start, L_end, L_width, L_strand, R_seqnames, R_start, R_end, R_width, R_strand, orientation, hybrid_selection)]
    setnames(unique.hybrid.dt, "read", "name")

    fwrite(unique.hybrid.dt , file = paste0("$sample_id", ".uniquehybrids.tsv.gz"), sep = "\t")
    """

}