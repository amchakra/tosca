#!/usr/bin/env nextflow

/*
========================================================================================
                                    amchakra/tosca
========================================================================================
hiCLIP analysis pipeline.
 #### Homepage / Documentation
 https://github.com/amchakra/tosca
----------------------------------------------------------------------------------------
*/

// Define DSL2
nextflow.preview.dsl=2

// Processes

workflow metadata {
    take: csv
    main:
        Channel
            .fromPath( csv )
            .splitCsv(header:true)
            .map { row -> [ row.sample_id, file(row.fastq, checkIfExists: true) ]  }
            .set { data }
    emit:
        data
}

process trim {
    tag "${sample_id}"
    publishDir 'results/trimmed', mode: 'copy', overwrite: false

    cpus 8

    input:
        tuple val(sample_id), path(reads)

    output:
        tuple val(sample_id), path("${sample_id}.trimmed.fastq.gz")

    shell:
    """
    cutadapt -j ${task.cpus} --minimum-length 16 -q 10 -a AGATCGGAAGAGC -o "${sample_id}.trimmed.fastq.gz" $reads > "${sample_id}_cutadapt.log"
    """
}

process premap {
    tag "${sample_id}"

    cpus 8

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
        --outFilterMismatchNoverReadLmax 0.04 \
        --outSAMattributes All --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 60000000000

    sambamba index -t ${task.cpus} ${sample_id}.Aligned.sortedByCoord.out.bam
    """

}

process filtersplicedreads {
    tag "${sample_id}"
    publishDir 'results/filtered', mode: 'copy', overwrite: false

    input:
        tuple val(sample_id), path(bam)

    output:
        tuple val(sample_id), path("${sample_id}.unspliced.fastq.gz")

    script:
    """
    #!/usr/bin/env python

    # Script to filter out reads that map to annotated splice junctions
    # A. M. Chakrabarti
    # 9th October 2018

    import sys
    import pysam
    from subprocess import run

    # ==========
    # Filtering function
    # ==========

    def FilterBam(bam_in, bam_out):

        readcount=0
        keptcount=0

        for read in bam_in.fetch(until_eof = True):
            readcount=readcount+1

            if read.is_unmapped:
                keptcount=keptcount+1
                bam_out.write(read)

            else:
                jM=read.get_tag("jM")
                if any(tag < 20 for tag in jM):
                    keptcount=keptcount+1
                    bam_out.write(read)

        bam_in.close()
        bam_out.close()

        # Print metrics
        print("Total reads:", readcount)
        print("Reads kept:", keptcount)
        print("Reads discarded:", readcount - keptcount)

    f_in = "$bam"
    f_out = f_in.replace('.Aligned.sortedByCoord.out.bam', '.unspliced.bam')

    bam_in = pysam.AlignmentFile(f_in, "rb")
    bam_out = pysam.AlignmentFile(f_out, 'wb', template = bam_in)

    FilterBam(bam_in, bam_out)

    # System call to bedtools to convert to fastq
    fastq_out = f_in.replace('.Aligned.sortedByCoord.out.bam', '.unspliced.fastq.gz')
    run(f'bedtools bamtofastq -i {f_out} -fq /dev/stdout | pigz > {fastq_out}', shell = True)
    """
}

process mapchimeras {

    tag "${sample_id}"
    publishDir 'results/mapped', mode: 'copy', overwrite: false

    cpus 8

    input:
        tuple val(sample_id), path(reads), path(star_transcript_index)

    output:
        tuple val(sample_id), path("${sample_id}.Aligned.sortedByCoord.out.bam"), path("${sample_id}.Aligned.sortedByCoord.out.bam.bai")

    shell:
    """
    STAR --runThreadN ${task.cpus} \
    --genomeDir $star_transcript_index --genomeLoad NoSharedMemory \
    --readFilesIn $reads --readFilesCommand zcat \
    --outFileNamePrefix ${sample_id}. \
    --outFilterMultimapNmax 1 \
    --alignIntronMin 10 --scoreGapNoncan 0 --scoreGapATAC 0 --scoreGapGCAG 0 --scoreGap 0 \
    --chimSegmentMin 12 --chimJunctionOverhangMin 12  --chimScoreJunctionNonGTAG 0 \
    --chimNonchimScoreDropMin 10 --chimOutType WithinBAM \
    --alignSJoverhangMin 12 --alignSJDBoverhangMin 12 --alignSJstitchMismatchNmax 5 -1 5 5 \
    --outSAMattributes All --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 9839723217

    sambamba index -t ${task.cpus} ${sample_id}.Aligned.sortedByCoord.out.bam
    """

}

process deduplicate {

    tag "${sample_id}"
    publishDir 'results/dedup', mode: 'copy', overwrite: false

    input:
        tuple val(sample_id), path(reads), path(bai)

    output:
        tuple val(sample_id), path("${sample_id}.dedup.bam"), path("${sample_id}.dedup.bam.bai")

    shell:
    """
    umi_tools dedup --umi-separator=":" -I $reads -S ${sample_id}.dedup.bam
    sambamba index ${sample_id}.dedup.bam
    """

}

process extracthybrids {

    tag "${sample_id}"
    publishDir 'results/hybrids', mode: 'copy', overwrite: false

    cpus 4
    memory '32G'
    time '12h'

    input:
        tuple val(sample_id), path(reads), path(bai)

    output:
        tuple val(sample_id), path("${sample_id}.hybrids.tsv.gz")

    script:
    """
    #!/usr/bin/env Rscript

    suppressPackageStartupMessages(library(data.table))
    suppressPackageStartupMessages(library(primavera))
    suppressPackageStartupMessages(library(tictoc))

    # Extract hybrids
    message("Extracting hybrids...")
    tic()
    hybrids.dt <- ExtractHybridsWithinBAM(aligned.bam = "$reads")
    hybrids.dt <- ReorientHybrids(hybrids.dt)
    toc()

    f_out <- gsub("dedup.bam", "hybrids.tsv.gz", "$reads")
    fwrite(hybrids.dt, f_out, sep = "\t")

    message("Completed!")
    """

}

process getbindingenergy {

    tag "${sample_id}"
    publishDir 'results/hybrids', mode: 'copy', overwrite: false

    time '24h'

    input:
        tuple val(sample_id), path(hybrids), path(transcript_fasta_csv)

    output:
        tuple val(sample_id), path("${sample_id}.hybrids.mfe.tsv.gz")

    script:
    """
    #!/usr/bin/env Rscript

    suppressPackageStartupMessages(library(data.table))
    suppressPackageStartupMessages(library(primavera))
    suppressPackageStartupMessages(library(rslurm))
    suppressPackageStartupMessages(library(tictoc))
    suppressPackageStartupMessages(library(parallel))

    # Load genome
    message("Loading genome...")
    tic()
    genome.dt <- fread("$transcript_fasta_csv")
    toc()

    hybrids.dt <- fread("$hybrids")

    # Get sequences
    message("Getting sequences...")
    tic()
    hybrids.dt <- GetSequence(hybrids.dt = hybrids.dt, genome.dt = genome.dt)
    toc()

    # Add ID
    hybrids.dt[, id := paste0("H", 1:.N)]

    # Getting MFE
    message("Calculating MFE...")
    tic()
    
    seq.df <- data.frame(hybrids.dt[, .(id, L_sequence, R_sequence)]) # Split out relevant part of hybrids.dt

    # Cluster jobs
    sjob <- slurm_apply(.slurm_GetMFE, seq.df, jobname = sapply(strsplit(basename("$hybrids"), "\\\\."), "[[", 1), nodes = 100, cpus_per_node = 1, slurm_options = list(time = "24:00:00"), submit = TRUE)
    Sys.sleep(60) # To give it enough time to submit before the first check

    status <- FALSE
    while(status == FALSE) {

        squeue.out <- system(paste("squeue -n", sjob\$jobname), intern = TRUE) # Get contents of squeue for this job
        if(length(squeue.out) == 1) status <- TRUE # i.e. only the header left
        Sys.sleep(60)

    }

    mfe <- get_slurm_out(sjob, outtype = 'raw')

    # Remove temporary files
    cleanup_files(sjob) 

    # Merge back
    names(mfe) <- NULL # Not sure why it is names with L_sequence - not anymore...?
    mfe.dt <- as.data.table(unlist(mfe), keep.rownames = TRUE)
    setnames(mfe.dt, c("id", "mfe"))
    setkey(mfe.dt, id)
    mfe.dt[mfe > 0, mfe := 0] # Get rid of positives
    mfe.dt[is.na(mfe), mfe := 0] # Fix for positives < 10?

    setkey(hybrids.dt, id)

    hybrids.dt <- merge(hybrids.dt, mfe.dt, by = "id")
    toc()

    f_out <- gsub("hybrids.tsv.gz", "hybrids.mfe.tsv.gz", "$hybrids")
    fwrite(hybrids.dt, f_out, sep = "\t")

    message("Completed!")
    """

}

process clusterhybrids {

    tag "${sample_id}"
    publishDir 'results/hybrids', mode: 'copy', overwrite: false

    cpus 8
    time '24h'

    input:
        tuple val(sample_id), path(hybrids)

    output:
        tuple val(sample_id), path("${sample_id}.intragenic_hybrids.mfe.clusters.tsv.gz")

    script:
    """
    #!/usr/bin/env Rscript

    suppressPackageStartupMessages(library(data.table))
    suppressPackageStartupMessages(library(primavera))
    suppressPackageStartupMessages(library(tictoc))

    setDTthreads(8)

    ptm <- proc.time()

    # Load hybrids
    hybrids.dt <- fread("$hybrids")

    # Get intragenic hybrids
    intragenic.hybrids.dt <- hybrids.dt[L_seqnames == R_seqnames][grep("ENSMUSG", L_seqnames)]
    fwrite(intragenic.hybrids.dt, gsub("\\\\.hybrids\\\\.", "\\\\.intragenic_hybrids\\\\.", "$hybrids"), sep = "\t")

    # Get Cluster
    message("Clustering...")
    tic()

    # Split out by gene
    intragenic.hybrids.list <- split(intragenic.hybrids.dt, intragenic.hybrids.dt\$L_seqnames)
    solo.intragenic.hybrids.list <- intragenic.hybrids.list[S4Vectors::elementNROWS(intragenic.hybrids.list) == 1] # Remove solos to add in later
    message(length(solo.intragenic.hybrids.list), " genes with one hybrid")
    toomany.intragenic.hybrids.list <- intragenic.hybrids.list[S4Vectors::elementNROWS(intragenic.hybrids.list) > 5000] # Remove too many
    message(length(toomany.intragenic.hybrids.list), " genes with >5000 hybrids")
    intragenic.hybrids.list <- intragenic.hybrids.list[S4Vectors::elementNROWS(intragenic.hybrids.list) > 1]
    intragenic.hybrids.list <- intragenic.hybrids.list[S4Vectors::elementNROWS(intragenic.hybrids.list) <= 5000]
    message(length(intragenic.hybrids.list), " genes to cluster")

    # TODO: add in check for length

    library(tictoc)
    tic()
    intragenic.hybrids.clusters.list <- lapply(1:length(intragenic.hybrids.list), function(i) {

    # message(i)
    ClusterHybrids(intragenic.hybrids.list[[i]], percent_overlap = 0.5)

    })

    # Name and id order flipped for genes without clusters, because of merging clusters back in, hence use.names = TRUE
    intragenic.hybrids.clusters.dt <- rbindlist(intragenic.hybrids.clusters.list, use.names = TRUE)
    solo.intragenic.hybrids.dt <- rbindlist(solo.intragenic.hybrids.list, use.names = TRUE) # Add solos back in
    toomany.intragenic.hybrids.dt <- ifelse(length(toomany.intragenic.hybrids.list) == 0, data.table(), rbindlist(toomany.intragenic.hybrids.list, use.names = TRUE)[, cluster := Inf]) # Not always have these then end up with a data frame with just Inf
    intragenic.hybrids.clusters.dt <- rbind(intragenic.hybrids.clusters.dt, solo.intragenic.hybrids.dt, toomany.intragenic.hybrids.dt, use.names = TRUE, fill = TRUE)
    toc()

    stopifnot(nrow(intragenic.hybrids.clusters.dt) == nrow(intragenic.hybrids.dt))

    f_out <- gsub("hybrids.mfe.tsv.gz", "intragenic_hybrids.mfe.clusters.tsv.gz", "$hybrids")
    fwrite(intragenic.hybrids.clusters.dt, f_out, sep = "\t")

    message("Completed!")
    """
}

process collapseclusters {

    tag "${sample_id}"
    publishDir 'results/clusters', mode: 'copy', overwrite: false

    cpus 8
    time '24h'

    input:
        tuple val(sample_id), path(hybrids), path(transcript_gtf)

    output:
        tuple val(sample_id), path("${sample_id}.clusters.tsv.gz"), path("${sample_id}.clusters.bed.gz")

    script:
    """
    #!/usr/bin/env Rscript

    suppressPackageStartupMessages(library(optparse))
    suppressPackageStartupMessages(library(data.table))
    suppressPackageStartupMessages(library(primavera))
    suppressPackageStartupMessages(library(tictoc))

    setDTthreads(8)

    hybrids.dt <- fread("$hybrids")

    clusters.dt <- CollapseClusters(hybrids.dt)
    f_out <- gsub(".intragenic_hybrids.mfe.clusters.tsv.gz", ".clusters.tsv.gz", "$hybrids")
    fwrite(clusters.dt, f_out, sep = "\t")

    message(nrow(clusters.dt[L_end >= R_start]), " clusters removed")
    clusters.dt <- clusters.dt[L_end < R_start]

    # ========== SPLIT OUT?

    # New method
    message("Converting coordinates...")
    genes.gr <- rtracklayer::import.gff2("$transcript_gtf")
    tic()
    f_out <- gsub("tsv.gz", "bed", f_out)
    ExportGenomicBED(seq.dt = clusters.dt, genes.gr = genes.gr, sam_tag = TRUE, filename = gsub("tsv.gz", "bed", f_out))
    system(paste("pigz", f_out))
    toc()

    message("Completed!")
    """
}

process convertcoordinates {

    tag "${sample_id}"
    publishDir 'results/bed', mode: 'copy', overwrite: false

    cpus 8
    time '24h'

    input:
        tuple val(sample_id), path(hybrids), path(transcript_gtf)

    output:
        tuple val(sample_id), path("${sample_id}.intragenic_hybrids.mfe.clusters.bed.gz")

    script:
    """
    #!/usr/bin/env Rscript

    suppressPackageStartupMessages(library(optparse))
    suppressPackageStartupMessages(library(data.table))
    suppressPackageStartupMessages(library(primavera))
    suppressPackageStartupMessages(library(tictoc))

    setDTthreads(8)

    genes.gr <- rtracklayer::import.gff2("$transcript_gtf")

    hybrids.dt <- fread("$hybrids")
    seq.dt <- hybrids.dt[overlapping_hybrid %in% c(NA, FALSE)] # Remove overlapping hybrids, NA is genomic orientation
    seq.dt <- seq.dt[grep("Mt", L_seqnames, invert = TRUE)] # Remove MT for now

    # New method
    message("Converting coordinates...")
    tic()
    f_out <- gsub("tsv.gz", "bed", "$hybrids")
    ExportGenomicBED(seq.dt = seq.dt, genes.gr = genes.gr, sam_tag = TRUE, filename = f_out)
    system(paste("pigz", f_out))
    toc()

    message("Completed!")
    """
}

process hybridbedtohybridbam {

    tag "${sample_id}"
    publishDir 'results/bam', mode: 'copy', overwrite: false

    cpus 8
    time '24h'

    input:
        tuple val(sample_id), path(bed), path(genome_fai)

    output:
        tuple val(sample_id), path("${sample_id}.intragenic_hybrids.mfe.clusters.bam")

    script:
    """
    #!/usr/bin/env python

    import sys
    import os
    import pysam
    from subprocess import run

    def AddCluster(bam_in, bam_out):

        readcount=0
        writecount=0

        for read in bam_in.fetch(until_eof = True):

            # print(type(read))
            readcount += 1

            # Get read name
            read_name = read.query_name
            tags = read_name.split("_")
            mfe = tags[-1]
            orientation = tags[-2]
            cluster = tags[-3]
            # print(cluster)
            
            read.tags += [("CL", cluster), ("BE", mfe), ("RO", orientation)]
            read.query_name = tags[0] # Remove extra tags from name

            bam_out.write(read)
            writecount += 1

        bam_in.close()
        bam_out.close()

        # Print metrics
        print("Read:", readcount)
        print("Written:", writecount)

    # ==========
    # Run
    # ==========

    # System call to bedtools to convert to BED12 to BAM

    f_in = "$bed"
    f_temp = f_in.replace('bed.gz', 'temp.bam')
    f_out = f_in.replace('bed.gz', 'bam')

    run(f'bedtools bedtobam -bed12 -i {f_in} -g $genome_fai > {f_temp}', shell = True)

    bam_in = pysam.AlignmentFile(f_temp, "rb")
    bam_out = pysam.AlignmentFile(f_out, 'wb', template = bam_in)

    AddCluster(bam_in, bam_out)

    pysam.sort("-o", f_out + '.tmp', f_out)
    os.rename(f_out + '.tmp', f_out)
    pysam.index(f_out)

    print("Completed")
    """

}

def hiclip_header() {

    return """
    -----------------------------------------------------------------
            __        _    ______  _____     _____  _______   
            [  |      (_) .' ___  ||_   _|   |_   _||_   __ \\  
            | |--.   __ / .'   \\_|  | |       | |    | |__) | 
            | .-. | [  || |         | |   _   | |    |  ___/  
            | | | |  | |\\ `.___.'\\ _| |__/ | _| |_  _| |_     
            [___]|__][___]`.____ .'|________||_____||_____|  

    -----------------------------------------------------------------
    """.stripIndent()

}

// Main workflow

params.input='metadata.csv'
params.star_genome_index = '/camp/lab/luscomben/home/users/chakraa2/projects/flora/mouse/ref/STAR_GRCm38_GencodeM24'
params.star_transcript_index = '/camp/lab/luscomben/home/users/chakraa2/projects/flora/mouse/ref/Mm_GencodeM24_rRNA_MT_genes'
params.transcript_fasta_csv = '/camp/lab/luscomben/home/users/chakraa2/projects/flora/mouse/ref/Mm_GencodeM24_rRNA_MT_genes.csv.gz'
params.genome_fai = '/camp/lab/luscomben/home/users/chakraa2/projects/flora/mouse/ref/GRCm38.primary_assembly.genome.fa.fai'
params.transcript_gtf = '/camp/lab/luscomben/home/users/chakraa2/projects/flora/mouse/ref/Mm_GencodeM24_rRNA_MT_genes.gtf.gz'

// Show banner
log.info hiclip_header()

// Create channels for static files
    ch_star_genome_index = Channel.fromPath(params.star_genome_index, checkIfExists: true)
    ch_star_transcript_index = Channel.fromPath(params.star_transcript_index, checkIfExists: true)
    ch_transcript_fasta_csv = Channel.fromPath(params.transcript_fasta_csv, checkIfExists: true)
    ch_genome_fai = Channel.fromPath(params.genome_fai, checkIfExists: true)
    ch_transcript_gtf = Channel.fromPath(params.transcript_gtf, checkIfExists: true)

workflow {

    // Get fastq paths 
    metadata(params.input)

    // Trim
    trim(metadata.out)

    // Filter spliced reads
    premap(trim.out.combine(ch_star_genome_index))
    filtersplicedreads(premap.out)

    // Map chimerias
    mapchimeras(filtersplicedreads.out.combine(ch_star_transcript_index))

    // Remove PCR duplicates
    deduplicate(mapchimeras.out)

    // Extract hybrids
    extracthybrids(deduplicate.out)

    // Get binding energies
    getbindingenergy(extracthybrids.out.combine(ch_transcript_fasta_csv))

    // Get clusters
    clusterhybrids(getbindingenergy.out)

    // Collapse clusters
    collapseclusters(clusterhybrids.out.combine(ch_transcript_gtf))

    // Convert coordinates
    convertcoordinates(clusterhybrids.out.combine(ch_transcript_gtf))

    // Write hybrid BAM
    hybridbedtohybridbam(convertcoordinates.out.combine(ch_genome_fai))

}