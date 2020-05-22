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

    input:
        tuple val(sample_id), path(reads)

    output:
        tuple val(sample_id), path("${sample_id}.trimmed.fastq.gz")

    shell:
    """
    cutadapt -j 8 --minimum-length 16 -q 10 -a AGATCGGAAGAGC -o "${sample_id}.trimmed.fastq.gz" $reads > "${sample_id}_cutadapt.log"
    """
}

process premap {
    tag "${sample_id}"

    input:
        tuple val(sample_id), path(reads), path(star_genome_index)
    
    output:
        tuple val(sample_id), path("${sample_id}.spliced.bam")

    shell:
    """
 
    """

}

process filtersplicedreads {

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

    f_in = "$crosslinks"
    f_out = f_in.replace('.spliced.bam', '.unspliced.bam')
    FilterBam(f_in, f_out)

    # System call to bedtools to convert to fastq
    fastq_out = f_in.replace('.spliced.bam', '.unspliced.fastq.gz')
    call(f'bedtools bamtofastq -i "{f_out}" -fq /dev/stdout | pigz > "{fastq_out}"')
    """
}

// Main workflow

workflow {

    // Get fastq paths 
    metadata(params.input)

    // Trim
    trim(metadata.out)




}