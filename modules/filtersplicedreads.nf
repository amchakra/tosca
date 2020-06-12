#!/usr/bin/env nextflow

// Specify DSL2
nextflow.preview.dsl = 2

process filtersplicedreads {

    tag "${sample_id}"
    publishDir "${params.outdir}/filtered", mode: 'copy', overwrite: true

    time '12h'

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