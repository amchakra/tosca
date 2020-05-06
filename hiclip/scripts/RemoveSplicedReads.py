# Script to filter our reads that map to annotated splice junctions
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

# ==========
# Run
# ==========

if len(sys.argv) == 3:
    bam_in = pysam.AlignmentFile(sys.argv[1], "rb")
    bam_out = pysam.AlignmentFile(sys.argv[2], 'wb', template = bam_in)
    FilterBam(bam_in, bam_out)
else:
    print("python RemoveSplicedReads.py <input_bam> <output_bam>")
