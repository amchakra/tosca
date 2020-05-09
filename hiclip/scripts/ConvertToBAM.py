# Script to convert hybrid BED12 to BAM with additional tags
# A. M. Chakrabarti
# 9th May 2020

import sys
import os
import pysam

# bam="name.bam"

bam_in = pysam.AlignmentFile(sys.argv[1], "rb")
bam_out = pysam.AlignmentFile(sys.argv[2], 'wb', template = bam_in)

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

if len(sys.argv) == 3:

    bam_in = pysam.AlignmentFile(sys.argv[1], "rb")
    bam_out = pysam.AlignmentFile(sys.argv[2], 'wb', template = bam_in)

    AddCluster(bam_in, bam_out)

    pysam.sort("-o", sys.argv[2] + '.tmp', sys.argv[2])
    os.rename(sys.argv[2] + '.tmp', sys.argv[2])
    pysam.index(sys.argv[2])

    print("Completed")

else:
    print("python ConvertToBAM.py <input_bam> <output_bam>")