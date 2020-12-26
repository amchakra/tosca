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

if len(sys.argv) == 3:

    f_in = sys.argv[1]
    f_temp = f_in.replace('bed.gz', 'temp.bam')
    f_out = f_in.replace('bed.gz', 'bam')
    genome_fai = sys.argv[2]

    print(f_in)
    print(f_temp)
    print(f_out)
    print(genome_fai)

    run(f'bedtools bedtobam -bed12 -i {f_in} -g {genome_fai} > {f_temp}', shell = True)

    bam_in = pysam.AlignmentFile(f_temp, "rb")
    bam_out = pysam.AlignmentFile(f_out, 'wb', template = bam_in)

    AddCluster(bam_in, bam_out)

    pysam.sort("-o", f_out + '.tmp', f_out)
    os.rename(f_out + '.tmp', f_out)
    pysam.index(f_out)

    print("Completed")

else:
    print("python convert_hybridbed_to_bam.py <input_bed> <genome_fai>")