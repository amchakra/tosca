#!/usr/bin/env python    

# Script to filter blat blast8 output
# A. M. Chakrabarti
# 27th March 2019

import sys
import gzip
import time

# ==========
# Functions
# ==========

# Function to count number of valid reads

def CountReads(blast_in, e_value, max_hits):

    counter = 0
    reads = {}

    with gzip.open(blast_in, mode = 'rt') as blast:
            for line in blast:
                
                counter += 1
                if counter % 1000000 == 0:
                    print(counter)
                
                mapping = line.rstrip('\n').rsplit('\t')

                read = mapping[0]
                evalue = float(mapping[10])
                s_start = int(mapping[8])
                s_end = int(mapping[9])

                if evalue <= e_value and s_start < s_end:
                    if read in reads:
                        reads[read] += 1
                    else:
                        reads[read] = 1
                    
    filtered_reads = {k: v for k,v in reads.items() if v > 1 and v <= max_hits}
    return(filtered_reads)

def FilterBlast(blast_in, blast_out, e_value, filtered_reads):

    counter = 0

    with gzip.open(blast_in, mode = 'rt') as blast:
        with gzip.open(blast_out, mode = 'wt') as blast_out:
            for line in blast:
                
                counter += 1
                if counter % 1000000 == 0:
                    print(counter)
                
                mapping = line.rstrip('\n').rsplit('\t')
                read = mapping[0]
                evalue = float(mapping[10])
                s_start = int(mapping[8])
                s_end = int(mapping[9])
                if evalue <= e_value and s_start < s_end and read in filtered_reads:
                    blast_out.write(line)

# ==========
# Run
# ==========

if len(sys.argv) == 5:
    blast_in = sys.argv[1]
    blast_out = sys.argv[2]
    e_value = float(sys.argv[3])
    max_hits = int(sys.argv[4])

    start = time.time()

    filtered_reads = CountReads(blast_in, e_value, max_hits)
    FilterBlast(blast_in, blast_out, e_value, filtered_reads)

    end = time.time()
    print((end-start)/60)

else:
    print("python FilterBlat.py <input_blast> <output_blast> <e_value> <max_hits>")