#!/usr/bin/env python

# Script to filter blat blast8 output
# A. M. Chakrabarti
# 27th March 2019 (Modified by I.A. Iosub: 3rd of March 2025)

import sys
import gzip
import time

# ==========
# Functions
# ==========

# Function to initialize log file
def CreateLogFile(log_file):
    with open(log_file, 'w') as log:
        log.write("Step\tReads\tBLAT mappings\tElapsed time (s)\n")

# Function to log information to the log file
def LogStep(step, reads, blat_mappings, start_time, log_file):
    elapsed_time = round(time.time() - start_time, 2)
    with open(log_file, 'a') as log:
        log.write(f"{step}\t{reads}\t{blat_mappings}\t{elapsed_time}\n")

# Function to count number of valid reads with granular logging
def CountReads(blast_in, e_value, max_hits, log_file):
    counter = 0
    reads = {}
    start_time = time.time()

    # Step 1: Count initial unique reads and total mappings
    initial_mappings = 0
    initial_reads = set()  # To count unique input reads

    with gzip.open(blast_in, mode='rt') as blast:
        for line in blast:
            initial_mappings += 1
            read = line.split('\t')[0]
            initial_reads.add(read)

    # Log initial counts before any filtering
    LogStep("Initial", len(initial_reads), initial_mappings, start_time, log_file)

    # Step 2: Apply e-value and orientation filters
    with gzip.open(blast_in, mode='rt') as blast:
        for line in blast:
            counter += 1
            if counter % 1000000 == 0:
                print(counter)

            mapping = line.rstrip('\n').rsplit('\t')
            read = mapping[0]
            evalue = float(mapping[10])
            s_start = int(mapping[8])
            s_end = int(mapping[9])

            # Filter based on e-value and orientation
            if evalue <= e_value and s_start < s_end:
                if read in reads:
                    reads[read] += 1
                else:
                    reads[read] = 1

    # Log after e-value and orientation filtering
    post_filter_unique_reads = len(reads)
    post_filter_mappings = sum(reads.values())
    LogStep("After e-value and orientation filtering", post_filter_unique_reads, post_filter_mappings, start_time, log_file)

    # Step 3: Apply max_hits filtering
    filtered_reads = {k: v for k, v in reads.items() if v > 1 and v <= max_hits}

    # Log after max hits filtering
    post_max_hits_unique_reads = len(filtered_reads)
    post_max_hits_mappings = sum(filtered_reads.values())
    LogStep("After max hits filtering", post_max_hits_unique_reads, post_max_hits_mappings, start_time, log_file)

    return filtered_reads


# Function to filter BLAT
def FilterBlast(blast_in, blast_out, filtered_reads, log_file):
    total_counter = 0        # To count all lines processed
    retained_mappings = 0    # To count only retained mappings
    start_time = time.time()
    retained_reads = set()

    with gzip.open(blast_in, mode='rt') as blast:
        with gzip.open(blast_out, mode='wt') as blast_out:
            for line in blast:
                total_counter += 1  # Count all lines processed
                if total_counter % 1000000 == 0:
                    print(total_counter)

                mapping = line.rstrip('\n').rsplit('\t')
                read = mapping[0]

                # Check if the read is in filtered_reads
                if read in filtered_reads:
                    blast_out.write(line)
                    retained_reads.add(read)
                    retained_mappings += 1  # Count only retained mappings

# ==========
# Run
# ==========

if len(sys.argv) == 6:
    blast_in = sys.argv[1]
    blast_out = sys.argv[2]
    e_value = float(sys.argv[3])
    max_hits = int(sys.argv[4])
    log_file = sys.argv[5]

    # Create and initialize log file
    CreateLogFile(log_file)

    start = time.time()

    filtered_reads = CountReads(blast_in, e_value, max_hits, log_file)
    FilterBlast(blast_in, blast_out, filtered_reads, log_file)

    end = time.time()
    print(f"Total time: {(end - start) / 60} minutes")

else:
    print("Usage: python filter_blat.py <input_blast> <output_blast> <e_value> <max_hits> <log_file>")