#!/bin/bash

# Preprocess PARIS FASTQ
# A. M. Chakrabarti
# 10th May 2020

# ==========
# PARIS
# ==========

# /5phos/DDDNNXXXXNNNNNNTACCCTTCGCTTCACACACAAG/iSp18/GGATCC/iSp18/TACTGAACCGCNNNNNN, where the DDD indicates non-cytosine bases, XXXX indicates the barcode, N indicates any of the 4 bases, iSp18 is spacer.
# The order of the barcodes are the same as described before (Spitale et al., 2015).
# The first random hexamer region is used to check for PCR duplication, whereas the second random hexamer is the primer region.

# DDDNNXXXXNNNNNN TACCCTTCGCTTCACACACAAG/iSp18/GGATCC/iSp18/TACTGAACCGCNNNNNN

# PARIS uses just the hexamer to collapse? But actually collapses unique sequences before mapping so includes DDDNN
# However iCLIP uses WWW as well, so here will use DDDNN-NNNNNN as the UMI
# So will be NNNNNNXXXXNNNNN

# UMI tools will extract the UMI, so need to remove 4 not 15

source activate ehiclip

if [ ! -d rbc ]; then
	mkdir rbc
fi

for i in *.fastq.gz; do
	sbatch -c 8 -t 12:00:00 --wrap="umi_tools extract -p NNNNNNXXXXNNNNN -I $i -S rbc/${i%%.*}.temp.fastq; \
	cutadapt -j 8 -u 4 -m 16 -o rbc/${i%%.*}.rbc.fastq.gz rbc/${i%%.*}.temp.fastq; \
	rm rbc/${i%%.*}.temp.fastq"
done

# ==========
# HuR
# ==========

# Spitale et al. (2015)

# iCLIP_ddRT_BC1: /5phos/DDDNNAACCNNNNAGATCGGAAGAGCGTCGTGAT/iSp18/GG ATCC/iSp18/TACTGAACCGC

# DDDNNAACCNNNNAGATCGGAAGAGCGTCGTGAT

# FAST-iCLIP just collapses unique sequences before mappint
# However iCLIP uses WWW as well, so here will use DDDNN-NNNN as the UMI
# So will be NNNNXXXXNNNNN

# UMI tools will extract the UMI, so need to remove 4 not 13

source activate ehiclip

if [ ! -d rbc ]; then
	mkdir rbc
fi

for i in *.fastq.gz; do
	sbatch -c 8 -t 12:00:00 --wrap="umi_tools extract -p NNNNXXXXNNNNN -I $i -S rbc/${i%%.*}.temp.fastq; \
	cutadapt -j 8 -u 4 -m 16 -o rbc/${i%%.*}.rbc.fastq.gz rbc/${i%%.*}.temp.fastq; \
	rm rbc/${i%%.*}.temp.fastq"
done