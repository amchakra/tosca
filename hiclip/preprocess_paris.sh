#!/bin/bash

# Preprocess PARIS FASTQ
# A. M. Chakrabarti
# 10th May 2020

# DDDNNXXXXNNNNNN TACCCTTCGCTTCACACACAAG/iSp18/GGATCC/iSp18/TACTGAACCGCNNNNNN
# Reverse oriented?

# UMI tools will extract the UMI, so need to remove 9 not 15

source activate ehiclip

mkdir rbc

for i in *.fastq.gz; do
	sbatch -c 8 -t 12:00:00 --wrap="umi_tools extract -p NNNNNNXXXX -I $i -S rbc/${i%%.*}.temp.fastq; \
	cutadapt -j 8 -u 9 -m 16 -o rbc/${i%%.*}.rbc.fastq.gz rbc/${i%%.*}.temp.fastq; \
	rm rbc/${i%%.*}.temp.fastq"
done