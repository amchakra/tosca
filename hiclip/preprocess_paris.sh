#!/bin/bash

# Preprocess PARIS FASTQ
# A. M. Chakrabarti
# 10th May 2020

# DDDNNXXXXNNNNNN TACCCTTCGCTTCACACACAAG/iSp18/GGATCC/iSp18/TACTGAACCGCNNNNNN

# PARIS uses just the hexamer to collapse?
# However iCLIP uses WWW as well, so here will use DDDNN-NNNNNN as the UMI
# So will be NNNNNNXXXXNNNNN

# UMI tools will extract the UMI, so need to remove 4 not 15

source activate ehiclip

mkdir rbc

for i in *.fastq.gz; do
	sbatch -c 8 -t 12:00:00 --wrap="umi_tools extract -p NNNNNNXXXXNNNNN -I $i -S rbc/${i%%.*}.temp.fastq; \
	cutadapt -j 8 -u 9 -m 16 -o rbc/${i%%.*}.rbc.fastq.gz rbc/${i%%.*}.temp.fastq; \
	rm rbc/${i%%.*}.temp.fastq"
done