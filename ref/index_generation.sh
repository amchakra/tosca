#!/bin/bash

# Script to create STAR indices
# A. M. Chakrabarti
# 11th June 2020

# ==========
# Rat
# ==========

# Genome
sbatch -c 8 --mem=40G -t 24:00:00 --wrap="STAR --runMode genomeGenerate --runThreadN 8 \
--genomeDir STAR_Rnor6_Ensembl100 \
--genomeFastaFiles Rattus_norvegicus.Rnor_6.0.dna.toplevel.fa \
--sjdbGTFfile Rattus_norvegicus.Rnor_6.0.100.gtf \
--sjdbOverhang 49"

# Transcriptome
sbatch -c 8 --mem=40G -t 24:00:00 --wrap="STAR --runMode genomeGenerate --runThreadN 8 \
--genomeDir STAR_Rn_Ens100Q_rRNA_MT_genes \
--genomeFastaFiles Rn_Ens100Q_rRNA_MT_genes.fa"

sbatch -c 8 --mem=40G -t 24:00:00 --wrap="STAR --runMode genomeGenerate --runThreadN 8 \
--genomeDir STAR_Rn_Ens100Q_rRNA_MT_genes_masked \
--genomeFastaFiles Rn_Ens100Q_rRNA_MT_genes_masked.fa"

# ==========
# Mouse
# ==========

# Genome
sbatch -c 8 --mem=40G -t 24:00:00 --wrap="STAR --runMode genomeGenerate --runThreadN 8 \
--genomeDir STAR_GRCm38_GencodeM24 \
--genomeFastaFiles GRCm38.primary_assembly.genome.fa \
--sjdbGTFfile gencode.vM24.annotation.gtf \
--sjdbOverhang 49"

# Transcriptome
sbatch -c 8 --mem=40G -t 24:00:00 --wrap="STAR --runMode genomeGenerate --runThreadN 8 \
--genomeDir STAR_Mm_GencodeM24_rRNA_MT_genes \
--genomeFastaFiles Mm_GencodeM24_rRNA_MT_genes.fa"

# ==========
# Human
# ==========

# Genome
sbatch -c 8 --mem=40G -t 24:00:00 --wrap="STAR --runMode genomeGenerate --runThreadN 8 \
--genomeDir STAR_GRCh38_GencodeV33 \
--genomeFastaFiles GRCh38.primary_assembly.genome.fa \
--sjdbGTFfile gencode.v33.annotation.gtf \
--sjdbOverhang 49"

# Transcriptome
sbatch -c 8 --mem=40G -t 24:00:00 --wrap="STAR --runMode genomeGenerate --runThreadN 8 \
--genomeDir STAR_Hs_GencodeV33_rRNA_MT_genes \
--genomeFastaFiles Hs_GencodeV33_rRNA_MT_genes.fa"