# Script to get rat custom fasta for ehiclipr
# A. M. Chakrabarti
# Last updated: 10th June 2020

# =========
# Get rRNA
# =========

library(rentrez)
library(ShortRead)
library(stringr)

rRNA.ids <- c("NR_033176.2","NR_046238.1","NR_046237.1","NR_046246.1") # NCBI rRNA 5S, 5.8S, 18S, 28S
rRNA.fa <- entrez_fetch(db="nucleotide", id=rRNA.ids, rettype="fasta")

fasta <- tempfile()
writeLines(rRNA.fa, fasta)

rRNA.fa <- readFasta(fasta)

id <- as.character(id(rRNA.fa))
id <- str_extract(id, "\\((.*)\\)")
id <- gsub("\\(|\\)", "", id)

rrna.seq <- sread(rRNA.fa)
names(rrna.seq) <- paste0("rRNA_", id)

# writeFasta(rrna.seq, "~/Ule/ref/Rn_rRNA.fa")
# min(14, (log2(sum(width(sread(readFasta("~/Ule/ref/Rn_rRNA.fa"))))/2) - 1))

file.remove(fasta)

# =========
# Get longest RNA (protein coding and non-coding)
# =========

# library(biomaRt)
library(data.table)
library(AnnotationDbi)
library(ShortRead)
library(GenomicFeatures)
library(BSgenome.Rnorvegicus.UCSC.rn6)
library(pbapply)
# library(ensembldb)

# =========
# Load annotation
# =========

genes.gr <- import.gff2("~/Dropbox (The Francis Crick)/rna_structure/ref/rat/Rn_Ens100_Q_v3.10k_0.01score.extended.gtf.gz")
genes.gr <- genes.gr[genes.gr$type == "gene"]
genes.gr <- keepStandardChromosomes(genes.gr, pruning.mode = "coarse")
seqlevelsStyle(genes.gr) <- "UCSC" # To match Ensembl

# findOverlaps(genes.gr, drop.self = TRUE, drop.redundant = TRUE)

# =========
# Get protein coding genes
# =========
pc.genes.gr <- genes.gr[genes.gr$gene_biotype == "protein_coding"] # needs to be gene_biotype as Ensembl
pc.genes.gr <- dropSeqlevels(pc.genes.gr, "chrM", pruning.mode = "coarse") # will split out mitochondrial separately

# ol <- findOverlaps(pc.genes.gr, drop.self = TRUE, drop.redundant = TRUE)
# length(unique(c(queryHits(ol), subjectHits(ol))))
# 
# unique(pc.genes.gr[sort(unique(c(queryHits(ol), subjectHits(ol))))]$gene_name)

# Get sequences
pc.regions.seq <- getSeq(Rnorvegicus, pc.genes.gr)
names(pc.regions.seq) <- paste0(pc.genes.gr$gene_name, "_", pc.genes.gr$gene_id)

# =========
# Get non coding genes
# =========
nc.genes.gr <- genes.gr[genes.gr$gene_biotype != "protein_coding"]
nc.genes.gr <- dropSeqlevels(nc.genes.gr, "chrM", pruning.mode = "coarse") # will split out mitochondrial separately

# Get sequences
nc.regions.seq <- getSeq(Rnorvegicus, nc.genes.gr)
names(nc.regions.seq) <- paste0(nc.genes.gr$gene_name, "_", nc.genes.gr$gene_id)

# =========
# Get mitochondrial genes
# =========
mt.genes.gr <- genes.gr[seqnames(genes.gr) == "chrM"]

# Get sequences
mt.regions.seq <- getSeq(Rnorvegicus, mt.genes.gr)
names(mt.regions.seq) <- paste0(mt.genes.gr$gene_name, "_", mt.genes.gr$gene_id)

# =========
# Write out fasta
# =========
all.seq <- c(rrna.seq, pc.regions.seq, nc.regions.seq, mt.regions.seq)
# writeFasta(all.seq, "~/Dropbox (The Francis Crick)/rna_structure/ref/rat/Rn_Ens100Q_rRNA_MT_genes.fa")
# system("pigz ~/Dropbox\\ \\(Lab\\)/rna_structure/ref/mouse/Rn_Ens100Q_rRNA_MT_genes.fa")

# Combine for coordinate conversion
seq.gr <- c(pc.genes.gr, nc.genes.gr, mt.genes.gr)
seq.gr$fasta_id <- paste0(seq.gr$gene_name, "_", seq.gr$gene_id)
seq.gr <- sort(seq.gr)
# export.gff2(seq.gr, "~/Dropbox (The Francis Crick)/rna_structure/ref/rat/Rn_Ens100Q_rRNA_MT_genes.gtf.gz")
# system("pigz ~/Dropbox\\ \\(Lab\\)/rna_structure/ref/rat/Rn_Ens100Q_rRNA_MT_genes.gtf")

# For new masked version
seq.gr <- c(pc.genes.gr, nc.genes.gr, mt.genes.gr)
seq.gr$name <- paste0(seq.gr$gene_name, "_", seq.gr$gene_id)
seq.gr$score[is.na(seq.gr$score)] <- -1
seq.gr <- sort(seq.gr)
seqlevelsStyle(seq.gr) <- "NCBI"
export.bed(seq.gr, "~/Dropbox (The Francis Crick)/rna_structure/ref/rat/Rn_Ens100Q_rRNA_MT_genes.bed.gz")

