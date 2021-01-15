# Script to get human custom fasta
# A. M. Chakrabarti
# 29th October 2020, updated 28th December 2020

library(rentrez)
library(stringr)
library(Biostrings)
library(rtracklayer)

setwd("/camp/lab/luscomben/home/users/chakraa2/projects/flora/ref/human")

# =========
# Get rRNA
# =========

get_ncbi_sequence <- function(accession) {
  
  fa <- entrez_fetch(db="nucleotide", id=accession, rettype="fasta")
  lines <- str_split(fa, "\n")[[1]]
  id <- gsub(">", "", lines[1])
  seq <- paste0(lines[-1], collapse = "")
  
  # Convert to DNAStringSet
  dss <- DNAStringSet(seq)
  names(dss) <- id
  
  return(dss)
  
}

rDNA.ds <- get_ncbi_sequence("U13369.1") # NCBI rDNA
rRNA_5S.ds <- get_ncbi_sequence("NR_023363.1") # NCBI rRNA 5S

rrna.ds <- c(rDNA.ds, rRNA_5S.ds)
names(rrna.ds) <- c("rDNA", "rRNA5S")

# =========
# Get genes
# =========

# Load annotation
genes.gr <- import.gff3("gencode.v33.annotation.gff3.gz")
genes.gr <- keepStandardChromosomes(genes.gr, pruning.mode = "coarse")
genes.gr <- dropSeqlevels(genes.gr, "chrY", pruning.mode = "coarse")
genes.gr <- genes.gr[genes.gr$type == "gene"]

# Get protein coding genes
biotypes <- unique(genes.gr$gene_type)
pc <- grep("protein_coding|IG_[A-Z]_gene|TR_[A-Z]_gene", biotypes, value = TRUE)
sel.biotypes <- c(pc, "lncRNA", "vault_RNA")
sel.genes.gr <- genes.gr[genes.gr$gene_type %in% sel.biotypes]
sel.genes.gr$ol <- countOverlaps(sel.genes.gr, drop.self = TRUE)

# Collapse overlapping
unique.sel.genes.gr <- sel.genes.gr[sel.genes.gr$ol == 0]
mcols(unique.sel.genes.gr) <- mcols(unique.sel.genes.gr)[, c("gene_name", "gene_id", "gene_type")]

multi.sel.genes.gr <- sel.genes.gr[sel.genes.gr$ol > 0]
multi.reduced.sel.genes.gr <- reduce(multi.sel.genes.gr, with.revmap = TRUE, min.gapwidth = 0) # gapwidth so doesn't merge immediately juxtaposed ranges

# Stitch together names
revmap <- mcols(multi.reduced.sel.genes.gr)$revmap
multi.sel.genes.gr.grl <- relist(multi.sel.genes.gr[unlist(revmap)], revmap)

multi.reduced.sel.genes.gr$gene_name <- sapply(multi.sel.genes.gr.grl, function(x) paste0(sort(x$gene_name), collapse = "_"))
multi.reduced.sel.genes.gr$gene_id <- sapply(multi.sel.genes.gr.grl, function(x) paste0(sort(x$gene_id), collapse = "_"))
multi.reduced.sel.genes.gr$gene_type <- sapply(multi.sel.genes.gr.grl, function(x) paste0(sort(unique(x$gene_type)), collapse = "-"))
multi.reduced.sel.genes.gr$revmap <- NULL

reduced.sel.genes.gr <- sort(c(unique.sel.genes.gr, multi.reduced.sel.genes.gr))
reduced.sel.genes.gr$name <- paste0(reduced.sel.genes.gr$gene_name, ":", reduced.sel.genes.gr$gene_id)

# =========
# Mask snRNAs and rRNA
# =========

snrna.gr <- genes.gr[genes.gr$gene_type %in% c("rRNA", "snRNA")]
mcols(snrna.gr) <- NULL
snrna.bed <- tempfile(tmpdir = ".", fileext = ".bed")
export.bed(snrna.gr, snrna.bed)

ref.fasta <- "GRCh38.primary_assembly.genome.fa"
masked.fasta <- tempfile(tmpdir = ".", fileext = ".fa")
cmd <- paste("bedtools maskfasta -fi", ref.fasta, "-bed", snrna.bed, "-fo", masked.fasta)
message(cmd)
system(cmd)

invisible(file.remove(snrna.bed))

reduced.sel.genes.bed <- tempfile(tmpdir = ".", fileext = ".bed")
export.bed(reduced.sel.genes.gr, reduced.sel.genes.bed)

cmd <- paste("bedtools getfasta -s -name -fi", masked.fasta, "-bed", reduced.sel.genes.bed, "| pigz > human.fa.gz")
message(cmd)
system(cmd)

invisible(file.remove(reduced.sel.genes.bed))
invisible(file.remove(masked.fasta))
invisible(file.remove(paste0(masked.fasta, ".fai")))

# Add rRNA
writeXStringSet(rrna.ds, "human.fa.gz", append = TRUE, compress = TRUE)

# =========
# Create transcript gtf for coordinate conversion
# =========

rrna.gr <- GRanges(seqnames = names(rrna.ds),
                   ranges = IRanges(start = 1, width = width(rrna.ds)),
                   strand = "+",
                   fasta_id = names(rrna.ds))

reduced.sel.genes.gr$fasta_id <- reduced.sel.genes.gr$name
reduced.sel.genes.gr$gene_name <- NULL
reduced.sel.genes.gr$gene_id <- NULL
reduced.sel.genes.gr$gene_type <- NULL
reduced.sel.genes.gr$name <- NULL

transcript.gtf <- c(rrna.gr, reduced.sel.genes.gr)
export.gff2(transcript.gtf, "human.gtf")
system("pigz human.gtf")