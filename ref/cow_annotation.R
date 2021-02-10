# Script to get cow custom fasta
# A. M. Chakrabarti
# 10th February

library(rentrez)
library(stringr)
library(Biostrings)
library(rtracklayer)

setwd("/camp/lab/luscomben/home/users/chakraa2/projects/virus/influenza/ref")

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

rDNA.ds <- get_ncbi_sequence("DQ222453.1") # NCBI rDNA
# NR_046257.1 # 45S
rRNA_5S.ds <- get_ncbi_sequence("XR_003031321.1") # NCBI rRNA 5S predicted (1 of 16)

rrna.ds <- c(rDNA.ds, rRNA_5S.ds)
names(rrna.ds) <- c("rDNA", "rRNA5S")

# =========
# Get genes
# =========


# Load annotation
genes.gr <- import.gff3("Bos_taurus.ARS-UCD1.2.102.gff3.gz")
genes.gr <- keepStandardChromosomes(genes.gr, pruning.mode = "coarse")
genes.gr <- dropSeqlevels(genes.gr, "chrY", pruning.mode = "coarse")
genes.gr <- genes.gr[genes.gr$type == "gene"]

# Get protein coding genes
biotypes <- unique(genes.gr$biotypes)
pc <- grep("protein_coding|IG_[A-Z]_gene|TR_[A-Z]_gene", biotypes, value = TRUE)
sel.biotypes <- c(pc)
sel.genes.gr <- genes.gr[genes.gr$biotype %in% sel.biotypes]
sel.genes.gr$ol <- countOverlaps(sel.genes.gr, drop.self = TRUE)

# Collapse overlapping
unique.sel.genes.gr <- sel.genes.gr[sel.genes.gr$ol == 0]
mcols(unique.sel.genes.gr) <- mcols(unique.sel.genes.gr)[, c("Name", "gene_id", "biotype")]

multi.sel.genes.gr <- sel.genes.gr[sel.genes.gr$ol > 0]
multi.reduced.sel.genes.gr <- reduce(multi.sel.genes.gr, with.revmap = TRUE, min.gapwidth = 0) # gapwidth so doesn't merge immediately juxtaposed ranges

# Stitch together names
revmap <- mcols(multi.reduced.sel.genes.gr)$revmap
multi.sel.genes.gr.grl <- relist(multi.sel.genes.gr[unlist(revmap)], revmap)

multi.reduced.sel.genes.gr$gene_name <- sapply(multi.sel.genes.gr.grl, function(x) paste0(sort(x$Name), collapse = "_"))
multi.reduced.sel.genes.gr$gene_id <- sapply(multi.sel.genes.gr.grl, function(x) paste0(sort(x$gene_id), collapse = "_"))
multi.reduced.sel.genes.gr$gene_type <- sapply(multi.sel.genes.gr.grl, function(x) paste0(sort(unique(x$biotype)), collapse = "-"))
multi.reduced.sel.genes.gr$revmap <- NULL

reduced.sel.genes.gr <- sort(c(unique.sel.genes.gr, multi.reduced.sel.genes.gr))
reduced.sel.genes.gr$name <- paste0(reduced.sel.genes.gr$Name, ":", reduced.sel.genes.gr$gene_id)

# =========
# Get fasta
# =========

ref.fasta <- "Bos_taurus.ARS-UCD1.2.dna_rm.toplevel.fa"
reduced.sel.genes.bed <- tempfile(tmpdir = ".", fileext = ".bed")
export.bed(reduced.sel.genes.gr, reduced.sel.genes.bed)

cmd <- paste("bedtools getfasta -s -name -fi", ref.fasta, "-bed", reduced.sel.genes.bed, "| pigz > cow.fa.gz")
message(cmd)
system(cmd)

invisible(file.remove(reduced.sel.genes.bed))

# Add rRNA
writeXStringSet(rrna.ds, "cow.fa.gz", append = TRUE, compress = TRUE)

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
reduced.sel.genes.gr$Name <- NULL
reduced.sel.genes.gr$biotype <- NULL

transcript.gtf <- c(rrna.gr, reduced.sel.genes.gr)
export.gff2(transcript.gtf, "cow.gtf")
system("pigz cow.gtf")