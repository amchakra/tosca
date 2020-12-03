# Script to get human custom fasta
# A. M. Chakrabarti
# 29th October 2020

# =========
# Get rRNA
# =========

library(rentrez)
library(stringr)
library(Biostrings)

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

library(rtracklayer)
library(BSgenome.Hsapiens.UCSC.hg38)

# Load annotation
genes.gr <- import.gff3("~/ref/gencode.v35.annotation.gff3.gz")
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

# Get sequences
# reduced.pc.genes.seq <- getSeq(Hsapiens, reduced.pc.genes.gr)
# names(reduced.pc.genes.seq) <- reduced.pc.genes.gr$id

# =========
# Mask snRNAs and rRNA
# =========
library(data.table)

snrna.gr <- genes.gr[genes.gr$gene_type %in% c("rRNA", "snRNA")]
mcols(snrna.gr) <- NULL
snrna.bed <- tempfile(tmpdir = "~/projects/hiclip", fileext = ".bed")
export.bed(snrna.gr, snrna.bed)

ref.fasta <- "~/ref/GRCh38.p13.genome.fa"
masked.fasta <- tempfile(tmpdir = "~/projects/hiclip", fileext = ".fa")
cmd <- paste("bedtools maskfasta -fi", ref.fasta, "-bed", snrna.bed, "-fo", masked.fasta)
system(cmd)

invisible(file.remove(snrna.bed))

reduced.sel.genes.bed <- tempfile(tmpdir = "~/projects/hiclip", fileext = ".bed")
export.bed(reduced.sel.genes.gr, reduced.sel.genes.bed)

cmd <- paste("bedtools getfasta -s -name -fi", masked.fasta, "-bed", reduced.sel.genes.bed, "| pigz > ~/projects/hiclip/human.fa.gz")
system(cmd)

invisible(file.remove(reduced.sel.genes.bed))
invisible(file.remove(masked.fasta))
invisible(file.remove(paste0(masked.fasta, ".fai")))
