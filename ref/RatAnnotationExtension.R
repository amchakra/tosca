library(rtracklayer)
library(GenomicFeatures)
library(data.table)
library(ggplot2)
library(scales)
library(cowplot)
library(parallel)
library(tictoc)

setwd("~/Ule/flora/rat/quantseq")

gtf <- import.gff2("~/Dropbox (The Francis Crick)/rna_structure/ref/rat/Rattus_norvegicus.Rnor_6.0.100.chr.gtf.gz")
genes.gr <- gtf[gtf$type == "gene"]
seqlevelsStyle(genes.gr) <- "UCSC"
pc.genes.gr <- genes.gr[genes.gr$gene_biotype == "protein_coding"]

# bedtools merge -s -c 4,5,6 -o distinct,sum,distinct -i polyAclusters.bed > /Users/chakraa2/Dropbox\ \(The\ Francis\ Crick\)/rna_structure/ref/rat/polyAclusters.merged.bed
pas.bed <- import.bed("~/Dropbox (The Francis Crick)/rna_structure/ref/rat/polyAclusters.merged.bed")
pas.gr <- resize(pas.bed, width = 1, fix = "end")

follow.ol <- follow(pc.genes.gr, pc.genes.gr)


# Get gene upstream of each PAS
follow.ol <- follow(pas.gr, pc.genes.gr)
table(is.na(follow.ol))
pas.gr <- pas.gr[!is.na(follow.ol)]

follow.ol <- follow(pas.gr, pc.genes.gr)
pas.genes.gr <- pc.genes.gr[follow.ol]
distance.ol <- distance(pas.gr, pas.genes.gr)

ggplot(data.table(dist = distance.ol), aes(x = dist)) +
  geom_density() +
  geom_vline(xintercept = 1e4, linetype = "dashed", colour = "red") +
  labs(x = "Distance from Quantseq site to annotated gene",
       y = "Density") +
  scale_x_log10(label = comma, breaks = c(10, 100, 1000, 1e4, 1e5, 1e6, 1e7)) +
  theme_cowplot()

# Get all PAS that are 10 kb downstream of the gene
pas.ol <- findOverlaps(resize(pc.genes.gr, fix = "start", width = width(pc.genes.gr) + 10000), pas.gr)
pas.genes.gr <- pc.genes.gr[queryHits(pas.ol)]
pas.genes.gr$pas <- end(pas.gr[subjectHits(pas.ol)])
pas.genes.gr$score <- pas.gr[subjectHits(pas.ol)]$score
pas.genes.gr$distance <- distance(pas.genes.gr, pas.gr[subjectHits(pas.ol)])

# Get distance to next gene
# Need to remove ones where it pA site is actually in next gene
follow.ol <- precede(pas.genes.gr, pc.genes.gr)
no_following.pas.genes.gr <- pas.genes.gr[is.na(follow.ol)]
no_following.pas.genes.gr$distance_to_next_gene <- Inf

following.pas.genes.gr <- pas.genes.gr[!is.na(follow.ol)]
follow.ol <- precede(following.pas.genes.gr, pc.genes.gr)
following.pas.genes.gr$distance_to_next_gene <- distance(following.pas.genes.gr, pc.genes.gr[follow.ol])

pas.genes.gr <- c(no_following.pas.genes.gr, following.pas.genes.gr)

# Remove
pas.genes.gr <- pas.genes.gr[pas.genes.gr$distance < pas.genes.gr$distance_to_next_gene]

# Adjust gene ends
# Only change if end is less than the PAS (+ve strand)
pos.pas.genes.gr <- pas.genes.gr[strand(pas.genes.gr) == "+" & end(pas.genes.gr) < pas.genes.gr$pas]
end(pos.pas.genes.gr) <- pos.pas.genes.gr$pas
neg.pas.genes.gr <- pas.genes.gr[strand(pas.genes.gr) == "-" & start(pas.genes.gr) > pas.genes.gr$pas]
start(neg.pas.genes.gr) <- neg.pas.genes.gr$pas

pas.genes.gr <- c(pos.pas.genes.gr, neg.pas.genes.gr)

pas.genes.grl <- split(pas.genes.gr, pas.genes.gr$gene_id)
tic()
cl <- makeForkCluster(4)
pas.genes.grl <- GRangesList(parLapply(cl = cl, 1:length(pas.genes.grl), function(i) {
  
  gr <- pas.genes.grl[[i]]
  gr <- gr[gr$score > 0.01 * max(gr$score)]
  gr <- gr[width(gr) == max(width(gr))]
  
  return(gr[1])
  
}))
stopCluster(cl)
toc()

stopifnot(all(elementNROWS(pas.genes.grl) == 1))
extended.genes.gr <- unlist(pas.genes.grl)
extended.genes.gr[extended.genes.gr$distance == 0]$score <- 0 # Assign scores of 0 if not actually extended

unextended.genes.gr <- genes.gr[!genes.gr$gene_id %in% extended.genes.gr$gene_id]
all.genes.gr <- sort(c(extended.genes.gr, unextended.genes.gr))

all.genes.bed <- all.genes.gr
all.genes.bed$score <- all.genes.bed$distance
all.genes.bed[is.na(all.genes.bed$score)]$score <- 0
all.genes.bed$name <- paste0(all.genes.bed$gene_name, "_", all.genes.bed$gene_id)

export.bed(all.genes.bed, "~/Dropbox (The Francis Crick)/rna_structure/ref/rat/Rn_Ens100_Q_v3.10k_0.01score.bed.gz")

ggplot(data.table(dist = all.genes.bed$score), aes(x = dist)) +
  # geom_histogram() +
  geom_density() +
  # geom_vline(xintercept = 1e4, linetype = "dashed", colour = "red") +
  labs(x = "Distance annotated gene has been extenced",
       y = "Density") +
  scale_x_log10(label = comma, breaks = c(10, 100, 1000, 1e4, 1e5, 1e6, 1e7)) +
  theme_cowplot()

# ==========
# Now create an additional transcript for extended genes
# ==========
# TxDb <- makeTxDbFromGFF("~/Dropbox (The Francis Crick)/rna_structure/ref/rat/Rattus_norvegicus.Rnor_6.0.100.chr.gtf.gz")
# saveDb(TxDb, "~/Dropbox (The Francis Crick)/rna_structure/ref/rat/Rattus_norvegicus.Rnor_6.0.100.chr.sqlite")

TxDb <- loadDb("~/Dropbox (The Francis Crick)/rna_structure/ref/rat/Rattus_norvegicus.Rnor_6.0.100.chr.sqlite")

utr3.grl <- threeUTRsByTranscript(TxDb, use.names = TRUE)

# First identify which transcript to extend
rosetta.dt <- as.data.table(mcols(gtf[gtf$type == "transcript"]))
pc.rosetta.dt <- rosetta.dt[gene_biotype == "protein_coding"][transcript_biotype == "protein_coding"]

seqlevelsStyle(gtf) <- "UCSC"
utr3 <- gtf[gtf$type == "three_prime_utr"]
tx <- gtf[gtf$type == "transcript"]
non.gene <- gtf[!is.na(gtf$transcript_id)]

tx.lengths.dt <- data.table(transcriptLengths(TxDb, with.utr5_len = TRUE, with.cds_len = TRUE, with.utr3_len = TRUE))
tx.lengths.dt <- tx.lengths.dt[tx_name %in% tx[tx$transcript_biotype == "protein_coding"]$transcript_id] # Select pc transcripts
setorder(tx.lengths.dt, gene_id, -tx_len, -utr3_len, -cds_len, -nexon, -utr5_len)
tx.longest.dt <- unique(tx.lengths.dt, by = "gene_id")

# i <- 1923
# i <- which(extended.genes.gr$gene_name == "Mapt")
# i <- which(extended.genes.gr$gene_name == "Glis1")
# i <- which(extended.genes.gr$gene_name == "Sobp")

tic()
cl <- makeForkCluster(4)
extended.tx.list <- parLapply(cl = cl, 1:length(extended.genes.gr), function(i) {
# extended.tx.list <- lapply(1:length(extended.genes.gr), function(i) {
    
  message(i)
  
  # Identify transcript to extend
  epas.gr <- resize(extended.genes.gr[i], width = 1, fix = "end")
  sel.utr3 <- utr3[utr3$gene_id == epas.gr$gene_id]
  
  if(length(sel.utr3) != 0) {
    
    nearest.utr3 <- sel.utr3[nearest(epas.gr, sel.utr3)] # Finds 3' UTR nearest to end of extended gene
    sel.tx <- non.gene[non.gene$transcript_id == nearest.utr3$transcript_id]
    
    # Extend tx
    if(as.character(strand(epas.gr)) == "+") end(sel.tx[sel.tx$type == "transcript"]) <- end(epas.gr)
    if(as.character(strand(epas.gr)) == "-") start(sel.tx[sel.tx$type == "transcript"]) <- start(epas.gr)
    
    # Extend utr3 and last exon
    # Need to pick closest 3' UTR segment in case spliced
    if(as.character(strand(epas.gr)) == "+") {
      
      end(sel.tx[sel.tx$type == "three_prime_utr"][nearest(epas.gr, sel.tx[sel.tx$type == "three_prime_utr"])]) <- end(epas.gr)
      end(sel.tx[sel.tx$type == "exon"][nearest(epas.gr, sel.tx[sel.tx$type == "exon"])]) <- end(epas.gr)
      
    }
    
    if(as.character(strand(epas.gr)) == "-") {
      
      start(sel.tx[sel.tx$type == "three_prime_utr"][nearest(epas.gr, sel.tx[sel.tx$type == "three_prime_utr"])]) <- start(epas.gr)
      start(sel.tx[sel.tx$type == "exon"][nearest(epas.gr, sel.tx[sel.tx$type == "exon"])]) <- start(epas.gr)
      
    } 
    
  }
  
  # Some extended protein coding genes don't have a 3' UTR
  if(length(sel.utr3) == 0) {
    
    sel.tx <- tx[tx$gene_id == epas.gr$gene_id]
    sel.tx <- sel.tx[sel.tx$transcript_id %in% tx.longest.dt$tx_name] # Pick longest pc
    nearest.tx <- sel.tx[nearest(epas.gr, sel.tx)]
    sel.tx <- non.gene[non.gene$transcript_id == nearest.tx$transcript_id]
    
    # Extend tx
    if(as.character(strand(epas.gr)) == "+") end(sel.tx[sel.tx$type == "transcript"]) <- end(epas.gr)
    if(as.character(strand(epas.gr)) == "-") start(sel.tx[sel.tx$type == "transcript"]) <- start(epas.gr)
    
    # Extend last exon
    if(as.character(strand(epas.gr)) == "+") {
      
      end.exon <- end(sel.tx[sel.tx$type == "exon"][nearest(epas.gr, sel.tx[sel.tx$type == "exon"])])
      end(sel.tx[sel.tx$type == "exon"][nearest(epas.gr, sel.tx[sel.tx$type == "exon"])]) <- end(epas.gr)
      
    }
    
    if(as.character(strand(epas.gr)) == "-") {
      
      start.exon <- start(sel.tx[sel.tx$type == "exon"][nearest(epas.gr, sel.tx[sel.tx$type == "exon"])])
      start(sel.tx[sel.tx$type == "exon"][nearest(epas.gr, sel.tx[sel.tx$type == "exon"])]) <- start(epas.gr)
      
    }    
    
    # Need to make 3' UTR entry
    new.utr3 <- sel.tx[sel.tx$type == "exon"]
    new.utr3 <- new.utr3[new.utr3$exon_number == max(new.utr3$exon_number)]
    new.utr3$type <- "three_prime_utr"
    new.utr3$exon_number <- NA
    new.utr3$exon_id <- NA
    new.utr3$exon_version <- NA
    new.utr3$protein_id <- NA
    new.utr3$protein_version <- NA
    new.utr3$tag <- NA
    
    if(as.character(strand(epas.gr)) == "+") {
      
      utr3.start <- end.exon + 1 # one after the exon end
      end(new.utr3) <- end(epas.gr)
      start(new.utr3) <- utr3.start
      sel.tx <- c(sel.tx, new.utr3) # Needs to be at the end
      
    }
    
    if(as.character(strand(epas.gr)) == "-") {
      
      utr3.end <- start.exon - 1 # one after the exon end
      start(new.utr3) <- start(epas.gr)
      end(new.utr3) <- utr3.end
      sel.tx <- c(sel.tx, new.utr3)  
      
    }    
    
  }
  
  # Amend tx id
  sel.tx$transcript_id <- gsub("^ENSRNOT0", "ENSRNOT1", sel.tx$transcript_id)
  sel.tx$transcript_name <-  paste0(sel.tx$transcript_name, "-E")
  sel.tx[!is.na(sel.tx$exon_id)]$exon_id <- gsub("^ENSRNOE0", "ENSRNOE1", sel.tx[!is.na(sel.tx$exon_id)]$exon_id)
  sel.tx[!is.na(sel.tx$protein_id)]$protein_id <- gsub("^ENSRNOP0", "ENSRNOP1", sel.tx[!is.na(sel.tx$protein_id)]$protein_id)
  
  # Add in other transcripts and structures
  sel.tx <- c(sel.tx, non.gene[non.gene$gene_id == epas.gr$gene_id])
  # return(sel.tx)
  return(c(extended.genes.gr[i], sel.tx))
  
})
stopCluster(cl)
extended.tx.gr <- unlist(GRangesList(extended.tx.list))
toc()

# seqlevelsStyle(extended.tx.gr) <- "NCBI"
# export(extended.tx.gr, "~/Dropbox (The Francis Crick)/ehiCLIP/Rn/extended_annotation/Rn_Ens99_Q_extended_only.gtf", format = "gtf")

# This will be a bit out of order, but should be alright as all still appropriately grouped by gene and transcript
other.gtf <- gtf[!gtf$gene_id %in% extended.tx.gr$gene_id]
stopifnot(all(start(extended.tx.gr) <= end(extended.tx.gr)))
# seqlevelsStyle(other.gtf) <- "NCBI"
export(c(extended.tx.gr, other.gtf), "~/Dropbox (The Francis Crick)/rna_structure/ref/rat/Rn_Ens100_Q_v3.10k_0.01score.extended.gtf", format = "gtf")
