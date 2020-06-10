# Script to get human custom fasta for ehiclipr
# A. M. Chakrabarti
# Last updated: 10th June 2020

# =========
# Get rRNA
# =========

library(rentrez)
library(ShortRead)
library(stringr)

rRNA.ids <- c("NR_023363.1","NR_003285.3","NR_003287.4","NR_003286.4") # NCBI rRNA 5S, 5.8S, 28S, 18S
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
library(BSgenome.Hsapiens.UCSC.hg38)
library(pbapply)
library(tictoc)
# library(ensembldb)

# =========
# Load annotation
# =========

genes.gr <- import.gff3("~/Dropbox (The Francis Crick)/rna_structure/ref/human/gencode.v34.annotation.gff3.gz")
genes.gr <- genes.gr[genes.gr$type == "gene"]

# findOverlaps(genes.gr, drop.self = TRUE, drop.redundant = TRUE)

# =========
# Get protein coding genes
# =========
pc.genes.gr <- genes.gr[genes.gr$gene_type == "protein_coding"]
pc.genes.gr <- dropSeqlevels(pc.genes.gr, "chrM", pruning.mode = "coarse") # will split out mitochondrial separately

ol <- findOverlaps(pc.genes.gr, drop.self = TRUE, drop.redundant = TRUE)
length(unique(c(queryHits(ol), subjectHits(ol))))

pc.genes.gr[sort(unique(c(queryHits(ol), subjectHits(ol))))]

# Get sequences
pc.regions.seq <- getSeq(Hsapiens, pc.genes.gr)
names(pc.regions.seq) <- paste0(pc.genes.gr$gene_name, "_", pc.genes.gr$gene_id)

# =========
# Get non coding genes
# =========
nc.genes.gr <- genes.gr[genes.gr$gene_type != "protein_coding"]
nc.genes.gr <- dropSeqlevels(nc.genes.gr, "chrM", pruning.mode = "coarse") # will split out mitochondrial separately

# Get sequences
nc.regions.seq <- getSeq(Hsapiens, nc.genes.gr)
names(nc.regions.seq) <- paste0(nc.genes.gr$gene_name, "_", nc.genes.gr$gene_id)

# =========
# Get mitochondrial genes
# =========
mt.genes.gr <- genes.gr[seqnames(genes.gr) == "chrM"]

# Get sequences
mt.regions.seq <- getSeq(Hsapiens, mt.genes.gr)
names(mt.regions.seq) <- paste0(mt.genes.gr$gene_name, "_", mt.genes.gr$gene_id)

# =========
# Write out fasta
# =========
all.seq <- c(rrna.seq, pc.regions.seq, nc.regions.seq, mt.regions.seq)
writeFasta(all.seq, "~/Dropbox (The Francis Crick)/rna_structure/ref/human/Hs_GencodeV34_rRNA_MT_genes.fa")
system("pigz ~/Dropbox\\ \\(Lab\\)/rna_structure/ref/human/Hs_GencodeV34_rRNA_MT_genes.fa")

# Combine for coordinate conversion
seq.gr <- c(pc.genes.gr, nc.genes.gr, mt.genes.gr)
seq.gr$fasta_id <- paste0(seq.gr$gene_name, "_", seq.gr$gene_id)
seq.gr <- sort(seq.gr)
export.gff2(seq.gr, "~/Dropbox (The Francis Crick)/rna_structure/ref/human/Hs_GencodeV34_rRNA_MT_genes.gtf")
system("pigz ~/Dropbox\\ \\(Lab\\)/rna_structure/ref/human/Hs_GencodeV34_rRNA_MT_genes.gtf")

# ======================================================
# Now annotate with regions of longest transcript
# ======================================================

# TxDb <- makeTxDbFromGFF("~/Dropbox (The Francis Crick)/rna_structure/ref/human/gencode.v34.annotation.gtf.gz")
# saveDb(TxDb, "~/Dropbox (The Francis Crick)/rna_structure/ref/human/gencode.v34.sqlite")
tic()
TxDb <- loadDb("~/Dropbox (The Francis Crick)/rna_structure/ref/human/gencode.v34.sqlite")
TxDb <- keepStandardChromosomes(TxDb, pruning.mode = "coarse")

# Get transcript lengths and add to annotation data.table
txlengths.dt <- data.table(transcriptLengths(TxDb, with.cds_len = TRUE, with.utr5_len = TRUE, with.utr3_len = TRUE))
setnames(txlengths.dt, c("tx_name", "gene_id"), c("transcript_id", "gene_id"))
txlengths.dt[, tx_id := NULL]

# Need to load GTF to get extra metadata columns
gtf <- import.gff2("~/Dropbox (The Francis Crick)/rna_structure/ref/human/gencode.v34.annotation.gtf.gz")
tx <- gtf[gtf$type == "transcript"]
gencode.dt <- as.data.table(mcols(tx))[, .(gene_id, gene_type, gene_name, level, transcript_id, transcript_type)]
stopifnot(nrow(txlengths.dt) == nrow(gencode.dt))
gencode.dt <- merge(txlengths.dt, gencode.dt, by = c("gene_id", "transcript_id"))
setnames(gencode.dt, c("gene_id", "transcript_id"), c("ensembl_gene_id", "ensembl_transcript_id"))

# Define protein coding and non-coding (according to transcript biotype)
tx_biotype <- unique(gencode.dt$transcript_type)
pc <- c("protein_coding", "IG_C_gene", "IG_D_gene", "IG_J_gene", "IG_LV_gene", "IG_V_gene", "TR_C_gene", "TR_D_gene", "TR_J_gene", "TR_V_gene")
rRNA <- c("rRNA", "Mt_rRNA")
tRNA <- c("Mt_tRNA")
miRNA <- "miRNA"
lncRNA <- c("lncRNA", "IG_C_pseudogene", "IG_D_pseudogene", "IG_pseudogene", "IG_V_pseudogene", "TR_J_pseudogene", "polymorphic_pseudogene", "processed_pseudogene", "processed_transcript", "pseudogene", "transcribed_processed_pseudogene", "transcribed_unitary_pseudogene", "transcribed_unprocessed_pseudogene", "translated_processed_pseudogene", "translated_unprocessed_pseudogene", "unitary_pseudogene", "unprocessed_pseudogene", "TR_V_pseudogene", "rRNA_pseudogene", "IG_J_pseudogene", "retained_intron")
ncRNA <- c("misc_RNA", "ribozyme", "scaRNA", "scRNA", "snoRNA", "snRNA", "sRNA", "vaultRNA")
other <- c("non_stop_decay", "nonsense_mediated_decay", "TEC")

stopifnot(all(tx_biotype %in% c(pc, rRNA, tRNA, miRNA, miRNA, lncRNA, ncRNA, other)))

# ==========
# Protein coding regions
# ==========

# Select longest pc
longest.pc.dt <- gencode.dt[transcript_type %in% pc] # selects pc
longest.pc.dt <- longest.pc.dt[, longest := max(tx_len), by = ensembl_gene_id][tx_len == longest] # selects longest

setorder(longest.pc.dt, -nexon, -utr3_len, -cds_len, -utr5_len)
setkey(longest.pc.dt, ensembl_gene_id)
longest.pc.dt <- unique(longest.pc.dt, by = "ensembl_gene_id") # removes duplicates prioritising most exons then longest utr3 etc.
stopifnot(anyDuplicated(longest.pc.dt[, ensembl_gene_id]) == 0)

# Create GRangesLists of transcript regions
utr5_tx.grl <- fiveUTRsByTranscript(TxDb, use.names = TRUE)
cds_tx.grl <- cdsBy(TxDb, by=c("tx"), use.names = TRUE)
intron_tx.grl <- intronsByTranscript(TxDb, use.names = TRUE)
utr3_tx.grl <- threeUTRsByTranscript(TxDb, use.names = TRUE)

# Function to create GRanges for regions of longest transcripts from GRangesLists above and add metadata
regionGR <- function(grl, region, coding = c(longest.nc.dt, longest.pc.dt)) {
  grl <- grl[names(grl) %in% coding$ensembl_transcript_id] # select transcripts (names = transcript names as GRangesLists created by "tx")
  gr <- unlist(grl) # unlist GRangesList to GRanges
  mcols(gr) <- NULL # reset metadata
  mcols(gr)$tx_name <- names(gr) # transfer transcript name from GRanges name
  mcols(gr)$ensembl_gene_id <- coding$ensembl_gene_id[match(gr$tx_name, coding$ensembl_transcript_id)] # match transcript names and add gene name
  mcols(gr)$gene_name <- coding$gene_name[match(gr$tx_name, coding$ensembl_transcript_id)] # match transcript names and add gene name
  mcols(gr)$biotype <- coding$transcript_type[match(gr$tx_name, coding$ensembl_transcript_id)]
  mcols(gr)$region <- region # add region name
  names(gr) <- NULL
  return(gr)
}

utr5.gr <- regionGR(utr5_tx.grl, region = "UTR5", coding = longest.pc.dt)
cds.gr <- regionGR(cds_tx.grl, region = "CDS", coding = longest.pc.dt)
introns.gr <- regionGR(intron_tx.grl, region = "INTRON", coding = longest.pc.dt)
utr3.gr <- regionGR(utr3_tx.grl, region = "UTR3", coding = longest.pc.dt)

pc.regions.gr <- c(utr5.gr, cds.gr, introns.gr, utr3.gr)


# Need to do it this way for regions that don't overlap the longest transcript
exons_tx.grl <- exonsBy(TxDb, by = "tx", use.names = TRUE)
# exons_g.grl <- exonsBy(TxDb, by = "gene")

cl <- makeForkCluster(4)
pc.genes.regions.grl <- pblapply(cl = cl, seq_along(pc.genes.gr), function(i) {
  
  # message(i)
  gr <- pc.genes.gr[i]
  gene <- gr$gene_id
  gene_name <- gr$gene_name
  biotype <- gr$gene_type
  
  gr <- unlist(tile(gr, width = 1))
  gr$annot <- as.character(NA) 
  annotated.gr <- GRanges()
  gr_length <- length(gr)
  
  # First annotate those that overlap longest transcript
  sel.pc.regions.gr <- pc.regions.gr[pc.regions.gr$ensembl_gene_id == gene] # Get matching gene
  
  ol <- findOverlaps(gr, sel.pc.regions.gr)
  
  # gr$ensembl_transcript_id <- as.character(NA)
  # gr$ensembl_gene_id <- as.character(NA)
  # gr$gene_name <- as.character(NA)
  # gr$biotype <- as.character(NA)
  # gr$region <- as.character(NA)
  # 
  # gr[queryHits(ol)]$ensembl_transcript_id <- sel.pc.regions.gr[subjectHits(ol)]$tx_name
  # gr[queryHits(ol)]$ensembl_gene_id <- sel.pc.regions.gr[subjectHits(ol)]$ensembl_gene_id
  # gr[queryHits(ol)]$gene_name <- sel.pc.regions.gr[subjectHits(ol)]$gene_name
  # gr[queryHits(ol)]$biotype <- sel.pc.regions.gr[subjectHits(ol)]$biotype
  # gr[queryHits(ol)]$region <- sel.pc.regions.gr[subjectHits(ol)]$region
  
  # Some protein coding genes don't have a protein coding transcript, e.g. ENSMUST00000155020.1
  if(length(ol) != 0) {
    
    gr[queryHits(ol)]$annot <- paste0(sel.pc.regions.gr[subjectHits(ol)]$tx_name, "|",
                                      sel.pc.regions.gr[subjectHits(ol)]$ensembl_gene_id, "|",
                                      sel.pc.regions.gr[subjectHits(ol)]$gene_name, "|",
                                      sel.pc.regions.gr[subjectHits(ol)]$biotype, "|",
                                      sel.pc.regions.gr[subjectHits(ol)]$region)
    
    # annotated.gr <- gr[!is.na(gr$ensembl_transcript_id)]
    # gr <- gr[is.na(gr$ensembl_transcript_id)]
    annotated.pc.gr <- gr[!is.na(gr$annot)]
    annotated.gr <- c(annotated.gr, annotated.pc.gr)
    gr <- gr[is.na(gr$annot)]
    
  }
  # First annotate the rest based on a heirarchy
  
  # UTR3
  sel.utr3.gr <- unlist(utr3_tx.grl[names(utr3_tx.grl) %in% gencode.dt[ensembl_gene_id == gene & transcript_type %in% pc]$ensembl_transcript_id])
  ol <- findOverlaps(gr, sel.utr3.gr)
  if(length(ol) != 0) {
    
    gr[queryHits(ol)]$annot <- paste0(names(sel.utr3.gr[subjectHits(ol)]), "|",
                                      gene, "|",
                                      gene_name, "|",
                                      biotype, "|",
                                      "OTHER_UTR3")
    
    annotated.utr3.gr <- gr[!is.na(gr$annot)]
    annotated.gr <- c(annotated.gr, annotated.utr3.gr)
    gr <- gr[is.na(gr$annot)]
    
  }
  
  # CDS
  sel.cds.gr <- unlist(cds_tx.grl[names(cds_tx.grl) %in% gencode.dt[ensembl_gene_id == gene & transcript_type %in% pc]$ensembl_transcript_id])
  ol <- findOverlaps(gr, sel.cds.gr) 
  
  if(length(ol) != 0) {
    
    gr[queryHits(ol)]$annot <- paste0(names(sel.cds.gr[subjectHits(ol)]), "|",
                                      gene, "|",
                                      gene_name, "|",
                                      biotype, "|",
                                      "OTHER_CDS")
    
    annotated.cds.gr <- gr[!is.na(gr$annot)]
    annotated.gr <- c(annotated.gr, annotated.cds.gr)
    gr <- gr[is.na(gr$annot)]
    
  }
  
  # UTR5
  sel.utr5.gr <- unlist(utr5_tx.grl[names(utr5_tx.grl) %in% gencode.dt[ensembl_gene_id == gene & transcript_type %in% pc]$ensembl_transcript_id])
  ol <- findOverlaps(gr, sel.utr5.gr) 
  
  if(length(ol) != 0) {
    
    gr[queryHits(ol)]$annot <- paste0(names(sel.utr5.gr[subjectHits(ol)]), "|",
                                      gene, "|",
                                      gene_name, "|",
                                      biotype, "|",
                                      "OTHER_UTR5")
    
    annotated.utr5.gr <- gr[!is.na(gr$annot)]
    annotated.gr <- c(annotated.gr, annotated.utr5.gr)    
    gr <- gr[is.na(gr$annot)]
    
  } 
  
  # Introns
  sel.intron.gr <- unlist(intron_tx.grl[names(intron_tx.grl) %in% gencode.dt[ensembl_gene_id == gene & transcript_type %in% pc]$ensembl_transcript_id])
  ol <- findOverlaps(gr, sel.intron.gr) 
  
  if(length(ol) != 0) {
    
    gr[queryHits(ol)]$annot <- paste0(names(sel.intron.gr[subjectHits(ol)]), "|",
                                      gene, "|",
                                      gene_name, "|",
                                      biotype, "|",
                                      "OTHER_INTRON")
    
    annotated.intron.gr <- gr[!is.na(gr$annot)]
    annotated.gr <- c(annotated.gr, annotated.intron.gr)    
    gr <- gr[is.na(gr$annot)]
    
  }   
  
  # Exons (i.e. non-protein coding)
  sel.exons.gr <- unlist(exons_tx.grl[names(exons_tx.grl) %in% gencode.dt[ensembl_gene_id == gene]$ensembl_transcript_id])
  ol <- findOverlaps(gr, sel.exons.gr) 
  
  if(length(ol) != 0) {
    
    gr[queryHits(ol)]$annot <- paste0(names(sel.exons.gr[subjectHits(ol)]), "|",
                                      gene, "|",
                                      gene_name, "|",
                                      biotype, "|",
                                      "OTHER_EXON")
    
    annotated.exon.gr <- gr[!is.na(gr$annot)]
    annotated.gr <- c(annotated.gr, annotated.exon.gr)   
    gr <- gr[is.na(gr$annot)]
    
  } 
  
  # Rest are introns
  
  if(length(gr) != 0) {
    
    gr$annot <- paste0("ENSMUST_ALL", "|",
                       gene, "|",
                       gene_name, "|",
                       biotype, "|",
                       "OTHER_INTRON")
    
    annotated.gr <- c(annotated.gr, gr)
    
  }
  
  annotated.gr <- sort(annotated.gr) # sort could be removed - just to check coord conversion
  
  stopifnot(length(annotated.gr) == gr_length)
  
  # Adjust coordinate depending on strand
  if(unique(strand(annotated.gr)) == "+") {
    
    start <- min(start(annotated.gr))
    tc.gr <- GRanges(seqnames = Rle(gene),
                     ranges = IRanges(start = start(annotated.gr) - start + 1,
                                      width = 1),
                     strand = Rle("+"),
                     annot = annotated.gr$annot)
    
    
  } else if(unique(strand(annotated.gr)) == "-") {
    
    start <- max(start(annotated.gr))
    tc.gr <- GRanges(seqnames = Rle(gene),
                     ranges = IRanges(start = abs(start(annotated.gr) - start - 1),
                                      width = 1),
                     strand = Rle("+"),
                     annot = annotated.gr$annot)
    
  }
  
  # Now need to collapse down (Not strictly necessary)
  annotated.grl <- split(tc.gr, tc.gr$annot)
  annotated.grl <- lapply(annotated.grl, reduce)
  reduced.annotated.gr <- sort(unlist(GRangesList(annotated.grl)))
  reduced.annotated.gr$annotation <- names(reduced.annotated.gr)
  names(reduced.annotated.gr) <- NULL
  
  stopifnot(length(reduce(reduced.annotated.gr)) == 1)
  stopifnot(width(reduce(reduced.annotated.gr)) == gr_length)
  
  return(reduced.annotated.gr)
  
})
stopCluster(cl)

# saveRDS(pc.genes.regions.grl, "~/Dropbox (The Francis Crick)/rna_structure/ref/mouse/gencode_vM24.pc.genes.regions.grl.rds")

pc.genes.regions.gr <- unlist(GRangesList(pc.genes.regions.grl))
# saveRDS(pc.genes.regions.gr, "~/Dropbox (The Francis Crick)/rna_structure/ref/mouse/gencode_vM24.pc.genes.regions.gr.rds")

# ==========
# Non coding regions
# ==========

# rRNA

rrna.gr <- GRanges(seqnames = names(rrna.seq),
                   ranges = IRanges(start = 1, width = width(rrna.seq)),
                   strand = "+",
                   annotation = paste0(names(rrna.seq), "|", names(rrna.seq), "|", sapply(strsplit(names(rrna.seq), "_"), "[[", 2), "|rRNA", "|rRNA"))

# mtRNA

cl <- makeForkCluster(4)
mt.genes.region.grl <- pblapply(cl = cl, seq_along(mt.genes.gr), function(i) {
  
  # message(i)
  gr <- mt.genes.gr[i]
  gene <- gr$gene_id
  gene_name <- gr$gene_name
  biotype <- gr$gene_type
  
  gr <- unlist(tile(gr, width = 1))
  gr$annot <- as.character(NA) 
  annotated.gr <- GRanges()
  gr_length <- length(gr)
  
  # Annotate exons
  sel.exons.gr <- unlist(exons_tx.grl[names(exons_tx.grl) %in% gencode.dt[ensembl_gene_id == gene]$ensembl_transcript_id])
  ol <- findOverlaps(gr, sel.exons.gr) 
  
  if(length(ol) != 0) {
    
    gr[queryHits(ol)]$annot <- paste0(names(sel.exons.gr[subjectHits(ol)]), "|",
                                      gene, "|",
                                      gene_name, "|",
                                      biotype, "|",
                                      "EXON")
    
    annotated.exon.gr <- gr[!is.na(gr$annot)]
    annotated.gr <- c(annotated.gr, annotated.exon.gr)   
    gr <- gr[is.na(gr$annot)]
    
  } 
  
  # Annotate introns
  sel.introns.gr <- unlist(intron_tx.grl[names(intron_tx.grl) %in% gencode.dt[ensembl_gene_id == gene]$ensembl_transcript_id])
  ol <- findOverlaps(gr, sel.introns.gr) 
  
  if(length(ol) != 0) {
    
    gr[queryHits(ol)]$annot <- paste0(names(sel.introns.gr[subjectHits(ol)]), "|",
                                      gene, "|",
                                      gene_name, "|",
                                      biotype, "|",
                                      "INTRON")
    
    annotated.exon.gr <- gr[!is.na(gr$annot)]
    annotated.gr <- c(annotated.gr, annotated.exon.gr)   
    gr <- gr[is.na(gr$annot)]
    
  } 
  
  # Rest are other introns
  if(length(gr) != 0) {
    
    gr$annot <- paste0("ENSMUST_ALL", "|",
                       gene, "|",
                       gene_name, "|",
                       biotype, "|",
                       "OTHER_INTRON")
    
    annotated.gr <- c(annotated.gr, gr)
    
  }
  
  annotated.gr <- sort(annotated.gr) # sort could be removed - just to check coord conversion
  
  stopifnot(length(annotated.gr) == gr_length)
  
  # Adjust coordinate depending on strand
  if(unique(strand(annotated.gr)) == "+") {
    
    start <- min(start(annotated.gr))
    tc.gr <- GRanges(seqnames = Rle(gene),
                     ranges = IRanges(start = start(annotated.gr) - start + 1,
                                      width = 1),
                     strand = Rle("+"),
                     annot = annotated.gr$annot)
    
    
  } else if(unique(strand(annotated.gr)) == "-") {
    
    start <- max(start(annotated.gr))
    tc.gr <- GRanges(seqnames = Rle(gene),
                     ranges = IRanges(start = abs(start(annotated.gr) - start - 1),
                                      width = 1),
                     strand = Rle("+"),
                     annot = annotated.gr$annot)
    
  }
  
  # Now need to collapse down (Not strictly necessary)
  annotated.grl <- split(tc.gr, tc.gr$annot)
  annotated.grl <- lapply(annotated.grl, reduce)
  reduced.annotated.gr <- sort(unlist(GRangesList(annotated.grl)))
  reduced.annotated.gr$annotation <- names(reduced.annotated.gr)
  names(reduced.annotated.gr) <- NULL
  
  stopifnot(length(reduce(reduced.annotated.gr)) == 1)
  stopifnot(width(reduce(reduced.annotated.gr)) == gr_length)
  
  return(reduced.annotated.gr)
  
})
stopCluster(cl)

mt.genes.regions.gr <- unlist(GRangesList(mt.genes.region.grl))
mt.genes.regions.gr$annotation <- gsub("\\|EXON$", "\\|MT_EXON", mt.genes.regions.gr$annotation)
mt.genes.regions.gr$annotation <- gsub("\\|INTRON$", "\\|MT_INTRON", mt.genes.regions.gr$annotation)
mt.genes.regions.gr$annotation <- gsub("\\|OTHER_INTRON$", "\\|MT_OTHER_INTRON", mt.genes.regions.gr$annotation)

# ncRNA

cl <- makeForkCluster(4)
nc.genes.region.grl <- pblapply(cl = cl, seq_along(nc.genes.gr), function(i) {
  
  # message(i)
  gr <- nc.genes.gr[i]
  gene <- gr$gene_id
  gene_name <- gr$gene_name
  biotype <- gr$gene_type
  
  gr <- unlist(tile(gr, width = 1))
  gr$annot <- as.character(NA) 
  annotated.gr <- GRanges()
  gr_length <- length(gr)
  
  # Annotate exons
  sel.exons.gr <- unlist(exons_tx.grl[names(exons_tx.grl) %in% gencode.dt[ensembl_gene_id == gene]$ensembl_transcript_id])
  ol <- findOverlaps(gr, sel.exons.gr) 
  
  if(length(ol) != 0) {
    
    gr[queryHits(ol)]$annot <- paste0(names(sel.exons.gr[subjectHits(ol)]), "|",
                                      gene, "|",
                                      gene_name, "|",
                                      biotype, "|",
                                      "EXON")
    
    annotated.exon.gr <- gr[!is.na(gr$annot)]
    annotated.gr <- c(annotated.gr, annotated.exon.gr)   
    gr <- gr[is.na(gr$annot)]
    
  } 
  
  # Annotate introns
  sel.introns.gr <- unlist(intron_tx.grl[names(intron_tx.grl) %in% gencode.dt[ensembl_gene_id == gene]$ensembl_transcript_id])
  ol <- findOverlaps(gr, sel.introns.gr) 
  
  if(length(ol) != 0) {
    
    gr[queryHits(ol)]$annot <- paste0(names(sel.introns.gr[subjectHits(ol)]), "|",
                                      gene, "|",
                                      gene_name, "|",
                                      biotype, "|",
                                      "INTRON")
    
    annotated.exon.gr <- gr[!is.na(gr$annot)]
    annotated.gr <- c(annotated.gr, annotated.exon.gr)   
    gr <- gr[is.na(gr$annot)]
    
  } 
  
  # Rest are introns
  if(length(gr) != 0) {
    
    gr$annot <- paste0("ENSMUST_ALL", "|",
                       gene, "|",
                       gene_name, "|",
                       biotype, "|",
                       "OTHER_INTRON")
    
    annotated.gr <- c(annotated.gr, gr)
    
  }
  
  annotated.gr <- sort(annotated.gr) # sort could be removed - just to check coord conversion
  
  stopifnot(length(annotated.gr) == gr_length)
  
  # Adjust coordinate depending on strand
  if(unique(strand(annotated.gr)) == "+") {
    
    start <- min(start(annotated.gr))
    tc.gr <- GRanges(seqnames = Rle(gene),
                     ranges = IRanges(start = start(annotated.gr) - start + 1,
                                      width = 1),
                     strand = Rle("+"),
                     annot = annotated.gr$annot)
    
    
  } else if(unique(strand(annotated.gr)) == "-") {
    
    start <- max(start(annotated.gr))
    tc.gr <- GRanges(seqnames = Rle(gene),
                     ranges = IRanges(start = abs(start(annotated.gr) - start - 1),
                                      width = 1),
                     strand = Rle("+"),
                     annot = annotated.gr$annot)
    
  }
  
  # Now need to collapse down (Not strictly necessary)
  annotated.grl <- split(tc.gr, tc.gr$annot)
  annotated.grl <- lapply(annotated.grl, reduce)
  reduced.annotated.gr <- sort(unlist(GRangesList(annotated.grl)))
  reduced.annotated.gr$annotation <- names(reduced.annotated.gr)
  names(reduced.annotated.gr) <- NULL
  
  stopifnot(length(reduce(reduced.annotated.gr)) == 1)
  stopifnot(width(reduce(reduced.annotated.gr)) == gr_length)
  
  return(reduced.annotated.gr)
  
})
stopCluster(cl)

nc.genes.region.gr <- unlist(GRangesList(nc.genes.region.grl))
nc.genes.region.gr$annotation <- gsub("\\|EXON$", "\\|NC_EXON", nc.genes.region.gr$annotation)
nc.genes.region.gr$annotation <- gsub("\\|INTRON$", "\\|NC_INTRON", nc.genes.region.gr$annotation)
nc.genes.region.gr$annotation <- gsub("\\|OTHER_INTRON$", "\\|NC_OTHER_INTRON", nc.genes.region.gr$annotation)

# Merge together
combined.seqlevels <- unique(c(seqlevels(rrna.gr), seqlevels(mt.genes.regions.gr), seqlevels(nc.genes.region.gr), seqlevels(pc.genes.regions.gr)))

seqlevels(rrna.gr) <- combined.seqlevels
seqlevels(mt.genes.regions.gr) <- combined.seqlevels
seqlevels(nc.genes.region.gr) <- combined.seqlevels
seqlevels(pc.genes.regions.gr) <- combined.seqlevels

all.genes.regions.gr <- c(rrna.gr,
                          mt.genes.regions.gr,
                          nc.genes.region.gr,
                          pc.genes.regions.gr)

saveRDS(all.genes.regions.gr, "~/Dropbox (The Francis Crick)/rna_structure/ref/human/gencode_v34.all.genes.regions.gr.rds")
export.gff2(all.genes.regions.gr, "~/Dropbox (The Francis Crick)/rna_structure/ref/human/gencode_v34.all.genes.regions.gff2")
all.genes.regions.gr$name <- all.genes.regions.gr$annotation
export.bed(all.genes.regions.gr, "~/Dropbox (The Francis Crick)/rna_structure/ref/human/gencode_v34.all.genes.regions.bed")
toc()

