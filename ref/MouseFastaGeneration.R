# Script to get mouse custom fasta for ehiclipr (Flora's Staufen 2)
# A. M. Chakrabarti
# Last updated: 28th January 2020

# =========
# Get rRNA
# =========

library(rentrez)
library(ShortRead)
library(stringr)

rRNA.ids <- c("NR_046144.1","NR_003280.2","NR_003278.3","NR_003279.1") # NCBI rRNA 5S, 5.8S, 18S, 28S - why are there lots of 5S version?
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
library(BSgenome.Mmusculus.UCSC.mm10)
library(pbapply)
# library(ensembldb)

# =========
# Load annotation
# =========

genes.gr <- import.gff3("~/Dropbox (The Francis Crick)/rna_structure/ref/mouse/gencode.vM24.annotation.gff3.gz")
genes.gr <- genes.gr[genes.gr$type == "gene"]

findOverlaps(genes.gr, drop.self = TRUE, drop.redundant = TRUE)

# =========
# Get protein coding genes
# =========
pc.genes.gr <- genes.gr[genes.gr$gene_type == "protein_coding"]
pc.genes.gr <- dropSeqlevels(pc.genes.gr, "chrM", pruning.mode = "coarse") # will split out mitochondrial separately

# Get sequences
pc.regions.seq <- getSeq(Mmusculus, pc.genes.gr)
names(pc.regions.seq) <- paste0(pc.genes.gr$gene_name, "_", pc.genes.gr$gene_id)

# =========
# Get non coding genes
# =========
nc.genes.gr <- genes.gr[genes.gr$gene_type != "protein_coding"]
nc.genes.gr <- dropSeqlevels(nc.genes.gr, "chrM", pruning.mode = "coarse") # will split out mitochondrial separately

# Get sequences
nc.regions.seq <- getSeq(Mmusculus, nc.genes.gr)
names(nc.regions.seq) <- paste0(nc.genes.gr$gene_name, "_", nc.genes.gr$gene_id)

# =========
# Get mitochondrial genes
# =========
mt.genes.gr <- genes.gr[seqnames(genes.gr) == "chrM"]

# Get sequences
mt.regions.seq <- getSeq(Mmusculus, mt.genes.gr)
names(mt.regions.seq) <- paste0(mt.genes.gr$gene_name, "_", mt.genes.gr$gene_id)

# =========
# Write out fasta
# =========
all.seq <- c(rrna.seq, pc.regions.seq, nc.regions.seq, mt.regions.seq)
writeFasta(all.seq, "~/Dropbox (The Francis Crick)/rna_structure/ref/mouse/Mm_GencodeM24_rRNA_MT_genes.fa")
system("pigz ~/Dropbox\\ \\(Lab\\)/rna_structure/ref/mouse/Mm_GencodeM24_rRNA_MT_genes.fa")

# Combine for coordinate conversion
seq.gr <- c(pc.genes.gr, nc.genes.gr, mt.genes.gr)
seq.gr$fasta_id <- paste0(seq.gr$gene_name, "_", seq.gr$gene_id)
seq.gr <- sort(seq.gr)
export.gff2(seq.gr, "~/Dropbox (The Francis Crick)/rna_structure/ref/mouse/Mm_GencodeM24_rRNA_MT_genes.gtf")
system("pigz ~/Dropbox\\ \\(Lab\\)/rna_structure/ref/mouse/Mm_GencodeM24_rRNA_MT_genes.gtf")







# ======================================================
# Now annotate with regions of longest transcript
# ======================================================

TxDb <- makeTxDbFromGFF("~/Dropbox (The Francis Crick)/ehiCLIP/Mm/ref/gencode.vM24.annotation.gtf.gz")
saveDb(TxDb, "~/Dropbox (The Francis Crick)/ehiCLIP/Mm/ref/gencode.vM24.sqlite")

TxDb <- keepStandardChromosomes(TxDb, pruning.mode = "coarse")

# Get transcript lengths and add to annotation data.table
txlengths.dt <- data.table(transcriptLengths(TxDb, with.cds_len = TRUE, with.utr5_len = TRUE, with.utr3_len = TRUE))
setnames(txlengths.dt, c("tx_name", "gene_id"), c("transcript_id", "gene_id"))
txlengths.dt[, tx_id := NULL]

gtf <- import.gff2("~/Dropbox (The Francis Crick)/ehiCLIP/Mm/ref/gencode.vM24.annotation.gtf.gz")
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
lncRNA <- c("lncRNA", "IG_C_pseudogene", "IG_D_pseudogene", "IG_pseudogene", "IG_V_pseudogene", "TR_J_pseudogene", "polymorphic_pseudogene", "processed_pseudogene", "processed_transcript", "pseudogene", "transcribed_processed_pseudogene", "transcribed_unitary_pseudogene", "transcribed_unprocessed_pseudogene", "translated_processed_pseudogene", "translated_unprocessed_pseudogene", "unitary_pseudogene", "unprocessed_pseudogene")
ncRNA <- c("misc_RNA", "ribozyme", "scaRNA", "scRNA", "snoRNA", "snRNA", "sRNA")
other <- c("non_stop_decay", "nonsense_mediated_decay", "TEC")

stopifnot(all(c(pc, rRNA, tRNA, miRNA, miRNA, lncRNA, ncRNA, other) %in% tx_biotype))

# ==========
# Protein coding regions
# ==========

# Select longest pc
gencode.dt[transcript_type %in% pc, longest := max(tx_len), by = ensembl_gene_id]
longest.pc.dt <- gencode.dt[transcript_type %in% pc & tx_len == longest] # selects longest
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
pc.regions.grl <- split(pc.regions.gr, pc.regions.gr$ensembl_gene_id)
pc.regions.grl <- sort(pc.regions.grl)

cl <- makeForkCluster(4)
pc.regions.tc.grl <- parLapply(cl = cl, 1:length(pc.regions.grl), function(i) {

  message(i)
  x <- pc.regions.grl[[i]]

  if(unique(strand(x)) == "+") {

    start <- min(start(x))
    tc.gr <- x
    tc.start <- start(tc.gr) - start + 1
    tc.end <- end(tc.gr) - start + 1

    tc.gr <- with(tc.gr, GRanges(seqnames = Rle(ensembl_gene_id),
                                 ranges = IRanges(start = tc.start, end = tc.end),
                                 strand = Rle("+"),
                                 ensembl_gene_id = ensembl_gene_id,
                                 gene_name = gene_name,
                                 biotype = biotype,
                                 region = region))

    stopifnot(all(width(tc.gr) == width(x)))
    stopifnot(length(reduce(tc.gr)) == 1)


  } else if(unique(strand(x)) == "-") {

    tc.gr <- x

    # Need to fix this bit so annotate those bits of the gene thar aren't overlapped by the longest transcript
    # gene.gr <- pc.genes.gr[pc.genes.gr$gene_id == unique(tc.gr$ensembl_gene_id)]
    # rest.gr <- setdiff(gene.gr, tc.gr)

    start <- max(end(x))


    flipped.start <- abs(end(tc.gr) - start + 1)
    flipped.end <- flipped.start + width(tc.gr) - 1

    tc.gr <- with(tc.gr, GRanges(seqnames = Rle(ensembl_gene_id), # Maybe add underscore name here
                                 ranges = IRanges(start = flipped.start, end = flipped.end),
                                 strand = Rle("+"),
                                 ensembl_gene_id = ensembl_gene_id,
                                 gene_name = gene_name,
                                 biotype = biotype,
                                 region = region))

    stopifnot(all(width(tc.gr) == width(x)))
    stopifnot(length(reduce(tc.gr)) == 1)

  }

  return(tc.gr)

})
stopCluster(cl)

pc.regions.tc.gr <- unlist(GRangesList(pc.regions.tc.grl))
saveRDS(pc.regions.tc.gr, "~/Dropbox (The Francis Crick)/ehiCLIP/Mm/ref/gencode_vM24_pc.regions.tc.gr.rds")
