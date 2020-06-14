tsv <- fread("~/Dropbox (The Francis Crick)/rna_structure/ref/rat/rn6.repeatmasker.tsv.gz")
mask.dt <- tsv[repClass %in% c("Simple_repeat", "Low_complexity")]
fwrite(mask.dt[, .(genoName, genoStart, genoEnd, repName, ".", strand)],
       "~/Dropbox (The Francis Crick)/rna_structure/ref/rat/rn6.repeatmasker.sr_lc.bed.gz",
       col.names = FALSE, sep = "\t")
bed <- import.bed("~/Dropbox (The Francis Crick)/rna_structure/ref/rat/rn6.repeatmasker.sr_lc.bed.gz")
bed <- keepStandardChromosomes(bed, pruning.mode = "coarse")
seqlevelsStyle(bed) <- "NCBI"
bed$score <- 0
export.bed(bed, "~/Dropbox (The Francis Crick)/rna_structure/ref/rat/rn6.repeatmasker.sr_lc.bed.gz")
