#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(toscatools))
suppressPackageStartupMessages(library(Rcpp))
suppressPackageStartupMessages(library(optparse))

option_list <- list(make_option(c("", "--hybrids"), action = "store", type = "character", help = "Hybrids file"),
            make_option(c("", "--genes"), action = "store", type = "character", help = "List of genes"),
            make_option(c("", "--fai"), action = "store", type = "character", help = "Transcript fasta index"),
            make_option(c("", "--bin_size"), action = "store", type = "integer", help = "Bin size"),
            make_option(c("", "--output"), action = "store", type = "character", help = "Output stem"))
opt_parser = OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# ==========
# Local functions, to be incorp. into toscatools when tested
# ==========

# Reorienting function adjusted for segments
reorient_hybrids_segments <- function(hybrids.dt) {

  # First do starts
  correct.dt <- hybrids.dt[L_start <= R_start]
  incorrect.dt <- hybrids.dt[L_start > R_start]

  renamed <- gsub("^L_", "X_", names(incorrect.dt))
  renamed <- gsub("^R_", "L_", renamed)
  renamed <- gsub("^X_", "R_", renamed)

  setnames(incorrect.dt, renamed)
  reoriented.dt <- rbindlist(list(correct.dt, incorrect.dt), use.names = TRUE)
  stopifnot(all(reoriented.dt$L_start <= reoriented.dt$R_start))

  # Then do segments (to make sure intergenics in same order)
  reoriented.dt[, `:=` (L_segment = tstrsplit(L_seqnames, "_")[[2]],
                        R_segment = tstrsplit(R_seqnames, "_")[[2]])]
  correct.dt <- reoriented.dt[L_segment <= R_segment]
  incorrect.dt <- reoriented.dt[L_segment > R_segment]

  renamed <- gsub("^L_", "X_", names(incorrect.dt))
  renamed <- gsub("^R_", "L_", renamed)
  renamed <- gsub("^X_", "R_", renamed)

  setnames(incorrect.dt, renamed)

  reoriented.dt <- rbindlist(list(correct.dt, incorrect.dt), use.names = TRUE)
  stopifnot(all(reoriented.dt$L_subject <= reoriented.dt$R_subject))
  stopifnot(nrow(reoriented.dt) == nrow(hybrids.dt))

  return(reoriented.dt)
}

src <- 
  '
  IntegerMatrix rcpp_get_contact_map(DataFrame hybrids, int gene_A_size, int gene_B_size) {
    
    IntegerVector L_start_v = hybrids["L_start"];   
    IntegerVector L_end_v = hybrids["L_end"];
    IntegerVector R_start_v = hybrids["R_start"];
    IntegerVector R_end_v = hybrids["R_end"];
    
    IntegerMatrix contact_map(gene_A_size, gene_B_size);

    int n_hybrids = hybrids.nrows();
    for (int i = 0; i < n_hybrids; i++) {
    
      int L_start = L_start_v[i] - 1;
      int L_end = L_end_v[i] - 1;
      int R_start = R_start_v[i] - 1;
      int R_end = R_end_v[i] - 1;
    
      for (int x = L_start; x <= L_end; x++) {
        for (int y = R_start; y <= R_end; y++) {
          contact_map(x, y) ++;
        }
      }
      
    }

    return(contact_map);
    
  }
  '
cppFunction(src)

# ==========

genes <- readLines(opt$genes)

# Intra-transcript maps

if(!any(grepl(",", genes))) {

    for(i in seq_along(genes)) {

        goi <- genes[i]
        message(goi)

        fai.dt <- fread(opt$fai, select = 1:2, col.names = c("gene", "length"))
        genome.size <- as.integer(fai.dt[gene == goi]$length)
        hybrid.dt <- fread(opt$hybrids)
        hybrid.dt <- hybrid.dt[type == "intragenic"][L_seqnames == goi]

        if(!nrow(hybrid.dt) == 0) {

            mat <- get_contact_map(hybrid.dt = hybrid.dt, genome.size = genome.size)
            saveRDS(mat, file = paste0(opt$output, ".", goi, ".mat.rds"))

            # Bin matrix and normalise
            if(opt$bin_size == 1) {
                binned.mat <- mat
            } else {
                binned.mat <- bin_matrix(mat, bin.size = opt$bin_size)
            }

            # binned.dt <- data.table(reshape2::melt(binned.mat))
            binned.dt <- melt(as.data.table(binned.mat)[, rn := 1:.N], id.vars = "rn")
            binned.dt[, variable := as.integer(gsub("^V", "", as.character(variable)))]
            binned.dt[, norm_value := value*1e6/nrow(hybrid.dt)]
            binned.dt <- binned.dt[value != 0]

            fwrite(binned.dt, file = paste0(opt$output, ".", goi, ".", opt$bin_size, "_binned_map.tsv.gz"), sep = "\t")

        } else {

            message("No hybrids")

        }

    }

} else if(all(grepl(",", genes))) {

    for(j in seq_along(genes)) {
        
        # Get two segments
        goi <- genes[j]
        # goi <- sort(unlist(tstrsplit(goi, ","))) # ensures alphabetical so matches L & R for reorient_hybrids
        goi <- unlist(tstrsplit(goi, ",")) # ensures alphabetical so matches L & R for reorient_hybrids
        gene_A <- goi[1]
        gene_B <- goi[2]
        
        # Get lengths for each segment
        fai.dt <- fread(opt$fai, select = 1:2, col.names = c("gene", "length"))
        gene_A.size <- as.integer(fai.dt[gene == gene_A]$length)
        gene_B.size <- as.integer(fai.dt[gene == gene_B]$length)
        
        hybrid.dt <- fread(opt$hybrids)
        # hybrid.dt <- reorient_hybrids_segments(hybrid.dt)
        hybrid.dt <- reorient_hybrids(hybrid.dt)
        hybrid.dt <- hybrid.dt[type == "intergenic"][L_seqnames == gene_A][R_seqnames == gene_B]

        if(!nrow(hybrid.dt) == 0) {

            # Now create matrix and fill
            # mat <- matrix(data = 0, nrow = gene_A.size, ncol = gene_B.size)
            
            # for (i in 1:nrow(hybrid.dt)) {
            # mat[
            #     hybrid.dt[i]$L_start:hybrid.dt[i]$L_end,
            #     hybrid.dt[i]$R_start:hybrid.dt[i]$R_end
            # ] <- mat[
            #     hybrid.dt[i]$L_start:hybrid.dt[i]$L_end,
            #     hybrid.dt[i]$R_start:hybrid.dt[i]$R_end
            # ] + 1
            # }
            
            mat <- rcpp_get_contact_map(hybrid.dt, gene_A.size, gene_B.size)
            saveRDS(mat, file = paste0(opt$output, ".", gene_A, "-", gene_B, ".mat.rds"))

            # Don't bin inter as current binning code assumes a square matrix
            binned.dt <- melt(as.data.table(mat)[, rn := 1:.N], id.vars = "rn")
            binned.dt[, variable := as.integer(gsub("^V", "", as.character(variable)))]
            binned.dt[, norm_value := value*1e6/nrow(hybrid.dt)]
            binned.dt <- binned.dt[value != 0]

            fwrite(binned.dt, file = paste0(opt$output, ".", gene_A, "-", gene_B, ".", 1, "_binned_map.tsv.gz"), sep = "\t")
        
        } else {

            message("No hybrids")

        }

    }

} else {

    message("--genes cannot be a mix of intra-transcript and inter-transcript yet.")

}