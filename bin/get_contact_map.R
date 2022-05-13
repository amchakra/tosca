#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(toscatools))
suppressPackageStartupMessages(library(optparse))

option_list <- list(make_option(c("", "--hybrids"), action = "store", type = "character", help = "Hybrids file"),
            make_option(c("", "--genes"), action = "store", type = "character", help = "List of genes"),
            make_option(c("", "--fai"), action = "store", type = "character", help = "Transcript fasta index"),
            make_option(c("", "--bin_size"), action = "store", type = "integer", help = "Bin size"),
            make_option(c("", "--output"), action = "store", type = "character", help = "Output stem"))
opt_parser = OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)


genes <- readLines(opt$genes) # Currently just one

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
        binned.dt[, `:=` variable := as.integer(gsub("^V", "", as.character(variable)))]
        binned.dt[, norm_value := value*1e6/nrow(hybrid.dt)]
        binned.dt <- binned.dt[value != 0]

        fwrite(binned.dt, file = paste0(opt$output, ".", goi, ".", opt$bin_size, "_binned_map.tsv.gz"), sep = "\t")

    } else {

        message("No hybrids")

    }

}