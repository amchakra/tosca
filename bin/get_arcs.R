#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(optparse))

option_list <- list(make_option(c("", "--clusters"), action = "store", type = "character", help = "Clusters file"),
            make_option(c("", "--genes"), action = "store", type = "character", help = "List of genes"),
            make_option(c("", "--breaks"), action = "store", type = "character", default = "0,0.3,0.8,1", help = "Comma separated string of breaks for colours [default: %default]"),
            make_option(c("", "--output"), action = "store", type = "character", help = "Output stem"))
opt_parser = OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

clusters.dt <- fread(opt$clusters)
genes <- readLines(opt$genes) # Currently just one

for(i in seq_along(genes)) {

    goi <- genes[i]
    message(goi)

    cluster.dt <- clusters.dt[L_seqnames == R_seqnames][L_seqnames == goi]

    if(opt$breaks == "none") {

        if(nrow(cluster.dt > 0)) {

            # Create arcs
            cluster.dt[, colour := 0]
            bp.dt <- cluster.dt[, .(L_genomic_seqnames, L_genomic_start, L_genomic_end, R_genomic_start, R_genomic_end, colour)]
            setorder(bp.dt, colour)   

            # Generate colour scheme
            cols <- "#08306b"
            cols.dt <- rbindlist(lapply(cols, function(x) data.table(paste(as.vector(col2rgb(x)), collapse = "\t"))))
            cols.dt[, `:=` (colour = "color:", annotation = "Count quantile 1")]
            setcolorder(cols.dt, c("colour", "V1", "annotation"))

            # Write out
            fname <- paste0(opt$output, ".", gsub(":", "_", goi), ".bp")
            fwrite(cols.dt, sep = "\t", col.names = FALSE, quote = FALSE, file = fname)
            fwrite(bp.dt, sep = "\t", col.names = FALSE, file = fname, append = TRUE)

        } else {

            message("No clusters")

        }

    } else {

        if(nrow(cluster.dt > 0)) {

            # Get breaks
            p <- as.numeric(unlist(strsplit(opt$breaks, ",")))

            # Create arcs
            cluster.dt[, colour := cut(count, breaks = quantile(count, probs = p), include.lowest = TRUE, labels = 0:(length(p) - 2))]
            bp.dt <- cluster.dt[, .(L_genomic_seqnames, L_genomic_start, L_genomic_end, R_genomic_start, R_genomic_end, colour)]
            setorder(bp.dt, colour)

            # Generate colour scheme
            cols <- ggthemes::tableau_seq_gradient_pal('Blue')(seq(0, 1, length = length(p) - 1))
            cols.dt <- rbindlist(lapply(cols, function(x) data.table(paste(as.vector(col2rgb(x)), collapse = "\t"))))
            cols.dt[, `:=` (colour = "color:", annotation = paste("Count quantile", 1:(length(p) - 1)))]
            setcolorder(cols.dt, c("colour", "V1", "annotation"))

            # Write out
            fname <- paste0(opt$output, ".", gsub(":", "_", goi), ".bp")
            fwrite(cols.dt, sep = "\t", col.names = FALSE, quote = FALSE, file = fname)
            fwrite(bp.dt, sep = "\t", col.names = FALSE, file = fname, append = TRUE)

        } else {

            message("No clusters")

        }

    }

}