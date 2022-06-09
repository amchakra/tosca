#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(toscatools))
suppressPackageStartupMessages(library(parallel))

option_list <- list(make_option(c("-b", "--blast8"), action = "store", type = "character", help = "Blat blast8"),
            make_option(c("-f", "--fasta"), action = "store", type = "character", help = "Reads fasta"),
            make_option(c("-o", "--output"), action = "store", type = "character", help = "Output file"),
            make_option(c("-t", "--threads"), action = "store", type = "integer", default = 8, help = "Number of threads [default: %default]"))
opt_parser = OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if(is.na(readLines(opt$blast8)[1])) {

    fwrite(data.table(), file = opt$output, sep = "\t", col.names = TRUE) # Accounts for empty filtered blast file

} else {

    # Load blast and get read lengths
    blast.dt <- load_blast8(opt$blast8)
    blast.dt <- add_read_lengths(blast.dt = blast.dt, fasta = opt$fasta)
    blast.dt <- calculate_blast8_metrics(blast.dt = blast.dt)

    blast.list <- split(blast.dt, blast.dt$query)
    cl <- makeForkCluster(opt$threads) # otherwise really slow...
    hybrids.list <- parLapply(cl = cl, blast.list, function(x) get_valid_hybrids(blast.query.dt = x))
    stopCluster(cl)

    # message(sum(S4Vectors::elementNROWS(hybrids.list) == 0), " out of ", length(hybrids.list), " reads did not have hybrids")
    # message(round(sum(S4Vectors::elementNROWS(hybrids.list) != 0)/length(hybrids.list), 4) * 100, "% of reads had hybrids")

    # Filter multi hits
    hybrids.dt <- rbindlist(hybrids.list)
    if(nrow(hybrids.dt != 0)) {
        valid.hybrids.dt <- filter_valid_hybrids(hybrids.dt)
    } else {
        valid.hybrids.dt <- data.table()
    }

    fwrite(valid.hybrids.dt, file = opt$output, sep = "\t", col.names = TRUE)

}