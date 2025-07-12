#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(toscatools))
suppressPackageStartupMessages(library(parallel))

option_list <- list(make_option(c("-b", "--blast8"), action = "store", type = "character", help = "Blat blast8"),
            make_option(c("-f", "--fasta"), action = "store", type = "character", help = "Reads fasta"),
            make_option(c("-o", "--output"), action = "store", type = "character", help = "Output file"),
            make_option(c("-l", "--log"), action = "store", type = "character", help = "Output file for logs"),
            make_option(c("-t", "--threads"), action = "store", type = "integer", default = 8, help = "Number of threads [default: %default]"))
opt_parser = OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Prepare data for logging
discarded_reasons <- c(
    "strong_contiguous_match_to_a_single_gene",
    "excessive_overlap_in_query_mappings",
    "excessive_gap_between_query_mappings",
    "excessive_overlap_in_subject_mappings"
)

if(is.na(readLines(opt$blast8)[1])) {

    message("Empty blast8 input — writing empty output and default log.")
    fwrite(data.table(), file = opt$output, sep = "\t", col.names = TRUE) # Accounts for empty filtered blast file

    log.dt <- data.table(
      type = c("initial_read_count", discarded_reasons, "remaining_read_count"),
      count = c(0L, rep(0L, length(discarded_reasons)), 0L)
    )
    fwrite(log.dt, file = opt$log, sep = "\t")

} else {

    # Load blast and get read lengths
    blast.dt <- toscatools::load_blast8(opt$blast8)
    blast.dt <- toscatools::add_read_lengths(blast.dt = blast.dt, fasta = opt$fasta)
    blast.dt <- toscatools::calculate_blast8_metrics(blast.dt = blast.dt)

    blast.list <- split(blast.dt, blast.dt$query)

    cl <- makeForkCluster(opt$threads) # otherwise really slow...
    hybrids.list <- parLapply(cl = cl, blast.list, function(x) get_valid_hybrids(blast.query.dt = x))
    stopCluster(cl)

    # get_valid_hybrids returns a data.table if valid, or a character string with a discard reason if not
    discarded_bools <- sapply(hybrids.list, is.character)

    # message(sum(S4Vectors::elementNROWS(hybrids.list) == 0), " out of ", length(hybrids.list), " reads did not have hybrids")
    # message(round(sum(S4Vectors::elementNROWS(hybrids.list) != 0)/length(hybrids.list), 4) * 100, "% of reads had hybrids")

    # Filter multi hits
    hybrids.dt <- rbindlist(hybrids.list[!discarded_bools])

    if(nrow(hybrids.dt != 0)) {
        valid.hybrids.dt <- filter_valid_hybrids(hybrids.dt)
    } else {
        valid.hybrids.dt <- data.table()
    }

    fwrite(valid.hybrids.dt, file = opt$output, sep = "\t", col.names = TRUE)

    # Prepare data for logging
    discarded_reads.list <- hybrids.list[discarded_bools]

    if (length(unlist(discarded_reads.list)) > 0) {
        message("Logging discard reasons...")
        discarded_vector <- factor(unlist(discarded_reads.list), levels = discarded_reasons) # force factor with full levels
        discarded_reads_log.dt <- as.data.table(table(discarded_vector)) # table() will include the unused levels with count = 0
        setnames(discarded_reads_log.dt, c("type", "count"))
    } else {
        message("No reads were discarded.")
        discarded_reads_log.dt <- data.table(type = discarded_reasons, count = rep(0L, length(discarded_reasons)))
    }

    log.dt <- rbind(
        data.table(type = "initial_read_count", count = length(blast.list)),
        discarded_reads_log.dt,
        data.table(type = "remaining_read_count", count = length(hybrids.list[!discarded_bools]))
    )

    fwrite(log.dt, file = opt$log, sep = "\t")

}