#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(toscatools))
suppressPackageStartupMessages(library(parallel))

# ============
# MONKEY PATCH
# ============

# This function is designed to be applied per query, returning either a data.table of valid hybrids or a character string with a discard reason.
get_valid_hybrids <- function(
    blast.query.dt,
    min_unmapped_length = 16,
    q_minoverlap = 4,
    q_maxgap = 4,
    s_minoverlap = 0,
    xlink_distance = 1000000) {

  # Stop if the query is not unique
  if (length(unique(blast.query.dt$query)) != 1) {
    stop("Error: get_valid_hybrids must be applied per unique query. Found multiple or no queries.")
    }

  # Keep best match for a given query region
  hybrids.dt <- blast.query.dt[evalue == min_evalue]

  # Match up with fasta read length and remove if enough of a continuous match for any hit
  # if(any(hybrids.dt$unmapped < (min_unmapped_length - q_minoverlap) & hybrids.dt$unmapped == 100)) {
  if (any(hybrids.dt$unmapped < (min_unmapped_length - q_minoverlap))) {
    return("strong_contiguous_match_to_a_single_gene") # the read has too little unmapped sequence
  } else {

    # Now get combinations
    hybrids.dt[, id := 1:.N]
    hybrids.dt <- merge(hybrids.dt, hybrids.dt, by = c("query"), allow.cartesian = TRUE)
    hybrids.dt <- hybrids.dt[id.y > id.x] # Remove duplicates
    hybrids.dt[, id := paste0(id.x, "_", id.y)]

    # Remove those with significant overlap in the query mappings
    hybrids.dt[, q_ol := min(q_end.x, q_end.y) - max(q_start.x, q_start.y) + 1, by = id]
    hybrids.dt <- hybrids.dt[q_ol <= q_minoverlap]
    if (nrow(hybrids.dt) == 0) {
        return("excessive_overlap_in_query_mappings")
    }
    # Remove those with too large a gap between the query mappings
    hybrids.dt <- hybrids.dt[q_ol >= -q_maxgap]
    if (nrow(hybrids.dt) == 0) {
        return("excessive_gap_between_query_mappings")
    }
    # Remove those with significant overlap in the subject mappings, if subjects are the same
    hybrids.dt <- hybrids.dt[, s_ol := ifelse(subject.x == subject.y,
      min(s_end.x, s_end.y) - max(s_start.x, s_start.y) + 1,
      0
    ), by = id]
    hybrids.dt <- hybrids.dt[s_ol <= s_minoverlap]
    if (nrow(hybrids.dt) == 0) {
        return("excessive_overlap_in_subject_mappings")
    }

    # # Remove those too far away from xlink position
    # hybrids.dt <- hybrids.dt[q_start.x < xlink_distance | q_start.y < xlink_distance]
    # if (nrow(hybrids.dt) == 0) {
    #     return("hybrids_too_far_from_xlink_site")
    # }

    # Rename columns
    n <- names(hybrids.dt)
    n[grep("\\.x$", n)] <- paste0("L_", gsub("\\.x$", "", n[grep("\\.x$", n)]))
    n[grep("\\.y$", n)] <- paste0("R_", gsub("\\.y$", "", n[grep("\\.y$", n)]))
    n <- gsub("_s_", "_", n)
    n <- gsub("_subject", "_seqnames", n)
    setnames(hybrids.dt, n)
    hybrids.dt[, `:=`(L_strand = "+", R_strand = "+")]

    return(hybrids.dt)
  }
}

# ============
# MONKEY PATCH END
# ============

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