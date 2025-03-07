#!/usr/bin/env Rscript

# Script to combine filter_blat.py and identify_hybrids.R logs from all read chunks
# I.A. Iosub
# 7th March 2025

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(parallel))

# ==========
# Functions
# ==========

# Function to read and aggregate filter_blat logs
aggregate_filter_blat_logs <- function(log_files, cores) {
  # Read all logs into a list of data.tables in parallel and remove "Elapsed time (s)" column
  log_list <- mclapply(log_files, function(log_file) {
    log.dt <- fread(log_file, sep = "\t", header = TRUE)
    log.dt[, `Elapsed time (s)` := NULL]  # Remove time column
    return(log.dt)
  }, mc.cores = cores)

  # Aggregate all logs by summing counts for each step
  aggregated.dt <- Reduce(function(x, y) {
    x[, .(`Reads remaining` = `Reads remaining` + y$`Reads remaining`,
        `Reads discarded` = `Reads discarded` + y$`Reads discarded`,
        `BLAT mappings remaining` = `BLAT mappings remaining` + y$`BLAT mappings remaining`),
        by = Step]
  }, log_list)

  return(aggregated.dt)
}

# Function to read and aggregate identify_hybrids logs
aggregate_identify_hybrids_logs <- function(log_files, cores) {
  # Read all logs into a list of data.tables in parallel
  log_list <- mclapply(log_files, function(log_file) {
    library(data.table)  # Needed inside parLapply
    log.dt <- fread(log_file, sep = "\t", header = TRUE)
    return(log.dt)
  }, mc.cores = cores)

  # Aggregate all logs by summing counts for each type
  aggregated.dt <- Reduce(function(x, y) {
    merge(x, y, by = "type", all = TRUE, suffixes = c(".x", ".y"))[
      , .(type, count = rowSums(.SD, na.rm = TRUE)), .SDcols = patterns("count")
    ]
  }, log_list)

  return(aggregated.dt)
}


# ==========
# Run
# ==========

option_list <- list(make_option(c("-l", "--logs"), action = "store", type = "character", help = "Logs directory"),
            make_option(c("-t", "--threads"), action = "store", type = "integer", default = 8, help = "Number of threads [default: %default]"),
            make_option(c("-o", "--output"), action = "store", type = "character", default = 8, help = "Output prefix"))
opt_parser = OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

filter_blat_logs.list <- list.files(opt$logs, pattern = ".filter_blat.log", full.names = TRUE)

identify_hybrids_logs.list <- list.files(opt$logs, pattern = ".identify_hybrids.log", full.names = TRUE)

# Ensure only logs for the specified sample are used
filter_blat_logs.list <- filter_blat_logs.list[str_detect(filter_blat_logs.list, opt$output)]
identify_hybrids_logs.list <- identify_hybrids_logs.list[str_detect(identify_hybrids_logs.list, opt$output)]

# Aggregate filter_blat logs
if (length(filter_blat_logs.list) > 0) {
  message("Aggregating filter_blat logs...\n")
  filter_blat_aggregated.dt <- aggregate_filter_blat_logs(filter_blat_logs.list, opt$threads)
  fwrite(filter_blat_aggregated.dt, file = paste0(opt$output,".filter_blat.log"), sep = "\t", col.names = TRUE)
} else {
  message("No filter_blat logs found.\n")
}

# Aggregate identify_hybrids logs
if (length(identify_hybrids_logs.list) > 0) {
  message("Aggregating identify_hybrids logs...\n")
  identify_hybrids_aggregated.dt <- aggregate_identify_hybrids_logs(identify_hybrids_logs.list, opt$threads)
  fwrite(identify_hybrids_aggregated.dt, file = paste0(opt$output,".identify_hybrids.log"), sep = "\t", col.names = TRUE)
} else {
  message("No identify_hybrids logs found.\n")
}

message("Aggregation complete.\n")