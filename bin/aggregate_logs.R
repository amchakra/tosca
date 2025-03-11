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
aggregate_filter_blat_logs <- function(log_files, num_cores) {
  # Read all logs into a list of data.tables in parallel and remove "Elapsed time (s)" column
  log_list <- mclapply(log_files, function(log_file) {
    log.dt <- fread(log_file, sep = "\t", header = TRUE)
    elapsed_time_col <- "Elapsed time (s)"
    if (elapsed_time_col %in% colnames(log.dt)) {
      log.dt[, (elapsed_time_col) := NULL]  # Remove time column
    }
    return(log.dt)
  }, mc.cores = num_cores)

  # Check if the list is empty
  if (length(log_list) == 0) {
    stop("The list of logs is empty.")
  }

  # If only one log file, return it directly
  if (length(log_list) == 1) {
    aggregated.dt <- log_list[[1]]
  } else {
    # Ensure all data.tables have the same structure
    col_names <- names(log_list[[1]])
    # Bind all tables together while keeping column names consistent
    merged.dt <- rbindlist(log_list, use.names = TRUE, fill = TRUE)
    # Sum across all numerical columns grouped by "Step"
    cols_to_sum <- setdiff(col_names, "Step")
    aggregated.dt <- merged.dt[, lapply(.SD, sum), .SDcols = cols_to_sum, by = Step]
    # Ensure column order remains the same as the input
    setcolorder(aggregated.dt, col_names)
  }

  return(aggregated.dt)
}

# Function to read and aggregate identify_hybrids logs
aggregate_identify_hybrids_logs <- function(log_files, num_cores) {
  # Read all logs into a list of data.tables in parallel
  log_list <- mclapply(log_files, function(log_file) {
    log.dt <- fread(log_file, sep = "\t", header = TRUE)
    return(log.dt)
  }, mc.cores = num_cores)

  # Check if the list is empty
  if (length(log_list) == 0) {
    stop("The list of logs is empty.")
  }

  # If only one log file, return it directly
  if (length(log_list) == 1) {
    aggregated.dt <- log_list[[1]]
  } else {
    # Ensure all data.tables have the same structure
    col_names <- names(log_list[[1]])
    # Bind all tables together while keeping column names consistent
    merged.dt <- rbindlist(log_list, use.names = TRUE, fill = TRUE)
    # Sum across all numerical columns grouped by "type"
    cols_to_sum <- setdiff(col_names, "type")
    aggregated.dt <- merged.dt[, lapply(.SD, sum), .SDcols = cols_to_sum, by = type]
    # Ensure column order remains the same as the input
    setcolorder(aggregated.dt, col_names)
  }

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