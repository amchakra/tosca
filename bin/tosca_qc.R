#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(data.table))

# ==========
# Get hybrid identification metrics
# ==========

message("Getting hybrid identification metrics...")

# Get input and spliced or trimmed reads
# filter_premap_reads.log (current) or filter_spliced_reads.log (older runs, no linear filter)
premap.logs <- list.files(".", pattern = "filter_(premap|spliced)_reads.log$", full.names = TRUE)

if (length(premap.logs) == 0) {
  message("No premap log files found (pre-mapping was skipped). Using cutadapt logs to get read counts.")

  # Extract sample names and trimmed read counts from cutadapt logs if available
  cutadapt.logs <- list.files(".", pattern = ".cutadapt.log$", full.names = TRUE)
  if (length(cutadapt.logs) > 0) {
    sample_names <- sapply(basename(cutadapt.logs), function(x) strsplit(x, "\\.")[[1]][1])

    # Extract the "Reads written (passing filters)" value from each log file
    read_counts <- sapply(cutadapt.logs, function(x) {
      # Find the line containing the number of trimmed reads
      target_line <- grep("Reads written \\(passing filters\\):", readLines(x), value = TRUE)
      as.numeric(sub(".*Reads written \\(passing filters\\):\\s*([0-9]+).*", "\\1", target_line))
  })
    premap.dt <- data.table(sample = sample_names,
                            total = read_counts,
                            unspliced = read_counts,
                            spliced = 0L,
                            linear = 0L)

  } else {
    stop("No premap or cutadapt log files found. Exiting.")
  }

} else {
  # Read counts by label, so the 4-line log (with Reads linear) and older 3-line logs both parse
  log_count <- function(lines, label) {
    hit <- grep(paste0("^", label, ":"), lines, value = TRUE)
    if (length(hit) == 0) return(0)
    as.numeric(sub(".*:\\s*", "", hit[1]))
  }
  premap.dt <- lapply(seq_along(premap.logs), function(i) {
    lines <- readLines(premap.logs[i])
    data.table(
      sample = tstrsplit(basename(premap.logs[i]), "\\.")[[1]],
      total = log_count(lines, "Total reads"),
      unspliced = log_count(lines, "Reads kept"),       # reads passed on to BLAT
      spliced = log_count(lines, "Reads discarded"),
      linear = log_count(lines, "Reads linear")
    )
  })
  premap.dt <- rbindlist(premap.dt)
}

# Get hybrid breakdown
hybrids.files <- list.files(".", pattern = ".hybrids.tsv.gz$", full.names = TRUE)
hybrid_identification.dt <- lapply(seq_along(hybrids.files), function(i) {

  dt <- fread(hybrids.files[i])
  hs.dt <- unique(dt[, .(name, hybrid_selection)])[, .N, by = hybrid_selection]
  setnames(hs.dt, "N", tstrsplit(basename(hybrids.files[i]), "\\.")[[1]])
  return(hs.dt)

})

# Merge the hybrid identification data and premap data
hybrid_identification.dt <- Reduce(function(x, y) merge(x, y, by = "hybrid_selection", all = TRUE),
                                   hybrid_identification.dt)
hybrid_identification.dt <- transpose(hybrid_identification.dt, keep.names = "sample", make.names = "hybrid_selection")
hybrid_identification.dt <- merge(hybrid_identification.dt, premap.dt, by = "sample", all.x = TRUE)
hybrid_identification.dt[is.na(hybrid_identification.dt)] <- 0  # Replace any NAs

# Ensure the required columns exist; if not, add them with default 0
hybrid_selection.cols <- c("ambiguous", "multi_overlap", "single")
for (col in hybrid_selection.cols) {
  if (!col %in% names(hybrid_identification.dt)) {
    hybrid_identification.dt[, (col) := 0]
  }
}

# Calculate nonhybrid as the reads passed to BLAT minus the sum of hybrid categories
hybrid_identification.dt[, nonhybrid := unspliced - ambiguous - multi_overlap - single]

# Select final columns
hybrid_identification.dt <- hybrid_identification.dt[, .(sample, spliced, linear, nonhybrid, ambiguous, multi_overlap, single)]

stopifnot(all(rowSums(hybrid_identification.dt[, -1]) %in% premap.dt$total)) # order might be different

fwrite(hybrid_identification.dt, "hybrid_identification.tsv", sep = "\t")

# ==========
# Get intra-inter, genomic-reverse
# ==========

message("Getting intra-inter, genomic-reverse metrics...")

hybrids.files <- list.files(".", pattern = ".hybrids.gc.annotated.tsv.gz$", full.names = TRUE)
if(length(hybrids.files) == 0) hybrids.files <- list.files(".", pattern = ".hybrids.gc.annotated.mfe.tsv.gz$", full.names = TRUE) # If analyse_structures
if(length(hybrids.files) == 0) hybrids.files <- list.files(".", pattern = ".hybrids.gc.annotated.mfe.shuffled.tsv.gz$", full.names = TRUE) # If shuffled_mfe

hybrids.list <- lapply(hybrids.files, fread)
hybrids.dt <- rbindlist(hybrids.list, use.names = TRUE)
hybrids.dt$sample <- rep(tstrsplit(basename(hybrids.files), "\\.")[[1]], S4Vectors::elementNROWS(hybrids.list))

intrainter.dt <- hybrids.dt[, .N, by = .(sample, type, orientation)]
intrainter.dt[type == "intragenic", type := paste0(type, "-", orientation)][, orientation := NULL]
intrainter.dt <- dcast.data.table(intrainter.dt, sample ~ type, value.var = "N")

fwrite(intrainter.dt, "intrainter.tsv", sep = "\t")

# ==========
# Get arm regions
# ==========

message("Getting arm region metrics...")

# Non-genomic reference entries (rRNA, tRNA, snRNA, 7SK, 7SL, Y RNA, whole mitochondrial genome) are labelled
# by class by toscatools::annotate_hybrids(); gene arms carry their regions.gtf region
nongenomic_regions <- c("rRNA", "tRNA", "snRNA", "7SK", "7SL", "Y_RNA", "Mt")
regions <- c(nongenomic_regions, "ncRNA", "UTR5", "CDS", "intron", "UTR3", "intergenic")

regions.dt <- rbindlist(list(hybrids.dt[, .(sample, L_region)],
                             hybrids.dt[, .(sample, R_region)]),
                        use.names = FALSE)
regions.dt <- regions.dt[, .N, by = .(sample, L_region)]

regions.dt <- dcast.data.table(regions.dt, sample ~ L_region, value.var = "N", fun.aggregate = sum, fill = 0L) # need fun.aggregate otherwise reverts to length for regions that are all 0 for some reason
missing.cols <- setdiff(regions, names(regions.dt))
for(col in missing.cols) regions.dt[, (col) := 0L] # 0 counts for regions absent from every sample (e.g. older references)
regions.dt <- regions.dt[, c("sample", regions), with = FALSE]

fwrite(regions.dt, "regions.tsv", sep = "\t")

# ==========
# Get links
# ==========

message("Getting link metrics...")

mRNA <- c("CDS", "intron", "UTR3", "UTR5")
links.dt <- hybrids.dt

links.dt[L_region == "rRNA" & R_region == "rRNA", link := "rRNA-rRNA"]
links.dt[L_region %in% mRNA & R_region %in% mRNA, link := "mRNA-mRNA"]
links.dt[L_region %in% mRNA & R_region == "rRNA", link := "rRNA-mRNA"]
links.dt[L_region == "rRNA" & R_region %in% mRNA, link := "rRNA-mRNA"]
links.dt[is.na(link), link := "other"]
links.dt <- links.dt[, .N, by = .(sample, link)]
# Reshaping to wide format
links.dt <-  dcast.data.table(links.dt, sample ~ link, value.var = "N")

# Set column order
links.cols <- c("sample", "rRNA-rRNA", "rRNA-mRNA", "mRNA-mRNA", "other") # all possible column names
missing.cols <- setdiff(links.cols, names(links.dt))
lapply(missing.cols, function(col) links.dt[, (col) := 0L]) # assign 0 counts to missing categories

links.dt <- links.dt[, .(sample, `rRNA-rRNA`, `rRNA-mRNA`, `mRNA-mRNA`, other)]

fwrite(links.dt, "links.tsv", sep = "\t")

# ==========
# Get intragenic mRNA regions
# ==========

message("Getting intragenic mRNA region metrics...")

regions <- c("ncRNA", "UTR5", "CDS", "intron", "UTR3")

intragenic_regions.dt <- hybrids.dt[type == "intragenic"][!L_region %in% c(nongenomic_regions, "intergenic")]

intragenic_regions.dt <- rbindlist(list(intragenic_regions.dt[, .(sample, L_region)],
                                        intragenic_regions.dt[, .(sample, R_region)]),
                        use.names = FALSE)
intragenic_regions.dt <- intragenic_regions.dt[, .N, by = .(sample, L_region)]

message(nrow(intragenic_regions.dt))

intragenic_regions.dt <-  dcast.data.table(intragenic_regions.dt, sample ~ L_region, value.var = "N", fun.aggregate = sum, fill = 0L)
missing.cols <- setdiff(regions, names(intragenic_regions.dt))
for(col in missing.cols) intragenic_regions.dt[, (col) := 0L] # 0 counts for regions absent from every sample
intragenic_regions.dt <- intragenic_regions.dt[, c("sample", regions), with = FALSE]

fwrite(intragenic_regions.dt, "intragenic_regions.tsv", sep = "\t")

# ==========
# Get clusters
# ==========

message("Getting cluster metrics...")

clusters.files <- list.files(".", pattern = ".clusters.gc.annotated.tsv.gz$", full.names = TRUE)
clusters.list <- lapply(clusters.files, fread)
clusters.dt <- rbindlist(clusters.list, use.names = TRUE)
clusters.dt$sample <- rep(tstrsplit(basename(clusters.files), "\\.")[[1]], S4Vectors::elementNROWS(clusters.list))
clusters.dt[, type := ifelse(L_seqnames == R_seqnames, "intragenic", "intergenic")]

duplexes.dt <- clusters.dt[, .N, by = .(sample, type)]
duplexes.dt <- dcast.data.table(duplexes.dt, sample ~ type, value.var = "N")

fwrite(duplexes.dt, "duplexes.tsv", sep = "\t")