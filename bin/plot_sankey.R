#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(networkD3))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(data.table))

# =========
# Options and paths
# =========

option_list <- list(make_option(c("-l", "--logs_dir"), action = "store", type = "character", default=NA, help = "Logs directory"),
make_option(c("-o", "--output"), action = "store", type = "character", default = 8, help = "Output name"))
opt_parser = OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# =========
# Load logs
# =========

logs.dir <- opt$logs_dir
all.logs <- list.files(logs.dir, full.names = TRUE)

# Load logs that are always present

cutadapt.log <- all.logs[str_detect(all.logs, ".cutadapt.log")] # trimming
# Extract sample name
sample_id <- str_split(basename(cutadapt.log), ".cutadapt.log")[[1]][1]
message("Analysing ", sample_id)

filter_blat.log <- all.logs[str_detect(all.logs, ".filter_blat.log")] # filter blat
identify_hybrids.log <- all.logs[str_detect(all.logs, "identify_hybrids.log")] # hybrid identification
dedup.log <- all.logs[str_detect(all.logs, ".dedup.log")] # dedup if umi dedup enabled; otherwise just logs for ambiguous removal

# Load optional logs

filter_spliced_reads.log <- all.logs[str_detect(all.logs, ".filter_spliced_reads.log")] # only if premapping enabled

# =========
# 1. Extract information from cutadapt and premapping logs
# =========

cutadapt.log <- readLines(cutadapt.log)
total.reads <- cutadapt.log[grep("^Total reads processed:", cutadapt.log)]
cutadapt_remaining.reads <- cutadapt.log[grep("^Reads written \\(passing filters\\):", cutadapt.log)]

min.length <- cutadapt.log[grep("^Command line parameters:", cutadapt.log)]
min.length <- parse_number(str_split(min.length, "--minimum-length ")[[1]][2])
total.reads <- parse_number(total.reads)
cutadapt_remaining.reads <- parse_number(cutadapt_remaining.reads)

cutadapt_discarded.reads <- total.reads - cutadapt_remaining.reads

if (length(cutadapt.log) != 0 & length(filter_spliced_reads.log) == 0) {

  message("premap log is empty")
  spliced.reads <- 0
  unspliced.reads <- cutadapt_remaining.reads

} else if (length(cutadapt.log) != 0 & length(filter_spliced_reads.log) == 1) {

  message("premap log is not empty")
  filter_spliced_reads.log <- readLines(filter_spliced_reads.log)
  trimmed.reads <- parse_number(filter_spliced_reads.log[grep("^Total reads:", filter_spliced_reads.log)])

  # Check that input reads for premapping is the same as the output reads from cutadapt
  stopifnot(trimmed.reads == cutadapt_remaining.reads)
  spliced.reads <- parse_number(filter_spliced_reads.log[grep("^Reads discarded:", filter_spliced_reads.log)])
  unspliced.reads <- parse_number(filter_spliced_reads.log[grep("^Reads kept:", filter_spliced_reads.log)])

}

# =========
# 2. Extract information from filter blat and identify hybrids logs
# =========

# pBLAT
blat_filter.log <- fread(filter_blat.log)
blat_mapped.reads <- blat_filter.log[blat_filter.log$Step == "initial"]$`Reads remaining`
blat_unmapped.reads <- cutadapt_remaining.reads - spliced.reads - blat_mapped.reads

too_high_evalue_antisense.reads <- blat_filter.log[str_detect(blat_filter.log$Step, "e-value|orientation")]$`Reads discarded`
too_many_blat_hits.reads <- blat_filter.log[str_detect(blat_filter.log$Step, "max hits")]$`Reads discarded`
blat_filter_remaining.reads <- blat_filter.log[str_detect(blat_filter.log$Step, "max hits")]$`Reads remaining`

# Extract `e-value` (using regex to match the number after ≤)
evalue_row <- blat_filter.log[str_detect(Step, "e-value|orientation")]
evalue <- str_extract(evalue_row$Step, "≤\\s*[0-9.]+") %>% str_remove("≤\\s*")
# Extract `max hits` (from e.g. "max hits (≤ 100)")
max_hits_row <- blat_filter.log[str_detect(Step, "max hits")]
max_hits <- str_extract(max_hits_row$Step, "≤\\s*[0-9]+") %>% str_remove("≤\\s*")

# Hybrid identification
identify_hybrids.log <- fread(identify_hybrids.log)
identify_hybrids_input.reads <- identify_hybrids.log[identify_hybrids.log$type =="initial_read_count"]$count
# Check that the input to hybrid identification is same as output of filtering blat
stopifnot(blat_filter_remaining.reads == identify_hybrids_input.reads)

strong_contiguous_match_to_a_single_gene.reads <- identify_hybrids.log[identify_hybrids.log$type =="strong_contiguous_match_to_a_single_gene"]$count
excessive_overlap_in_query_mappings.reads <- identify_hybrids.log[identify_hybrids.log$type =="excessive_overlap_in_query_mappings"]$count
excessive_gap_between_query_mappings.reads <- identify_hybrids.log[identify_hybrids.log$type =="excessive_gap_between_query_mappings"]$count
excessive_overlap_in_subject_mappings.reads <- identify_hybrids.log[identify_hybrids.log$type =="excessive_overlap_in_subject_mappings"]$count
identify_hybrids_remaining.reads <- identify_hybrids.log[identify_hybrids.log$type =="remaining_read_count"]$count

# =========
# 3. Extract information from dedup logs
# =========

dedup.log <- readLines(dedup.log)
ambiguous.reads <- parse_number(dedup.log[grep("ambiguous hybrids$", dedup.log)])
dedup_input.reads <- parse_number(dedup.log[grep("^Hybrids in:", dedup.log)])

# Check input reads for dedup + ambiguous are same as output reads from hybrid identification
stopifnot(ambiguous.reads + dedup_input.reads == identify_hybrids_remaining.reads)

dedup_remaining.reads <- parse_number(dedup.log[grep("^Hybrids out:", dedup.log)])
duplicated.reads <- dedup_input.reads - dedup_remaining.reads

# =========
# Build the data for the actual Sankey plot
# =========

read_list <- list()

# Cutadapt filter
read_list['total_reads'] <- total.reads

read_list['cutadapt_discarded'] <- cutadapt_discarded.reads
read_list['cutadapt_remaining'] <- cutadapt_remaining.reads

# Remove spliced reads
read_list['spliced'] <- spliced.reads
read_list['unspliced'] <- unspliced.reads

# Map
read_list['blat_unmapped'] <- blat_unmapped.reads
read_list['blat_mapped'] <- blat_mapped.reads

# Filter BLAT
read_list['too_high_evalue_antisense'] <- too_high_evalue_antisense.reads
read_list['too_many_blat_hits'] <- too_many_blat_hits.reads
read_list['blat_remaining'] <- blat_filter_remaining.reads

# Identify hybrids
read_list['strong_contiguous_match_to_a_single_gene'] <- strong_contiguous_match_to_a_single_gene.reads
read_list['excessive_overlap_in_query_mappings'] <- excessive_overlap_in_query_mappings.reads
read_list['excessive_gap_between_query_mappings'] <- excessive_gap_between_query_mappings.reads
read_list['excessive_overlap_in_subject_mappings']<- excessive_overlap_in_subject_mappings.reads
read_list['identify_hybrids_remaining'] <- identify_hybrids_remaining.reads

# Dedup and/or ambiguous
read_list['ambiguous'] <- ambiguous.reads
read_list['duplicated'] <- duplicated.reads
read_list['final'] <- dedup_remaining.reads

# Define nodes based on the new structure of read_list
nodes <- data.frame(
  name = c(
    'Total',
    'Cutadapt Discarded',
    paste0('Trimmed ≥ ', min.length) ,
    'Spliced',
    'Unspliced',
    'BLAT Unmapped',
    'BLAT Mapped',
    paste0('Too High (> ', evalue, ') E-value or Antisense'),
    paste0('Too Many (> ', max_hits, ') BLAT Hits'),
    'BLAT Remaining',
    'Strong Match to Single Gene',
    'Excessive Overlap in Query',
    'Excessive Gap in Query',
    'Excessive Overlap in Subject',
    'Valid Hybrids',
    'Ambiguous',
    'Duplicated',
    'Single or Multi-overlap',
    'Final'
  ),
  group = c(
    'a', 'b', 'b', 'c', 'c', 'd', 'd',
    'e', 'e', 'e', 'f', 'f', 'f', 'f', 'f',
    'g', 'g', 'g',
    'h'
  )
)

# Define links based on the updated read_list keys
links <- as.data.frame(rbind(
  c(0, 1, read_list[['cutadapt_discarded']]),
  c(0, 2, read_list[['cutadapt_remaining']]),
  c(2, 3, read_list[['spliced']]),
  c(2, 4, read_list[['unspliced']]),
  c(4, 5, read_list[['blat_unmapped']]),
  c(4, 6, read_list[['blat_mapped']]),
  c(6, 7, read_list[['too_high_evalue_antisense']]),
  c(6, 8, read_list[['too_many_blat_hits']]),
  c(6, 9, read_list[['blat_remaining']]),
  c(9, 10, read_list[['strong_contiguous_match_to_a_single_gene']]),
  c(9, 11, read_list[['excessive_overlap_in_query_mappings']]),
  c(9, 12, read_list[['excessive_gap_between_query_mappings']]),
  c(9, 13, read_list[['excessive_overlap_in_subject_mappings']]),
  c(9, 14, read_list[['identify_hybrids_remaining']]),
  c(14, 15, read_list[['ambiguous']]),
  c(14, 16, read_list[['duplicated']]),
  c(14, 17, read_list[['final']]),
  c(17, 18, read_list[['final']])
))

# Set column names for links
colnames(links) <- c('source', 'target', 'value')
links$link_source <- nodes$name[links$source + 1]

# Define color scheme
my_color <- 'd3.scaleOrdinal().domain(["a", "b", "c", "d", "e", "f", "g", "h"]).range(["#F3ECD9", "#F0F1E3", "#D7E0D8", "#C6D5D0", "#889C9B", "#7D7A70", "#5C625C"])'

# Plotting Sankey diagram
p <- sankeyNetwork(
  Links = links,
  Nodes = nodes,
  Source = 'source',
  Target = 'target',
  Value = 'value',
  LinkGroup = 'link_source',
  NodeID = 'name',
  NodeGroup = 'group',
  units = 'reads',
  fontSize = 16,
  nodeWidth = 50,
  fontFamily = "sans-serif",
  width = 1800,
  height = 1000,
  margin = c(top = 5, right = 1, bottom = 5, left = 1),
  nodePadding = 20,
  iterations = 10,
  colourScale = my_color,
  sinksRight = FALSE
)

saveNetwork(p, opt$output, selfcontained = FALSE)