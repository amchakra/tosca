

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(Biostrings))
suppressPackageStartupMessages(library(pbapply))
suppressPackageStartupMessages(library(GenomicInteractions))

option_list <- list(make_option(c("-b", "--blast"), action = "store", type = "character", help = "Filtered blast output"),
                    make_option(c("-f", "--fasta"), action = "store", type = "character", help = "Fasta of reads"),
                    make_option(c("-o", "--output"), action = "store", type = "character", help = "Output file"),
                    make_option(c("-t", "--threads"), action = "store", type = "character", help = "Number of threads"))
opt_parser = OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# opt <- list(blast="results/filtered/demux_NNNNATCTGGGANNNNN.00.filtered.blast8.gz",
#   output = "test.tsv",
#   threads = 8)

# ==========
# Functions
# ==========

# GetValidHybrids <- function(blast8.query.dt, fa.dt, min_unmapped_length = 21, q_minoverlap = 4, q_maxgap = 4, s_minoverlap = 0, xlink_distance = 5) {
  
#   hybrids.dt <- blast8.query.dt

#   # Keep best match for a given query region
#   hybrids.dt[, min_evalue := min(evalue), by = .(q_start, q_end)]
#   hybrids.dt <- hybrids.dt[evalue == min_evalue]  

#   # Match up with fasta read length and remove if enough of a continuous match
#   stopifnot(unique(hybrids.dt$query) %in% fa.dt$query)
#   hybrids.dt <- merge(hybrids.dt, fa.dt, by = "query", all.x = TRUE)
#   hybrids.dt[, unmapped := readlength - alignment_length]
#   hybrids.dt[, mapped := alignment_length/readlength]  
  
#   if(any(hybrids.dt$unmapped < min_unmapped_length)) return(hybrids.dt)

#   # Now get combinations
#   hybrids.dt[, id := 1:.N]
#   hybrids.dt <- merge(hybrids.dt, hybrids.dt, by = c("query"), allow.cartesian = TRUE)
#   hybrids.dt <- hybrids.dt[id.y > id.x] # Remove duplicates
#   hybrids.dt[, id := paste0(id.x, "_", id.y)]
  
#   # Remove those with significant overlap in the query mappings and too large a gap between the query mappings
#   hybrids.dt[, q_ol := min(q_end.x, q_end.y) - max(q_start.x, q_start.y) + 1, by = id]
#   hybrids.dt <- hybrids.dt[q_ol <= q_minoverlap][q_ol > -q_maxgap]
  
#   # Remove those with significant overlap in the subject mappings, if subjects are the same
#   hybrids.dt <- hybrids.dt[, s_ol := ifelse(subject.x == subject.y,
#                                 min(s_end.x, s_end.y) - max(s_start.x, s_start.y) + 1,
#                                 0), by = id]
#   hybrids.dt <- hybrids.dt[s_ol <= s_minoverlap]
  
#   # Remove those too far away from xlink position
#   hybrids.dt[q_start.x < xlink_distance | q_start.y < xlink_distance]

#   return(hybrids.dt)
  
# }

GetValidHybrids <- function(blast8.query.dt, min_unmapped_length = 21, q_minoverlap = 4, q_maxgap = 4, s_minoverlap = 0, xlink_distance = 5) {
  
  hybrids.dt <- blast8.query.dt
  
  # Keep best match for a given query region
  hybrids.dt <- hybrids.dt[evalue == min_evalue]  
  
  # Match up with fasta read length and remove if enough of a continuous match
  if(any(hybrids.dt$unmapped < (min_unmapped_length - q_minoverlap) & hybrids.dt$unmapped == 100)) {
    
    return(data.table())
    
  } else {
  
    # Now get combinations
    hybrids.dt[, id := 1:.N]
    hybrids.dt <- merge(hybrids.dt, hybrids.dt, by = c("query"), allow.cartesian = TRUE)
    hybrids.dt <- hybrids.dt[id.y > id.x] # Remove duplicates
    hybrids.dt[, id := paste0(id.x, "_", id.y)]
    
    # Remove those with significant overlap in the query mappings and too large a gap between the query mappings
    hybrids.dt[, q_ol := min(q_end.x, q_end.y) - max(q_start.x, q_start.y) + 1, by = id]
    hybrids.dt <- hybrids.dt[q_ol <= q_minoverlap][q_ol >= -q_maxgap]
    
    # Remove those with significant overlap in the subject mappings, if subjects are the same
    hybrids.dt <- hybrids.dt[, s_ol := ifelse(subject.x == subject.y,
                                              min(s_end.x, s_end.y) - max(s_start.x, s_start.y) + 1,
                                              0), by = id]
    hybrids.dt <- hybrids.dt[s_ol <= s_minoverlap]
    
    # Remove those too far away from xlink position
    hybrids.dt <- hybrids.dt[q_start.x < xlink_distance | q_start.y < xlink_distance]
    
    return(hybrids.dt)
  
  }
  
}

FilterValidHybrids <- function(valid.hybrids.dt) {

  valid.hybrids.dt[, multi := .N, by = query]
  single.hybrids.dt <- valid.hybrids.dt[multi == 1]
  multi.hybrids.dt <- valid.hybrids.dt[multi > 1]
  
  if(nrow(multi.hybrids.dt) > 0) {

  # ==========
  # 1. First calculate all the metrics
  # ==========
    
  multi.hybrids.dt[, id := 1:.N]
  multi.hybrids.dt[, q_length := max(q_end.x, q_end.y) - min(q_start.x, q_start.y) + 1, by = id]
  multi.hybrids.dt[, q_length_exc_gap := max(q_end.x, q_end.y) - min(q_start.x, q_start.y) + 1 + q_ol, by = id]
  multi.hybrids.dt[, total_alignment_length := alignment_length.x + alignment_length.y]
  
  # ==========
  # 2. Select ones that overlap single hits
  # ==========
  
  L.gr <- with(single.hybrids.dt, GRanges(seqnames = subject.x,
                                          ranges = IRanges(start = s_start.x, end = s_end.x),
                                          strand = "+",
                                          evalue = evalue.x,
                                          bitscore = bit_score.x,
                                          read = query))
  
  R.gr <- with(single.hybrids.dt, GRanges(seqnames = subject.y,
                                          ranges = IRanges(start = s_start.y, end = s_end.y),
                                          strand = "+",
                                          evalue = evalue.y,
                                          bitscore = bit_score.y,
                                          read = query))
  
  seqlevels(L.gr) <- unique(c(seqlevels(L.gr), seqlevels(R.gr)))
  seqlevels(R.gr) <- unique(c(seqlevels(L.gr), seqlevels(R.gr)))
  single.gi <- GInteractions(L.gr, R.gr)
  single.gi <- swapAnchors(single.gi)
  
  L.gr <- with(multi.hybrids.dt, GRanges(seqnames = subject.x,
                                         ranges = IRanges(start = s_start.x, end = s_end.x),
                                         strand = "+",
                                         evalue = evalue.x,
                                         bitscore = bit_score.x,
                                         read = query,
                                         id = id))
  
  R.gr <- with(multi.hybrids.dt, GRanges(seqnames = subject.y,
                                         ranges = IRanges(start = s_start.y, end = s_end.y),
                                         strand = "+",
                                         evalue = evalue.y,
                                         bitscore = bit_score.y,
                                         read = query,
                                         id = id))
  
  seqlevels(L.gr) <- unique(c(seqlevels(L.gr), seqlevels(R.gr)))
  seqlevels(R.gr) <- unique(c(seqlevels(L.gr), seqlevels(R.gr)))
  multi.gi <- GInteractions(L.gr, R.gr)
  multi.gi <- swapAnchors(multi.gi)

  seqlevels(single.gi) <- unique(c(seqlevels(single.gi), seqlevels(multi.gi)))
  seqlevels(multi.gi) <- unique(c(seqlevels(single.gi), seqlevels(multi.gi)))  
  multi.gi$ol <- countOverlaps(multi.gi, single.gi, ignore.strand = FALSE, use.region = "both")
  
  multi.ol.dt <- data.table(query = multi.gi$anchor1.read,
                            id = multi.gi$anchor1.id,
                            ol = multi.gi$ol)
  multi.ol.dt <- multi.ol.dt[ol > 0]
  multi.ol.dt[, max_ol := max(ol), by = query]
  multi.ol.dt <- multi.ol.dt[ol == max_ol] # Select one with most overlaps
  
  ol.multi.hybrids.dt <- multi.hybrids.dt[id %in% multi.ol.dt$id]

  ol.multi.hybrids.dt[, max_q := max(q_length_exc_gap), by = query] # select one with most aligned
  ol.multi.hybrids.dt <- ol.multi.hybrids.dt[q_length_exc_gap == max_q]

  ol.multi.hybrids.dt[, max_q_ol := max(q_ol), by = query] # select one with shortest gap
  ol.multi.hybrids.dt <- ol.multi.hybrids.dt[q_ol == max_q_ol]

  ol.multi.hybrids.dt[, multi := .N, by = query]
  ol.multi.hybrids.dt <- ol.multi.hybrids.dt[multi == 1] # only keep ones with now one solution
  ol.multi.hybrids.dt[, h_sel := "multi_overlap"]
  
  # ==========
  # 2. Select ones based on quality of hybrid
  # ==========
  
  qual.multi.hybrids.dt <- multi.hybrids.dt[!query %in% ol.multi.hybrids.dt$query]
  
  qual.multi.hybrids.dt[, max_q := max(q_length_exc_gap), by = query] # select one with most aligned 
  qual.multi.hybrids.dt <- qual.multi.hybrids.dt[q_length_exc_gap == max_q]

  qual.multi.hybrids.dt[, max_q_ol := max(q_ol), by = query] # select one with shortest gap  
  qual.multi.hybrids.dt <- qual.multi.hybrids.dt[q_ol == max_q_ol]

  qual.multi.hybrids.dt[, multi := .N, by = query]
  qual.multi.hybrids.dt <- qual.multi.hybrids.dt[multi == 1]  
  qual.multi.hybrids.dt[, h_sel := "qual"]
  
  # ==========
  # 3. Select ones if intragenic hit
  # ==========

  intra.multi.hybrids.dt <- multi.hybrids.dt[!query %in% c(ol.multi.hybrids.dt$query, qual.multi.hybrids.dt$query)]
  intra.multi.hybrids.dt <- intra.multi.hybrids.dt[subject.x == subject.y]
  intra.multi.hybrids.dt[, multi := .N, by = query]
  intra.multi.hybrids.dt <- intra.multi.hybrids.dt[multi == 1]
  intra.multi.hybrids.dt[, h_sel := "intra"]
  
  } else {

    ol.multi.hybrids.dt <- data.table()
    qual.multi.hybrids.dt <- data.table()
    intra.multi.hybrids.dt <- data.table()

  }

  # ==========
  # 4. Collect together
  # ==========  
  
  single.hybrids.dt[, h_sel := "single"]
  assigned.multi.hybrid.dt <- rbind(single.hybrids.dt, ol.multi.hybrids.dt, qual.multi.hybrids.dt, intra.multi.hybrids.dt, fill = TRUE)
  
  message(nrow(assigned.multi.hybrid.dt), " out of ", length(unique(valid.hybrids.dt$query)), " reads assigned - ",
          round(nrow(assigned.multi.hybrid.dt)/length(unique(valid.hybrids.dt$query)), 2) * 100, "%")

  # ==========
  # 4. Tidy up
  # ==========  
  
  assigned.multi.hybrid.dt <- assigned.multi.hybrid.dt[, .(query, h_sel,
                                                           subject.x, q_start.x, q_end.x, s_start.x, s_end.x, evalue.x, bit_score.x,
                                                           subject.y, q_start.y, q_end.y, s_start.y, s_end.y, evalue.y, bit_score.y)]
  setnames(assigned.multi.hybrid.dt, c("read", "hybrid_selection",
                                       "L_seqnames", "L_read_start", "L_read_end", "L_start", "L_end", "L_eval", "L_bitscore",
                                       "R_seqnames", "R_read_start", "R_read_end", "R_start", "R_end", "R_eval", "R_bitscore"))
  
  return(assigned.multi.hybrid.dt)
  
}

# ==========
# Run
# ==========

ptm <- proc.time()

# Load fasta for read lengths
fasta <- readDNAStringSet(opt$fasta)
fasta.dt <- data.table(query = names(fasta),
                       readlength = width(fasta))
setkey(fasta.dt, query)

# Load blast and get read lengths
blast.dt <- fread(input = opt$blast, col.names = c("query", "subject", "identity", "alignment_length", "mismatches", "gap_openings","q_start","q_end","s_start","s_end","evalue","bit_score"))
setorder(blast.dt, query)
setkey(blast.dt, query)
stopifnot(all(blast.dt$query %in% fasta.dt$query))
blast.dt <- fasta.dt[blast.dt]

# Add calculations
blast.dt[, min_evalue := min(evalue), by = .(query, q_start, q_end)]
blast.dt[, q_length := q_end - q_start + 1]
blast.dt[, `:=` (unmapped = readlength - q_length,
                 mapped = q_length/readlength)]

# Split for parallel processing
blast.list <- split(blast.dt, blast.dt$query)
rm(fasta)
rm(fasta.dt)
rm(blast.dt)
invisible(gc())

# Get valid hybrids
# cl <- makeCluster(as.integer(opt$threads))
# clusterExport(cl = cl, varlist = c("GetValidHybrids", ":="))
cl <- makeForkCluster(opt$threads) # otherwise really slow...
hybrids.list <- pblapply(cl = cl, blast.list, function(x) GetValidHybrids(blast8.query.dt = x))
stopCluster(cl)

message(sum(elementNROWS(hybrids.list) == 0), " out of ", length(hybrids.list), " reads did not have hybrids")
message(round(sum(elementNROWS(hybrids.list) != 0)/length(hybrids.list), 4) * 100, "% of reads had hybrids")
message()

# Filter multi hits
hybrids.dt <- rbindlist(hybrids.list)
filtered.hybrids.dt <- FilterValidHybrids(valid.hybrids.dt = hybrids.dt)
filtered.hybrids.dt[, .N, by = hybrid_selection]

# Write out
fwrite(filtered.hybrids.dt, file = opt$output, sep = "\t", col.names = TRUE)

message("Completed in: ", round((proc.time() - ptm)[3]/60), " minute(s)")