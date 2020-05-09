# Plot hybrid QC metrics
# A. M. Chakrabarti
# 9th May 2020

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(scales))
suppressPackageStartupMessages(library(ggthemes))
suppressPackageStartupMessages(library(cowplot))

# ==========
# Orientation
# ==========

message("Plotting hybrid read orientation...")

hybrid.files <- snakemake@input$hybrids

orientation.list <- lapply(1:length(hybrid.files), function(i) {
  
  hybrid.dt <- fread(hybrid.files[i])
  orientation.dt <- hybrid.dt[, .N, by = orientation]
  orientation.dt$exp <- sapply(strsplit(basename(hybrid.files[i]), "\\."), "[[", 1)
  
  return(orientation.dt)
  
})

orientation.dt <- rbindlist(orientation.list)

p1 <- ggplot(orientation.dt, aes(x = exp, y = N, fill = orientation)) +
  geom_bar(stat = "identity") +
  scale_fill_tableau() +
  scale_y_continuous(labels = comma) +
  labs(title = "Hybrid read orientation",
       x = "",
       y = "Count",
       fill = "") +
  coord_flip() +
  theme_cowplot() + theme(legend.position = "bottom")

p2 <- ggplot(orientation.dt, aes(x = exp, y = N, fill = orientation)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_tableau() +
  scale_y_continuous(labels = comma) +
  labs(title = "Hybrid read orientation",
       x = "",
       y = "Percentage",
       fill = "") +
  coord_flip() +
  theme_cowplot() + theme(legend.position = "bottom")

ggsave(p1, filename = snakemake@output[[1]], width = 297, height = 210, units = "mm")
ggsave(p2, filename = snakemake@output[[2]], width = 297, height = 210, units = "mm")

# ==========
# Read counts
# ==========

message("Plotting hybrid recovery...")

total.list <- lapply(1:length(hybrid.files), function(i) {
  
  star.out <- gsub("hybrids/", "preprocessed/", hybrid.files[i])
  star.out <- gsub("hybrids.tsv.gz$", "Log.final.out", star.out)
  
  reads.dt <- data.table(exp = sapply(strsplit(basename(hybrid.files[i]), "\\."), "[[", 1),
                         total = as.integer(sapply(strsplit(readLines(star.out)[6], "\t"), "[[", 2)))
  
  return(reads.dt)
  
})

total.dt <- rbindlist(total.list)
total.dt <- merge(total.dt, orientation.dt[, sum(N), by = exp], by = "exp")
setnames(total.dt, "V1", "hybrids")
total.dt[, lost := total- hybrids, by = exp]
total.dt <- melt.data.table(total.dt[, .(exp, hybrids, lost)], id.vars = "exp")

p1 <- ggplot(total.dt, aes(x = exp, y = value, fill = variable)) +
  geom_bar(stat = "identity") +
  scale_fill_tableau() +
  scale_y_continuous(labels = comma) +
  labs(title = "Hybrid recovery",
       x = "",
       y = "Count",
       fill = "") +
  coord_flip() +
  theme_cowplot() + theme(legend.position = "bottom")

p2 <- ggplot(total.dt, aes(x = exp, y = value, fill = variable)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_tableau() +
  scale_y_continuous(labels = comma) +
  labs(title = "Hybrid recovery",
       x = "",
       y = "Percentage",
       fill = "") +
  coord_flip() +
  theme_cowplot() + theme(legend.position = "bottom")

ggsave(p1, filename = snakemake@output[[3]], width = 297, height = 210, units = "mm")
ggsave(p2, filename = snakemake@output[[4]], width = 297, height = 210, units = "mm")

# ==========
# SJ motifs
# ==========

message("Plotting hybrid SJ motifs...")

# All mRNA

sj.list <- lapply(1:length(hybrid.files), function(i) {
  
  hybrid.dt <- fread(hybrid.files[i])
  # hybrid.dt <- hybrid.dt[L_seqnames == R_seqnames][grep("ENSMUSG", L_seqnames)]
  sj.dt <- hybrid.dt[, .N, by = sj]
  setorder(sj.dt, -N)
  sj.dt$exp <- sapply(strsplit(basename(hybrid.files[i]), "\\."), "[[", 1)
  
  return(sj.dt)
  
})

sj.dt <- rbindlist(sj.list)
sj.dt <- sj.dt[nchar(sj) == 4]
sj.dt[, z := scale(N), by = exp]

p1 <- ggplot(sj.dt, aes(x = z, colour = exp)) +
  geom_density() +
  scale_colour_tableau(palette = "Tableau 20") +
  labs(title = "Hybrid SJ motifs",
       x = "Motif z-score",
       y = "Density",
       colour = "") +
  facet_wrap(~ exp, scales = "free") +
  theme_cowplot() + theme(legend.position = "none")

# Intragenic mRNA

sj.list <- lapply(1:length(hybrid.files), function(i) {
  
  hybrid.dt <- fread(hybrid.files[i])
  hybrid.dt <- hybrid.dt[L_seqnames == R_seqnames][grep("ENSMUSG", L_seqnames)]
  sj.dt <- hybrid.dt[, .N, by = sj]
  setorder(sj.dt, -N)
  sj.dt$exp <- sapply(strsplit(basename(hybrid.files[i]), "\\."), "[[", 1)
  
  return(sj.dt)
  
})

sj.dt <- rbindlist(sj.list)
sj.dt <- sj.dt[nchar(sj) == 4]
sj.dt[, z := scale(N), by = exp]

p2 <- ggplot(sj.dt, aes(x = z, colour = exp)) +
  geom_density() +
  scale_colour_tableau(palette = "Tableau 20") +
  labs(title = "Intragenic mRNA hybrid SJ motifs",
       x = "Motif z-score",
       y = "Density",
       colour = "") +
  facet_wrap(~ exp, scales = "free") +
  theme_cowplot() + theme(legend.position = "none")

ggsave(p1, filename = snakemake@output[[5]], width = 297, height = 210, units = "mm")
ggsave(p2, filename = snakemake@output[[6]], width = 297, height = 210, units = "mm")

# ==========
# Binding energy
# ==========

message("Plotting cluster binding energy...")

cluster.files <- snakemake@input$clusters

mfe.list <- lapply(1:length(cluster.files), function(i) {
  
  cluster.dt <- fread(cluster.files[i])
  cluster.dt$exp <- sapply(strsplit(basename(cluster.files[i]), "\\."), "[[", 1)
  
  return(cluster.dt[, .(exp, mfe)])
  
})

mfe.dt <- rbindlist(mfe.list)

p <- ggplot(mfe.dt, aes(x = mfe, colour = exp)) +
  geom_density() +
  scale_colour_tableau(palette = "Tableau 20") +
  labs(title = "Cluster binding energy",
       x = "Mean binding energy",
       y = "Density",
       colour = "") +
  # facet_wrap(~ exp, scales = "free") +
  theme_cowplot() + theme(legend.position = "right")

ggsave(p, filename = snakemake@output[[7]], width = 297, height = 210, units = "mm")
