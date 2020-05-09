# Genome FASTA to csv script
# A. M. Chakrabarti
# 9th May 2020

library(data.table)

ensg <- Biostrings::readDNAStringSet("/home/camp/chakraa2/working/nobby/projects/flora/mouse/ref/Mm_GencodeM24_rRNA_MT_genes.fa")
ensg.dt <- data.table(gene_id = names(ensg),
                      sequence = as.character(ensg))

fwrite(ensg.dt, "/home/camp/chakraa2/working/nobby/projects/flora/mouse/ref/Mm_GencodeM24_rRNA_MT_genes.csv.gz")