#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl=2

process getcontactmap {

    tag "${sample_id}"
    publishDir "${params.outdir}/maps", mode: 'copy', overwrite: true

    time '24h'
    memory '32 G'

    input:
        tuple val(sample_id), path(hybrids)

    output:
        tuple val(sample_id), path("${sample_id}.hybrid.mat.rds"), emit: rds
        tuple val(sample_id), path("${sample_id}.binned.mat.tsv"), emit: tsv

    script:
    
    genomesize = "$params.genomesize"

    """
    #!/usr/bin/env Rscript

    suppressPackageStartupMessages(library(data.table))

    # Functions

    bin_matrix <- function(mat, bin.size) {
    
    mat.size <- nrow(mat)/bin.size
    bin.mat <- matrix(data = 0, nrow = mat.size, ncol = mat.size)
    
    for(i in seq_len(mat.size)) {
        
        for(j in seq_len(mat.size)) {
        
        bin.mat[i, j] <- sum(mat[((i*bin.size) - bin.size + 1):(i*bin.size), ((j*bin.size) - bin.size + 1):(j*bin.size)]) 
        
        }
        
    }
    
    return(bin.mat)
    
    }

    # Get matrix

    genome.size <- as.numeric("$genomesize")

    message(basename(hybrid.file))
    hybrid.dt <- fread(hybrid.file)
    print(hybrid.dt[, .N, by = hybrid_selection])
    hybrid.dt <- hybrid.dt[hybrid_selection == "single"]
    
    hybrid.dt <- ReorientHybrids(hybrid.dt)
    
    mat <- matrix(data = 0, nrow = genome.size, ncol = genome.size)
    
    for(i in 1:nrow(hybrid.dt)) {
        
        if(i%%1e5 == 0) message(i)
        mat[hybrid.dt[i]\$L_start:hybrid.dt[i]\$L_end,
            hybrid.dt[i]\$R_start:hybrid.dt[i]\$R_end] <- mat[hybrid.dt[i]\$L_start:hybrid.dt[i]\$L_end,
                                                            hybrid.dt[i]\$R_start:hybrid.dt[i]\$R_end] + 1
        
    }
    
    saveRDS(mat, file = paste0("$sample_id", ".hybrid.mat.rds"))

    # Bin matrix
    binned.mat <- bin_matrix(rep1.mat, bin.size = 100)
    binned.dt <- data.table(reshape2::melt(binned.mat))

    binned.dt[, norm_value := value*1e6/nrow(hybrid.dt)]
    binned.dt <- binned.dt[value != 0]
    fwrite(binned.dt, file = paste0("$sample_id", ".binned.mat.tsv"), sep = "\t")

    """

}