#!bin/bash
# Script to create conda environment and install R packages

conda create -n clip -c bioconda -c conda-forge umi_tools star fastqc r-base bowtie2 cutadapt pysam samtools snakemake bedtools

conda activate clip

R -e 'install.packages(c("devtools", "data.table", "R.utils", "ggplot2", "scales", "tictoc", "optparse", "rslurm", "cowplot", "ggthemes", "BiocManager"), repos = "https://cloud.r-project.org/")'
R -e 'BiocManager::install(c("GenomicAlignments", "rtracklayer", "graph", "RBGL"))'
R -e 'devtools::install_github("amchakra/primavera")'
