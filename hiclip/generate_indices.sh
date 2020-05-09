sbatch -c 8 -t 12:00:00 --mem=64G --wrap="STAR --runMode genomeGenerate --runThreadN 8 --genomeDir Mm_GencodeM24_rRNA_MT_genes --genomeFastaFiles Mm_GencodeM24_rRNA_MT_genes.fa --limitGenomeGenerateRAM 39550915285"

sbatch -c 8 -t 12:00:00 --wrap="bowtie2-build --threads 8 Mm_GencodeM24_rRNA_MT_genes.fa Mm_GencodeM24_rRNA_MT_genes"