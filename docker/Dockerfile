FROM nfcore/base:1.12.1
LABEL authors="Anob Chakrabarti" \
      description="Docker image containing all software requirements for Tosca"

# Install the conda environment
COPY environment.yml /
RUN conda env create --quiet -f /environment.yml && conda clean -a

# Add conda installation dir to PATH (instead of doing 'conda activate')
ENV PATH /opt/conda/envs/tosca-0.2/bin:$PATH

# Dump the details of the installed packages to a file for posterity
RUN conda env export --name tosca-0.2 > tosca-0.2.yml

# Install R packages
RUN R -e 'install.packages(c("devtools", "data.table", "R.utils", "ggplot2", "scales", "tictoc", "optparse", "parallel", "rslurm", "cowplot", "ggthemes", "BiocManager"), repos = "https://cloud.r-project.org/")'
RUN R -e 'BiocManager::install(c("GenomicAlignments", "rtracklayer", "graph", "RBGL"))'
RUN R -e 'devtools::install_github("amchakra/primavera")'

# Instruct R processes to use these empty files instead of clashing with a local version
RUN touch .Rprofile
RUN touch .Renviron