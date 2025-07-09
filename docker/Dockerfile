FROM condaforge/miniforge3:25.3.0-3
LABEL authors="Anob Chakrabarti" \
      description="Docker image containing all software requirements for Tosca"

# Install the conda environment
COPY environment.yml /
RUN conda env create --quiet -f /environment.yml && conda clean -a

# Add conda installation dir to PATH (instead of doing 'conda activate')
ENV PATH=/opt/conda/envs/tosca-2.0.0dev/bin:$PATH

# Dump the details of the installed packages to a file for posterity
RUN conda env export --name tosca-2.0.0dev > tosca-2.0.0dev.yml

# Install R packages
RUN R -e 'remotes::install_github("amchakra/toscatools", ref = "dev")'

# Instruct R processes to use these empty files instead of clashing with a local version
RUN touch .Rprofile
RUN touch .Renviron