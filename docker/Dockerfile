FROM nfcore/base:1.12.1
LABEL authors="Anob Chakrabarti" \
      description="Docker image containing all software requirements for Tosca"

# Install the conda environment
COPY environment.yml /
RUN conda env create --quiet -f /environment.yml && conda clean -a

# Add conda installation dir to PATH (instead of doing 'conda activate')
ENV PATH /opt/conda/envs/tosca-1.0.0/bin:$PATH

# Dump the details of the installed packages to a file for posterity
RUN conda env export --name tosca-1.0.0 > tosca-1.0.0.yml

# Install R packages
RUN R -e 'devtools::install_github("amchakra/toscatools")'

# Instruct R processes to use these empty files instead of clashing with a local version
RUN touch .Rprofile
RUN touch .Renviron