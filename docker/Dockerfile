FROM continuumio/miniconda3:4.8.2
LABEL authors="anob.chakrabarti@crick.ac.uk"

# Install procps so that Nextflow can poll CPU usage
RUN apt-get update && apt-get install -y procps && apt-get clean -y

# Install the conda environment
COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a

# Add conda installation dir to PATH (instead of doing 'conda activate')
ENV PATH /opt/conda/envs/ac-ehiclip-0.1/bin:$PATH

# Dump the details of the installed packages to a file for posterity
RUN conda env export --name ac-ehiclip-0.1 > ac-ehiclip-0.1.yml