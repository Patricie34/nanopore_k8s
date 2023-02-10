# Set the base image to Ubuntu 18.04 and NVIDIA GPU
FROM genomicpariscentre/guppy-gpu:latest
ENV PATH "/opt/ont/guppy/bin:$PATH"

FROM mambaorg/micromamba:1.3.0
COPY --chown=$MAMBA_USER:$MAMBA_USER environment.yml /tmp/env.yaml
#overcome "  Command 'ps' required by nextflow to collect task metrics cannot be found"
RUN apt-get update && apt-get install -y procps && rm -rf /var/lib/apt/lists/* 
RUN micromamba install -y -n base -f /tmp/env.yaml && \
    micromamba clean --all --yes
# # Install miniconda
# ENV CONDA_DIR /opt/conda
# RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh && \
#      /bin/bash ~/miniconda.sh -b -p /opt/conda

# # Put conda in path so we can use conda activate
# ENV PATH=$CONDA_DIR/bin:$PATH
# # Install the conda environment
# COPY environment.yml /
# RUN conda env create -f /environment.yml && conda clean --all -f -y

# # Add conda installation dir to PATH (instead of doing 'conda activate')
# ENV PATH /opt/conda/envs/nanopore_mapping/bin:$PATH

USER 1000