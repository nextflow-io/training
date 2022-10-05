FROM gitpod/workspace-base

USER root

# Install util tools.
RUN apt-get update \
 && apt-get install -y \
  apt-utils \
  sudo \
  git \
  less \
  wget \
  tree \
  graphviz 

RUN mkdir -p /workspace/data \
    && chown -R gitpod:gitpod /workspace/data
  
# Install Conda
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh && \
    bash Miniconda3-latest-Linux-x86_64.sh -b -p /opt/conda && \
    rm Miniconda3-latest-Linux-x86_64.sh

ENV PATH="/opt/conda/bin:$PATH"

RUN chown -R gitpod:gitpod /opt/conda

USER gitpod

# Install nextflow, nf-core, Mamba, and pytest-workflow
RUN conda update -n base -c defaults conda && \
    conda install nextflow nf-core pytest-workflow mamba -n base -c conda-forge -c bioconda && \
    nextflow self-update && \
    conda clean --all -f -y

RUN export NXF_EDGE=1

ENV NXF_VER=22.09.7-edge

RUN nextflow self-update

RUN unset JAVA_TOOL_OPTIONS

RUN export PS1='\t -> '


