FROM debian:stretch-slim

USER root

# Install util tools.
RUN apt-get update \
 && apt-get install -y \
  apt-transport-https \
  apt-utils \
  sudo \
  git \
  less \
  wget \
  curl \
  tree \
  graphviz


RUN useradd -ms /bin/bash gitpod

RUN mkdir -p /workspace/data \
    && chown -R gitpod:gitpod /workspace/data


# Install Conda
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh && \
    bash Miniconda3-latest-Linux-x86_64.sh -b -p /opt/conda && \
    rm Miniconda3-latest-Linux-x86_64.sh

ENV PATH="/opt/conda/bin:$PATH"

RUN chown -R gitpod:gitpod /opt/conda

USER gitpod

#WORKDIR /home/gitpod

RUN export NXF_EDGE=1

ENV NXF_VER=22.09.7-edge

# Install nextflow, nf-core, Mamba, and pytest-workflow
RUN conda update -n base -c defaults conda && \
    conda install nextflow nf-core pytest-workflow mamba -n base -c conda-forge -c bioconda && \
    conda clean --all -f -y

RUN chown -R gitpod:gitpod /opt/conda

RUN unset JAVA_TOOL_OPTIONS

RUN export PS1='\t -> '
