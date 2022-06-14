FROM gitpod/workspace-full:latest

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
  
RUN mkdir /home/gitpod/.conda
# Install conda
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh && \
    /bin/bash ~/miniconda.sh -b -p /opt/conda && \
    rm ~/miniconda.sh && \
    ln -s /opt/conda/etc/profile.d/conda.sh /etc/profile.d/conda.sh
    
RUN chown -R gitpod:gitpod /opt/conda \
    && chmod -R 777 /opt/conda \
    && chown -R gitpod:gitpod /home/gitpod/.conda \
    && chmod -R 777 /home/gitpod/.conda

# unset JAVA_TOOL_OPTIONS
# cd nf-training
# curl -s https://get.nextflow.io | bash
# chmod +x nextflow
# sudo mv nextflow /usr/local/bin/
# docker pull nextflow/rnaseq-nf
# alias conda_activate=". /opt/conda/etc/profile.d/conda.sh; conda activate base"

# Give back control
USER root

# Cleaning
RUN apt-get clean
