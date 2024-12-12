FROM gitpod/workspace-base

USER root

# Install util tools.
# software-properties-common is needed to add ppa support for Apptainer installation
RUN apt-get update --quiet && \
    apt-get install --quiet --yes \
        apt-transport-https \
        apt-utils \
        sudo \
        git \
        less \
        wget \
        curl \
        tree \
        graphviz \
        software-properties-common && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Taken from: https://github.com/nf-core/tools/blob/master/nf_core/gitpod/gitpod.Dockerfile
# Install Apptainer (Singularity)
RUN add-apt-repository -y ppa:apptainer/ppa && \
    apt-get update --quiet && \
    apt install -y apptainer && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Install Conda
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh && \
    bash Miniconda3-latest-Linux-x86_64.sh -b -p /opt/conda && \
    rm Miniconda3-latest-Linux-x86_64.sh

ENV PATH="/opt/conda/bin:$PATH"

# User permissions
RUN mkdir -p /workspace/data \
    && chown -R gitpod:gitpod /opt/conda /workspace/data

# Install Tower Agent
RUN curl -fSL https://github.com/seqeralabs/tower-agent/releases/latest/download/tw-agent-linux-x86_64 > tw-agent && \
    chmod +x tw-agent && \
    mv tw-agent /usr/local/bin/tw-agent

# Change user to gitpod
USER gitpod

# Uncomment if we need to pin the Nextflow version
ENV NXF_EDGE=0
ENV NXF_VER=24.10.0

# Install nextflow, nf-core, Mamba, and pytest-workflow
RUN conda config --add channels defaults && \
    conda config --add channels bioconda && \
    conda config --add channels conda-forge && \
    conda config --set channel_priority strict && \
    conda update --quiet --yes --all && \
    conda install --quiet --yes --name base \
        mamba \
        nextflow \
        nf-core \
        nf-test \
        black \
        prettier \
        pre-commit \
        linkify-it-py \
        pytest-workflow && \
    conda clean --all --force-pkgs-dirs --yes

# Update Nextflow
RUN nextflow self-update && nextflow -version

RUN unset JAVA_TOOL_OPTIONS

RUN export PS1='\t -> '
