FROM mcr.microsoft.com/devcontainers/base:ubuntu

USER root

# Replace shell with bash so we can source files
SHELL [ "/bin/bash", "-c" ]
ENV SHELL=/bin/bash
RUN rm /bin/sh && ln -s /bin/bash /bin/sh
ENV CONDA_DIR=/opt/conda
ENV PATH=${CONDA_DIR}/bin:${PATH}

# Install util tools.
# software-properties-common is needed to add ppa support for Apptainer installation
RUN apt-get update --quiet && \
    export DEBIAN_FRONTEND=noninteractive && \
    apt-get install --quiet --yes \
        apt-transport-https \
        apt-utils \
        bzip2 \
        ca-certificates \
        curl \
        default-jre \
        git \
        graphviz \
        less \
        python3 \
        sudo \
        tini\
        tree \
        wget \
        software-properties-common && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Install Apptainer (Singularity)
RUN add-apt-repository -y ppa:apptainer/ppa && \
    apt-get update --quiet && \
    apt install -y apptainer && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Install Seqera Platform "Tower Agent"
RUN curl -fSL https://github.com/seqeralabs/tower-agent/releases/latest/download/tw-agent-linux-x86_64 > tw-agent && \
    chmod +x tw-agent && \
    mv tw-agent /usr/local/bin/tw-agent

# Set up directories
RUN mkdir -p /workspaces/.nextflow && \
    mkdir -p /workspaces/training/

# Uncomment if we need to pin the Nextflow version
ENV NXF_EDGE=0
ENV NXF_VER=24.10.4
ENV NXF_HOME=/workspaces/.nextflow

# Install Conda
RUN wget --no-hsts -q -O Miniforge3.sh "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh" && \
    bash Miniforge3.sh -b -p "${CONDA_DIR}"
#     && \
#    rm Miniforge3.sh && \
#    conda clean --tarballs --index-cache --packages --yes && \
#    find "${CONDA_DIR}" -follow -type f -name '*.a' -delete && \
#    find "${CONDA_DIR}" -follow -type f -name '*.pyc' -delete && \
#    conda clean --force-pkgs-dirs --all --yes  && \
#    echo ". ${CONDA_DIR}/etc/profile.d/conda.sh && conda activate base" >> /etc/skel/.bashrc && \
#    echo ". ${CONDA_DIR}/etc/profile.d/conda.sh && conda activate base" >> ~/.bashrc && \
#    conda init && \
#    conda activate

# Install nextflow, nf-core, Mamba, and pytest-workflow
RUN conda config --add channels bioconda && \
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

# Cleanup
RUN printf 'unset JAVA_TOOL_OPTIONS\n' >> $HOME/.bashrc

# Update Nextflow
RUN nextflow self-update && nextflow -version

RUN unset JAVA_TOOL_OPTIONS

# RUN export PS1='\t -> '
RUN export PS1='\[\e[3;36m\]${PWD/*//} ->\[\e[0m\] '
