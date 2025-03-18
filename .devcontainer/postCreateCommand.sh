#! /usr/bin/env bash

# Set safe working directory
git config --global --add safe.directory /workspaces/training

# Set conda channels
conda config --remove channels defaults && \
    conda config --add channels bioconda && \
    conda config --add channels conda-forge && \
    conda config --set channel_priority strict && \
    conda install --quiet --yes --name base \
        nf-core \
        nf-test \
        black \
        prettier \
        pre-commit \
        linkify-it-py \
        pytest-workflow && \
    conda clean --all --force-pkgs-dirs --yes

# Install tw agent
curl -fSL https://github.com/seqeralabs/tower-agent/releases/latest/download/tw-agent-linux-x86_64 > tw-agent && \
    chmod +x tw-agent && \
    mv tw-agent /usr/local/bin/tw-agent

# Set up directories
mkdir -p /workspaces/.nextflow && \
    mkdir -p /workspaces/training/

# Install Apptainer (Singularity)
add-apt-repository -y ppa:apptainer/ppa && \
    apt-get update --quiet && \
    apt install -y apptainer && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Output Env Details
nextflow -version; if [ -z \"$CODESPACES\" ]; then echo \"Devcontainers Development\"; else echo \"Codespaces Development\";  fi
