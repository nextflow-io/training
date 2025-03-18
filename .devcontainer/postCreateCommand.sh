#! /usr/bin/env bash

# Set up directories
mkdir -p /workspaces/.nextflow && \
    mkdir -p /workspaces/training/

# Set safe working directory
git config --global --add safe.directory /workspaces/training

# Install packages with uv
cat .devcontainer/requirements.txt | xargs -I {} bash -c "uv tool install {}"

# Set conda channels
conda config --add channels bioconda && \
    conda config --set channel_priority strict

# Install tw agent
curl -fSL https://github.com/seqeralabs/tower-agent/releases/latest/download/tw-agent-linux-x86_64 > tw-agent && \
    chmod +x tw-agent && \
    mv tw-agent /usr/local/bin/tw-agent

# Install pre-commit hooks
pre-commit install --install-hooks

# Output Env Details
nextflow -version; if [ -z \"$CODESPACES\" ]; then echo \"Devcontainers Development\"; else echo \"Codespaces Development\";  fi
