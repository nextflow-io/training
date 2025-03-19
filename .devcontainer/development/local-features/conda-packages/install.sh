#!/usr/bin/env bash

# Install conda packages needed for Nextflow Training

conda config --add channels bioconda
conda config --add channels conda-forge
conda config --set channel_priority strict
conda update --quiet --yes --all

conda install --quiet --yes --name base \
    nextflow \
    nf-core \
    nf-test \
    pre-commit \
    linkify-it-py \
    pytest-workflow

conda clean --all --force-pkgs-dirs --yes
