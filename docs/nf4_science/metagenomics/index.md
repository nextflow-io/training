---
title: Nextflow for Metagenomics
hide:
  - toc
---

# Nextflow for Metagenomics

!!! note

    This course is a community effort led by [Jeferyd Yepes](https://jeferydyepes.com/) and edited by the Nextflow Training team. It will be re-visited every now and then for updates and suggestions from the community members.

The course is designed for researchers on focused metagenomics (WGS/shotgun) data analysis who are interested in developing or customizing taxonomic annotation pipelines.
It builds on the [Hello Nextflow](../../hello_nextflow/) and [Nextflow for RNAseq](../rnaseq/) beginner training and demonstrates how to use Nextflow in the specific context of metagenomics data analysis.

Specifically, this course demonstrates how to implement a simple read taxonomic annotation starting from removing host sequences, passing through re-estimating species abundance with Bayesian statistics, until generating complete reports.

Let's get started! Click on the "Open in GitHub Codespaces" button below to launch the training environment (preferably in a separate tab), then read on while it loads.

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

## Learning objectives

By the end of this course, you will have learnt how to apply foundational Nextflow concepts and tooling to a typical metagenomics use case.

Concretely, you will be able to:

- Write a linear workflow to perform host removal, taxonomic annotation and species abundance re-estimation
- Handle domain-specific files such as Kraken2 and Bracken reports resources appropriately
- Run analysis for a single sample or leverage on Nextflow's dataflow paradigm to parallelize multi-sample analysis
- Separate the processes and workflow in a more structured manner attempting to a first step in following [nf-core](https://nf-co.re/) guidelines
- Use conditionals and operators to control workflow execution
- Include custom scripts to be run within a given process

## Prerequisites

The course assumes some minimal familiarity with the following:

- Tools and file formats commonly used in this scientific domain
- Experience with the command line
- Foundational Nextflow concepts and tooling covered in the [Hello Nextflow](../../hello_nextflow/) and [Nextflow for RNAseq](../rnaseq/) beginner training.

For technical requirements and environment setup, see the [Environment Setup](../../envsetup/) mini-course.
