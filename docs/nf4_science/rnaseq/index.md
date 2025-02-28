---
title: Nextflow for RNAseq
hide:
  - toc
---

# Nextflow for RNAseq

This training course is intended for researchers in transcriptomics and related fields who are interested in developing or customizing data analysis pipelines.
It builds on the [Hello Nextflow](../../hello_nextflow/) beginner training and demonstrates how to use Nextflow in the specific context of bulk RNAseq analysis.

Specifically, this course demonstrates how to implement a simple bulk RNAseq processing pipeline to trim adapter sequences, align the reads to a genome reference and performs quality control (QC) at several stages.

Let's get started! Click on the "Open in GitHub Codespaces" button below to launch the training environment (preferably in a separate tab), then read on while it loads.

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

## Learning objectives

By working through this course, you will learn how to apply foundational Nextflow concepts and tooling to a typical RNAseq use case.

By the end of this workshop you will be able to:

- Write a linear workflow to apply basic RNAseq processing and QC methods
- Handle domain-specific files such as FASTQ and reference genome resources appropriately
- Handle single-end and paired-end sequencing data
- Leverage Nextflow's dataflow paradigm to parallelize per-sample RNAseq processing
- Aggregate QC reports across multiple steps and samples using relevant channel operators

<!-- TODO
- Configure pipeline execution and manage and optimize resource allocations
- Implement per-step and end-to-end pipeline tests that handle RNAseq-specific idiosyncrasies appropriately
-->
<!-- TODO for future expansion: add metadata/samplesheet handling -->

## Prerequisites

The course assumes some minimal familiarity with the following:

- Tools and file formats commonly used in this scientific domain
- Experience with the command line
- Foundational Nextflow concepts and tooling covered in the [Hello Nextflow](../../hello_nextflow/) beginner training.

For technical requirements and environment setup, see the [Environment Setup](../../envsetup/) mini-course.
