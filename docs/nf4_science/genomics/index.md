---
title: Nextflow for Genomics
hide:
  - toc
---

# Nextflow for Genomics

This training course is intended for researchers in genomics and related fields who are interested in developing or customizing data analysis pipelines.
It builds on the [Hello Nextflow](../hello_nextflow/) beginner training and demonstrates how to use Nextflow in the specific context of the genomics domain.

Specifically, this course demonstrates how to implement a simple variant calling pipeline with [GATK](https://gatk.broadinstitute.org/) (Genome Analysis Toolkit), a widely used software package for analyzing high-throughput sequencing data.

!!! note

    Don't worry if you're not familiar with GATK specifically.
    We'll summarize the necessary concepts as we go, and the workflow implementation principles we demonstrate here apply broadly to any command line tool that processes genomics data.

## Learning objectives

By working through this course, you will learn how to apply foundational Nextflow concepts and tooling to a typical genomics use case.

By the end of this workshop you will be able to:

- Write a linear workflow to apply variant calling to a single sample
- Handle accessory files such as index files and reference genome resources appropriately
- Leverage Nextflow's dataflow paradigm to parallelize per-sample variant calling
- Implement multi-sample variant calling using relevant channel operators
- Configure pipeline execution and manage and optimize resource allocations
- Implement per-step and end-to-end pipeline tests that handle genomics-specific idiosyncrasies appropriately

<!-- TODO for future expansion: add metadata/samplesheet handling -->

## Prerequisites

The course assumes some minimal familiarity with the following:

- Tools and file formats commonly used in this scientific domain
- Experience with the command line
- Foundational Nextflow concepts and tooling covered in the [Hello Nextflow](../hello_nextflow/) beginner training.

For technical requirements and environment setup, see the [Environment Setup](../envsetup/) mini-course.

## Get started

To get started, open the training environment by clicking the 'Open in GitHub Codespaces' button below.

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)
