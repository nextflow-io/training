---
title: Nextflow for Genomics
hide:
    - toc
---

# Nextflow for Genomics

This training course is intended for researchers in genomics and related fields who are interested in developing or customizing data analysis pipelines.
It builds on the [Hello Nextflow](../hello_nextflow/index.md) beginner training and demonstrates how to use Nextflow in the specific context of the genomics domain.

## Learning objectives

By working through this course, you will learn how to apply foundational Nextflow concepts and tooling to a typical genomics use case.

By the end of this workshop you will be able to:

-   Write a linear workflow to apply variant calling to a single sample
-   Handle accessory files such as index files and reference genome resources appropriately
-   Leverage Nextflow's dataflow paradigm to parallelize per-sample variant calling
-   Implement multi-sample variant calling using relevant channel operators
-   Configure pipeline execution and manage and optimize resource allocations
-   Implement per-step and end-to-end pipeline tests that handle genomics-specific idiosyncrasies appropriately

<!-- TODO for future expansion: add metadata/samplesheet handling -->

## Prerequisites

The course assumes some minimal familiarity with the following:

-   Tools and file formats commonly used in this scientific domain
-   Experience with the command line
-   Foundational Nextflow concepts and tooling covered in the [Hello Nextflow](../hello_nextflow/index.md) beginner training.

For technical requirements and environment setup, see the [Environment Setup](../envsetup/index.md) mini-course.

## Get started

To get started, open the training environment by clicking the 'Open in Gitpod' button below.

[![Open in Gitpod](https://img.shields.io/badge/Gitpod-%20Open%20in%20Gitpod-908a85?logo=gitpod)](https://gitpod.io/#https://github.com/nextflow-io/training)
