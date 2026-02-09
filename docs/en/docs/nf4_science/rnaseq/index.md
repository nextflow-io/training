---
title: Nextflow for RNAseq
hide:
  - toc
page_type: index_page
index_type: course
additional_information:
  technical_requirements: true
  learning_objectives:
    - Write a linear workflow to apply basic RNAseq processing and QC methods
    - Handle domain-specific files such as FASTQ and reference genome resources appropriately
    - Handle single-end and paired-end sequencing data
    - Leverage Nextflow's dataflow paradigm to parallelize per-sample RNAseq processing
    - Aggregate QC reports across multiple steps and samples using relevant channel operators
  audience_prerequisites:
    - "**Audience:** This course is designed for researchers in transcriptomics and related fields who want to develop or customize data analysis pipelines."
    - "**Skills:** Some familiarity with the command line, basic scripting concepts, and common RNAseq file formats is assumed."
    - "**Prerequisites:** Foundational Nextflow concepts and tooling covered in [Hello Nextflow](../../hello_nextflow/)."
---

# Nextflow for RNAseq

**A hands-on course applying Nextflow to a real-world transcriptomics use case: bulk RNAseq processing with Trim Galore, HISAT2 and FastQC.**

This course builds on the [Hello Nextflow](../../hello_nextflow/) beginner training and demonstrates how to use Nextflow in the specific context of bulk RNAseq analysis.
You will implement a processing pipeline that trims adapter sequences, aligns reads to a genome reference and performs quality control (QC) at several stages.

<!-- additional_information -->

## Course overview

This course is hands-on, with goal-oriented exercises structured to introduce information gradually.

You will start by running the processing tools manually in the terminal to understand the methodology, then progressively build up a Nextflow pipeline that automates and scales the analysis.

### Lesson plan

We've broken this down into three parts that each focus on specific aspects of applying Nextflow to an RNAseq use case.

| Course chapter                                                         | Summary                                                                                                 | Estimated duration |
| ---------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------- | ------------------ |
| [Part 1: Method overview](./01_method.md)                              | Understanding the RNAseq processing methodology and running the tools manually                          | 30 mins            |
| [Part 2: Single-sample implementation](./02_single-sample.md)          | Building a pipeline that trims, aligns and QCs a single sample, then scaling to handle multiple samples | 60 mins            |
| [Part 3: Multi-sample paired-end implementation](./03_multi-sample.md) | Extending the pipeline to handle paired-end data and aggregate QC reports across samples                | 45 mins            |

By the end of this course, you will be able to apply foundational Nextflow concepts and tooling to a typical RNAseq use case.

Ready to take the course?

[Get started :material-arrow-right:](00_orientation.md){ .md-button .md-button--primary }

<!-- Clearfix for float -->
<div style="content: ''; clear: both; display: table;"></div>
