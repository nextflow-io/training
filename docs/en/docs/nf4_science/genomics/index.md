---
title: Nextflow for Genomics
hide:
  - toc
page_type: index_page
index_type: course
additional_information:
  technical_requirements: true
  learning_objectives:
    - Write a linear workflow to apply variant calling to a single sample
    - Handle accessory files such as index files and reference genome resources appropriately
    - Leverage Nextflow's dataflow paradigm to parallelize per-sample variant calling
    - Implement multi-sample joint calling using relevant channel operators
  audience_prerequisites:
    - "**Audience:** This course is designed for researchers in genomics and related fields who want to develop or customize data analysis pipelines."
    - "**Skills:** Some familiarity with the command line, basic scripting concepts, and common genomics file formats is assumed."
    - "**Prerequisites:** Foundational Nextflow concepts and tooling covered in [Hello Nextflow](../../hello_nextflow/)."
---

# Nextflow for Genomics

**A hands-on course applying Nextflow to a real-world genomics use case: variant calling with GATK.**

This course builds on the [Hello Nextflow](../../hello_nextflow/) beginner training and demonstrates how to use Nextflow in the specific context of the genomics domain.
You will implement a variant calling pipeline with [GATK](https://gatk.broadinstitute.org/) (Genome Analysis Toolkit), a widely used software package for analyzing high-throughput sequencing data.

<!-- additional_information -->

## Course overview

This course is hands-on, with goal-oriented exercises structured to introduce information gradually.

You will start by running the variant calling tools manually in the terminal to understand the methodology, then progressively build up a Nextflow pipeline that automates and scales the analysis.

### Lesson plan

We've broken this down into three parts that each focus on specific aspects of applying Nextflow to a genomics use case.

| Course chapter                                                           | Summary                                                                                         | Estimated duration |
| ------------------------------------------------------------------------ | ----------------------------------------------------------------------------------------------- | ------------------ |
| [Part 1: Method overview](./01_method.md)                                | Understanding the variant calling methodology and running the tools manually                    | 30 mins            |
| [Part 2: Per-sample variant calling](./02_per_sample_variant_calling.md) | Building a pipeline that indexes BAM files and calls variants, then scaling to multiple samples | 60 mins            |
| [Part 3: Joint calling on a cohort](./03_joint_calling.md)               | Adding multi-sample joint genotyping using channel operators to aggregate per-sample outputs    | 45 mins            |

By the end of this course, you will be able to apply foundational Nextflow concepts and tooling to a typical genomics use case.

Ready to take the course?

[Get started :material-arrow-right:](00_orientation.md){ .md-button .md-button--primary }

<!-- Clearfix for float -->
<div style="content: ''; clear: both; display: table;"></div>
