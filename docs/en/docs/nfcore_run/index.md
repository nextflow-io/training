---
title: nf-core Run
hide:
  - toc
page_type: index_page
index_type: course
additional_information:
  technical_requirements: true
  learning_objectives:
    - Find, retrieve, and run nf-core community pipelines
    - Configure pipeline execution using parameters and configuration files
    - Understand how nf-core pipelines validate parameters and input data
    - Run a production-scale pipeline (nf-core/rnaseq) and override its default resource allocations
  audience_prerequisites:
    - "**Audience:** This course is designed for learners who are new to Nextflow and nf-core and want to run existing community pipelines."
    - "**Skills:** Some familiarity with the command line, basic scripting concepts and common file formats is assumed."
    - "**Domain:** The exercises use bioinformatics pipelines, but no prior scientific domain knowledge is required."
---

# nf-core Run

**nf-core Run is a hands-on introduction to finding, running, and configuring nf-core community pipelines.**

Working through practical examples and guided exercises, you will learn to find and retrieve nf-core pipelines, run them using their built-in test profiles, and customize their execution through parameters and configuration files.

You will take away the skills and confidence to start running nf-core pipelines for your own analyses.

<!-- additional_information -->

## Course overview

This course is hands-on, with goal-oriented exercises structured to introduce information gradually.

You will start with `nf-core/demo`, a minimal pipeline maintained by the nf-core project for training purposes, then apply what you've learned to `nf-core/rnaseq`, a widely-used production pipeline for bulk RNA sequencing analysis.

This course focuses on running pipelines.
If you're looking for an intro to developing nf-core-compatible pipelines, see [Hello nf-core](../hello_nf-core/index.md).

### Lesson plan

| Course chapter                                  | Summary                                                                                                                      | Estimated duration |
| ----------------------------------------------- | ---------------------------------------------------------------------------------------------------------------------------- | ------------------ |
| [Part 1: Run a demo pipeline](./01_run_demo.md) | Find and retrieve an nf-core pipeline, run it using its test profile, configure its execution, and run a production pipeline | 60 mins            |

By the end of this course, you will be able to take advantage of the wealth of community pipelines offered by the nf-core project.

Ready to take the course?

[Start learning :material-arrow-right:](00_orientation.md){ .md-button .md-button--primary }
