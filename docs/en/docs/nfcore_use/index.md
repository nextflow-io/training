---
title: Use nf-core
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
    - Find and evaluate nf-core pipelines for a specific scientific question
    - Prepare and troubleshoot input samplesheets for your own data
    - Tune resource allocations based on execution reports
  audience_prerequisites:
    - "**Audience:** This course is designed for learners who already know how to run local Nextflow pipelines and are new to nf-core, and want to run existing community pipelines."
    - "**Skills:** Some familiarity with the command line, basic scripting concepts and common file formats is assumed."
    - "**Courses:** Must have completed [Nextflow Run](../nextflow_run/index.md) or otherwise be comfortable running a local pipeline with `nextflow run`."
    - "**Domain:** The exercises use bioinformatics pipelines, but no prior scientific domain knowledge is required."
---

# Use nf-core

**Use nf-core is a hands-on introduction to finding, running, and configuring nf-core community pipelines.**

Working through practical examples and guided exercises, you will learn to find and retrieve nf-core pipelines, run them using their built-in test profiles, and customize their execution through parameters and configuration files.

You will take away the skills and confidence to start running nf-core pipelines for your own analyses.

<!-- additional_information -->

## Course overview

This course is hands-on, with goal-oriented exercises structured to introduce information gradually.

You will start with `nf-core/demo`, a minimal pipeline maintained by the nf-core project for training purposes, then apply what you've learned to `nf-core/rnaseq`, a widely-used production pipeline for bulk RNA sequencing analysis.

This course focuses on running pipelines.
If you're looking for an intro to developing nf-core-compatible pipelines, see [Build with nf-core](../nfcore_build/index.md).

### Lesson plan

| Course chapter                                                       | Summary                                                                                      | Estimated duration |
| -------------------------------------------------------------------- | -------------------------------------------------------------------------------------------- | ------------------ |
| [Part 1: Run a demo pipeline](./01_run_demo.md)                      | Find and retrieve an nf-core pipeline and run it using its test profile                      | 20 mins            |
| [Part 2: Configure pipeline execution](./02_configure_execution.md)  | Set parameters, understand validation, and customize resource allocation and tool arguments  | 20 mins            |
| [Part 3: Run a production pipeline](./03_run_production_pipeline.md) | Pin and run nf-core/rnaseq, fit it to a small machine, and cancel a run                      | 25 mins            |
| [Part 4: Choose a pipeline](./04_choose_pipeline.md)                 | Find candidate pipelines for a scientific question and evaluate whether they fit             | 40 mins            |
| [Part 5: Prepare your inputs](./05_inputs.md)                        | Write a samplesheet for your own data or for public data, and troubleshoot validation errors | 45 mins            |
| [Part 6: Scale up and tune compute resources](./06_resources.md)     | Run on a bigger machine, then use execution reports to right-size resource allocations       | 40 mins            |

By the end of this course, you will be able to take advantage of the wealth of community pipelines offered by the nf-core project.

Ready to take the course?

[Start learning :material-arrow-right:](00_orientation.md){ .md-button .md-button--primary }
