---
title: Nextflow Run
hide:
  - toc
page_type: index_page
index_type: course
additional_information:
  technical_requirements: true
  learning_objectives:
    - Launch and manage execution of Nextflow workflows
    - Find and interpret outputs (results) and log files
    - Troubleshoot basic issues
    - Recognize core Nextflow components in a simple multi-step workflow
    - Configure pipeline execution to run on common computing platforms including HPC and cloud
    - Summarize best practices for reproducibility, portability and code re-use that make pipelines FAIR, including code modularity and software containers
  audience_prerequisites:
    - "**Audience:** This course is designed for learners who are completely new to Nextflow and want to execute existing pipelines."
    - "**Skills:** Some familiarity with the command line, basic scripting concepts and common file formats is assumed."
    - "**Domain:** The exercises are all domain-agnostic, so no prior scientific knowledge is required."
---

# Nextflow Run

**Nextflow Run is a hands-on introduction to running reproducible and scalable data analysis workflows.**

Working through practical examples and guided exercises, you will learn the fundamentals of using Nextflow, including how to execute pipelines, manage files and software dependencies, parallelize execution effortlessly, and run workflows across different computing environments. You will take away the skills and confidence to start running workflows with Nextflow.

<!-- additional_information -->

## Course overview

The entire course is designed to be hands-on, with goal-oriented exercises structured to introduce information gradually.

You will execute several versions of a simple Nextflow pipeline, starting with one consisting of a single step and progressively building up to a multi-step version that takes some text inputs, runs a few transformation steps, and outputs a single text file containing an ASCII picture of a character saying the transformed text.

### Lesson plan

We've broken this down into three parts that will each focus on specific aspects of running and managing pipelines written in Nextflow.

| Course chapter                                 | Summary                                                                                          | Estimated duration |
| ---------------------------------------------- | ------------------------------------------------------------------------------------------------ | ------------------ |
| [Part 1: Run basic operations](./01_basics.md) | Launching and managing execution of a simple workflow                                            | 30 mins            |
| [Part 2: Run real pipelines](./02_pipeline.md) | Processing complex inputs, running multi-step workflows and parallelizing execution effortlessly | 60 mins            |
| [Part 3: Run configuration](./03_config.md)    | Customizing pipeline behavior and optimizing usage in different computational environments       | 60 mins            |

By the end of this course, you will be well-prepared for tackling the next steps in your journey to run reproducible workflows for your scientific computing needs.

Ready to take the course?

[Start learning :material-arrow-right:](00_orientation.md){ .md-button .md-button--primary }
