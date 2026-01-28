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
    - Run pipelines from remote repositories like GitHub
    - Find and interpret outputs (results) and log files
    - Troubleshoot common issues when tasks fail
    - Recognize core Nextflow components in a simple multi-step workflow
    - Configure pipeline execution to run on common computing platforms including HPC and cloud
    - Summarize best practices for reproducibility, portability and code re-use that make pipelines FAIR, including code modularity and software containers
  audience_prerequisites:
    - "**Audience:** This course is designed for learners who want to **run existing pipelines** created by others (such as nf-core pipelines). You will learn to recognize code patterns but won't write pipelines from scratch. If you want to **develop your own pipelines**, take [Hello Nextflow](../hello_nextflow/index.md) instead."
    - "**Skills:** Some familiarity with the command line, basic scripting concepts and common file formats is assumed."
    - "**Domain:** The exercises are all domain-agnostic, so no prior scientific knowledge is required."
---

# Nextflow Run

**Nextflow Run is a hands-on introduction to running reproducible and scalable data analysis workflows.**

You will learn how to execute existing pipelines, interpret their outputs, configure them for different computing environments, and troubleshoot common issues.
By the end, you'll have the skills and confidence to run Nextflow pipelines created by others, such as those from [nf-core](https://nf-co.re).

<!-- additional_information -->

## Course overview

### Which course should I take?

We offer two introductory Nextflow courses with different goals:

| Course                                           | Goal                       | You will...                                                                      |
| ------------------------------------------------ | -------------------------- | -------------------------------------------------------------------------------- |
| **Nextflow Run** (this course)                   | Run existing pipelines     | Learn to execute, configure, and troubleshoot pipelines created by others        |
| **[Hello Nextflow](../hello_nextflow/index.md)** | Develop your own pipelines | Build pipelines from scratch, learning channels, operators, and modules in depth |

Both courses use the same training environment and cover some overlapping ground (running pipelines, configuration basics), so you don't need to take both unless you want comprehensive coverage.

### What you'll do

This course is hands-on, with goal-oriented exercises structured to introduce information gradually.
It's named after the core Nextflow command:

```bash
nextflow run <pipeline>
```

You will execute several versions of a Nextflow pipeline that processes text inputs.
You'll start with a simple version that consists of a single step, and eventually progress to a multi-step version that takes a CSV file of tabular text inputs, runs a few transformation steps, and outputs a single text file containing an ASCII picture of a character saying the transformed text.

### Lesson plan

We've broken this down into four parts that will each focus on specific aspects of running and managing pipelines written in Nextflow.

| Course chapter                                     | Summary                                                                                                            | Estimated duration |
| -------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------ | ------------------ |
| [Part 1: Run basic operations](./01_basics.md)     | Launching and managing execution of a simple workflow, running remote pipelines                                    | 45 mins            |
| [Part 2: Run real pipelines](./02_pipeline.md)     | Processing complex inputs, running multi-step workflows, using containers and parallelizing execution effortlessly | 60 mins            |
| [Part 3: Run configuration](./03_config.md)        | Customizing pipeline behavior and optimizing usage in different computational environments                         | 60 mins            |
| [Part 4: Troubleshooting](./04_troubleshooting.md) | Diagnosing task failures, common error messages, debugging techniques                                              | 30 mins            |

By the end of this course, you will be well-prepared for tackling the next steps in your journey to run reproducible workflows for your scientific computing needs.

Ready to take the course?

[Start learning :material-arrow-right:](00_orientation.md){ .md-button .md-button--primary }
