---
title: Nextflow Run
hide:
  - toc
page_type: index_page
index_type: course
additional_information:
  technical_requirements: true
  learning_objectives:
    - Retrieve, launch and manage execution of nf-core pipelines
    - Describe the code structure and project organization of nf-core pipelines
    - Create a basic nf-core compatible pipeline from a template
    - Upgrade a plain Nextflow workflow to fit nf-core standards
    - Add nf-core modules to an nf-core compatible pipeline
    - Contribute your own modules to nf-core
    - Validate inputs and parameters using nf-core tooling
  audience_prerequisites:
    - "**Audience:** This course is designed for learners who are already familiar with basic Nextflow and want to learn to use nf-core resources and best practices."
    - "**Skills:** Familiarity with the command line, basic scripting concepts and common file formats is assumed."
    - "**Courses:** Must have completed the [Hello Nextflow](../hello_nextflow/index.md) course or equivalent."
    - "**Domain:** The exercises are all domain-agnostic, so no prior scientific knowledge is required."
---

# Hello nf-core

**Hello nf-core is a hands-on introduction to using nf-core resources and best practices for developing Nextflow pipelines.**

Working through practical examples and guided exercises, you will learn to use and develop nf-core compatible modules and pipelines, and to utilize nf-core tooling effectively.

You will take away the skills and confidence to start developing pipelines according to nf-core best practices.

![nf-core logo](./img/nf-core-logo.png)

<!-- additional_information -->

## Course overview

This course is designed to be hands-on, with goal-oriented exercises structured to introduce information gradually.

You will be introduced to [**nf-core**](https://nf-co.re/), a community effort to develop and maintain a curated set of scientific pipelines built using Nextflow, as well as relevant tooling and guidelines that promote open development, testing, and peer review ([Nat Biotechnol 38, 276–278 (2020)](https://www.nature.com/articles/s41587-020-0439-x), [Genome Biol 26, 228 (2025)](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-025-03673-9)).

The pipelines developed by the nf-core community are designed to be modular, scalable, and portable, allowing researchers to easily adapt and execute them using their own data and compute resources.
The best practices guidelines enforced by the project further ensure that the pipelines are robust, well-documented, and validated against real-world datasets.
This helps to increase the reliability and reproducibility of scientific analyses and ultimately enables researchers to accelerate their scientific discoveries.

We won't cover everything there is to know about nf-core pipelines in this course, because nf-core encompasses many features and conventions developed by the community over years.
Instead, we will focus on the essential concepts that will help you get started and understand how nf-core works.

### Lesson plan

We've broken this down into five parts that will each focus on specific aspects of using nf-core resources.

| Course chapter                                             | Summary                                                                                                                                                            | Estimated duration |
| ---------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------ | ------------------ |
| [Part 1: Run a demo pipeline](./01_run_demo.md)            | Run an existing nf-core pipeline and examine its code structure to get a sense of what makes these pipelines different from basic Nextflow workflows               | 30 mins            |
| [Part 2: Rewrite Hello for nf-core](./02_rewrite_hello.md) | Adapt an existing workflow to the nf-core template scaffold, starting from the simple workflow produced in the [Hello Nextflow](../hello_nextflow/index.md) course | 60 mins            |
| [Part 3: Use an nf-core module](./03_use_module.md)        | Explore the community modules library and learn to integrate pre-built, tested modules that wrap common bioinformatics tools                                       | 30 mins            |
| [Part 4: Make an nf-core module](./04_make_module.md)      | Create your own nf-core-style module using the specific structure, naming conventions, and metadata requirements established by nf-core                            | 30 mins            |
| [Part 5: Add input validation](./05_input_validation.md)   | Implement input validation for both command-line parameters and input data files using nf-schema                                                                   | 30 mins            |

By the end of this course, you will be able to take advantage of the enormous wealth of resources offered by the nf-core project.

Ready to take the course?

[Start learning :material-arrow-right:](00_orientation.md){ .md-button .md-button--primary }
