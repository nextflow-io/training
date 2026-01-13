---
title: Hello nf-core
hide:
  - toc
---

# Hello nf-core

**nf-core** is a community effort to develop and maintain a curated set of scientific pipelines built using Nextflow, as well as relevant tooling and guidelines that promote open development, testing, and peer review.

![nf-core logo](./img/nf-core-logo.png)

The pipelines developed by the nf-core community are designed to be modular, scalable, and portable, allowing researchers to easily adapt and execute them using their own data and compute resources.
The best practices guidelines enforced by the project further ensure that the pipelines are robust, well-documented, and validated against real-world datasets.
This helps to increase the reliability and reproducibility of scientific analyses and ultimately enables researchers to accelerate their scientific discoveries.

During this training, you will be introduced to nf-core in a series of hands-on exercises as described further below.

**Additional information:** You can learn more about the project's origins and governance at [https://nf-co.re/about](https://nf-co.re/about).

**Reference publication:** nf-core is published in Nature Biotechnology: [Nat Biotechnol 38, 276–278 (2020). Nature Biotechnology](https://www.nature.com/articles/s41587-020-0439-x).
There is another recent publication in Genome Biology [Genome Biol 26, 228 (2025). Genome Biology](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-025-03673-9).

## Audience & prerequisites

This training is intended for learners who have at least basic Nextflow skills and wish to level up to using nf-core resources and best practices in their work.

**Prerequisites**

- A GitHub account OR a local installation as described [here](../envsetup/02_local).
- Experience with command line and basic scripting.
- Completed the [Hello Nextflow](../hello_nextflow/index.md) course or equivalent.

## Learning objectives

You will learn to use and develop nf-core compatible modules and pipelines, and to utilize nf-core tooling effectively.

By the end of this training, you will be able to:

- Find and run nf-core pipelines
- Describe the code structure and project organization of nf-core pipelines
- Create a basic nf-core compatible pipeline from a template
- Convert basic Nextflow modules to nf-core compatible modules
- Add nf-core modules to an nf-core compatible pipeline
- Validate inputs and parameters using nf-core tooling

## Detailed lesson plan

This training course aims to teach you the core concepts for running nf-core-style pipelines.
We won't cover everything there is to know about nf-core pipelines, because nf-core encompasses many features and conventions developed by the community over years.
Instead, we will focus on the essential concepts that will help you get started and understand how nf-core works.

#### Part 1: Run a demo pipeline

First, you'll **run an existing nf-core pipeline** and examine its code structure to get a sense of what makes these pipelines different from basic Nextflow workflows.
The elaborate directory structure, configuration system, and standardized conventions might seem like a lot at first, but the benefits will become clear as you learn to decode and utilize these resources effectively.

#### Part 2: Rewrite Hello for nf-core

Next, you'll **adapt an existing workflow to the nf-core template scaffold**, starting from the simple workflow produced in the [Hello Nextflow](../hello_nextflow/index.md) course.
Many pipeline development efforts start from existing code, so learning how to restructure an existing workflow to leverage nf-core's nested workflow system is a practical skill you're likely to use repeatedly in your work.

#### Part 3: Use an nf-core module

Then you'll discover one of nf-core's biggest advantages: the **community modules library**.
Instead of writing every process from scratch, you'll learn to integrate pre-built, tested modules that wrap common bioinformatics tools.
This approach saves time and ensures consistency across pipelines.

#### Part 4: Make an nf-core module

Of course, the modules library doesn't have everything, so you'll also learn to **create your own nf-core-style module**.
You'll learn to work with the specific structure, naming conventions, and metadata requirements that make modules shareable and maintainable by the community.

#### Part 5: Add input validation

Finally, you'll implement **input validation** for both command-line parameters and input data files using nf-schema.
This catches errors before pipelines start to run, providing fast feedback and clear error messages. This type of upfront validation makes pipelines more robust and easier to use.

**By the end of the course, you'll have transformed a basic Nextflow workflow into an nf-core-style pipeline with standardized structure, reusable components, and robust validation.**

Ready to take the course?

[Start learning :material-arrow-right:](00_orientation.md){ .md-button .md-button--primary }
