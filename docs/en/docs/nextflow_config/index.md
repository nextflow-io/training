---
title: Nextflow Dot Config
hide:
  - toc
page_type: index_page
index_type: course
additional_information:
  technical_requirements: true
  learning_objectives:
    - Manage workflow input parameters with run-specific configuration and parameter files
    - Control where and how workflow outputs get published
    - Switch software packaging technology between Docker and Conda
    - Select an execution platform and control compute resource allocations
    - Define and combine profiles to switch between preset configurations
  audience_prerequisites:
    - "**Audience:** This course is designed for learners who already know how to launch local Nextflow pipelines and want to configure execution in more depth."
    - "**Skills:** Some familiarity with the command line is assumed."
    - "**Courses:** Must have completed [Nextflow Run](../nextflow_run/index.md) or otherwise be comfortable running a local pipeline with `nextflow run`."
---

# Nextflow Dot Config

**Nextflow Dot Config is a hands-on introduction to configuring Nextflow pipeline execution.**

Working through goal-oriented exercises, you will learn how to manage inputs and outputs, adapt a pipeline to different compute environments, and bundle configuration into switchable profiles.

You will take away the skills and confidence to configure Nextflow pipeline execution like a pro.

<!-- additional_information -->

## Course overview

This course is hands-on, and builds on the skills covered in [Nextflow Run](../nextflow_run/index.md).

You will take the same multi-step pipeline from that course and progressively extend its configuration: first managing inputs and outputs, then adapting execution to different compute environments, and finally bundling everything into profiles you can switch between at runtime.

### Lesson plan

| Course chapter                                                                     | Summary                                                                                              | Estimated duration |
| ---------------------------------------------------------------------------------- | ---------------------------------------------------------------------------------------------------- | ------------------ |
| [Part 1: Manage inputs and outputs](./01_parameters_and_outputs.md)                | Manage input parameters with configuration and parameter files, and control output publishing        | 20 mins            |
| [Part 2: Adapt to your compute environment](./02_packaging_execution_resources.md) | Switch software packaging technology, select an execution platform, and control resource allocations | 25 mins            |
| [Part 3: Use profiles to switch configurations](./03_profiles.md)                  | Define and combine profiles, and inspect the fully resolved configuration                            | 15 mins            |

By the end of this course, you will be comfortable configuring Nextflow pipelines for a range of compute environments, and switching between them with minimal hassle.

Ready to take the course?

[Start learning :material-arrow-right:](00_orientation.md){ .md-button .md-button--primary }
