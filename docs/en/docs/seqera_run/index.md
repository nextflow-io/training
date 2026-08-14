---
title: Scale with Seqera
hide:
  - toc
page_type: index_page
index_type: course
additional_information:
  technical_requirements: true
  learning_objectives:
    - Sign up for Seqera Platform and explore the Community Showcase
    - Add a pipeline to a workspace and launch it from the web interface
    - Authenticate and launch pipelines from the command line with the `tw` CLI
    - Register a GitHub-hosted pipeline and launch it both ways
  audience_prerequisites:
    - "**Audience:** This course is designed for learners who want to run Nextflow pipelines at scale using Seqera Platform."
    - "**Skills:** Familiarity with running nf-core pipelines from the command line is assumed."
    - "**Courses:** Must have completed [Use nf-core](../nfcore_run/index.md) or otherwise be comfortable running `nf-core/demo` and `nf-core/rnaseq`."
---

# Scale with Seqera

**Scale with Seqera is a hands-on introduction to launching and monitoring Nextflow pipelines with Seqera Platform.**

Working through practical examples, you will set up access to Seqera Platform, launch a production-scale pipeline from both the web interface and the command line, and add a new pipeline to your workspace.

You will take away the skills and confidence to run and monitor your own pipelines on Seqera Platform.

<!-- additional_information -->

## Course overview

This course is hands-on, and builds on the pipelines you already ran in [Use nf-core](../nfcore_run/index.md).

You will start by signing up for Seqera Platform and launching `nf-core/rnaseq`, a production-scale pipeline, from the web interface.
Then you'll switch to the `tw` command-line tool to do the same from a terminal, and finally register a new pipeline, `nf-core/demo`, and launch it both ways.

### Lesson plan

| Course chapter                                                             | Summary                                                                                      | Estimated duration |
| -------------------------------------------------------------------------- | -------------------------------------------------------------------------------------------- | ------------------ |
| [Part 1: Launch pipelines from the web interface](./01_run_with_seqera.md) | Set up Seqera Platform access and launch a production-scale pipeline from the web interface  | 20 mins            |
| [Part 2: Launch pipelines from the command line](./02_launch_from_cli.md)  | Authenticate the `tw` CLI, launch a saved pipeline, and register a new pipeline from the CLI | 25 mins            |

By the end of this course, you will be comfortable launching and monitoring Nextflow pipelines on Seqera Platform, whether you prefer working from the web interface or the command line.

Ready to take the course?

[Start learning :material-arrow-right:](00_orientation.md){ .md-button .md-button--primary }
