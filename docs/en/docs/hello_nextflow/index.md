---
title: Hello Nextflow
hide:
    - toc
page_type: index_page
index_type: course
additional_information:
    technical_requirements: true
    learning_objectives:
        - Write processes and workflows, and connect them into a multi-step pipeline
        - Use channel factories and operators to feed inputs into a workflow and process them in parallel
        - Publish workflow outputs and organize them in the results directory
        - Add command-line parameters, with their default and test values in the configuration
        - Organize pipeline code into reusable modules
        - Run processes in software containers, with a Conda package as an alternative, for reproducibility and portability
        - Use task directories and error messages to troubleshoot basic issues
    audience_prerequisites:
        - "**Audience:** This course is designed for learners who are new to developing Nextflow pipelines and want to write their own."
        - "**Skills:** Some familiarity with the command line, basic scripting concepts and common file formats is assumed."
        - "**Courses:** Must have completed [Nextflow Run](../nextflow_run/index.md), which this course builds on."
        - "**Domain:** The exercises are all domain-agnostic, so no prior scientific knowledge is required."
    videos_playlist: https://youtube.com/playlist?list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&si=eF7cLR62goy-lc6n
---

# Hello Nextflow

**Hello Nextflow is a hands-on introduction to building reproducible and scalable data analysis workflows.**

Working through practical examples and guided exercises, you will learn the fundamentals of developing pipelines with Nextflow, including how to define processes, connect them into pipelines, manage files and software dependencies, and parallelize execution effortlessly.

You will take away the skills and confidence to start developing and running your own workflows with Nextflow.

<!-- additional_information -->

## Course overview

This course is designed to be hands-on, with goal-oriented exercises structured to introduce information gradually.

In [Nextflow Run](../nextflow_run/index.md), you ran the finished Hello pipeline: it reads names from a CSV file, greets each one, translates the greetings into another language, collects them into a single file, and turns the result into ASCII art.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello_pipeline_complete.svg"
</figure>

In this course, you write that pipeline yourself.
You start from a script that does a single thing, writing 'Hello World!' to a file, and add one piece at a time.

### Lesson plan

Each part adds a piece of the pipeline, and focuses on the Nextflow concepts needed to build it.

| Course chapter | What you add | Estimated duration |
| --- | --- | --- |
| [Part 1: Hello World](./01_hello_world.md) | The `sayHello` process, greeting a name provided on the command line: processes, workflows, publishing outputs and parameters | 30 mins |
| [Part 2: Hello Channels](./02_hello_channels.md) | Reading the names from a CSV file: channels, operators and parallel execution | 45 mins |
| [Part 3: Hello Workflow](./03_hello_workflow.md) | The `translateGreeting` and `collectGreetings` steps: connecting processes with channels, writing a process in Python, and handling multiple inputs and outputs | 60 mins |
| [Part 4: Hello Modules](./04_hello_modules.md) | Modules for each process, and parameter values in the configuration: organizing code for reuse and maintenance | 30 mins |
| [Part 5: Hello Containers](./05_hello_containers.md) | The `cowpy` step: running a process in a container | 60 mins |

After Part 5, the optional [Challenges](./challenges.md) apply the Hello pipeline to Nextflow features that the course does not cover.

By the end of this course, you will be well-prepared for tackling the next steps in your journey to develop reproducible workflows for your scientific computing needs.

Ready to take the course?

[Get started :material-arrow-right:](00_orientation.md){ .md-button .md-button--primary }

<!-- Clearfix for float -->
<div style="content: ''; clear: both; display: table;"></div>
