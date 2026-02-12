---
title: Nextflow for {DOMAIN}
hide:
  - toc
page_type: index_page
index_type: course
additional_information:
  technical_requirements: true
  learning_objectives:
    - Write a linear workflow to apply {METHOD} to a single sample
    - Handle accessory files such as {ACCESSORY_FILES} appropriately
    - Leverage Nextflow's dataflow paradigm to parallelize per-sample processing
    - Implement multi-sample aggregation using relevant channel operators
  audience_prerequisites:
    - "**Audience:** This course is designed for researchers in {DOMAIN} and related fields who want to develop or customize data analysis pipelines."
    - "**Skills:** Some familiarity with the command line, basic scripting concepts, and common {DOMAIN} file formats is assumed."
    - "**Prerequisites:** Foundational Nextflow concepts and tooling covered in [Hello Nextflow](../../hello_nextflow/)."
---

# Nextflow for {DOMAIN}

**A hands-on course applying Nextflow to a real-world {DOMAIN} use case: {METHOD_SHORT_DESCRIPTION}.**

This course builds on the [Hello Nextflow](../../hello_nextflow/) beginner training and demonstrates how to use Nextflow in the specific context of the {DOMAIN} domain.
You will implement a {METHOD} pipeline with [{TOOL_A}]({TOOL_A_URL}) and [{TOOL_B}]({TOOL_B_URL}).

<!-- additional_information -->

## Course overview

This course is hands-on, with goal-oriented exercises structured to introduce information gradually.

You will start by running the analysis tools manually in the terminal to understand the methodology, then progressively build up a Nextflow pipeline that automates and scales the analysis.

### Lesson plan

We've broken this down into three parts that each focus on specific aspects of applying Nextflow to a {DOMAIN} use case.

| Course chapter                                            | Summary                                                                                           | Estimated duration |
| --------------------------------------------------------- | ------------------------------------------------------------------------------------------------- | ------------------ |
| [Part 1: Method overview](./01_method.md)                 | Understanding the {METHOD} methodology and running the tools manually                             | 30 mins            |
| [Part 2: Single-sample processing](./02_single_sample.md) | Building a pipeline that {PART2_SUMMARY}, then scaling to multiple samples                        | 60 mins            |
| [Part 3: Multi-sample aggregation](./03_multi_sample.md)  | Adding multi-sample {AGGREGATION_SUMMARY} using channel operators to aggregate per-sample outputs | 45 mins            |

By the end of this course, you will be able to apply foundational Nextflow concepts and tooling to a typical {DOMAIN} use case.

Ready to take the course?

[Get started :material-arrow-right:](00_orientation.md){ .md-button .md-button--primary }

<!-- Clearfix for float -->
<div style="content: ''; clear: both; display: table;"></div>
