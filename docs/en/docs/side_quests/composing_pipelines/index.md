---
title: Composing Pipelines
hide:
  - toc
---

# Composing Pipelines

Every Nextflow pipeline is built from smaller pieces.
Modules hold the processes that do the work, workflows group processes into reusable stages, and a pipeline wires those stages together behind its own parameters and outputs.
Composition happens at each of those levels, and at one more: a whole pipeline can be included and called by another pipeline.

This course follows one small project up that ladder, from modules to workflows, from workflows to a pipeline, and from a pipeline to a pipeline of pipelines.
You start with a greetings pipeline that has two stages: one creates timestamped greetings, and the other transforms them into uppercase and reversed text.
In Parts 1 and 2, you compose its modules into named workflows and give those workflows typed interfaces.
In Parts 3 to 6, the whole greetings pipeline becomes a building block: you make it includable, then connect it to a separately developed report pipeline that summarizes its outputs, in the same way that a differential abundance analysis builds on the outputs of an RNA-seq pipeline.

That last step has traditionally needed glue outside Nextflow: a script that runs one pipeline, waits for it, and builds the next pipeline's inputs from the first one's output paths.
By the end of this course, you'll replace that glue with a single pipeline that includes both and connects them through typed interfaces, while each of them still runs on its own.

## Audience & prerequisites

Parts 1 and 2 are about composition within a single pipeline and are relevant to anyone writing Nextflow pipelines with more than a handful of processes.
Parts 3-6 are about composition between pipelines and are aimed at developers who build on the outputs of other pipelines, or whose pipelines others build on.

**Prerequisites**

- A GitHub account OR a local installation as described [here](../../envsetup/02_local.md).
- Completed the [Hello Nextflow](../../hello_nextflow/index.md) course or equivalent.
- Completed the [Types and Records](../types_and_records/index.md) side quest before starting Part 2.

**Working directory:** `side-quests/composing_pipelines`

!!! note "Requires Nextflow 26.10"

    This course requires Nextflow 26.10 or later.
    Pipeline composition, used in Parts 3-6, is an experimental feature.

#### Open the training codespace

If you haven't yet done so, make sure to open the training environment as described in the [Environment Setup](../../envsetup/index.md).

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

## Learning objectives

By the end of this training, you will be able to:

**Composing within a pipeline (Parts 1-2):**

- Break a pipeline into named workflows with `take:` and `emit:` interfaces
- Compose named workflows from an entry workflow and pass data between them
- Give workflow interfaces static types, with records at the boundary between workflows
- Check the wiring between workflows with `nextflow lint`

**Composing between pipelines (Parts 3-6):**

- Give a pipeline a typed `params {}` block and an `output {}` block that together form its interface
- Include a whole pipeline and call it like a workflow
- Feed one pipeline's outputs into another pipeline in a single run
- Pass parameters and configuration to included pipelines, and avoid the common pitfalls
- Write pipelines that work both standalone and when included by another pipeline

## Lesson plan

#### Part 1: Workflows of Workflows

Turn the greeting and transform stages into named workflows with clear interfaces, and compose them into one pipeline.

#### Part 2: Typed Interfaces

Write the contract between the two workflows down with static types and records, so that `nextflow lint` reports a broken connection at the call.

#### Part 3: Pipelines as Workflows

Give the greetings pipeline a typed `params {}` block and a complete `output {}` block, run it on its own, then include it in a new pipeline and call it like a workflow.

#### Part 4: Chaining Pipelines

Connect the greetings pipeline to a downstream report pipeline, first the traditional way with two runs and a hand-written samplesheet, then as one composed run with a single DAG and one `-resume`.

#### Part 5: Parameters and Configuration

Pass parameters and configuration to included pipelines from the command line, from config and from the calling code, and control which configuration applies across pipeline boundaries.

#### Part 6: Making Pipelines Composable

Track down a bug that only appears when the report pipeline is included by another pipeline, fix it, and finish with a checklist for writing pipelines that work both standalone and composed.

Ready to take the course?

[Start learning :material-arrow-right:](01_workflows_of_workflows.md){ .md-button .md-button--primary }
