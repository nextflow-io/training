---
title: Building with Seqera AI
hide:
  - toc
---

# Building with Seqera AI

CoScientist is Seqera's AI assistant for bioinformatics, integrated directly into the Seqera Platform.
This side quest drives it hands-on against the [`nextflow-io/rnaseq-nf`](https://github.com/nextflow-io/rnaseq-nf) demo pipeline.
It takes you from a first conversation through development, debugging, opening a pull request, and writing a test.

## Audience & prerequisites

This side quest is designed for learners who have completed the [Hello Nextflow](../../hello_nextflow/index.md) course and are comfortable working on the command line with basic Nextflow concepts.

**Prerequisites**

- Completed [Hello Nextflow](../../hello_nextflow/index.md) or equivalent.
- Comfortable with the CLI and familiar with basic Nextflow concepts.

**Provided sandbox:** a training Platform workspace with CoScientist enabled, a working compute environment, and a GitHub account to fork into.

## Learning objectives

By the end of this side quest, you will be able to:

- Converse with CoScientist in the web chat UI to explore and plan pipeline work.
- Understand how MCP connects CoScientist to the Seqera Platform and GitHub.
- Register a pipeline on the Launchpad and inspect compute and data assets.
- Fork and run `rnaseq-nf` from the Launchpad.
- Debug a failed run using CoScientist.
- Open a pull request from a CoScientist conversation.
- Write an nf-test and run it to green from the CLI.

## The example pipeline

[`nextflow-io/rnaseq-nf`](https://github.com/nextflow-io/rnaseq-nf) is a minimal RNA-seq pipeline that runs INDEX and FASTQC in parallel, followed by QUANT (Salmon) and MULTIQC.
It ships test data in `data/ggal/` and has no nf-test suite yet. That makes it a good target for adding one.

**Time estimate:** ~2 hours

[Get started :material-arrow-right:](00_orientation.md){ .md-button .md-button--primary }
