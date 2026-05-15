---
title: Troubleshooting Workflows
hide:
  - toc
---

# Troubleshooting Workflows

Debugging is a critical skill that can save you hours of frustration and help you become a more effective Nextflow developer.
Throughout your career, especially when you're starting out, you'll encounter bugs while building and maintaining your workflows.

This mini-course covers two complementary skills.
The first is recognising the shape of common errors, so that the error message in your terminal becomes a signpost rather than a wall.
The second is the toolkit of techniques you reach for when an error isn't immediately obvious from its message.

## Audience & prerequisites

These lessons are aimed at Nextflow users who have completed the basics and want to debug their own workflows confidently.

**Prerequisites**

- Completed the [Hello Nextflow](../../hello_nextflow/index.md) tutorial or equivalent.
- Comfortable with basic Nextflow concepts (processes, channels, operators).
- Docker installed (the examples use containerised processes).

**Optional:** We recommend completing the [IDE Features for Nextflow Development](../dev_environment/index.md) side quest first.
That covers comprehensive coverage of IDE features that support debugging (syntax highlighting, error detection, etc.), which we use heavily here.

**Working directory:** `side-quests/debugging`

## Learning objectives

By the end of this mini-course, you will be able to:

**Recognising common errors (Part 1):**

- Read Nextflow error messages and locate the relevant code
- Recognise syntax errors, channel structure errors, and process errors
- Apply targeted fixes for each error category

**The Nextflow debugging toolkit (Part 2):**

- Inspect a process work directory to find out what actually ran
- Validate workflow logic before execution with `-preview`
- Stream process output in real time with `debug true`
- Iterate on workflow logic without running real commands using `-stub-run`
- Diagnose cache invalidation problems with `-dump-hashes`
- Apply a systematic four-phase debugging method

## Lesson plan

#### Part 1: Common errors and how to fix them

A catalogue of the most common errors you'll meet, with an example for each.
Read this one through, or use it as a reference when a specific error message lands in your terminal.

#### Part 2: The Nextflow debugging toolkit

Every tool is applied to the same small pipeline, so you can see how each tool fits into a real development workflow.
The lesson culminates with a `-dump-hashes` walkthrough of three cache-invalidation experiments, and a practical debugging exercise on an unfamiliar pipeline.

Ready to start?

[Start with Part 1 :material-arrow-right:](01_common_errors.md){ .md-button .md-button--primary }
