---
title: Hello Plugins
hide:
  - toc
---

# Hello Plugins

Nextflow's plugin system allows you to extend the language with custom functions, executors, trace observers, and more.
Plugins enable the community to add features to Nextflow without modifying its core, making them ideal for sharing reusable functionality across pipelines.

During this training, you will learn how to use existing plugins and optionally create your own.

## Audience & prerequisites

This training is designed for two audiences:

1. **All Nextflow users** (Part 1): Discover, install, and use existing plugins to extend your pipelines.
2. **Developers** (Parts 2-6): Create your own plugins with custom functions, observers, and configuration.

**Prerequisites**

- A GitHub account OR a local installation as described [here](../envsetup/02_local).
- Completed the [Hello Nextflow](../hello_nextflow/index.md) course or equivalent.
- For plugin development (Parts 2-6): Java 21+ installed and basic familiarity with object-oriented programming.

## Learning objectives

By the end of this training, you will be able to:

**Using plugins (Part 1):**

- Install and configure existing plugins in your workflows
- Import and use plugin functions

**Developing plugins (Parts 2-6):**

- Create a new plugin project using Nextflow's scaffolding command
- Implement custom functions callable from workflows
- Build, test, and install your plugin locally
- Create trace observers to hook into workflow lifecycle events
- Add configuration options to make plugins customizable
- Distribute your plugin

## Lesson plan

#### Part 1: Plugin basics

Use an existing plugin from the registry in a Nextflow workflow.

#### Part 2: Create a plugin project

Scaffold a new plugin project and examine its structure.

#### Part 3: Custom functions

Implement custom functions, build your plugin, and run it in a workflow.

#### Part 4: Testing

Write and run unit tests using the Spock framework.

#### Part 5: Trace observers

Hook into workflow lifecycle events to build a task counter.

#### Part 6: Configuration & Distribution

Read settings from `nextflow.config` to make your plugin customizable, then learn how to share it.

Ready to take the course?

[Start learning :material-arrow-right:](00_orientation.md){ .md-button .md-button--primary }
