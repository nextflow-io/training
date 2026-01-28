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

1. **All Nextflow users** (Parts 1-2): Learn to discover, install, and use existing plugins to extend your pipelines.
2. **Developers** (Parts 3-7): Learn to create your own plugins with custom functions, observers, and configuration.

**Prerequisites**

- A GitHub account OR a local installation as described [here](../envsetup/02_local).
- Completed the [Hello Nextflow](../hello_nextflow/index.md) course or equivalent.
- For plugin development (Parts 3-7): Java 21+ installed and basic familiarity with object-oriented programming.

## Learning objectives

By the end of this training, you will be able to:

**Using plugins (Parts 1-2):**

- Understand what plugins are and how they extend Nextflow
- Install and configure existing plugins in your workflows
- Import and use plugin functions

**Developing plugins (Parts 3-7):**

- Create a new plugin project using Nextflow's scaffolding command
- Implement custom functions callable from workflows
- Create trace observers to hook into workflow lifecycle events
- Add configuration options to make plugins customizable
- Build, test, and distribute your plugin

## Detailed lesson plan

This training course teaches you both how to use plugins and how to create them.

#### Part 1: Plugin basics

You'll start by understanding **what plugins are and how they work**.
You'll learn about the different types of extensions plugins can provide (functions, observers, executors, filesystems) and how Nextflow's plugin system is built on PF4J.
Then you'll discover and use existing plugins from the Nextflow Plugin Registry.

#### Part 2: Create a plugin project

Next, you'll **scaffold a new plugin project** using Nextflow's built-in command.
You'll examine the generated project structure and understand how the different components (Plugin, Extension, Factory, Observer) work together.

#### Part 3: Custom functions

You'll **implement custom functions** that can be called from Nextflow workflows.
Using the `@Function` annotation, you'll create functions that transform data and make them available via the familiar `include` syntax.

#### Part 4: Build and test

You'll learn the **plugin development cycle**: building with Gradle, writing unit tests with Spock, and installing locally for testing.
You'll see how to iterate quickly and debug common issues.

#### Part 5: Trace observers

You'll explore **trace observers**, which let your plugin hook into workflow lifecycle events.
You'll build an observer that responds to task completion and workflow events, enabling features like custom logging, metrics collection, or notifications.

#### Part 6: Configuration

You'll make your plugin **configurable** by reading settings from `nextflow.config`.
Users will be able to customize plugin behavior without modifying code.

#### Part 7: Distribution

Finally, you'll learn **how to share your plugin** with others.
You can publish to the public Nextflow Plugin Registry or set up internal distribution for proprietary plugins.

**By the end of the course, you'll have created a fully functional Nextflow plugin with custom functions, lifecycle observers, and configuration support.**

Ready to take the course?

[Start learning :material-arrow-right:](00_orientation.md){ .md-button .md-button--primary }
