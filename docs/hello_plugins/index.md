---
title: Hello Plugins
hide:
  - toc
---

# Hello Plugins

Nextflow's plugin system allows you to extend the language with custom functions, monitoring hooks, execution backends, and more.
Plugins enable the community to add features to Nextflow without modifying its core, making them ideal for sharing reusable functionality across pipelines.

During this training, you will learn how to use existing plugins and optionally create your own.

## Audience & prerequisites

Using existing plugins (Part 1) is straightforward and valuable for all Nextflow users.

Parts 2-6 cover building your own plugins, which many Nextflow developers may never need to do.
However, creating a plugin can be a good way to share reusable functions across pipelines or within a team.
These lessons involve Groovy code, build tools, and concepts that may be new if your background is in scripting or data analysis rather than software engineering.
No prior Java or Groovy experience is required; the lessons explain these concepts as they come up.

If the development sections feel challenging, focus on Part 1 and return to the rest when you have a specific plugin idea in mind.

**Prerequisites**

- A GitHub account OR a local installation as described [here](../envsetup/02_local).
- Completed the [Hello Nextflow](../hello_nextflow/index.md) course or equivalent.
- For plugin development (Parts 2-6): Java 21+ installed.

## Learning objectives

By the end of this training, you will be able to:

**Using plugins (Part 1):**

- Install and configure existing plugins in your workflows
- Import and use plugin functions

**Developing plugins (Parts 2-6):**

- Create a new plugin project using Nextflow's built-in project generator
- Implement custom functions callable from workflows
- Build, test, and install your plugin locally
- Monitor workflow events (e.g., task completion, pipeline start/end) for custom logging or notifications
- Add configuration options to make plugins customizable
- Distribute your plugin

## Lesson plan

#### Part 1: Plugin basics

Use an existing plugin from the registry in a Nextflow workflow.

#### Part 2: Create a plugin project

Generate a new plugin project and examine its structure.

#### Part 3: Custom functions

Implement custom functions, build your plugin, and run it in a workflow.

#### Part 4: Testing

Write and run unit tests using the Spock framework.

#### Part 5: Workflow monitoring

Respond to events like task completion to build a task counter.

#### Part 6: Configuration & Distribution

Read settings from `nextflow.config` to make your plugin customizable, then learn how to share it.

Ready to take the course?

[Start learning :material-arrow-right:](00_orientation.md){ .md-button .md-button--primary }
