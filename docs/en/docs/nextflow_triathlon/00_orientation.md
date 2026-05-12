# Getting started

## Start a training environment

To use the pre-built environment we provide on GitHub Codespaces, click the "Open in GitHub Codespaces" button below. For other options, see [Environment options](../envsetup/index.md).

We recommend opening the training environment in a new browser tab or window (use right-click, ctrl-click or cmd-click depending on your equipment) so that you can read on while the environment loads.
You will need to keep these instructions open in parallel to work through the course.

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

### Environment basics

This training environment contains all the software, code and data necessary to work through the training course, so you don't need to install anything yourself.

The codespace is set up with a VSCode interface, which includes a filesystem explorer, a code editor and a terminal shell.
All instructions given during the course (e.g. 'open the file', 'edit the code' or 'run this command') refer to those three parts of the VSCode interface unless otherwise specified.

If you are working through this course by yourself, please acquaint yourself with the [environment basics](../envsetup/01_setup.md) for further details.

### Machine size

This course demonstrates parallel execution of multiple pipeline tasks.
To see parallelism in action, you need a machine with multiple cores.

!!! tip "Use an 8-core machine"

    The default Codespaces machine (2 cores) will serialize task execution and prevent you from observing parallel runs.
    Before starting, upgrade your codespace to an **8-core** machine:

    1. At [github.com/codespaces](https://github.com/codespaces), find your codespace and click the **...** menu.
    2. Select **Change machine type**.
    3. Choose **8-core** and click **Update codespace**.
    4. Restart the codespace for the change to take effect.

### Version requirements

This course requires Nextflow 25.10.2 or later and nf-core tools 3.5.2 or later.

- **Part 1** runs with the **v2 syntax parser enabled** (the default in 25.10+)
- **Part 2** requires the **v1 syntax parser** — the lesson will tell you when to switch

If you are using a local or custom environment, please make sure you are using the correct settings as documented [here](../info/nxf_versions.md).

## Course structure

This course covers three disciplines in sequence.
Part 3 uses the Seqera Platform web interface.
Parts 1 and 2 are run from the command line, each from its own subdirectory:

- **Part 1** uses `nextflow-triathlon/basics/`
- **Parts 2 and 3** use `nextflow-triathlon/nf-core/`

Each part begins with a tip reminding you which directory to use.

## Readiness checklist

Think you're ready to dive in?

- [ ] I understand the goal of this course and its prerequisites
- [ ] My environment is up and running

If you can check all the boxes, you're good to go.

**To continue to [Part 1: Run Nextflow](./01_run_basics.md), click on the arrow in the bottom right corner of this page.**
