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

### Version requirements

This course requires Nextflow 25.10.2 or later and nf-core tools 3.5.2 or later.

- **Part 1** runs with the **v2 syntax parser enabled** (the default in 25.10+)
- **Part 2** requires the **v1 syntax parser** — the lesson will tell you when to switch

If you are using a local or custom environment, please make sure you are using the correct settings as documented [here](../info/nxf_versions.md).

## Get ready to work

Once your codespace is running, there are two things to do before diving in: set your working directory, and take a look at the materials provided.

### Set the working directory

By default, the codespace opens at the root of all training courses.
For this course, change to the `triathlon/` directory:

```bash
cd /workspaces/training/triathlon
```

Then set VSCode to focus on this directory, so only the relevant files appear in the file explorer sidebar:

```bash
code .
```

!!! tip

    Each part of the course runs from its own subdirectory.
    The first step of each part will tell you which one to use.

### Explore the materials provided

You can explore the course materials using the file explorer on the left, or with the `tree` command.
Run the following from the terminal to see the full structure:

```bash
tree . -L 3
```

??? abstract "Directory contents"

    ```console
    triathlon
    ├── basics
    │   ├── 1-hello.nf
    │   ├── 2-inputs.nf
    │   ├── data
    │   │   ├── greetings-extended.csv
    │   │   └── greetings.csv
    │   ├── main.nf
    │   ├── modules
    │   │   ├── collectGreetings.nf
    │   │   ├── convertToUpper.nf
    │   │   ├── cowpy.nf
    │   │   └── sayHello.nf
    │   └── nextflow.config
    └── nf-core
        └── laptop.config
    ```

The **`basics/`** directory contains everything for Part 1: three workflow scripts of increasing complexity, input data, process modules, and a configuration file.

The **`nf-core/`** directory contains a `laptop.config` file used in Parts 2 and 3 to cap resource usage when running nf-core pipelines locally.

In Part 3 you will install the `seqera` CLI and log in with `seqera login`; no config file is needed in the working directory.

## Readiness checklist

Think you're ready to dive in?

- [ ] I understand the goal of this course and its prerequisites
- [ ] My environment is up and running

If you can check all the boxes, you're good to go.

**To continue to [Part 1: Run Nextflow](./01_run_nextflow.md), click on the arrow in the bottom right corner of this page.**
