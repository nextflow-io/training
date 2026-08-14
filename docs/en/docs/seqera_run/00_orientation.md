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

## Get ready to work

Once your codespace is running, there are two things to do before diving in: set your working directory, and take a look at the materials provided.

### Set the working directory

By default, the codespace opens at the root of all training courses.
For this course, change to the `seqera-run/` directory:

```bash
cd seqera-run/
```

Then set VSCode to focus on this directory, so only the relevant files appear in the file explorer sidebar:

```bash
code .
```

!!! tip

    If for whatever reason you move out of this directory (e.g. your codespace goes to sleep), you can always use the full path to return to it, assuming you're running this within the Github Codespaces training environment:

    ```bash
    cd /workspaces/training/seqera-run
    ```

### Explore the materials provided

You can explore the course materials using the file explorer on the left, or with the `tree` command.
Run the following from the terminal to see the full structure:

```bash
tree .
```

??? abstract "Directory contents"

    ```console
    .
    └── .seqera_config
    ```

The **`.seqera_config`** file is a stub you will fill in during section 3 to configure the `tw` CLI with your Seqera access token and workspace.

## Readiness checklist

Think you're ready to dive in?

- [ ] I understand the goal of this course and its prerequisites
- [ ] My environment is up and running
- [ ] I've set my working directory appropriately

If you can check all the boxes, you're good to go.

**To continue to [Part 1: Launch pipelines from the web interface](./01_run_with_seqera.md), click on the arrow in the bottom right corner of this page.**
