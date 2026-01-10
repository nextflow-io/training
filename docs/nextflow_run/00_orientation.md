# Getting started

## Start a training environment

To use the pre-built environment we provide on GitHub Codespaces, click the "Open in GitHub Codespaces" button below. For other options, see [Environment options](../envsetup/index.md).

We recommend opening the training environment in a new browser tab or window (use right-click, ctrl-click or cmd-click depending on your equipment) so that you can read on while the environment loads.
You will need to keep these instructions open in parallel to work through the course.

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

### Environment basics

This training environment contains all the software, code and data necessary to work through the training course, so you don't need to install anything yourself.

The codespace is set up with a VSCode interface, which includes a filesystem explorer, a code editor and a terminal shell.
All instructions given during the course (e.g. 'open the file', 'edit the code' or 'run this command') refer to those three parts of the VScode interface unless otherwise specified.

If you are working through this course by yourself, please acquaint yourself with the [environment basics](../envsetup/01_setup.md) for further details.

### Version requirements

This training is designed for Nextflow 25.10.2 or later **with the v2 syntax parser ENABLED**.
If you are using a local or custom environment, please make sure you are using the correct settings as documented [here](../nxf_versions.md).

## Get ready to work

Once your codespace is running, there are two things you need to do before diving into the training: set your working directory for this specific course, and take a look at the materials provided.

### Set the working directory

By default, the codespace opens with the work directory set at the root of all training courses, but for this course, we'll be working in the `nextflow-run/` directory.

Change directory now by running this command in the terminal:

```bash
cd nextflow-run/
```

!!! tip

    If for whatever reason you move out of this directory (e.g. your codespace goes to sleep), you can always use the full path to return to it, assuming you're running this within the Github Codespaces training environment:

    ```bash
    cd /workspaces/training/nextflow-run
    ```

Now let's have a look at the contents.

### Explore the materials provided

You can explore the contents of this directory by using the file explorer on the left-hand side of the training workspace.
Alternatively, you can use the `tree` command.

Throughout the course, we use the output of `tree` to represent directory structure and contents in a readable form, sometimes with minor modifications for clarity.

Here we generate a table of contents to the second level down:

```bash
tree . -L 2
```

??? abstract "Directory contents"

    ```console
    .
    ├── 1-hello.nf
    ├── 2a-inputs.nf
    ├── 2b-multistep.nf
    ├── 2c-modules.nf
    ├── 2d-container.nf
    ├── 3-main.nf
    ├── modules
    │   ├── collectGreetings.nf
    │   ├── convertToUpper.nf
    │   ├── cowpy.nf
    │   └── sayHello.nf
    ├── nextflow.config
    └── test-params.yaml

    1 directory, 12 files
    ```

Click on the colored box to expand the section and view its contents.
We use collapsible sections like this to include expected command output in a concise way.

### Content guide

- **The `.nf` files** are workflow scripts that are numbered based on what part of the course they're used in.

- **The file `nextflow.config`** is a configuration file that sets minimal environment properties.
  You can ignore it for now.

- **The file `greetings.csv`** contains input data we'll use in most of the course. It is described in Part 2, when we introduce it for the first time.

- **The file `test-params.yaml`** is a file we'll use in Part 3. You can ignore it for now.

## Readiness checklist

Think you're ready to dive in?

- [ ] I understand the goal of this course and its prerequisites
- [ ] My environment is up and running
- [ ] I've set my working directory appropriately

If you can check all the boxes, you're good to go.

**To continue to Part 1, click on the arrow in the bottom right corner of this page.**
