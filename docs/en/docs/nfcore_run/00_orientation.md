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

This training works with **Nextflow 25.10.2** or later **with the v2 syntax parser DISABLED**.

#### If you are using our training environment:

You MUST run the following command before going any further:

```bash
export NXF_SYNTAX_PARSER=v1
```

#### If you are using a local or custom environment:

Please make sure you are using the correct settings as documented [here](../info/nxf_versions.md).

The training additionally requires **nf-core tools 3.5.2**.
If you use a different version of nf-core tooling, you may have difficulty following along.

You can check what version is installed in your environment using the command `nf-core --version`.

## Get ready to work

Once your codespace is running, there are two things you need to do before diving into the training: set your working directory for this specific course, and take a look at the materials provided.

### Set the working directory

By default, the codespace opens with the work directory set at the root of all training courses, but for this course, we'll be working in the `nfcore-run/` directory.

Change directory now by running this command in the terminal:

```bash
cd nfcore-run/
```

!!! tip

    If for whatever reason you move out of this directory (e.g. your codespace goes to sleep), you can always use the full path to return to it, assuming you're running this within the Github Codespaces training environment:

    ```bash
    cd /workspaces/training/nfcore-run
    ```

Next, explore the contents of this directory.

### Explore the materials provided

You can explore the contents of this directory by using the file explorer on the left-hand side of the training workspace.
Alternatively, you can use the `tree` command.

```bash
tree .
```

??? abstract "Directory contents"

    ```console
    .
    └── laptop.config
    ```

- **The `laptop.config` file** is a configuration file we'll use in section 4 to cap resource usage when running a production-scale pipeline locally.
  You can ignore it until then.

## Readiness checklist

Think you're ready to dive in?

- [ ] I understand the goal of this course and its prerequisites
- [ ] My environment is up and running
- [ ] I've made certain that the syntax parser is set to **v1**
- [ ] I've set my working directory appropriately

If you can check all the boxes, you're good to go.

**To continue to Part 1, click on the arrow in the bottom right corner of this page.**
