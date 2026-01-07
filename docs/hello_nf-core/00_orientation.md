# Getting started

To start the course, launch the training environment by clicking the "Open in GitHub Codespaces" button below.
We recommend opening the training environment in a new browser tab (use right-click, ctrl-click or cmd-click depending on your equipment) so that you can read on while the environment loads.
You will need to keep these instructions open in parallel.

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

## Training environment

Our training environment runs on GitHub Codespaces (free Github account required) and contains all the software, code and data necessary to work through this training course, so you don't need to install anything yourself.

The codespace is set up with a VSCode interface, which includes a filesystem explorer, a code editor and a terminal shell.
All instructions given during the course (e.g. 'open the file', 'edit the code' or 'run this command') refer to those three parts of the VScode interface unless otherwise specified.

If you are working through this course by yourself, please go through the [Environment Setup](../../envsetup/) mini-course for further details before going any further.

!!! warning

    This training is designed for nf-core tools version 3.4.1, which should be the version installed in the codespace we provide.
    If you use a different version of nf-core tooling, you may have difficulty following along.

    You can check what version is installed using the command`nf-core --version`.

## Get ready to work

Once your codespace is running, there are two things you need to do before diving into the training: set your working directory for this specific course, and take a look at the materials provided.

### Set the working directory

By default, the codespace opens with the work directory set at the root of all training courses, but for this course, we'll be working in the `hello-nf-core/` directory.

Change directory now by running this command in the terminal:

```bash
cd hello-nf-core/
```

!!! tip

    If for whatever reason you move out of this directory (e.g. your codespace goes to sleep), you can always use the full path to return to it, assuming you're running this within the Github Codespaces training environment:

    ```bash
    cd /workspaces/training/hello-nf-core
    ```

Now let's have a look at the contents of this directory.

### Check out the materials provided

You can explore the contents of this directory by using the file explorer on the left-hand side of the training workspace.
Alternatively, you can use the `tree` command.

Throughout the course, we use the output of `tree` to represent directory structure and contents in a readable form, sometimes with minor modifications for clarity.

Here we generate a table of contents to the second level down:

```bash
tree . -L 2
```

If you run this inside `hello-nf-core`, you should see the following output.

??? abstract "Directory contents"

    ```console
    .
    ├── greetings.csv
    ├── original-hello
    │   ├── hello.nf
    │   ├── modules
    │   └── nextflow.config
    └── solutions
        ├── composable-hello
        ├── core-hello-part2
        ├── core-hello-part3
        ├── core-hello-part4
        ├── core-hello-part5
        └── core-hello-start

    9 directories, 3 files
    ```

!!! note

    We use collapsible sections like this to include expected command output in a concise way.
    Click on the colored box to expand the section and view its contents.

**Content guide:**

- **The `greetings.csv` file** is a CSV containing some minimal columnar data we use for testing purposes.

- **The `original-hello` directory** contains a copy of the source code produced by working through the complete Hello Nextflow training series (with Docker enabled).

- **The `solutions` directory** contains the completed workflow scripts that result from each step of the course.
  They are intended to be used as a reference to check your work and troubleshoot any issues.

## Readiness checklist

Think you're ready to dive in?

- [ ] I understand the goal of this course and its prerequisites
- [ ] My codespace is up and running
- [ ] I've set my working directory appropriately

If you can check all the boxes, you're good to go.

**To continue to Part 1, click on the arrow in the bottom right corner of this page.**
