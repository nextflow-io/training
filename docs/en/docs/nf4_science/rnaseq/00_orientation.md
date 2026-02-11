# Getting started

## Start a training environment

To use the pre-built environment we provide on GitHub Codespaces, click the "Open in GitHub Codespaces" button below. For other options, see [Environment options](../../envsetup/index.md).

We recommend opening the training environment in a new browser tab or window (use right-click, ctrl-click or cmd-click depending on your equipment) so that you can read on while the environment loads.
You will need to keep these instructions open in parallel to work through the course.

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

### Environment basics

This training environment contains all the software, code and data necessary to work through the training course, so you don't need to install anything yourself.

The codespace is set up with a VSCode interface, which includes a filesystem explorer, a code editor and a terminal shell.
All instructions given during the course (e.g. 'open the file', 'edit the code' or 'run this command') refer to those three parts of the VScode interface unless otherwise specified.

If you are working through this course by yourself, please acquaint yourself with the [environment basics](../../envsetup/01_setup.md) for further details.

### Version requirements

This training is designed for Nextflow 25.10.2 or later **with the v2 syntax parser ENABLED**.
If you are using a local or custom environment, please make sure you are using the correct settings as documented [here](../../info/nxf_versions.md).

## Get ready to work

Once your codespace is running, there are two things you need to do before diving into the training: set your working directory for this specific course, and take a look at the materials provided.

### Set the working directory

By default, the codespace opens with the work directory set at the root of all training courses, but for this course, we'll be working in the `nf4-science/rnaseq/` directory.

Change directory now by running this command in the terminal:

```bash
cd nf4-science/rnaseq/
```

You can set VSCode to focus on this directory, so that only the relevant files show in the file explorer sidebar:

```bash
code .
```

!!! tip

    If for whatever reason you move out of this directory (e.g. your codespace goes to sleep), you can always use the full path to return to it, assuming you're running this within the Github Codespaces training environment:

    ```bash
    cd /workspaces/training/nf4-science/rnaseq
    ```

Now let's have a look at the contents.

### Explore the materials provided

You can explore the contents of this directory by using the file explorer on the left-hand side of the training workspace.
Alternatively, you can use the `tree` command.

Throughout the course, we use the output of `tree` to represent directory structure and contents in a readable form, sometimes with minor modifications for clarity.

Here we generate a table of contents to the third level down:

```bash
tree . -L 3
```

??? abstract "Directory contents"

    ```console
    .
    ├── data
    │   ├── genome.fa
    │   ├── paired-end.csv
    │   ├── reads
    │   │   ├── ENCSR000COQ1_1.fastq.gz
    │   │   ├── ENCSR000COQ1_2.fastq.gz
    │   │   ├── ENCSR000COQ2_1.fastq.gz
    │   │   ├── ENCSR000COQ2_2.fastq.gz
    │   │   ├── ENCSR000COR1_1.fastq.gz
    │   │   ├── ENCSR000COR1_2.fastq.gz
    │   │   ├── ENCSR000COR2_1.fastq.gz
    │   │   ├── ENCSR000COR2_2.fastq.gz
    │   │   ├── ENCSR000CPO1_1.fastq.gz
    │   │   ├── ENCSR000CPO1_2.fastq.gz
    │   │   ├── ENCSR000CPO2_1.fastq.gz
    │   │   └── ENCSR000CPO2_2.fastq.gz
    │   └── single-end.csv
    ├── nextflow.config
    ├── rnaseq.nf
    └── solutions
        ├── modules
        │   ├── fastqc.nf
        │   ├── fastqc_pe.nf
        │   ├── hisat2_align.nf
        │   ├── hisat2_align_pe.nf
        │   ├── multiqc.nf
        │   ├── trim_galore.nf
        │   └── trim_galore_pe.nf
        ├── rnaseq-2.1.nf
        ├── rnaseq-2.2.nf
        ├── rnaseq-2.3.nf
        ├── rnaseq-3.1.nf
        ├── rnaseq-3.2.nf
        └── rnaseq_pe-3.3.nf
    ```

Click on the colored box to expand the section and view its contents.
We use collapsible sections like this to display expected command output as well as directory and file contents in a concise way.

- **The `rnaseq.nf` file** is a scaffold for a workflow script that you'll build up as you go through the course.

- **The `modules` directory** contains scaffolds for process module that you'll fill in during the course.

- **The file `nextflow.config`** is a configuration file that sets minimal environment properties.
  You can ignore it for now.

- **The `data` directory** contains input data and related resources, described later in the course.

- **The `solutions` directory** contains completed workflow scripts and modules that result from each step of the course.
  They are intended to be used as a reference to check your work and troubleshoot any issues.
  The Part 2 solution can be used as a starting point for Part 3.

## Readiness checklist

Think you're ready to dive in?

- [ ] I understand the goal of this course and its prerequisites
- [ ] My environment is up and running
- [ ] I've set my working directory appropriately

If you can check all the boxes, you're good to go.

**To continue to [Part 1: Method overview](./01_method.md), click on the arrow in the bottom right corner of this page.**
