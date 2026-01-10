# Getting started

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/G3CV-FcV-rc?si=nyLvwhrSB2m1NPc5&amp;list=PLPZ8WHdZGxmXiHf8B26oB_fTfoKQdhlik" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

/// caption
:fontawesome-brands-youtube:{ .youtube } See [the whole playlist](https://www.youtube.com/playlist?list=PLPZ8WHdZGxmXiHf8B26oB_fTfoKQdhlik) on the Nextflow YouTube channel.

:green_book: The video transcript is available [here](./transcripts/00_orientation.md).
///

!!! tip

    The YouTube videos have some super powers!

    - :fontawesome-solid-closed-captioning: High quality (manually curated) captions / subtitles. Switch them on with the :material-subtitles: icon
    - :material-bookmark: Video chapters in the timeline that correspond to page headings.

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

By default, the codespace opens with the work directory set at the root of all training courses, but for this course, we'll be working in the `hello-nextflow/` directory.

Change directory now by running this command in the terminal:

```bash
cd hello-nextflow/
```

!!! tip

    If for whatever reason you move out of this directory (e.g. your codespace goes to sleep), you can always use the full path to return to it, assuming you're running this within the Github Codespaces training environment:

    ```bash
    cd /workspaces/training/hello-nextflow
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

If you run this inside `hello-nextflow`, you should see the following output.

??? abstract "Directory contents"

    ```console
    .
    ├── greetings.csv
    ├── hello-channels.nf
    ├── hello-config.nf
    ├── hello-containers.nf
    ├── hello-modules.nf
    ├── hello-workflow.nf
    ├── hello-world.nf
    ├── nextflow.config
    ├── solutions
    │   ├── 1-hello-world
    │   ├── 2-hello-channels
    │   ├── 3-hello-workflow
    │   ├── 4-hello-modules
    │   ├── 5-hello-containers
    │   └── 6-hello-config
    └── test-params.json

    7 directories, 9 files
    ```

!!! tip

    We use collapsible sections like this to include expected command output in a concise way.
    Click on the colored box to expand the section and view its contents.

### Content guide

- **The `.nf` files** are workflow scripts that are named based on what part of the course they're used in.

- **The file `nextflow.config`** is a configuration file that sets minimal environment properties.
  You can ignore it for now.

- **The file `greetings.csv`** contains input data we'll use in most of the course. It is described in Part 1, when we introduce it for the first time.

- **The file `test-params.json`** is a file we'll use in Part 6. You can ignore it for now.

- **The `solutions` directory** contains the completed workflow scripts that result from each step of the course.
  They are intended to be used as a reference to check your work and troubleshoot any issues.
  The name and number in the filename correspond to the step of the relevant part of the course.
  For example, the file `hello-world-4.nf` is the expected result of completing steps 1 through 4 of Part 1: Hello World.

## Readiness checklist

Think you're ready to dive in?

- [ ] I understand the goal of this course and its prerequisites
- [ ] My environment is up and running
- [ ] I've set my working directory appropriately

If you can check all the boxes, you're good to go.

**To continue to Part 1, click on the arrow in the bottom right corner of this page.**
