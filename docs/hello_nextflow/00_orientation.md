# Orientation

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

## GitHub Codespaces

The GitHub Codespaces environment contains all the software, code and data necessary to work through this training course, so you don't need to install anything yourself.
However, you do need a (free) GitHub account to log in, and you should take a few minutes to familiarize yourself with the interface.

If you have not yet done so, please go through the [Environment Setup](../../envsetup/) mini-course before going any further.

## Working directory

Throughout this training course, we'll be working in the `hello-nextflow/` directory.

Change directory now by running this command in the terminal:

```bash
cd hello-nextflow/
```

!!!tip

    If for whatever reason you move out of this directory, you can always use the full path to return to it, assuming you're running this within the Github Codespaces training environment:

    ```bash
    cd /workspaces/training/hello-nextflow
    ```

Now let's have a look at the contents of this directory.

## Materials provided

You can explore the contents of this directory by using the file explorer on the left-hand side of the training workspace.
Alternatively, you can use the `tree` command.

Throughout the course, we use the output of `tree` to represent directory structure and contents in a readable form, sometimes with minor modifications for clarity.

Here we generate a table of contents to the second level down:

```bash
tree . -L 2
```

If you run this inside `hello-nextflow`, you should see the following output:

```console title="Directory contents"
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

**Here's a summary of what you should know to get started:**

- **The `.nf` files** are workflow scripts that are named based on what part of the course they're used in.

- **The file `nextflow.config`** is a configuration file that sets minimal environment properties.
  You can ignore it for now.

- **The file `greetings.csv`** contains input data we'll use in most of the course. It is described in Part 1, when we introduce it for the first time.

- **The file `test-params.json`** is a file we'll use in Part 6. You can ignore it for now.

- **The `solutions` directory** contains the completed workflow scripts that result from each step of the course.
  They are intended to be used as a reference to check your work and troubleshoot any issues.
  The name and number in the filename correspond to the step of the relevant part of the course.
  For example, the file `hello-world-4.nf` is the expected result of completing steps 1 through 4 of Part 1: Hello World.

**Now, to begin the course, click on the arrow in the bottom right corner of this page.**
