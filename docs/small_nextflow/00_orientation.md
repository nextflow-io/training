# Orientation

This orientation assumes you have already opened the training environment by clicking on the "Open in GitHub Codespaces" button.
If not, please do so now, ideally in a second browser window or tab so you can refer back to these instructions.

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

## GitHub Codespaces

The GitHub Codespaces environment contains all the software, code and data necessary to work through this training course, so you don't need to install anything yourself.
However, you do need a (free) GitHub account to log in, and you should take a few minutes to familiarize yourself with the interface.

If you have not yet done so, please go through the [Environment Setup](../../envsetup/) mini-course before going any further.

## Working directory

Throughout this training course, we'll be working in the `small_nextflow/` directory.

Change directory now by running this command in the terminal:

```bash
cd small_nextflow/
```

!!!tip

    If for whatever reason you move out of this directory, you can always use the full path to return to it, assuming you're running this within the GitHub Codespaces training environment:

    ```bash
    cd /workspaces/training/small_nextflow
    ```

Now let's have a look at the contents of this directory.

## Materials provided

You can explore the contents of this directory by using the file explorer on the left-hand side of the training workspace.
Alternatively, you can use the `ls` command.

Here we list the contents of the directory:

```bash
ls -la
```

If you run this inside `small_nextflow`, you should see a minimal directory structure:

```console title="Directory contents"
.
├── .stuff/
│   ├── cat_me.sh
│   ├── classify.py
│   └── pyproject.toml
└── main.nf
```

**Here's a summary of what you should know to get started:**

- **The `.stuff/` directory** contains helper scripts and configuration files we'll use throughout the workshop.
  You can think of this as a toolbox we'll pull from as we build our workflow.

- **The file `main.nf`** is where we'll write our Nextflow workflow.
  It starts nearly empty, and we'll build it up step by step.

- **The `cat_me.sh` script** fetches random cat images from an API for our workflow to process.

- **The `classify.py` script** is a Python program that uses machine learning to classify images.

- **The `pyproject.toml` file** describes the Python dependencies needed for the classification script.

Throughout this workshop, we'll start with this minimal setup and progressively build a complete image classification workflow.

**Now, to begin the course, click on the arrow in the bottom right corner of this page.**
