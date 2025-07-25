# Orientation

This orientation assumes you have already opened the training environment by clicking on the "Open in GitHub Codespaces" button.
If not, please do so now, ideally in a second browser window or tab so you can refer back to these instructions.

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

## GitHub Codespaces

The GitHub Codespaces environment contains all the software, code and data necessary to work through this training course, so you don't need to install anything yourself.
However, you do need a (free) GitHub account to log in, and you should take a few minutes to familiarize yourself with the interface.

If you have not yet done so, please go through the [Environment Setup](../../envsetup/) mini-course before going any further.

## Working directory

Throughout this training course, we'll be working in the `nextflow-run/` directory.

Change directory now by running this command in the terminal:

```bash
cd nextflow-run/
```

!!!tip

    If for whatever reason you move out of this directory, you can always use the full path to return to it, assuming you're running this within the GitHub Codespaces training environment:

    ```bash
    cd /workspaces/training/nextflow-run
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

If you run this inside `nextflow-run`, you should see the following output:

```console title="Directory contents"
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

**Here's a summary of what you should know to get started:**

- **The `.nf` files** are workflow scripts that are numbered based on what part of the course they're used in.

- **The file `nextflow.config`** is a configuration file that sets minimal environment properties.
  You can ignore it for now.

- **The file `greetings.csv`** contains input data we'll use in most of the course. It is described in Part 2, when we introduce it for the first time.

- **The file `test-params.yaml`** is a file we'll use in Part 3. You can ignore it for now.

**Now, to begin the course, click on the arrow in the bottom right corner of this page.**
