# Orientation

## GitHub Codespaces

The GitHub Codespaces environment contains all the software, code and data necessary to work through this training course, so you don't need to install anything yourself.
However, you do need a (free) GitHub account to log in, and you should take a few minutes to familiarize yourself with the interface.

If you have not yet done so, please go through the [Environment Setup](../../envsetup/) mini-course before going any further.

!!! warning

    This training is designed for nf-core tools version 3.4.1, which should be the version installed in the codespace. If you use a different version of nf-core tooling you may have difficulty following along.

    You can check what version is installed using the command`nf-core --version`.

## Working directory

Throughout this training course, we'll be working in the `hello-nf-core/` directory.

Change directory now by running this command in the terminal:

```bash
cd hello-nf-core/
```

!!! tip

    If for whatever reason you move out of this directory, you can always use the full path to return to it, assuming you're running this within the Github Codespaces training environment:

    ```bash
    cd /workspaces/training/hello-nf-core
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

If you run this inside `hello-nf-core`, you should see the following output:

```console title="Directory contents"
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
    └── core-hello-start

8 directories, 3 files
```

**Here's a summary of what you should know to get started:**

- **The `greetings.csv` file** is a CSV containing some minimal columnar data we use for testing purposes.

- **The `original-hello` directory** contains a copy of the source code produced by working through the complete Hello Nextflow training series (with Docker enabled).

- **The `solutions` directory** contains the completed workflow scripts that result from each step of the course.
  They are intended to be used as a reference to check your work and troubleshoot any issues.

**Now, to begin the course, click on the arrow in the bottom right corner of this page.**
