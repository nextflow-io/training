# Orientation

The GitHub Codespaces environment contains all the software, code and data necessary to work through this training course, so you don't need to install anything yourself.
However, you do need a (free) account to log in, and you should take a few minutes to familiarize yourself with the interface.

If you have not yet done so, please follow [this link](../../envsetup/) before going any further.

## Materials provided

Throughout this training course, we'll be working in the `side-quests/` directory.
This directory contains all the code files, test data and accessory files you will need.

Feel free to explore the contents of this directory; the easiest way to do so is to use the file explorer on the left-hand side of the GitHub Codespaces workspace.
Alternatively, you can use the `tree` command.
Throughout the course, we use the output of `tree` to represent directory structure and contents in a readable form, sometimes with minor modifications for clarity.

Here we generate a table of contents to the second level down:

```bash
tree . -L 2
```

If you run this inside `side-quests`, you should see the following output:

```console title="Directory contents"
.
├── metadata
├── nf-core
├── nf-test
├── solutions
├── splitting_and_grouping
└── workflows_of_workflows
```

**Here's a summary of what you should know to get started:**

- **Each directory corresponds to an individual side quest.**
  Their contents are detailed on the corresponding side quest's page.

- **The `solutions` directory** contains the completed workflow and/or module scripts that result from running through various steps of each side quest.
  They are intended to be used as a reference to check your work and troubleshoot any issues.

!!!tip

    If for whatever reason you move out of this directory, you can always run this command to return to it:

    ```bash
    cd /workspaces/training/side-quests
    ```

Now, to begin the course, click on the arrow in the bottom right corner of this page.
