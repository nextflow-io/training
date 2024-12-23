# Part 1: Hello World

A "Hello World!" is a minimalist example that is meant to demonstrate the basic syntax and structure of a programming language or software framework. The example typically consists of printing the phrase "Hello, World!" to the output device, such as the console or terminal, or writing it to a file.

In this first part of the Hello Nextflow training course, we ease into the topic with a very simple domain-agnostic Hello World example, which we'll progressively build up to demonstrate the usage of foundational Nextflow logic and components.

---

## 0. Warmup: Run Hello World directly

Let's demonstrate this with a simple command that we run directly in the terminal, to show what it does before we wrap it in Nextflow.

### 0.1. Make the terminal say hello

```bash
echo 'Hello World!'
```

### 0.2. Now make it write the text output to a file

```bash
echo 'Hello World!' > output.txt
```

### 0.3. Verify that the output file is there using the `ls` command

```bash
ls
```

### 0.4. Show the file contents

```bash
cat output.txt
```

!!! tip

    In the Gitpod environment, you can also find the output file in the file explorer, and view its contents by clicking on it. Alternatively, you can use the `code` command to open the file for viewing.

    ```bash
    code output.txt
    ```

### Takeaway

You now know how to run a simple command in the terminal that outputs some text, and optionally, how to make it write the output to a file.

### What's next?

Discover what that would look like written as a Nextflow workflow.

---

## 1. Examine the Hello World workflow starter script

As mentioned in the orientation, we provide you with a fully functional if minimalist workflow script named `hello-world.nf` that does the same thing as before (write out 'Hello World!') but with Nextflow.

To get you started, we'll first open up the workflow script so you can get a sense of how it's structured

### 1.1. Examine the overall code structure

Let's open the `hello-world.nf` script in the editor pane.

!!! note

    The file is in the `hello-nextflow` directory, which should be your current working directory.
    You can either click on the file in the file explorer, or type `ls` in the terminal and Cmd+Click (MacOS) or Ctrl+Click (PC) on the file to open it.

```groovy title="hello-world.nf" linenums="1"
#!/usr/bin/env nextflow

/*
 * Use echo to print 'Hello World!' to a file
 */
process sayHello {

    output:
        path 'output.txt'

    script:
    """
    echo 'Hello World!' > output.txt
    """
}

workflow {

    // emit a greeting
    sayHello()
}
```

As you can see, a Nextflow script involves two main types of core components: one or more **processes**, and the **workflow** itself.
Each **process** describes what operation(s) the corresponding step in the pipeline should accomplish, while the **workflow** describes the dataflow logic that connects the various steps.

Let's take a closer look at the **process** block first, then we'll look at the **workflow** block.

### 1.2 The `process` definition

The first block of code describes a **process**.
The process definition starts with the keyword `process`, followed by the process name and finally the process body delimited by curly braces.
The process body must contain a script block which specifies the command to run, which can be anything you would be able to run in a command line terminal.

Here we have a **process** called `sayHello` that writes its **output** to a file named `output.txt`.

```groovy title="hello-world.nf" linenums="3"
/*
 * Use echo to print 'Hello World!' to a file
 */
process sayHello {

    output:
        path 'output.txt'

    script:
    """
    echo 'Hello World!' > output.txt
    """
}
```

This a very minimal process definition that just contains an `output` definition and the `script` to execute.

The `output` definition includes the `path` qualifier, which tells Nextflow this should be handled as a path (includes both directory paths and files).
Another common qualifier is `val`.

!!! note

    The output definition does not _determine_ what output will be created.
    It simply _declares_ what is the expected output, so that Nextflow can look for it once execution is complete.
    This is necessary for verifying that the command was executed successfully and for passing the output to downstream processes if needed.

!!! warning

    This example is brittle because we hardcoded the output filename in two separate places (the script and the output blocks).
    If we change one but not the other, the script will break.
    Later, you'll learn how to use variables to avoid this problem.

In a real-world pipeline, a process usually contains additional blocks such as directives, inputs, and conditional clauses, which we'll introduce later in this training course.

### 1.3 The `workflow` definition

The second block of code describes the **workflow** itself.
The workflow definition starts with the keyword `workflow`, followed by an optional name, then the workflow body delimited by curly braces.

Here we have a **workflow** that consists of one call to the `sayHello` process.

```groovy title="hello-world.nf" linenums="16"
workflow {

    // emit a greeting
    sayHello()
}
```

This a very minimal **workflow** definition.
In a real-world pipeline, the workflow typically contains multiple calls to **processes** connected by **channels**.
You'll learn how to add more processes and connect them by channels in a little bit.

### Takeaway

You now know how a simple Nextflow workflow is structured.

### What's next?

Learn to launch the workflow, monitor execution and find your outputs.

---

## 2. Run the workflow

Looking at code is not nearly as fun as running it, so let's try this out in practice.

### 2.1. Launch the workflow and monitor execution

In the terminal, run the following command.

```bash
nextflow run hello-world.nf
```

You console output should look something like this:

```console title="Output"
 N E X T F L O W   ~  version 24.10.0

 ┃ Launching `hello-world.nf` [reverent_carson] DSL2 - revision: 463b611a35

executor >  local (1)
[1c/7d08e6] sayHello [100%] 1 of 1 ✔
```

Congratulations, you just ran your first Nextflow workflow!

The most important output here is the last line (line 6), which reports that the `sayHello` process was successfully executed once.

Okay, that's great, but where do we find the output?
The `sayHello` process definition said that the output would be sent to a file, but where is it?

### 2.2. Find the output and logs in the `work` directory

When you run Nextflow for the first time in a given directory, it creates a directory called `work` where it will write all files (and symlinks) generated in the course of execution.
Have a look inside; you'll find a subdirectory named with a hash (in order to make it unique; we'll discuss why in a bit), nested two levels deep and containing a handful of log files.

!!! tip

    If you browse the contents of the task subdirectory in the Gitpod's VSCode file explorer, you'll see all these files right away.
    However, these files are set to be invisible in the terminal, so if you want to use `ls` or `tree` to view them, you'll need to set the relevant option for displaying invisible files.

    ```bash
    tree -a work
    ```

    You should see something like this, though the exact subdirectory names will be different on your system.

    ```console title="Directory contents"
    work
    └── 1c
        └── 7d08e685a7aa7060b9c21667924824
            ├── .command.begin
            ├── .command.err
            ├── .command.log
            ├── .command.out
            ├── .command.run
            ├── .command.sh
            └── .exitcode
    ```

You may have noticed that the subdirectory names appeared (in truncated form) in the output from the workflow run, in the line that says:

```console title="Output"
[1c/7d08e6] sayHello [100%] 1 of 1 ✔
```

This tells you what is the subdirectory path for that specific process call (sometimes called task).

!!! note

    Nextflow creates a separate unique subdirectory for each process call.
    It stages the relevant input files, script, and other helper files there, and writes any output files and logs there as well.

If we look inside the subdirectory, we find the following log files:

-   **`.command.begin`**: Metadata related to the beginning of the execution of the process task
-   **`.command.err`**: Error messages (stderr) emitted by the process task
-   **`.command.log`**: Complete log output emitted by the process task
-   **`.command.out`**: Regular output (stdout) by the process task
-   **`.command.sh`**: The command that was run by the process task call
-   **`.exitcode`**: The exit code resulting from the command

[TODO] UPDATE DESCRIPTION TO LOOK AT ACTUAL OUTPUT FILE INSTEAD OF STDOUT

In this case, you can look for your output in the `output.txt` file.
Open it and you will find the `Hello World!` greeting, which was the expected result of our minimalist workflow.

It's also worth having a look at the `.command.sh` file, which tells you what command Nextflow actually executed. In this case it's very straightforward, but later in the course you'll see commands that involve some interpolation of variables. When you're dealing with that, you need to be able to check exactly what was run, especially when troubleshooting an issue.

### Takeaway

You know how to decipher a simple Nextflow script, run it and find the output and logs in the work directory.

### What's next?

Learn how to manage your workflow executions conveniently.

---

## 3. Manage workflow executions

[TODO] A FEW USEFUL TIPS

### 3.1. Re-launch a workflow with `-resume`

Nextflow has an option called `-resume` that allows you to re-run a pipeline you've already launched previously.
When launched with `-resume` any processes that have already been run with the exact same code, settings and inputs will be skipped.
Using this mode means Nextflow will only run processes that are either new, have been modified or are being provided new settings or inputs.

There are two key advantages to doing this:

-   If you're in the middle of developing your pipeline, you can iterate more rapidly since you only effectively have to run the process(es) you're actively working on in order to test your changes.
-   If you're running a pipeline in production and something goes wrong, in many cases you can fix the issue and relaunch the pipeline, and it will resume running from the point of failure, which can save you a lot of time and compute.

To use it, simply add `-resume` to your command:

```bash
nextflow run hello-world.nf -resume
```

The console output should look similar.

```console title="Output"
 N E X T F L O W   ~  version 24.10.0

 ┃ Launching `hello-world.nf` [thirsty_gautier] DSL2 - revision: 6654bc1327

[10/15498d] sayHello [100%] 1 of 1, cached: 1 ✔
```

Notice the additional `cached:` bit in the process status line, which means that Nextflow has recognized that it has already done this work and simply re-used the result from the last run.

!!! note

    When your re-run a pipeline with `resume`, Nextflow does not overwrite any files written to a publishDir directory by any process call that was previously run successfully.

### 3.2. Delete older work directories

[TODO] SHOW HOW TO USE `nextflow clean -before <run_name> -n` (dry run) THEN `nextflow clean -before <run_name>` (do it). Using dry run first because deletion is permanent.

[TODO] Using `-before` bc it feels like the most likely to be commonly useful. Link to other options here: https://www.nextflow.io/docs/latest/reference/cli.html#clean

[TODO] Note on how deleting work directories breaks ability to resume from those directories and deletes outputs so save your outputs!

### Takeaway

You know how to to relaunch a pipeline without repeating steps that were already run in an identical way, and how to use the `nextflow clean` command to clean up old work directories.

### What's next?

Learn to 'publish' outputs to a location outside the work directory for convenience.

---

## 4. Publish outputs

You'll have noticed that the output is buried in a working directory several layers deep.
Nextflow is in control of this directory and we are not supposed to interact with it.

Let's look at how to use the `publishDir` directive for managing this more conveniently.

!!! note

    A newer syntax option had been proposed to make it possible to declare and publish workflow-level outputs, documented [here](https://www.nextflow.io/docs/latest/workflow.html#publishing-outputs).
    This will eventually make using `publishDir` at the process level redundant for completed pipelines.
    However, we expect that `publishDir` will still remain very useful during pipeline development.

### 4.1. Add a `publishDir` directive to the process

To make the output file more accessible, we can utilize the `publishDir` directive.
This directive tells Nextflow to automatically copy the output file to a designated output directory.
It allows us to get easy access to the desired output file without having to drill down into the work directory.

_Before:_

```groovy title="hello-world.nf" linenums="6"
process sayHello {

    output:
        path 'output.txt'
```

_After:_

```groovy title="hello-world.nf" linenums="6"
process sayHello {

    publishDir 'results', mode: 'copy'

    output:
        path 'output.txt'
```

### 4.2. Run the workflow again

```bash
nextflow run hello-world.nf
```

The log output should start looking very familiar:

```console title="Output"
 N E X T F L O W   ~  version 24.10.0

 ┃ Launching `hello-world.nf` [mighty_lovelace] DSL2 - revision: 6654bc1327

executor >  local (1)
[10/15498d] sayHello [100%] 1 of 1 ✔
```

This time, Nextflow will have created a new directory called `results/`.
Our `output.txt` file is in this directory.
If you check the contents it should match the output in our work/task directory.
This is how we move results files outside of the working directories.

!!! note

    It is also possible to set the `publishDir` directive to make a symbolic link to the file instead of actually copying it.
    This is useful when you're dealing with very large files.
    However, if you delete the work directory as part of a cleanup operation, you will lost access to the file, so always make sure you have actual copies of everything you care about before deleting anything.

### Takeaway

You know how to use the `publishDir` directive to move files outside of the Nextflow working directory.

More generally, you know how to interpret a simple Nextflow workflow, manage its execution, and retrieve outputs.

### What's next?

[TODO] LEARN HOW TO USE CHANNELS TO PROVIDE INPUTS TO A PROCESS
