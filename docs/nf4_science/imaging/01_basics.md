# Part 1: Run basic operations

In this first part of the Nextflow for Bioimaging training course, we'll use a very basic domain-agnostic Hello World example to demonstrate essential operations and point out the corresponding Nextflow code components.

## 1. Run the workflow

We provide you with a workflow script named `hello-world.nf` that takes an input via a command-line argument named `--greeting` and produces a text file containing that greeting.
We're not going to look at the code yet; first let's see what it looks like to run it.

### 1.1. Launch the workflow and monitor execution

In the terminal, run the following command:

```bash
nextflow run hello-world.nf --greeting 'Hello World!'
```

Your console output should look something like this:

```console title="Output" linenums="1"
 N E X T F L O W   ~  version 25.04.3

Launching `hello-world.nf` [goofy_torvalds] DSL2 - revision: c33d41f479

executor >  local (1)
[a3/7be2fa] sayHello | 1 of 1 ✔
```

Congratulations, you just ran your first Nextflow workflow!

The most important output here is the last line (line 6):

```console title="Output" linenums="6"
[a3/7be2fa] sayHello | 1 of 1 ✔
```

This tells us that the `sayHello` process was successfully executed once (`1 of 1 ✔`).

That's great, but you may be wondering: where is the output?

### 1.2. Find the output file in the `results` directory

This workflow is configured to publish its output to a directory called `results`.
If you look at your current directory, you will see that when you ran the workflow, Nextflow created a new directory called `results`, which contains a file called `output.txt`.

```console title="results/" linenums="1"
results
└── output.txt
```

Open the file; the contents should match the greeting you specified on the command line.

<details>
  <summary>File contents</summary>

```console title="results/output.txt" linenums="1"
Hello World!
```

</details>

That's great, our workflow did what it was supposed to do!

However, be aware that the 'published' result is a copy (or in some cases a symlink) of the actual output produced by Nextflow when it executed the workflow.

So now, we are going to peek under the hood to see where Nextflow actually executed the work.

!!! warning

    Not all workflows will be set up to publish outputs to a results directory, and/or the directory name may be different.
    A little further in this section, we will show you how to find out where this behavior is specified.

### 1.3. Find the original output and logs in the `work/` directory

When you run a workflow, Nextflow creates a distinct 'task directory' for every single invocation of each process in the workflow (=every step in the pipeline).
For each one, it will stage the necessary inputs, execute the relevant instruction(s) and write outputs and log files within that one directory, which is named automatically using a hash in order to make it unique.

All of these task directories will live under a directory called `work` within your current directory (where you're running the command).

That may sound confusing, so let's see what that looks like in practice.

Going back to the console output for the workflow we ran earlier, we had this line:

```console title="Excerpt of command output" linenums="6"
[a3/7be2fa] sayHello | 1 of 1 ✔
```

See how the line starts with `[a3/7be2fa]`?
That is a truncated form of the task directory path for that one process call, and tells you where to find the output of the `sayHello` process call within the `work/` directory path.

You can find the full path by typing the following command (replacing `a3/7be2fa` with what you see in your own terminal) and pressing the tab key to autocomplete the path or adding an asterisk:

```bash
tree work/a3/7be2fa*
```

This should yield the full path directory path: `work/a3/7be2fa7be2fad5e71e5f49998f795677fd68`

Let's take a look at what's in there.

!!! Tip

    If you browse the contents of the task subdirectory in the VSCode file explorer, you'll see all the files right away.
    However, the log files are set to be invisible in the terminal, so if you want to use `ls` or `tree` to view them, you'll need to set the relevant option for displaying invisible files.

    ```bash
    tree -a work
    ```

The exact subdirectory names will be different on your system.

<details>
  <summary>Directory contents</summary>

```console title="work/"
work
└── a3
    └── 7be2fad5e71e5f49998f795677fd68
        ├── .command.begin
        ├── .command.err
        ├── .command.log
        ├── .command.out
        ├── .command.run
        ├── .command.sh
        ├── .exitcode
        └── output.txt
```

</details>

You should immediately recognize the `output.txt` file, which is in fact the original output of the `sayHello` process that got published to the `results` directory.
If you open it, you will find the `Hello World!` greeting again.

<details>
  <summary>File contents of output.txt</summary>

```console title="work/a3/7be2fa7be2fad5e71e5f49998f795677fd68/output.txt" linenums="1"
Hello World!
```

</details>

So what about all those other files?

These are the helper and log files that Nextflow wrote as part of the task execution:

- **`.command.begin`**: Sentinel file created as soon as the task is launched.
- **`.command.err`**: Error messages (`stderr`) emitted by the process call
- **`.command.log`**: Complete log output emitted by the process call
- **`.command.out`**: Regular output (`stdout`) by the process call
- **`.command.run`**: Full script run by Nextflow to execute the process call
- **`.command.sh`**: The command that was actually run by the process call
- **`.exitcode`**: The exit code resulting from the command

The `.command.sh` file is especially useful because it shows you the main command Nextflow executed not including all the bookkeeping and task/environment setup.

<details>
  <summary>File contents</summary>

```console title="work/a3/7be2fa7be2fad5e71e5f49998f795677fd68/.command.sh" linenums="1"
#!/bin/bash -ue
echo 'Hello World!' > output.txt

```

</details>

!!! Tip

    When something goes wrong and you need to troubleshoot what happened, it can be useful to look at the `command.sh` script to check exactly what command Nextflow composed based on the workflow instructions, variable interpolation and so on.

### 1.4. Optional exercise: re-run with different greetings

Try re-running the workflow a few times with different values for the `--greeting` argument, then look at both the contents of the `results/` directory and the task directories.

Observe how the outputs and logs of isolated task directories are preserved, whereas the contents of the `results` directory are overwritten by the output of subsequent executions.

### Takeaway

You know how to run a simple Nextflow script, monitor its execution and find its outputs.

### What's next?

Learn how to read a basic Nextflow script and identify how its components relate to its functionality.

---

## 2. Examine the Hello World workflow starter script

What we did there was basically treating the workflow script like a black box.
Now that we've seen what it does, let's open the box and look inside.

_The goal here is not to memorize the syntax of Nextflow code, but to form some basic intuition of what are the main components and how they are organized._

### 2.1. Examine the overall code structure

Let's open the `hello-world.nf` script in the editor pane.

<details>
  <summary>Code</summary>

```groovy title="hello-world.nf" linenums="1"
#!/usr/bin/env nextflow

/*
 * Use echo to print a greeting to a file
 */
process sayHello {

    publishDir 'results', mode: 'copy'

    input:
        val greeting

    output:
        path 'output.txt'

    script:
    """
    echo '$greeting' > output.txt
    """
}

workflow {

    // emit a greeting
    sayHello(params.greeting)
}
```

</details>

A Nextflow script involves two main types of core components: one or more **processes**, and the **workflow** itself.
Each **process** describes what operation(s) the corresponding step in the pipeline should accomplish, while the **workflow** describes the dataflow logic that connects the various steps.

Let's take a closer look at the **process** block first, then we'll look at the **workflow** block.

### 2.2. The `process` definition

The first block of code describes a **process**.
The process definition starts with the keyword `process`, followed by the process name and finally the process body delimited by curly braces.
The process body must contain a script block which specifies the command to run, which can be anything you would be able to run in a command line terminal.

Here we have a **process** called `sayHello` that takes an **input** variable called `greeting` and writes its **output** to a file named `output.txt`.

<details>
  <summary>Code</summary>

```groovy title="hello-world.nf" linenums="3"
/*
 * Use echo to print a greeting to a file
 */
process sayHello {

    publishDir 'results', mode: 'copy'

    input:
        val greeting

    output:
        path 'output.txt'

    script:
    """
    echo '$greeting' > output.txt
    """
}
```

</details>

This is a very minimal process definition that just contains an `input` definition, an `output` definition and the `script` to execute.

The `input` definition includes the `val` qualifier, which tells Nextflow to expect a value of some kind (can be a string, a number, whatever).

The `output` definition includes the `path` qualifier, which tells Nextflow this should be handled as a path (includes both directory paths and files).

!!! Tip

    The output definition does not _determine_ what output will be created.
    It simply _declares_ where to find the expected output file(s), so that Nextflow can look for it once execution is complete.

    This is necessary for verifying that the command was executed successfully and for passing the output to downstream processes if needed.
    Output produced that doesn't match what is declared in the output block will not be passed to downstream processes.

In a real-world pipeline, a process usually contains additional information such as process directives, which we'll introduce in a little bit.

### 2.3. The `workflow` definition

The second block of code describes the **workflow** itself.
The workflow definition starts with the keyword `workflow`, followed by an optional name, then the workflow body delimited by curly braces.

Here we have a **workflow** that consists of one call to the `sayHello` process, which takes an input, `params.greeting`, which holds the value we gave to the `--greeting` parameter.

```groovy title="hello-world.nf" linenums="22"
workflow {

    // emit a greeting
    sayHello(params.greeting)
}
```

This is a very minimal **workflow** definition.
In a real-world pipeline, the workflow typically contains multiple calls to **processes** connected by **channels**, and there may be default values set up for the variable inputs.

We'll see this in action when we run nf-core/molkart in Part 2 of the course.

### 2.4. The `params` system of command-line parameters

The `params.greeting` we provide to the `sayHello()` process call is a neat bit of Nextflow code and is worth spending an extra minute on.

As mentioned above, that's how we pass the value of the `--greeting` command-line parameter to the `sayHello()` process call.
In fact, simply declaring `params.someParameterName` will enable us to give the workflow a parameter named `--someParameterName` from the command-line.

!!! Tip

    These workflow parameters declared using the `params` system always take two dashes (`--`).
    This distinguishes them from Nextflow-level parameters, which only take one dash (`-`).

### Takeaway

You now know how a simple Nextflow workflow is structured, and how the basic components relate to its functionality.

### What's next?

Learn to manage your workflow executions conveniently.

---

## 3. Manage workflow executions

Knowing how to launch workflows and retrieve outputs is great, but you'll quickly find there are a few other aspects of workflow management that will make your life easier.

Here we show you how to take advantage of the `resume` feature for when you need to re-launch the same workflow, how to inspect the execution logs with `nextflow log`, and how to delete older work directories with `nextflow clean`.

### 3.1. Re-launch a workflow with `-resume`

Sometimes, you're going to want to re-run a pipeline that you've already launched previously without redoing any work that was already completed successfully.

Nextflow has an option called `-resume` that allows you to do this.
Specifically, in this mode, any processes that have already been run with the exact same code, settings and inputs will be skipped.
This means Nextflow will only run processes that you've added or modified since the last run, or to which you're providing new settings or inputs.

There are two key advantages to doing this:

- If you're in the middle of developing a pipeline, you can iterate more rapidly since you only have to run the process(es) you're actively working on in order to test your changes.
- If you're running a pipeline in production and something goes wrong, in many cases you can fix the issue and relaunch the pipeline, and it will resume running from the point of failure, which can save you a lot of time and compute.

To use it, simply add `-resume` to your command and run it:

```bash
nextflow run hello-world.nf --greeting 'Hello World!' -resume
```

The console output should look similar.

<details>
  <summary>Command output</summary>

```console linenums="1"
 N E X T F L O W   ~  version 25.04.3

Launching `hello-world.nf` [tiny_noyce] DSL2 - revision: c33d41f479

[a3/7be2fa] process > sayHello [100%] 1 of 1, cached: 1 ✔
```

</details>

Look for the `cached:` bit that has been added in the process status line (line 5), which means that Nextflow has recognized that it has already done this work and simply reused the result from the previous successful run.

You can also see that the work subdirectory hash is the same as in the previous run.
Nextflow is literally pointing you to the previous execution and saying "I already did that over there."

!!! Tip

    When your re-run a pipeline with `resume`, Nextflow does not overwrite any files written to a `publishDir` directory by any process call that was previously run successfully.

### 3.2. Inspect the log of past executions

Whenever you launch a nextflow workflow, a line gets written to a log file called `history`, under a hidden directory called `.nextflow` in the current working directory.

A more convenient way to access this information is to use the `nextflow log` command.

```bash
nextflow log
```

This will output the contents of the log file to the terminal, showing you the timestamp, run name, status, and full command line for every Nextflow run that has been launched from within the current working directory.

### 3.3. Delete older work directories

During the development process, you'll typically run your draft pipelines a large number of times, which can lead to an accumulation of very many files across many subdirectories.
Since the subdirectories are named randomly, it is difficult to tell from their names what are older vs. more recent runs.

Nextflow includes a convenient `clean` subcommand that can automatically delete the work subdirectories for past runs that you no longer care about, with several [options](https://www.nextflow.io/docs/latest/reference/cli.html#clean) to control what will be deleted.

You can use the Nextflow log to look up a run based on its timestamp and/or command line, then use `nextflow clean -before <run_name> -f` to delete work directories from earlier runs.

!!! Warning

    Deleting work subdirectories from past runs removes them from Nextflow's cache and deletes any outputs that were stored in those directories.
    That means it breaks Nextflow's ability to resume execution without re-running the corresponding processes.

    You are responsible for saving any outputs that you care about or plan to rely on! If you're using the `publishDir` directive for that purpose, make sure to use the `copy` mode, not the `symlink` mode.

### Takeaway

You know how to relaunch a pipeline without repeating steps that were already run in an identical way, inspect the execution log, and use the `nextflow clean` command to clean up old work directories.

### What's next?

Now that you understand basic Nextflow operations, you're ready to run a real bioimaging pipeline with nf-core/molkart.
