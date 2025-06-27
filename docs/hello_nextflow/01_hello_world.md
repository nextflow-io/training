# Part 1: Hello World

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/8X2hHI-9vms?si=F0t9LFYLjAWoyRXj&amp;list=PLPZ8WHdZGxmXiHf8B26oB_fTfoKQdhlik" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

/// caption
:fontawesome-brands-youtube:{ .youtube } See [the whole playlist](https://www.youtube.com/playlist?list=PLPZ8WHdZGxmXiHf8B26oB_fTfoKQdhlik) on the Nextflow YouTube channel.

:green_book: The video transcript is available [here](./transcripts/01_hello_world.md).
///

In this first part of the Hello Nextflow training course, we ease into the topic with a very basic domain-agnostic Hello World example, which we'll progressively build up to demonstrate the usage of foundational Nextflow logic and components.

!!! note

    A "Hello World!" is a minimalist example that is meant to demonstrate the basic syntax and structure of a programming language or software framework. The example typically consists of printing the phrase "Hello, World!" to the output device, such as the console or terminal, or writing it to a file.

---

## 0. Warmup: Run Hello World directly

Let's demonstrate this with a simple command that we run directly in the terminal, to show what it does before we wrap it in Nextflow.

!!! tip

    Remember that you should now be inside the `hello-nextflow/` directory as described in the Orientation.

### 0.1. Make the terminal say hello

```bash
echo 'Hello World!'
```

This outputs the text 'Hello World' to the terminal.

```console title="Output"
Hello World!
```

### 0.2. Now make it write the text output to a file

```bash
echo 'Hello World!' > output.txt
```

This does not output anything to the terminal.

```console title="Output"

```

### 0.3. Show the file contents

```bash
cat output.txt
```

The text 'Hello World' is now in the output file we specified.

```console title="output.txt" linenums="1"
Hello World!
```

!!! tip

    In the training environment, you can also find the output file in the file explorer, and view its contents by clicking on it. Alternatively, you can use the `code` command to open the file for viewing.

    ```bash
    code output.txt
    ```

### Takeaway

You now know how to run a simple command in the terminal that outputs some text, and optionally, how to make it write the output to a file.

### What's next?

Find out what that would look like written as a Nextflow workflow.

---

## 1. Examine the Hello World workflow starter script

As mentioned in the orientation, we provide you with a fully functional if minimalist workflow script named `hello-world.nf` that does the same thing as before (write out 'Hello World!') but with Nextflow.

To get you started, we'll first open up the workflow script so you can get a sense of how it's structured.

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

### 1.2. The `process` definition

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

This is a very minimal process definition that just contains an `output` definition and the `script` to execute.

The `output` definition includes the `path` qualifier, which tells Nextflow this should be handled as a path (includes both directory paths and files).
Another common qualifier is `val`.

!!! note

    The output definition does not _determine_ what output will be created.
    It simply _declares_ what is the expected output, so that Nextflow can look for it once execution is complete.
    This is necessary for verifying that the command was executed successfully and for passing the output to downstream processes if needed. Output produced that doesn't match what is declared in the output block will not be passed to downstream processes.

!!! warning

    This example is brittle because we hardcoded the output filename in two separate places (the script and the output blocks).
    If we change one but not the other, the script will break.
    Later, you'll learn how to use variables to avoid this problem.

In a real-world pipeline, a process usually contains additional blocks such as directives and inputs, which we'll introduce in a little bit.

### 1.3. The `workflow` definition

The second block of code describes the **workflow** itself.
The workflow definition starts with the keyword `workflow`, followed by an optional name, then the workflow body delimited by curly braces.

Here we have a **workflow** that consists of one call to the `sayHello` process.

```groovy title="hello-world.nf" linenums="17"
workflow {

    // emit a greeting
    sayHello()
}
```

This is a very minimal **workflow** definition.
In a real-world pipeline, the workflow typically contains multiple calls to **processes** connected by **channels**, and the processes expect one or more variable **input(s)**.

You'll learn how to add variable inputs later in this training module; and you'll learn how to add more processes and connect them by channels in Part 3 of this course.

### Takeaway

You now know how a simple Nextflow workflow is structured.

### What's next?

Learn to launch the workflow, monitor execution and find your outputs.

---

## 2. Run the workflow

Looking at code is not nearly as fun as running it, so let's try this out in practice.

### 2.1. Launch the workflow and monitor execution

In the terminal, run the following command:

```bash
nextflow run hello-world.nf
```

You console output should look something like this:

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

Importantly, this line also tells you where to find the output of the `sayHello` process call.
Let's look at that now.

### 2.2. Find the output and logs in the `work` directory

When you run Nextflow for the first time in a given directory, it creates a directory called `work` where it will write all files (and any symlinks) generated in the course of execution.

Within the `work` directory, Nextflow organizes outputs and logs per process call.
For each process call, Nextflow creates a nested subdirectory, named with a hash in order to make it unique, where it will stage all necessary inputs (using symlinks by default), write helper files, and write out logs and any outputs of the process.

The path to that subdirectory is shown in truncated form in square brackets in the console output.
Looking at what we got for the run shown above, the console log line for the sayHello process starts with `[a3/7be2fa]`. That corresponds to the following directory path: `work/a3/7be2fa7be2fad5e71e5f49998f795677fd68`

Let's take a look at what's in there.

!!! tip

    If you browse the contents of the task subdirectory in the VSCode file explorer, you'll see all the files right away.
    However, the log files are set to be invisible in the terminal, so if you want to use `ls` or `tree` to view them, you'll need to set the relevant option for displaying invisible files.

    ```bash
    tree -a work
    ```

You should see something like this, though the exact subdirectory names will be different on your system:

```console title="Directory contents"
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

These are the helper and log files:

- **`.command.begin`**: Metadata related to the beginning of the execution of the process call
- **`.command.err`**: Error messages (`stderr`) emitted by the process call
- **`.command.log`**: Complete log output emitted by the process call
- **`.command.out`**: Regular output (`stdout`) by the process call
- **`.command.run`**: Full script run by Nextflow to execute the process call
- **`.command.sh`**: The command that was actually run by the process call
- **`.exitcode`**: The exit code resulting from the command

The `.command.sh` file is especially useful because it tells you what command Nextflow actually executed.
In this case it's very straightforward, but later in the course you'll see commands that involve some interpolation of variables.
When you're dealing with that, you need to be able to check exactly what was run, especially when troubleshooting an issue.

The actual output of the `sayHello` process is `output.txt`.
Open it and you will find the `Hello World!` greeting, which was the expected result of our minimalist workflow.

```console title="output.txt" linenums="1"
Hello World!
```

### Takeaway

You know how to decipher a simple Nextflow script, run it and find the output and relevant log files in the work directory.

### What's next?

Learn how to manage your workflow executions conveniently.

---

## 3. Manage workflow executions

Knowing how to launch workflows and retrieve outputs is great, but you'll quickly find there are a few other aspects of workflow management that will make your life easier, especially if you're developing your own workflows.

Here we show you how to use the `publishDir` directive to store in an output folder all the main results from your pipeline run, the `resume` feature for when you need to re-launch the same workflow, and how to delete older work directories with `nextflow clean`.

### 3.1. Publish outputs

As you have just learned, the output produced by our pipeline is buried in a working directory several layers deep.
This is done on purpose; Nextflow is in control of this directory and we are not supposed to interact with it.

However, that makes it inconvenient to retrieve outputs that we care about.

Fortunately, Nextflow provides a way to manage this more conveniently, called the `publishDir` directive, which acts at the process level.
This directive tells Nextflow to publish the output(s) of the process to a designated output directory. By default, the outputs are published as symbolic links from the `work` directory.
It allows us to retrieve the desired output file without having to dig down into the work directory.

#### 3.1.1. Add a `publishDir` directive to the `sayHello` process

In the workflow script file `hello-world.nf`, make the following code modification:

=== "After"

    ```groovy title="hello-world.nf" linenums="6" hl_lines="3"
    process sayHello {

        publishDir 'results', mode: 'copy'

        output:
            path 'output.txt'
    ```

=== "Before"

    ```groovy title="hello-world.nf" linenums="6"
    process sayHello {

        output:
            path 'output.txt'
    ```

#### 3.1.2. Run the workflow again

Now run the modified workflow script:

```bash
nextflow run hello-world.nf
```

The log output should look very familiar:

```console title="Output" linenums="1"
 N E X T F L O W   ~  version 25.04.3

Launching `hello-world.nf` [jovial_mayer] DSL2 - revision: 35bd3425e5

executor >  local (1)
[62/49a1f8] sayHello | 1 of 1 ✔
```

This time, Nextflow has created a new directory called `results/`.
Our `output.txt` file is in this directory.
If you check the contents it should match the output in the work subdirectory.
This is how we publish results files outside of the working directories conveniently.

When you're dealing with very large files that you don't need to retain for long, you may prefer to set the `publishDir` directive to make a symbolic link to the file instead of copying it.
However, if you delete the work directory as part of a cleanup operation, you will lose access to the file, so always make sure you have actual copies of everything you care about before deleting anything.

!!! note

    A newer syntax option documented [here](https://www.nextflow.io/docs/latest/workflow.html#publishing-outputs) has been proposed to make it possible to declare and publish workflow-level outputs.
    This will eventually make using `publishDir` at the process level redundant for completed pipelines.
    However, we expect that `publishDir` will still remain very useful during pipeline development.

### 3.2. Re-launch a workflow with `-resume`

Sometimes, you're going to want to re-run a pipeline that you've already launched previously without redoing any steps that already completed successfully.

Nextflow has an option called `-resume` that allows you to do this.
Specifically, in this mode, any processes that have already been run with the exact same code, settings and inputs will be skipped.
This means Nextflow will only run processes that you've added or modified since the last run, or to which you're providing new settings or inputs.

There are two key advantages to doing this:

- If you're in the middle of developing your pipeline, you can iterate more rapidly since you only have to run the process(es) you're actively working on in order to test your changes.
- If you're running a pipeline in production and something goes wrong, in many cases you can fix the issue and relaunch the pipeline, and it will resume running from the point of failure, which can save you a lot of time and compute.

To use it, simply add `-resume` to your command and run it:

```bash
nextflow run hello-world.nf -resume
```

The console output should look similar.

```console title="Output" linenums="1"
 N E X T F L O W   ~  version 25.04.3

Launching `hello-world.nf` [golden_cantor] DSL2 - revision: 35bd3425e5

[62/49a1f8] sayHello | 1 of 1, cached: 1 ✔
```

Look for the `cached:` bit that has been added in the process status line (line 5), which means that Nextflow has recognized that it has already done this work and simply re-used the result from the previous successful run.

You can also see that the work subdirectory hash is the same as in the previous run.
Nextflow is literally pointing you to the previous execution and saying "I already did that over there."

!!! note

    When your re-run a pipeline with `resume`, Nextflow does not overwrite any files written to a `publishDir` directory by any process call that was previously run successfully.

### 3.3. Delete older work directories

During the development process, you'll typically run your draft pipelines a large number of times, which can lead to an accumulation of very many files across many subdirectories.
Since the subdirectories are named randomly, it is difficult to tell from their names what are older vs. more recent runs.

Nextflow includes a convenient `clean` subcommand that can automatically delete the work subdirectories for past runs that you no longer care about, with several [options](https://www.nextflow.io/docs/latest/reference/cli.html#clean) to control what will be deleted.

Here we show you an example that deletes all subdirectories from runs before a given run, specified using its run name.
The run name is the machine-generated two-part string shown in square brackets in the `Launching (...)` console output line.

First we use the dry run flag `-n` to check what will be deleted given the command:

```bash
nextflow clean -before golden_cantor -n
```

The output should look like this:

```console title="Output"
Would remove /workspaces/training/hello-nextflow/work/a3/7be2fad5e71e5f49998f795677fd68
```

If you don't see any lines output, you either did not provide a valid run name or there are no past runs to delete.

If the output looks as expected and you want to proceed with the deletion, re-run the command with the `-f` flag instead of `-n`:

```bash
nextflow clean -before golden_cantor -f
```

You should now see the following:

```console title="Output"
Removed /workspaces/training/hello-nextflow/work/a3/7be2fad5e71e5f49998f795677fd68
```

!!! Warning

    Deleting work subdirectories from past runs removes them from Nextflow's cache and deletes any outputs that were stored in those directories.
    That means it breaks Nextflow's ability to resume execution without re-running the corresponding processes.

    You are responsible for saving any outputs that you care about or plan to rely on! If you're using the `publishDir` directive for that purpose, make sure to use the `copy` mode, not the `symlink` mode.

### Takeaway

You know how to publish outputs to a specific directory, relaunch a pipeline without repeating steps that were already run in an identical way, and use the `nextflow clean` command to clean up old work directories.

### What's next?

Learn to provide a variable input via a command-line parameter and utilize default values effectively.

---

## 4. Use a variable input passed on the command line

In its current state, our workflow uses a greeting hardcoded into the process command.
We want to add some flexibility by using an input variable, so that we can more easily change the greeting at runtime.

### 4.1. Modify the workflow to take and use a variable input

This requires us to make three changes to our script:

1. Tell the process to expect a variable input by adding an `input:` block
2. Edit the process to use the input
3. Set up a command-line parameter and provide its value as an input to the process call

Let's make these changes one at a time.

#### 4.1.1. Add an input block to the process definition

First we need to adapt the process definition to accept an input called `greeting`.

In the process block, make the following code change:

=== "After"

    ```groovy title="hello-world.nf" linenums="6" hl_lines="5 6"
    process sayHello {

        publishDir 'results', mode: 'copy'

        input:
            val greeting

        output:
            path 'output.txt'
    ```

=== "Before"

    ```groovy title="hello-world.nf" linenums="6"
    process sayHello {

        publishDir 'results', mode: 'copy'

        output:
            path 'output.txt'
    ```

The `greeting` variable is prefixed by `val` to tell Nextflow it's a value (not a path).

#### 4.1.2. Edit the process command to use the input variable

Now we swap the original hardcoded value for the value of the input variable we expect to receive.

In the process block, make the following code change:

=== "After"

    ```groovy title="hello-channels.nf" linenums="16" hl_lines="3"
    script:
    """
    echo '$greeting' > output.txt
    """
    ```

=== "Before"

    ```groovy title="hello-channels.nf" linenums="16"
    script:
    """
    echo 'Hello World!' > output.txt
    """
    ```

Make sure to prepend the `$` symbol to tell Nextflow this is a variable name that needs to be replaced with the actual value (=interpolated).

#### 4.1.3. Set up a CLI parameter and provide it as input to the process call

Now we need to actually set up a way to provide an input value to the `sayHello()` process call.

We could simply hardcode it directly by writing `sayHello('Hello World!')`.
However, when we're doing real work with our workflow, we're often going to want to be able to control its inputs from the command line.

Good news: Nextflow has a built-in workflow parameter system called `params`, which makes it easy to declare and use CLI parameters. The general syntax is to declare `params.<parameter_name>` to tell Nextflow to expect a `--<parameter_name>` parameter on the command line.

Here, we want to create a parameter called `--greeting`, so we need to declare `params.greeting` somewhere in the workflow.
In principle we can write it anywhere; but since we're going to want to give it to the `sayHello()` process call, we can plug it in there directly by writing `sayHello(params.greeting)`.

!!! note

    The parameter name (at the workflow level) does not have to match the input variable name (at the process level).
    We're just using the same word because that's what makes sense and keeps the code readable.

In the workflow block, make the following code change:

=== "After"

    ```groovy title="hello-world.nf" linenums="24" hl_lines="2"
    // emit a greeting
    sayHello(params.greeting)
    ```

=== "Before"

    ```groovy title="hello-world.nf" linenums="24"
    // emit a greeting
    sayHello()
    ```

This tells Nextflow to run the `sayHello` process on the value provided through the `--greeting` parameter.

#### 4.1.4. Run the workflow command again

Let's run it!

```bash
nextflow run hello-world.nf --greeting 'Bonjour le monde!'
```

If you made all three edits correctly, you should get another successful execution:

```console title="Output" linenums="1"
 N E X T F L O W   ~  version 25.04.3

Launching `hello-world.nf` [elated_lavoisier] DSL2 - revision: 7c031b42ea

executor >  local (1)
[4b/654319] sayHello | 1 of 1 ✔
```

Be sure to open up the output file to check that you now have the new version of the greeting.

```console title="results/output.txt" linenums="1"
Bonjour le monde!
```

Voilà!

!!! tip

    You can readily distinguish Nextflow-level parameters from pipeline-level parameters.

    - Parameters that apply to a pipeline always take a double hyphen (`--`).
    - Parameters that modify a Nextflow setting, _e.g._ the `-resume` feature we used earlier, take a single hyphen (`-`).

### 4.2. Use default values for command line parameters

In many cases, it makes sense to supply a default value for a given parameter so that you don't have to specify it for every run.

#### 4.2.1. Set a default value for the CLI parameter

Let's give the `greeting` parameter with a default value by declaring it before the workflow definition.

```groovy title="hello-world.nf" linenums="22"
/*
 * Pipeline parameters
 */
params.greeting = 'Holà mundo!'
```

!!! tip

    You can put the parameter declaration inside the workflow block if you prefer. Whatever you choose, try to group similar things in the same place so you don't end up with declarations all over the place.

#### 4.2.2. Run the workflow again without specifying the parameter

Now that you have a default value set, you can run the workflow again without having to specify a value in the command line.

```bash
nextflow run hello-world.nf
```

The console output should look the same.

```console title="Output" linenums="1"
 N E X T F L O W   ~  version 25.04.3

Launching `hello-world.nf` [determined_edison] DSL2 - revision: 3539118582

executor >  local (1)
[72/394147] sayHello | 1 of 1 ✔
```

Check the output in the results directory:

```console title="results/output.txt" linenums="1"
Holà mundo!
```

Nextflow used the default value of the greeting parameter to create the output.

#### 4.2.3. Run the workflow again with the parameter to override the default value

If you provide the parameter on the command line, the CLI value will override the default value.

Try it out:

```bash
nextflow run hello-world.nf --greeting 'Konnichiwa!'
```

The console output should look the same.

```console title="Output" linenums="1"
 N E X T F L O W   ~  version 25.04.3

Launching `hello-world.nf` [elegant_faraday] DSL2 - revision: 3539118582

executor >  local (1)
[6f/a12a91] sayHello | 1 of 1 ✔
```

Now you will have the corresponding new output in your results directory.

```console title="results/output.txt" linenums="1"
Konnichiwa!
```

!!! note

    In Nextflow, there are multiple places where you can specify values for parameters.
    If the same parameter is set to different values in multiple places, Nexflow will determine what value to use based on the order of precedence that is described [here](https://www.nextflow.io/docs/latest/config.html).

### Takeaway

You know how to use a simple variable input provided at runtime via a command-line parameter, as well as set up, use and override default values.

More generally, you know how to interpret a simple Nextflow workflow, manage its execution, and retrieve outputs.

### What's next?

Take a little break, you've earned it!
When you're ready, move on to Part 2 to learn how to use channels to feed inputs into your workflow, which will allow you to take advantage of Nextflow's built-in dataflow parallelism and other powerful features.
