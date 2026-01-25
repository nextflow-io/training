# Part 1: Run basic operations

In this first part of the Nextflow Run training course, we ease into the topic with a very basic domain-agnostic Hello World example, which we'll use to demonstrate essential operations and point out the corresponding Nextflow code components.

??? info "What is a Hello World example?"

    A "Hello World!" is a minimalist example that is meant to demonstrate the basic syntax and structure of a programming language or software framework.
    The example typically consists of printing the phrase "Hello, World!" to the output device, such as the console or terminal, or writing it to a file.

---

## 1. Run a Hello World directly

Let's demonstrate this concept with a simple command that we run directly in the terminal, to show what it does before we wrap it in Nextflow.

!!! tip

    Remember that you should now be inside the `nextflow-run/` directory as described on the [Getting Started](00_orientation.md) page.

### 1.1. Make the terminal say hello

Run the following command in your terminal.

```bash
echo 'Hello World!'
```

??? success "Command output"

    ```console
    Hello World!
    ```

This outputs the text 'Hello World' right there in the terminal.

### 1.2. Write the output to a file

Running pipelines mostly involves reading data from files and writing results to other files, so let's modify the command to write the text output to a file to make the example a bit more relevant.

```bash
echo 'Hello World!' > output.txt
```

??? success "Command output"

    ```console

    ```

This does not output anything to the terminal.

### 1.3. Find the output

The text 'Hello World' should now be in the output file we specified, named `output.txt`.
You can open it in the file explorer or from the command line using the `cat` utility, for example.

??? abstract "File contents"

    ```console title="output.txt" linenums="1"
    Hello World!
    ```

This is what we're going to try to replicate with our very first Nextflow workflow.

### Takeaway

You now know how to run a simple command in the terminal that outputs some text, and optionally, how to make it write the output to a file.

### What's next?

Find out what it takes to run a Nextflow workflow that achieves the same result.

---

## 2. Run the workflow

We provide you with a workflow script named `1-hello.nf` that takes an input greeting via a command-line argument named `--input` and produces a text file containing that greeting.

We're not going to look at the code yet; first let's see what it looks like to run it.

### 2.1. Launch the workflow and monitor execution

In the terminal, run the following command:

```bash
nextflow run 1-hello.nf --input 'Hello World!'
```

??? success "Command output"

    ```console hl_lines="7"
    N E X T F L O W   ~  version 25.10.2

    Launching `1-hello.nf` [goofy_torvalds] DSL2 - revision: c33d41f479

    executor >  local (1)
    [a3/7be2fa] sayHello | 1 of 1 ✔
    ```

If your console output looks something like that, then congratulations, you just ran your first Nextflow workflow!

The most important output here is the last line, which is highlighted in the output above:

```console
[a3/7be2fa] sayHello | 1 of 1 ✔
```

This tells us that the `sayHello` process was successfully executed once (`1 of 1 ✔`).

That's great, but you may be wondering: where is the output?

### 2.2. Find the output file in the `results` directory

This workflow is configured to publish its output to a results directory.
If you look at your current directory, you will see that when you ran the workflow, Nextflow created a new directory called `results`, as well as a subdirectory called `1-hello` under that, containing a file called `output.txt`.

```console title="results/"
results
└── 1-hello
    └── output.txt
```

Open the file; the contents should match the string you specified on the command line.

```console title="results/1-hello/output.txt" linenums="1"
Hello World!
```

That's great, our workflow did what it was supposed to do!

However, be aware that the 'published' result is a copy (or in some cases a symbolic link) of the actual output produced by Nextflow when it executed the workflow.

So now, we are going to peek under the hood to see where Nextflow actually executed the work.

!!! Warning

    Not all workflows will be set up to publish outputs to a results directory, and/or the directory names and structure may be different.
    A little further in this section, we will show you how to find out where this behavior is specified.

### 2.3. Find the original output and logs in the `work/` directory

When you run a workflow, Nextflow creates a distinct 'task directory' for every single invocation of each process in the workflow (=every step in the pipeline).
For each one, it will stage the necessary inputs, execute the relevant instruction(s) and write outputs and log files within that one directory, which is named automatically using a hash in order to make it unique.

All of these task directories will live under a directory called `work` within your current directory (where you're running the command).

That may sound confusing, so let's see what that looks like in practice.

Going back to the console output for the workflow we ran earlier, we had this line:

```console
[a3/7be2fa] sayHello | 1 of 1 ✔
```

See how the line starts with `[a3/7be2fa]`?
That is a truncated form of the task directory path for that one process call, and tells you where to find the output of the `sayHello` process call within the `work/` directory path.

You can find the full path by typing the following command (replacing `a3/7be2fa` with what you see in your own terminal) and pressing the tab key to autocomplete the path or adding an asterisk:

```bash
ls work/a3/7be2fa*
```

This should yield the full path directory path: `work/a3/7be2fa7be2fad5e71e5f49998f795677fd68`

Let's take a look at what's in there.

??? abstract "Directory contents"

    ```console
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

??? question "Don't see the same thing?"

    The exact subdirectory names will be different on your system.

    If you browse the contents of the task subdirectory in the VSCode file explorer, you'll see all the files right away.
    However, the log files are set to be invisible in the terminal, so if you want to use `ls` or `tree` to view them, you'll need to set the relevant option for displaying invisible files.

    ```bash
    tree -a work
    ```

You should immediately recognize the `output.txt` file, which is in fact the original output of the `sayHello` process that got published to the `results` directory.
If you open it, you will find the `Hello World!` greeting again.

```console title="work/a3/7be2fa7be2fad5e71e5f49998f795677fd68/output.txt"
Hello World!
```

So what about all those other files?

These are the helper and log files that Nextflow wrote as part of the task execution:

- **`.command.begin`**: Sentinel file created as soon as the task is launched.
- **`.command.err`**: Error messages (`stderr`) emitted by the process call
- **`.command.log`**: Complete log output emitted by the process call
- **`.command.out`**: Regular output (`stdout`) by the process call
- **`.command.run`**: Full script run by Nextflow to execute the process call
- **`.command.sh`**: The command that was actually run by the process call
- **`.exitcode`**: The exit code resulting from the command

The `.command.sh` file is especially useful because it shows you the main command Nextflow executed, not including all the bookkeeping and task/environment setup.

```console title="work/a3/7be2fa7be2fad5e71e5f49998f795677fd68/command.sh"
#!/bin/bash -ue
echo 'Hello World!' > output.txt

```

So this confirms that the workflow composed the same command we ran directly on the command-line earlier.

When something goes wrong and you need to troubleshoot what happened, it can be useful to look at the `command.sh` script to check exactly what command Nextflow composed based on the workflow instructions, variable interpolation and so on.

### 2.4. Re-run the workflow with different greetings

Try re-running the workflow a few times with different values for the `--input` argument, then look at the task directories.

??? abstract "Directory contents"

    ```console
    work
    ├── 0f
    │   └── 52b7e07b0e274a80843fca48ed21b8
    │       ├── .command.begin
    │       ├── .command.err
    │       ├── .command.log
    │       ├── .command.out
    │       ├── .command.run
    │       ├── .command.sh
    │       ├── .exitcode
    │       └── output.txt
    ├── 67
    │   ├── 134e6317f90726c6c17ad53234a32b
    │   │   ├── .command.begin
    │   │   ├── .command.err
    │   │   ├── .command.log
    │   │   ├── .command.out
    │   │   ├── .command.run
    │   │   ├── .command.sh
    │   │   ├── .exitcode
    │   │   └── output.txt
    │   └── e029f2e75305874a9ab263d21ebc2c
    │       ├── .command.begin
    │       ├── .command.err
    │       ├── .command.log
    │       ├── .command.out
    │       ├── .command.run
    │       ├── .command.sh
    │       ├── .exitcode
    │       └── output.txt
    ├── 6c
    │   └── d4fd787e0b01b3c82e85696c297500
    │       ├── .command.begin
    │       ├── .command.err
    │       ├── .command.log
    │       ├── .command.out
    │       ├── .command.run
    │       ├── .command.sh
    │       ├── .exitcode
    │       └── output.txt
    └── e8
        └── ab99fad46ade52905ec973ff39bb80
            ├── .command.begin
            ├── .command.err
            ├── .command.log
            ├── .command.out
            ├── .command.run
            ├── .command.sh
            ├── .exitcode
            └── output.txt
    ```

You see that a new subdirectory with a complete set of output and log files has been created for each run.

In contrast, if you look at the `results` directory, there is still only one set of results, and the content of the output file corresponds to whatever you ran last.

??? abstract "Directory contents"

    ```console title="results/"
    results
    └── 1-hello
        └── output.txt
    ```

This shows you that the published results will get overwritten by subsequent executions, whereas the task directories under `work/` are preserved.

### Takeaway

You know how to run a simple Nextflow script, monitor its execution and find its outputs.

### What's next?

Learn how to read a basic Nextflow script and identify how its components relate to its functionality.

---

## 3. Examine the Hello World workflow starter script

What we did there was basically treating the workflow script like a black box.
Now that we've seen what it does, let's open the box and look inside.

Our goal here is not to memorize the syntax of Nextflow code, but to form some basic intuition of what are the main components and how they are organized.

### 3.1. Examine the overall code structure

You'll find the `1-hello.nf` script in your current directory, which should be `nextflow-run`. Open it in the editor pane.

??? full-code "Full code file"

    ```groovy title="1-hello.nf" linenums="1"
    #!/usr/bin/env nextflow

    /*
    * Use echo to print 'Hello World!' to a file
    */
    process sayHello {

        input:
        val greeting

        output:
        path 'output.txt'

        script:
        """
        echo '${greeting}' > output.txt
        """
    }

    /*
    * Pipeline parameters
    */
    params {
        input: String
    }

    workflow {

        main:
        // emit a greeting
        sayHello(params.input)

        publish:
        first_output = sayHello.out
    }

    output {
        first_output {
            path '1-hello'
            mode 'copy'
        }
    }
    ```

A Nextflow workflow script typically includes one or more **process** definitions, the **workflow** itself, and a few optional blocks such as **params** and **output**.

Each **process** describes what operation(s) the corresponding step in the pipeline should accomplish, while the **workflow** describes the dataflow logic that connects the various steps.

Let's take a closer look at the **process** block first, then we'll look at the **workflow** block.

### 3.2. The `process` definition

The first block of code describes a **process**.
The process definition starts with the keyword `process`, followed by the process name and finally the process body delimited by curly braces.
The process body must contain a script block which specifies the command to run, which can be anything you would be able to run in a command line terminal.

```groovy title="1-hello.nf" linenums="3"
/*
* Use echo to print a greeting to a file
*/
process sayHello {

    input:
    val greeting

    output:
    path 'output.txt'

    script:
    """
    echo '${greeting}' > output.txt
    """
}
```

Here we have a **process** called `sayHello` that takes an **input** variable called `greeting` and writes its **output** to a file named `output.txt`.

<figure class="excalidraw">
--8<-- "docs/nextflow_run/img/sayhello_with_input.svg"
</figure>

This is a very minimal process definition that just contains an `input` definition, an `output` definition and the `script` to execute.

The `input` definition includes the `val` qualifier, which tells Nextflow to expect a value of some kind (can be a string, a number, whatever).

The `output` definition includes the `path` qualifier, which tells Nextflow this should be handled as a path (includes both directory paths and files).

### 3.3. The `workflow` definition

The second block of code describes the **workflow** itself.
The workflow definition starts with the keyword `workflow`, followed by an optional name, then the workflow body delimited by curly braces.

Here we have a **workflow** that consists of a `main:` block and a `publish:` block.
The `main:` block is the main body of the workflow and the `publish:` block lists the outputs that should be published to the `results` directory.

```groovy title="1-hello.nf" linenums="27"
workflow {

    main:
    // emit a greeting
    sayHello(params.input)

    publish:
    first_output = sayHello.out
}
```

In this case the `main:` block contains a call to the `sayHello` process and gives it an input called `params.input` to use as the greeting.

As we'll discuss in more detail in a moment, `params.input` holds the value we gave to the `--input` parameter in our command line.

The `publish:` block lists the output of the `sayHello()` process call, which it refers to as `sayHello.out` and gives the name `first_output` (this can be anything the workflow author wants).

This is a very minimal **workflow** definition.
In a real-world pipeline, the workflow typically contains multiple calls to **processes** connected by **channels**, and there may be default values set up for the variable inputs.

We'll get into that in Part 2 of the course.
For now, let's take a closer look at how our workflow is handling inputs and outputs.

### 3.4. The `params` system of command-line parameters

The `params.input` we provide to the `sayHello()` process call is a neat bit of Nextflow code and is worth spending an extra minute on.

As mentioned above, that's how we pass the value of the `--input` command-line parameter to the `sayHello()` process call.
In fact, simply declaring `params.someParameterName` is enough to give the workflow a parameter named `--someParameterName` from the command-line.

Here we've formalized that parameter declaration by setting up a `params` block that specifies the type of input the workflow expects (Nextflow 25.10.2 and later).

```groovy title="1-hello.nf" linenums="20"
/*
 * Pipeline parameters
 */
params {
    input: String
}
```

Supported types include `String`, `Integer`, `Float`, `Boolean`, and `Path`.

!!! tip

    Workflow parameters declared using the `params` system always take two dashes on the command line (`--`).
    This distinguishes them from Nextflow-level parameters, which only take one dash (`-`).

### 3.5. The `publish` directive

On the other end of the workflow, we've already glanced at the `publish:` block.
That's one half of the output handling system; the other half is the `output` block located below.

```groovy title="1-hello.nf" linenums="37"
output {
    first_output {
        path '1-hello'
        mode 'copy'
    }
}
```

This specifies that the `first_output` output listed in the `publish:` block should be copied to a subdirectory called `1-hello` under the default `results` output directory.

The `mode 'copy'` line overrides the system's default behavior, which is to make a symbolic link (or symlink) to the original file in the `work/` directory instead of a proper copy.

There are more options than displayed here for controlling the publishing behavior; we'll cover a few later on.
You'll also see that when a workflow generates multiple outputs, each one gets listed this way in the `output` block.

??? info "Older syntax for publishing outputs using `publishDir`"

    Until very recently, the established way to publish outputs was to do it at the level of each individual process using a `publishDir` directive.

    You will still find this code pattern all over the place in older Nextflow pipelines and process modules, so it's important to be aware of it.

    Instead of having a `publish:` block in the workflow and an `output` block at the top level, you would see a `publishDir` line in the sayHello` process definition:

    ```groovy title="Syntax example" linenums="1" hl_lines="3"
    process sayHello {

        publishDir 'results/1-hello', mode: 'copy'

        output:
        path 'output.txt'

        script:
        """
        echo 'Hello World!' > output.txt
        """
    }
    ```

    However, we do not recommend using this in any new work as it will eventually be disallowed in future versions of the Nextflow language.

### Takeaway

You now know how a simple Nextflow workflow is structured, and how the basic components relate to its functionality.

### What's next?

Learn to manage your workflow executions conveniently.

---

## 4. Manage workflow executions

Knowing how to launch workflows and retrieve outputs is great, but you'll quickly find there are a few other aspects of workflow management that will make your life easier.

Here we show you how to take advantage of the `resume` feature for when you need to re-launch the same workflow, how to inspect the execution logs with `nextflow log`, and how to delete older work directories with `nextflow clean`.

### 4.1. Re-launch a workflow with `-resume`

Sometimes, you're going to want to re-run a pipeline that you've already launched previously without redoing any work that was already completed successfully.

Nextflow has an option called `-resume` that allows you to do this.
Specifically, in this mode, any processes that have already been run with the exact same code, settings and inputs will be skipped.
This means Nextflow will only run processes that you've added or modified since the last run, or to which you're providing new settings or inputs.

There are two key advantages to doing this:

- If you're in the middle of developing a pipeline, you can iterate more rapidly since you only have to run the process(es) you're actively working on in order to test your changes.
- If you're running a pipeline in production and something goes wrong, in many cases you can fix the issue and relaunch the pipeline, and it will resume running from the point of failure, which can save you a lot of time and compute.

To use it, simply add `-resume` to your command and run it:

```bash
nextflow run 1-hello.nf --greeting 'Hello World!' -resume
```

??? success "Command output"

    ```console linenums="1"
    N E X T F L O W   ~  version 25.04.3

    Launching `1-hello.nf` [tiny_noyce] DSL2 - revision: c33d41f479

    [a3/7be2fa] process > sayHello [100%] 1 of 1, cached: 1 ✔
    ```

The console output should look familiar, but there's one thing that's a little different compared to before.

Look for the `cached:` bit that has been added in the process status line (line 5), which means that Nextflow has recognized that it has already done this work and simply reused the result from the previous successful run.

You can also see that the work subdirectory hash is the same as in the previous run.
Nextflow is literally pointing you to the previous execution and saying "I already did that over there."

!!! tip

    When your re-run a pipeline with `resume`, Nextflow does not overwrite any files published outside of the work directory by any executions that were run successfully previously.

### 4.2. Inspect the log of past executions

Whenever you launch a nextflow workflow, a line gets written to a log file called `history`, under a hidden directory called `.nextflow` in the current working directory.

??? abstract "File contents"

    ```txt title=".nextflow/history" linenums="1"
    2025-07-04 19:27:09	1.8s	wise_watson	OK	3539118582ccde68dde471cc2c66295c	a02c9c46-c3c7-4085-9139-d1b9b5b194c8	nextflow run 1-hello.nf --input 'Hello World'
    2025-07-04 19:27:20	2.9s	spontaneous_blackwell	OK	3539118582ccde68dde471cc2c66295c	59a5db23-d83c-4c02-a54e-37ddb73a337e	nextflow run 1-hello.nf --input Bonjour
    2025-07-04 19:27:31	1.8s	gigantic_yonath	OK	3539118582ccde68dde471cc2c66295c	5acaa83a-6ad6-4509-bebc-cb25d5d7ddd0	nextflow run 1-hello.nf --input 'Dobry den'
    2025-07-04 19:27:45	2.4s	backstabbing_swartz	OK	3539118582ccde68dde471cc2c66295c	5f4b3269-5b53-404a-956c-cac915fbb74e	nextflow run 1-hello.nf --input Konnichiwa
    2025-07-04 19:27:57	2.1s	goofy_wilson	OK	3539118582ccde68dde471cc2c66295c	5f4b3269-5b53-404a-956c-cac915fbb74e	nextflow run 1-hello.nf --input Konnichiwa -resume
    ```

This file gives you the timestamp, run name, status, revision ID, session ID and full command line for every Nextflow run that has been launched from within the current working directory.

A more convenient way to access this information is to use the `nextflow log` command.

```bash
nextflow log
```

??? success "Command output"

    ```console linenums="1"
    TIMESTAMP               DURATION        RUN NAME                STATUS  REVISION ID     SESSION ID                              COMMAND
    2025-07-04 19:27:09     1.8s            wise_watson             OK       3539118582     a02c9c46-c3c7-4085-9139-d1b9b5b194c8    nextflow run 1-hello.nf --input 'Hello World'
    2025-07-04 19:27:20     2.9s            spontaneous_blackwell   OK       3539118582     59a5db23-d83c-4c02-a54e-37ddb73a337e    nextflow run 1-hello.nf --input Bonjour
    2025-07-04 19:27:31     1.8s            gigantic_yonath         OK       3539118582     5acaa83a-6ad6-4509-bebc-cb25d5d7ddd0    nextflow run 1-hello.nf --input 'Dobry den'
    2025-07-04 19:27:45     2.4s            backstabbing_swartz     OK       3539118582     5f4b3269-5b53-404a-956c-cac915fbb74e    nextflow run 1-hello.nf --input Konnichiwa
    2025-07-04 19:27:57     2.1s            goofy_wilson            OK       3539118582     5f4b3269-5b53-404a-956c-cac915fbb74e    nextflow run 1-hello.nf --input Konnichiwa -resume
    ```

This will output the contents of the log file to the terminal, augmented with a header line.

You'll notice that the session ID changes whenever you run a new `nextflow run` command, EXCEPT if you're using the `-resume` option.
In that case, the session ID stays the same.

Nextflow uses the session ID to group run caching information under the `cache` directory, also located under `.nextflow`.

### 4.3. Delete older work directories

If you run a lot of pipelines, you may end up accumulating very many files across many subdirectories.
Since the subdirectories are named randomly, it is difficult to tell from their names what are older vs. more recent runs.

Fortunately Nextflow includes a helpful `clean` subcommand that can automatically delete the work subdirectories for past runs that you no longer care about.

#### 4.3.1. Determine deletion criteria

There are multiple [options](https://www.nextflow.io/docs/latest/reference/cli.html#clean) to determine what to delete.

Here we show you an example that deletes all subdirectories from runs before a given run, specified using its run name.

Look up the most recent successful run where you didn't use `-resume`; in our case the run name was `backstabbing_swartz`.

The run name is the machine-generated two-part string shown in square brackets in the `Launching (...)` console output line.
You can also use the Nextflow log to look up a run based on its timestamp and/or command line.

#### 4.3.2. Do a dry run

First we use the dry run flag `-n` to check what will be deleted given the command:

```bash
nextflow clean -before backstabbing_swartz -n
```

??? success "Command output"

    ```console
    Would remove /workspaces/training/hello-nextflow/work/eb/1a5de36637b475afd88fca7f79e024
    Would remove /workspaces/training/hello-nextflow/work/6b/19b0e002ea13486d3a0344c336c1d0
    Would remove /workspaces/training/hello-nextflow/work/45/9a6dd7ab771f93003d040956282883
    ```

Your output will have different task directory names and may have a different number of lines, but it should look similar to the example.

If you don't see any lines output, you either did not provide a valid run name or there are no past runs to delete. Make sure to change `backstabbing_swartz` in the example command to whatever is the corresponding latest run name in your log.

#### 4.3.3. Proceed with deletion

If the output looks as expected and you want to proceed with the deletion, re-run the command with the `-f` flag instead of `-n`:

```bash
nextflow clean -before backstabbing_swartz -f
```

??? success "Command output"

    ```console
    Removed /workspaces/training/hello-nextflow/work/eb/1a5de36637b475afd88fca7f79e024
    Removed /workspaces/training/hello-nextflow/work/6b/19b0e002ea13486d3a0344c336c1d0
    Removed /workspaces/training/hello-nextflow/work/45/9a6dd7ab771f93003d040956282883
    ```

The output should be similar to before, but now saying 'Removed' instead of 'Would remove'.
Note that this does not remove the two-character subdirectories (like `eb/` above) but it does empty their contents.

!!! Warning

    Deleting work subdirectories from past runs removes them from Nextflow's cache and deletes any outputs that were stored in those directories.
    That means it breaks Nextflow's ability to resume execution without re-running the corresponding processes.

    You are responsible for saving any outputs that you care about! That is the main reason we prefer to use the `copy` mode rather than the `symlink` mode for the `publish` directive.

### Takeaway

You know how to relaunch a pipeline without repeating steps that were already run in an identical way, inspect the execution log, and use the `nextflow clean` command to clean up old work directories.

### What's next?

Take a little break! You've just absorbed the building blocks of Nextflow syntax and basic usage instructions.

In the next section of this training, we're going to look at four successively more realistic versions of the Hello World pipeline that will demonstrate how Nextflow allows you to process multiple inputs efficiently, run workflows composed of multiple steps connected together, leverage modular code components, and utilize containers for greater reproducibility and portability.
