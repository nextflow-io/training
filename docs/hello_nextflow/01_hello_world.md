# Part 1: Hello World

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/8X2hHI-9vms?si=F0t9LFYLjAWoyRXj&amp;list=PLPZ8WHdZGxmXiHf8B26oB_fTfoKQdhlik" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

/// caption
:fontawesome-brands-youtube:{ .youtube } See [the whole playlist](https://www.youtube.com/playlist?list=PLPZ8WHdZGxmXiHf8B26oB_fTfoKQdhlik) on the Nextflow YouTube channel.

:green_book: The video transcript is available [here](./transcripts/01_hello_world.md).
///

In this first part of the Hello Nextflow training course, we ease into the topic with a very basic domain-agnostic Hello World example, which we'll progressively build up to demonstrate the usage of foundational Nextflow logic and components.

??? info "What is a Hello World example?"

    A "Hello World!" is a minimalist example that is meant to demonstrate the basic syntax and structure of a programming language or software framework. The example typically consists of printing the phrase "Hello, World!" to the output device, such as the console or terminal, or writing it to a file.

---

## 0. Warmup: Run Hello World directly

Let's demonstrate this with a simple command that we run directly in the terminal, to show what it does before we wrap it in Nextflow.

!!! tip

    Remember that you should now be inside the `hello-nextflow/` directory as described on the [Getting Started](00_orientation.md) page.

### 0.1. Make the terminal say hello

Run the following command in your terminal.

```bash
echo 'Hello World!'
```

??? success "Command output"

    ```console
    Hello World!
    ```

This outputs the text 'Hello World' right there in the terminal.

### 0.2. Write the output to a file

Running pipelines mostly involves reading data from files and writing results to other files, so let's modify the command to write the text output to a file to make the example a bit more relevant.

```bash
echo 'Hello World!' > output.txt
```

??? success "Command output"

    ```console

    ```

This does not output anything to the terminal.

### 0.3. Find the output

The text 'Hello World' should now be in the output file we specified, named `output.txt`.
You can open it in the file explorer or from the command line using the `cat` utility, for example.

??? abstract "File contents"

    ```console title="output.txt" linenums="1"
    Hello World!
    ```

### Takeaway

You now know how to run a simple command in the terminal that outputs some text, and optionally, how to make it write the output to a file.

### What's next?

Find out what that would look like written as a Nextflow workflow.

---

## 1. Examine the script and run it

We provide you with a fully functional if minimalist workflow script named `hello-world.nf` that does the same thing as before (write out 'Hello World!') but with Nextflow.

To get you started, let's open up the workflow script so you can get a sense of how it's structured.
Then we'll run it and look for its outputs.

### 1.1. Examine the code

You'll find the `hello-world.nf` script in your current director, which should be `hello-nextflow`. Open it in the editor pane.

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

    main:
    // emit a greeting
    sayHello()
}
```

As you can see, a Nextflow script involves two main types of core components: one or more **processes**, and the **workflow** itself.
Each **process** describes what operation(s) the corresponding step in the pipeline should accomplish, while the **workflow** describes the dataflow logic that connects the various steps.

Let's take a closer look at the **process** block first, then we'll look at the **workflow** block.

#### 1.1.1. The `process` definition

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

Importantly, the output definition does not _determine_ what output will be created.
It simply _declares_ what is the expected output, so that Nextflow can look for it once execution is complete.
This is necessary for verifying that the command was executed successfully and for passing the output to downstream processes if needed. Output produced that doesn't match what is declared in the output block will not be passed to downstream processes.

!!! warning

    This example is brittle because we hardcoded the output filename in two separate places (the script and the output blocks).
    If we change one but not the other, the script will break.
    Later on, you'll learn ways to use variables to mitigate this problem.

In a real-world pipeline, a process usually contains additional blocks such as directives and inputs, which we'll introduce in a little bit.

#### 1.1.2. The `workflow` definition

The second block of code describes the **workflow** itself.
The workflow definition starts with the keyword `workflow`, followed by an optional name, then the workflow body delimited by curly braces.

Here we have a **workflow** that consists of a `main:` block (which says 'this is the main body of the workflow') containing a call to the `sayHello` process.

```groovy title="hello-world.nf" linenums="17"
workflow {

    main:
    // emit a greeting
    sayHello()
}
```

This is a very minimal **workflow** definition.
In a real-world pipeline, the workflow typically contains multiple calls to **processes** connected by **channels**, and the processes expect one or more variable **input(s)**.

You'll learn how to add variable inputs later in this training module; and you'll learn how to add more processes and connect them by channels in Part 3 of this course.

!!! tip

    Technically the `main:` line is not required for simple workflows like this, so you may encounter workflows that don't have it.
    But we'll need it for taking advantage of workflow-level outputs, so we might as well include it from the start.

### 1.2. Run the workflow

Looking at code is not nearly as fun as running it, so let's try this out in practice.

#### 1.2.1. Launch the workflow and monitor execution

In the terminal, run the following command:

```bash
nextflow run hello-world.nf
```

??? success "Command output" hl_lines="6"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `hello-world.nf` [goofy_torvalds] DSL2 - revision: c33d41f479

    executor >  local (1)
    [65/7be2fa] sayHello | 1 of 1 ✔
    ```

If your console output looks something like that, then congratulations, you just ran your first Nextflow workflow!

The most important output here is the last line, which is highlighted in the output above:

```console
[65/7be2fa] sayHello | 1 of 1 ✔
```

This tells us that the `sayHello` process was successfully executed once (`1 of 1 ✔`).

Importantly, this line also tells you where to find the output of the `sayHello` process call.
Let's look at that now.

#### 1.2.2. Find the output and logs in the `work` directory

When you run Nextflow for the first time in a given directory, it creates a directory called `work` where it will write all files (and any symlinks) generated in the course of execution.

Within the `work` directory, Nextflow organizes outputs and logs per process call.
For each process call, Nextflow creates a nested subdirectory, named with a hash in order to make it unique, where it will stage all necessary inputs (using symlinks by default), write helper files, and write out logs and any outputs of the process.

The path to that subdirectory is shown in truncated form in square brackets in the console output.
Looking at what we got for the run shown above, the console log line for the sayHello process starts with `[65/7be2fa]`. That corresponds to the following directory path: `work/65/7be2fa7be2fad5e71e5f49998f795677fd68`

Let's take a look at what's in there.

!!! tip

    If you browse the contents of the task subdirectory in the VSCode file explorer, you'll see all the files right away.
    However, the log files are set to be invisible in the terminal, so if you want to use `ls` or `tree` to view them, you'll need to set the relevant option for displaying invisible files.

    ```bash
    tree -a work
    ```

??? abstract "Directory contents"

    ```console
    work
    └── 65
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

The exact subdirectory names will be different on your system, but you should see the following helper and log files:

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

Finally the actual output of the `sayHello` process is `output.txt`.
Open it and you will find the `Hello World!` greeting, which was the expected result of our minimalist workflow.

??? abstract "File contents"

    ```console title="output.txt"
    Hello World!
    ```

All that for that!

### Takeaway

You know how to decipher a simple Nextflow script, run it and find the output and relevant log files in the work directory.

### What's next?

Learn how to 'publish' the workflow outputs to a more convenient location.

---

## 2. Publish outputs

As you have just learned, the output produced by our pipeline is buried in a working directory several layers deep.
This is done on purpose; Nextflow is in control of this directory and we are not supposed to interact with it.
However, that makes it inconvenient to retrieve outputs that we care about.

Fortunately, Nextflow provides a way to 'publish' outputs to a designated directory using [workflow-level output definitions](https://www.nextflow.io/docs/latest/workflow.html#workflow-outputs).

### 2.1. Basic usage

This is going to involve two new pieces of code:

1. A `publish:` block inside the `workflow` body, declaring process outputs.
2. An `output` block to the script specifying output options such as mode and location.

#### 2.1.1. Declare the output of the `sayHello` process

We need to add a `publish:` block to the workflow body (same kind of code element as the `main:` block) and list the output of the `sayHello()` process.

In the workflow script file `hello-world.nf`, add the following lines of code:

=== "After"

    ```groovy title="hello-world.nf" linenums="17" hl_lines="23-24"
    workflow {

        main:
        // emit a greeting
        sayHello()

        publish:
        first_output = sayHello.out
    }
    ```

=== "Before"

    ```groovy title="hello-world.nf" linenums="17"
    workflow {

        main:
        // emit a greeting
        sayHello()
    }
    ```

You see that we can refer to the output of the process simply by doing `sayHello().out`, and assign it an arbitrary name, `first_output`.

#### 2.1.2. Add an `output:` block to the script

Now we just need to add the `output:` block where the output directory path will be specified. Note that this new block sits **outside** and **below** the `workflow` block within the script.

In the workflow script file `hello-world.nf`, add the following lines of code:

=== "After"

    ```groovy title="hello-world.nf" linenums="17" hl_lines="27-31"
    workflow {

        main:
        // emit a greeting
        sayHello()

        publish:
        first_output = sayHello.out
    }

    output {
        first_output {
            path '.'
        }
    }
    ```

=== "Before"

    ```groovy title="hello-world.nf" linenums="17"
    workflow {

        main:
        // emit a greeting
        sayHello()

        publish:
        first_output = sayHello.out
    }
    ```

We can use this to assign specific paths to any process outputs declared in the `workflow` block.
Later, you'll learn about ways to generate sophisticated output directory structures, but for now, we're just hardcoding a minimal path for simplicity.

#### 2.1.3. Run the workflow

Now run the modified workflow script:

```bash
nextflow run hello-world.nf
```

??? success "Command output"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `hello-world.nf` [jovial_mayer] DSL2 - revision: 35bd3425e5

    executor >  local (1)
    [9f/48ef97] sayHello | 1 of 1 ✔
    ```

The terminal output should look familiar. Externally, nothing has changed.

However, check your file explorer: this time, Nextflow has created a new directory called `results/`.

??? abstract "Directory contents"

    ```console hl_lines="10-11 22"
    .
    ├── greetings.csv
    ├── hello-channels.nf
    ├── hello-config.nf
    ├── hello-containers.nf
    ├── hello-modules.nf
    ├── hello-workflow.nf
    ├── hello-world.nf
    ├── nextflow.config
    ├── results
    │   └── output.txt -> /workspaces/training/hello-nextflow/work/9f/48ef97f110b0dbd83635d7cbe288d2/output.txt
    ├── solutions
    │   ├── 1-hello-world
    │   ├── 2-hello-channels
    │   ├── 3-hello-workflow
    │   ├── 4-hello-modules
    │   ├── 5-hello-containers
    │   └── 6-hello-config
    ├── test-params.json
    └── work
        ├── 65
        └── 9f
    ```

Inside the `results` directory, we find a symbolic link to the `output.txt` produced in the work directory by the command we just ran.

This allows us to easily retrieve output files without having to dig through the work subdirectory.

### 2.2. Set a custom location

Having a default location is great, but you might want to customize where the results are saved and how they are organized.

For example, you may want to organize your outputs into subdirectories.
The simplest way to do that is to assign specific output path per output.

#### 2.2.1. Modify the output path

Once again, modifying the publish behavior for a specific output is really straightforward.
To set a custom location, just edit the `path` accordingly:

=== "After"

    ```groovy title="hello-world.nf" linenums="27" hl_lines="3"
    output {
        first_output {
            path 'hello_world'
        }
    }
    ```

=== "Before"

    ```groovy title="hello-world.nf" linenums="27"
    output {
        first_output {
            path '.'
        }
    }
    ```

Since this is set at the level of the individual output, you can specify different locations and subdirectories to suit your needs.

#### 2.2.2. Run the workflow again

Let's try it out.

```bash
nextflow run hello-world.nf
```

??? success "Command output"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `hello-world.nf` [tiny_shaw] DSL2 - revision: 757723adc1

    executor >  local (1)
    [8c/79499c] process > sayHello [100%] 1 of 1 ✔
    ```

This time the result gets written under the specified subdirectory.

??? abstract "Directory contents"

    ```console hl_lines="2-3"
    results/
    ├── hello_world
    │   └── output.txt -> /workspaces/training/hello-nextflow/work/8c/79499c2e506b79e2e01acb808d9d12/output.txt
    └── output.txt -> /workspaces/training/hello-nextflow/work/65/f56f2cd75df1352e106fcdd084b97b/output.txt
    ```

You can use as many levels of nesting as you'd like.

It is also possible to use the process name or other variables to name the directories used to organize results, and it is possible to change the default name of the top-level output directory (which is controlled by the special variable `outputDir`).
We will cover these options in later trainings.

### 2.3. Set the publish mode to copy

By default, the outputs are published as symbolic links from the `work` directory.
That means there is only a single file on the filesystem.

This is great when you're dealing with very large files, for which you don't want to store multiple copies.
However, if you delete the work directory at any point (we'll cover cleanup operations shortly), you will lose access to the file.
So you need to have a plan for saving copies of any important files to a secure place.

One easy option is to switch the publish mode to copy for the outputs you care about.

#### 2.3.1. Add the mode directive

This bit is really straightforward.
Just add `mode 'copy'` to the relevant workflow-level output definition:

=== "After"

    ```groovy title="hello-world.nf" linenums="27" hl_lines="4"
    output {
        first_output {
            path 'hello_world'
            mode 'copy'
        }
    }
    ```

=== "Before"

    ```groovy title="hello-world.nf" linenums="27"
    output {
        first_output {
            path 'hello_world'
        }
    }
    ```

This sets the publish mode for that specific output.

#### 2.3.2. Run the workflow again

Let's try it out.

```bash
nextflow run hello-world.nf
```

??? success "Command output"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `hello-world.nf` [tiny_shaw] DSL2 - revision: 757723adc1

    executor >  local (1)
    [df/521638] process > sayHello [100%] 1 of 1 ✔
    ```

This time, if you look at the results, the file is a proper copy instead of just a symlink.

??? abstract "Directory contents"

    ```console hl_lines="3"
    results/
    ├── hello_world
    │   └── output.txt
    └── output.txt -> /workspaces/training/hello-nextflow/work/65/f56f2cd75df1352e106fcdd084b97b/output.txt
    ```

Since this too is set at the level of the individual output, it allows you to set the publish mode in a granular way.
This is come in especially handy later when we move on to multi-step pipelines, where you might want to only copy final outputs and leave intermediate outputs as symlinks, for example.

As noted earlier, there are other, more sophisticated options for controlling how outputs are published.
We'll show you how to use them in due time in your Nextflow journey.

### 2.4. Process-level `publishDir` directive (FYI)

Until very recently, the established way to publish outputs was to do it at the level of each individual process using a `publishDir` directive.

To achieve what we just did for the outputs of the `sayHello` process, we would have instead added the following line to the process definition:

```groovy title="hello-world.nf" linenums="6" hl_lines="2"
process sayHello {

    publishDir 'results/hello_world', mode: 'copy'

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

You will still find this code pattern all over the place in older Nextflow pipelines and process modules, so it's important to be aware of it.
However, we do not recommend using it in any new work as it will eventually be disallowed in future versions of the Nextflow language.

### Takeaway

You know how to publish workflow outputs to a more convenient location.

### What's next?

Learn to provide a variable input via a command-line parameter and utilize default values effectively.

---

## 3. Use a variable input passed on the command line

In its current state, our workflow uses a greeting hardcoded into the process command.
We want to add some flexibility by using an input variable, so that we can more easily change the greeting at runtime.

This requires us to make three sets of changes to our script:

1. Change the process to expect a variable input
2. Set up a command-line parameter to capture user input
3. Pass the input to the process in the workflow body

Let's make these changes one at a time.

### 3.1. Change the `sayHello` process to expect a variable input

We need to edit the process definition to (1) accept an input variable and (2) use that variable in the command line.

#### 3.1.1. Add an input block to the process definition

First, let's adapt the process definition to accept an input called `greeting`.

In the process block, make the following code change:

=== "After"

    ```groovy title="hello-world.nf" linenums="6" hl_lines="3-4"
    process sayHello {

        input:
        val greeting

        output:
        path 'output.txt'
    ```

=== "Before"

    ```groovy title="hello-world.nf" linenums="6"
    process sayHello {

        output:
        path 'output.txt'
    ```

The `greeting` variable is prefixed by `val` to tell Nextflow it's a value (not a path).

#### 3.1.2. Edit the process command to use the input variable

Now we swap the original hardcoded value for the value of the input variable we expect to receive.

In the process block, make the following code change:

=== "After"

    ```groovy title="hello-channels.nf" linenums="14" hl_lines="3"
    script:
    """
    echo '${greeting}' > output.txt
    """
    ```

=== "Before"

    ```groovy title="hello-channels.nf" linenums="14"
    script:
    """
    echo 'Hello World!' > output.txt
    """
    ```

The `$` symbol and curly braces (`{ }`) tell Nextflow this is a variable name that needs to be replaced with the actual input value (=interpolated).

!!! tip

    The curly braces (`{ }`) were technically optional in previous versions of Nextflow, so you might see older workflows where this is written as `echo '$greeting' > output.txt`.

Now that the `sayHello()` process is ready to accept a variable input, we need a way to provide an input value to the process call at the workflow level.

### 3.2. Set up a command-line parameter to capture user input

We could simply hardcode an input directly by making the process call `sayHello('Hello World!')`.
However, when we're doing real work with our workflow, we're going to want to be able to control its inputs from the command line.

Good news: Nextflow has a built-in workflow parameter system called `params`, which makes it easy to declare and use CLI parameters.

The general syntax is to declare `params.<parameter_name>` to tell Nextflow to expect a `--<parameter_name>` parameter on the command line.

Here, we want to create a parameter called `--input`, so we need to declare `params.input` somewhere in the workflow.
In principle we can write it anywhere; but since we're going to want to give it to the `sayHello()` process call, we can plug it in there directly by writing `sayHello(params.input)`.

!!! tip

    The parameter name (at the workflow level) does not have to match the input variable name (at the process level).
    We're just using the same word because that's what makes sense and keeps the code readable.

In the workflow block, make the following code change:

=== "After"

    ```groovy title="hello-world.nf" linenums="23" hl_lines="2"
    // emit a greeting
    sayHello(params.input)
    ```

=== "Before"

    ```groovy title="hello-world.nf" linenums="23"
    // emit a greeting
    sayHello()
    ```

This tells Nextflow to run the `sayHello` process on the value provided through the `--input` parameter.

### 3.3. Run the workflow command

Let's run it!

```bash
nextflow run hello-world.nf --input 'Bonjour le monde!'
```

??? success "Command output"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `hello-world.nf` [elated_lavoisier] DSL2 - revision: 7c031b42ea

    executor >  local (1)
    [4b/654319] sayHello | 1 of 1 ✔
    ```

If you made all these edits correctly, you should get another successful execution.

Be sure to open up the output file to check that you now have the new version of the greeting.

```console title="results/output.txt" linenums="1"
Bonjour le monde!
```

Voilà!

!!! tip

    You can readily distinguish Nextflow-level parameters from pipeline-level parameters.

    - Parameters that apply to a pipeline always take a double hyphen (`--`).
    - Parameters that modify a Nextflow setting, _e.g._ the `-resume` feature we used earlier, take a single hyphen (`-`).

### 3.4. Use default values for command line parameters

In many cases, it makes sense to supply a default value for a given parameter so that you don't have to specify it for every run.

#### 3.4.1. Set a default value for the CLI parameter

Let's give the `input` parameter a default value by declaring it before the workflow definition.

```groovy title="hello-world.nf" linenums="20"
/*
 * Pipeline parameters
 */
params {
    input: String = 'Holà mundo!'
}
```

As you see, we can specify the type of input that the workflow expects (Nextflow 25.10.2 and later).
The syntax is `name: Type = default_value`.
Supported types include `String`, `Integer`, `Float`, `Boolean`, and `Path`.

!!! info

    In older workflows, you may see that whole `params` block written as just `input = 'Holà mundo!'`.

As you add more parameters to your pipeline, you should add them all to this block, whether or not you need to give them a default value.
This will make it easy to find all configurable parameters at a glance.

#### 3.4.2. Run the workflow again without specifying the parameter

Now that you have a default value set, you can run the workflow again without having to specify a value in the command line.

```bash
nextflow run hello-world.nf
```

??? success "Command output"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `hello-world.nf` [determined_edison] DSL2 - revision: 3539118582

    executor >  local (1)
    [72/394147] sayHello | 1 of 1 ✔
    ```

Check the output in the results directory:

```console title="results/hello_world/output.txt"
Holà mundo!
```

Nextflow used the default value of the greeting parameter to create the output.

#### 3.4.3. Override the default value

If you provide the parameter on the command line, the CLI value will override the default value.

Try it out:

```bash
nextflow run hello-world.nf --input 'Konnichiwa!'
```

??? success "Command output"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `hello-world.nf` [elegant_faraday] DSL2 - revision: 3539118582

    executor >  local (1)
    [6f/a12a91] sayHello | 1 of 1 ✔
    ```

Now you will have the corresponding new output in your results directory.

```console title="results/output.txt"
Konnichiwa!
```

!!! note

    In Nextflow, there are multiple places where you can specify values for parameters.
    If the same parameter is set to different values in multiple places, Nexflow will determine what value to use based on the order of precedence that is described [here](https://www.nextflow.io/docs/latest/config.html).

### Takeaway

You know how to use a simple variable input provided at runtime via a command-line parameter, as well as set up, use and override default values.

### What's next?

Learn how to manage executions more conveniently.

---

## 4. Manage workflow executions

Knowing how to launch workflows and retrieve outputs is great, but you'll quickly find there are a few other aspects of workflow management that will make your life easier, especially if you're developing your own workflows.

Here we show you how to use the `resume` feature for when you need to re-launch the same workflow, how to inspect the log of past executions, and how to delete older work directories with `nextflow clean`.

<!-- Any other cool options we should include? Added log -->

### 4.1. Re-launch a workflow with `-resume`

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

??? success "Command output"

    ```console hl_lines="5"
    N E X T F L O W   ~  version 25.10.2

    Launching `hello-world.nf` [golden_cantor] DSL2 - revision: 35bd3425e5

    [62/49a1f8] sayHello | 1 of 1, cached: 1 ✔
    ```

The console output should look familiar, but there's one thing that's a little different compared to before.

Look for the `cached:` bit that has been added in the process status line (line 5), which means that Nextflow has recognized that it has already done this work and simply re-used the result from the previous successful run.

You can also see that the work subdirectory hash is the same as in the previous run.
Nextflow is literally pointing you to the previous execution and saying "I already did that over there."

!!! note

    When your re-run a pipeline with `resume`, Nextflow does not overwrite any files published outside of the work directory by any executions that were run successfully previously.

### 4.2. Inspect the log of past executions

Whether you're developing a new pipeline or running pipelines in production, at some point you'll probably need to look up information about past runs.
Here's how to do that.

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

During the development process, you'll typically run your draft pipeline a large number of times, which can lead to an accumulation of many files across many subdirectories.

Nextflow includes a convenient `clean` subcommand that can automatically delete the work subdirectories for past runs that you no longer care about, with several [options](https://www.nextflow.io/docs/latest/reference/cli.html#clean) to control what will be deleted.

Here we show you an example that deletes all subdirectories from runs before a given run, specified using its run name.
The run name is the machine-generated two-part string shown in square brackets in the `Launching (...)` console output line.

First we use the dry run flag `-n` to check what will be deleted given the command:

```bash
nextflow clean -before golden_cantor -n
```

??? success "Command output"

```console
Would remove /workspaces/training/hello-nextflow/work/a3/7be2fad5e71e5f49998f795677fd68
```

If you don't see any lines output, you either did not provide a valid run name or there are no past runs to delete. (Make sure to change `golden_cantor` in the example command to whatever is the corresponding latest run name in your log.)

If the output looks as expected and you want to proceed with the deletion, re-run the command with the `-f` flag instead of `-n`:

```bash
nextflow clean -before golden_cantor -f
```

??? success "Command output"

```console
Removed /workspaces/training/hello-nextflow/work/a3/7be2fad5e71e5f49998f795677fd68
```

This will not remove the two-character subdirectories (like `a3/` above) but it will empty their contents.

!!! Warning

    Deleting work subdirectories from past runs removes them from Nextflow's cache and deletes any outputs that were stored in those directories.
    That means it breaks Nextflow's ability to resume execution without re-running the corresponding processes.

    You are responsible for saving any outputs that you care about or plan to rely on! One option is to use the `copy` mode rather than the `symlink` mode for the `publish` directive.

### Takeaway

You know how to publish outputs to a specific directory, relaunch a pipeline without repeating steps that were already run in an identical way, and use the `nextflow clean` command to clean up old work directories.

More generally, you know how to interpret a simple Nextflow workflow, manage its execution, and retrieve outputs.

### What's next?

Take a little break, you've earned it!
When you're ready, move on to Part 2 to learn how to use channels to feed inputs into your workflow, which will allow you to take advantage of Nextflow's built-in dataflow parallelism and other powerful features.
