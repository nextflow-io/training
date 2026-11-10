# Part 1: Hello World

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/tOukLxWCHiA?si=y8lAedhEHWaTV4zd&amp;list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&amp;cc_load_policy=1" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

/// caption
:fontawesome-brands-youtube:{ .youtube } See [the whole playlist](https://youtube.com/playlist?list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&si=eF7cLR62goy-lc6n) on the Nextflow YouTube channel.

:green_book: The video transcript is available [here](./transcripts/01_hello_world.md).
///

In this course, you write the Hello pipeline you ran in [Nextflow Run](../nextflow_run/index.md), one piece at a time; see the [course overview](./index.md) for the plan.

In this first part, you start from a script that does a single thing, writing 'Hello World!' to a file.
You examine and run it, publish its output, and make it greet a name provided on the command line.

---

## 1. Examine the script and run it

We provide you with a minimal workflow script named `hello-world.nf`.
It does the same thing as running `echo 'Hello World!' > output.txt` in the terminal, but with Nextflow.

!!! tip

    Remember that you should now be inside the `hello-nextflow/` directory as described on the [Getting Started](00_orientation.md) page.

### 1.1. Examine the code

Open the `hello-world.nf` script in the editor pane.

??? full-code "Full code file"

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

A Nextflow workflow script typically includes one or more [**process**](https://nextflow.io/docs/latest/process.html) definitions and the [**workflow**](https://nextflow.io/docs/latest/workflow.html) itself, plus a few optional blocks (not present here) which we'll introduce later on.

Each **process** describes what operation(s) the corresponding step in the pipeline should accomplish, while the **workflow** describes the dataflow logic that connects the various steps.

We're going to take a closer look at the **process** block first, then we'll look at the **workflow** block.

#### 1.1.1. The `process` definition

The first block of code describes a **process**.

The process definition starts with the keyword `process`, followed by the process name and finally the process body delimited by curly braces.
The process body must contain a script block which specifies the command to run, which can be anything you would be able to run in a command line terminal.

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

Here we have a **process** called `sayHello` that writes its **output** to a file named `output.txt`.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello_world.svg"
</figure>

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

In the terminal, run the following command:

```bash
nextflow run hello-world.nf
```

??? success "Command output"

    ```console hl_lines="6"
    N E X T F L O W   ~  version 26.04.4

    Launching `hello-world.nf` [grave_wiles] revision: dc38bebf45

    executor >  local (1)
    [11/2ab905] sayHello | 1 of 1 ✔
    ```

The `sayHello` process ran once (`1 of 1 ✔`).
As you saw in Nextflow Run, the prefix `[11/2ab905]` is the truncated path to the task directory under `work/`, where Nextflow ran the process.

Open that directory, for example `work/11/2ab905b8f0bf001eaefd9112486273` (yours will have a different name), and look at two of its files:

- **`output.txt`** contains `Hello World!`, the output declared in the process.
- **`.command.sh`** contains the command Nextflow ran:

  ```console title=".command.sh"
  #!/bin/bash -ue
  echo 'Hello World!' > output.txt
  ```

The `script` block of the process became that command.
Here the command is hardcoded, but later on the script will include variables, and `.command.sh` is where you check how Nextflow filled them in.

??? question "Don't see the `.command.sh` file?"

    Files whose names start with a dot are hidden in the terminal.
    Use `ls -a` to list them, or browse the task directory in the VS Code file explorer.

### Takeaway

You know the basic structure of a Nextflow script, with a process that defines a step and a workflow that calls it, and how the script block of a process becomes the command that Nextflow runs.

### What's next?

Make the output easy to find by publishing it.

---

## 2. Publish outputs

The output of the `sayHello` process is buried in a task directory under `work/`.
This is on purpose: Nextflow is in control of the `work/` directory, and you are not supposed to interact with it.

In Nextflow Run, the pipeline published its outputs to a `results/` directory.
The pipeline author sets that up in the workflow script using [workflow output definitions](https://nextflow.io/docs/latest/workflow.html#workflow-outputs).
This involves two new pieces of code:

1. A `publish:` section inside the `workflow` body, listing the process outputs to publish.
2. An `output` block, outside the `workflow` body, specifying how to publish them.

### 2.1. Default settings

Start by publishing the output of `sayHello` with Nextflow's default settings.

#### 2.1.1. Declare the output of the `sayHello` process

Add a `publish:` section to the workflow body (the same kind of code element as the `main:` section) and list the output of the `sayHello` process.

In the workflow script file `hello-world.nf`, add the following lines of code:

=== "After"

    ```groovy title="hello-world.nf" linenums="17" hl_lines="7-8"
    workflow {

        main:
        // emit a greeting
        sayHello()

        publish:
        greetings = sayHello.out
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

You refer to the output of the process with `sayHello.out`, and assign it a name of your choice, here `greetings`.

#### 2.1.2. Add an `output` block to the script

Now add the `output` block, which describes how to publish each output listed in the `publish:` section.
This block sits **outside** and **below** the `workflow` block.

In the workflow script file `hello-world.nf`, add the following lines of code:

=== "After"

    ```groovy title="hello-world.nf" linenums="17" hl_lines="11-14"
    workflow {

        main:
        // emit a greeting
        sayHello()

        publish:
        greetings = sayHello.out
    }

    output {
        greetings {
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
        greetings = sayHello.out
    }
    ```

The entry `greetings` matches the name in the `publish:` section.
It is empty for now, so Nextflow publishes the output with its default settings.

#### 2.1.3. Run the workflow

Run the modified workflow script:

```bash
nextflow run hello-world.nf
```

??? success "Command output"

    ```console
    N E X T F L O W   ~  version 26.04.4

    Launching `hello-world.nf` [evil_carson] revision: 43c662561c

    executor >  local (1)
    [99/bf7b2a] sayHello | 1 of 1 ✔

    Outputs:

      /workspaces/training/hello-nextflow/results

      greetings: output.txt
    ```

The console output now ends with the `Outputs:` summary you saw in Nextflow Run, and Nextflow has created a `results/` directory:

??? abstract "Directory contents"

    ```console
    results
    └── output.txt -> /workspaces/training/hello-nextflow/work/99/bf7b2a1e183d53bad9c4dc7553ec96/output.txt
    ```

With the default settings, the output is published at the top level of `results/`, as a symbolic link to the file in the task directory.
As you saw in Nextflow Run, the person running the pipeline can change the output directory with `-output-dir` and choose a different default publish mode.

### 2.2. Customize the output settings

As the pipeline author, you can set where each output is published within the output directory, and how.

#### 2.2.1. Publish to a subdirectory

In the workflow script file `hello-world.nf`, add a `path` to the `greetings` output:

=== "After"

    ```groovy title="hello-world.nf" linenums="27" hl_lines="3"
    output {
        greetings {
            path 'hello'
        }
    }
    ```

=== "Before"

    ```groovy title="hello-world.nf" linenums="27"
    output {
        greetings {
        }
    }
    ```

The `path` sets where the files are published within the output directory, here a `hello` subdirectory.

#### 2.2.2. Set the publish mode

In the workflow script file `hello-world.nf`, add a `mode` to the `greetings` output:

=== "After"

    ```groovy title="hello-world.nf" linenums="27" hl_lines="4"
    output {
        greetings {
            path 'hello'
            mode 'copy'
        }
    }
    ```

=== "Before"

    ```groovy title="hello-world.nf" linenums="27"
    output {
        greetings {
            path 'hello'
        }
    }
    ```

The `mode` sets how the files are published.
With `'copy'`, Nextflow publishes a copy of the file instead of a symbolic link, so the results stay available even if you delete the `work/` directory.
A mode set in the `output` block takes precedence over the default that users can set in their configuration.

#### 2.2.3. Run the workflow

Run the workflow again:

```bash
nextflow run hello-world.nf
```

??? success "Command output"

    ```console
    N E X T F L O W   ~  version 26.04.4

    Launching `hello-world.nf` [cranky_kalam] revision: afe6b9ece6

    executor >  local (1)
    [0d/1b76d4] sayHello | 1 of 1 ✔

    Outputs:

      /workspaces/training/hello-nextflow/results

      greetings: hello/output.txt
    ```

The `Outputs:` summary now lists the file under `hello/`.

??? abstract "Directory contents"

    ```console
    results
    ├── hello
    │   └── output.txt
    └── output.txt -> /workspaces/training/hello-nextflow/work/99/bf7b2a1e183d53bad9c4dc7553ec96/output.txt
    ```

The new `output.txt` in `results/hello/` is a regular file, not a symbolic link.
The symbolic link at the top level of `results/` is left over from the previous run: Nextflow does not delete files published by earlier runs.

### 2.3. Note on process-level `publishDir` directives

Until very recently, the established way to publish outputs was to do it at the level of each individual process using a `publishDir` directive.

To achieve what we just did for the outputs of the `sayHello` process, we would have instead added the following line to the process definition:

```groovy title="hello-world.nf" linenums="6" hl_lines="3"
process sayHello {

    publishDir 'results/hello', mode: 'copy'

    output:
    path 'output.txt'

    script:
    """
    echo 'Hello World!' > output.txt
    """
}
```

You will still find this code pattern all over the place in older Nextflow pipelines and process modules, so it's important to be aware of it.
However, we do not recommend using it in any new work as it will eventually be disallowed in future versions of the Nextflow language.

### Takeaway

You know how to publish the outputs of a workflow with a `publish:` section and an `output` block, and how to set the subdirectory and publish mode of each output.

### What's next?

Learn to greet a name provided via a command-line parameter, and give it a default value.

---

## 3. Use a variable input passed on the command line

In its current state, our workflow writes the same hardcoded greeting every time.
We want to add some flexibility by using an input variable, so that it can greet a name we provide at runtime.

This requires us to make three sets of changes to our script:

1. Change the process to expect a variable input
2. Set up a command-line parameter to capture user input
3. Pass the input to the process in the workflow body

Let's make these changes one at a time.

### 3.1. Change the `sayHello` process to expect a variable input

We need to edit the process definition to (1) accept an input variable and (2) use that variable in the command line.

#### 3.1.1. Add an input block to the process definition

First, let's adapt the process definition to accept an input called `name`.

In the process block, make the following code change:

=== "After"

    ```groovy title="hello-world.nf" linenums="3" hl_lines="2 6-7"
    /*
     * Use echo to print a greeting to a file
     */
    process sayHello {

        input:
        val name

        output:
        path 'output.txt'
    ```

=== "Before"

    ```groovy title="hello-world.nf" linenums="3" hl_lines="2"
    /*
     * Use echo to print 'Hello World!' to a file
     */
    process sayHello {

        output:
        path 'output.txt'
    ```

The `name` variable is prefixed by `val` to tell Nextflow it's a value (not a path).
We also update the comment above the process, since it will no longer always print 'Hello World!'.

#### 3.1.2. Edit the process command to use the input variable

Now we insert the input variable into the greeting, in place of the hardcoded `World`.

In the process block, make the following code change:

=== "After"

    ```groovy title="hello-world.nf" linenums="14" hl_lines="3"
        script:
        """
        echo 'Hello, ${name}!' > output.txt
        """
    ```

=== "Before"

    ```groovy title="hello-world.nf" linenums="14" hl_lines="3"
        script:
        """
        echo 'Hello World!' > output.txt
        """
    ```

The `$` symbol and curly braces (`{ }`) tell Nextflow this is a variable name that needs to be replaced with the actual input value (=interpolated).

!!! tip

    The curly braces (`{ }`) were technically optional in previous versions of Nextflow, so you might see older workflows where this is written as `echo 'Hello, $name!' > output.txt`.

Now that the `sayHello()` process is ready to accept a variable input, we need a way to provide an input value to the process call at the workflow level.

### 3.2. Set up a command-line parameter to capture user input

We could simply hardcode an input directly by making the process call `sayHello('World')`.
However, when we're doing real work with our workflow, we're going to want to be able to control its inputs from the command line, so we can do something like this:

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello_world_input.svg"
</figure>

Fortunately, Nextflow has a built-in workflow parameter system called [`params`](https://nextflow.io/docs/latest/config.html#params) that makes it easy to declare and use CLI parameters.

The general syntax is to declare `params.<parameter_name>` to tell Nextflow to expect a `--<parameter_name>` parameter on the command line.

Here, we want to create a parameter called `--name`, so we need to declare `params.name` somewhere in the workflow.
In principle we can write it anywhere; but since we're going to want to give it to the `sayHello()` process call, we can plug it in there directly by writing `sayHello(params.name)`.

In the workflow block, make the following code change:

=== "After"

    ```groovy title="hello-world.nf" linenums="23" hl_lines="2"
        // emit a greeting
        sayHello(params.name)
    ```

=== "Before"

    ```groovy title="hello-world.nf" linenums="23" hl_lines="2"
        // emit a greeting
        sayHello()
    ```

This tells Nextflow to run the `sayHello` process on the value provided through the `--name` parameter.

In effect, we've accomplished steps (2) and (3) outlined at the start of the section in a single go.

### 3.3. Run the workflow command

Let's run it!

```bash
nextflow run hello-world.nf --name 'Mateo'
```

??? success "Command output"

    ```console
    N E X T F L O W   ~  version 26.04.4

    Launching `hello-world.nf` [hungry_boyd] revision: c7b2d15844

    executor >  local (1)
    [c2/128cf0] sayHello | 1 of 1 ✔

    Outputs:

      /workspaces/training/hello-nextflow/results

      greetings: hello/output.txt
    ```

If you made all these edits correctly, you should get another successful execution.

Be sure to open up the output file to check that the greeting now includes the name you provided.

??? abstract "File contents"

    ```console title="results/hello/output.txt"
    Hello, Mateo!
    ```

Note how the new execution has overwritten the output file published to the `results/hello/` directory.

### 3.4. Use default values for command line parameters

Ok, that was convenient, but in many cases, it makes sense to supply a default value for a given parameter so that you don't have to specify it for every run.

#### 3.4.1. Set a default value for the CLI parameter

Let's give the `name` parameter a default value by declaring it before the workflow definition.

In the workflow script, make the following code change:

=== "After"

    ```groovy title="hello-world.nf" linenums="18" hl_lines="3-8"
    }

    /*
     * Pipeline parameters
     */
    params {
        name: String = 'Friend'
    }

    workflow {
    ```

=== "Before"

    ```groovy title="hello-world.nf" linenums="18"
    }

    workflow {
    ```

As you see, we can specify the type of input that the workflow expects (Nextflow 25.10.2 and later).
The syntax is `#!groovy name: Type = default_value`.
Supported types include `String`, `Integer`, `Float`, `Boolean`, and `Path`.

!!! info

    In older workflows, you may see that whole `params` block written as just `name = 'Friend'`.

As you add more parameters to your pipeline, you should add them all to this block, whether or not you need to give them a default value.
This will make it easy to find all configurable parameters at a glance.

#### 3.4.2. Run the workflow again without specifying the parameter

Now that you have a default value set, you can run the workflow again without having to specify a value in the command line.

```bash
nextflow run hello-world.nf
```

??? success "Command output"

    ```console
    N E X T F L O W   ~  version 26.04.4

    Launching `hello-world.nf` [high_lamarck] revision: 07609cbd83

    executor >  local (1)
    [23/6db169] sayHello | 1 of 1 ✔

    Outputs:

      /workspaces/training/hello-nextflow/results

      greetings: hello/output.txt
    ```

??? question "If it didn't work"

    If that failed with an error that looks like this:

    ```
    ERROR ~ Script compilation error
    - file : /workspaces/training/hello-nextflow/hello-world.nf
    - cause: you tried to assign a value to the class 'java.lang.String'
    @ line 24, column 11.
          name: String = 'Friend'
                 ^

    1 error


    -- Check '.nextflow.log' file for details
    ```

    Then you're probably using the older v1 Nextflow language parser.
    This was mentioned at the start of the course, but maybe you missed it.
    Check the [Nextflow versions](../info/nxf_versions.md) help material.

    The v2 parser is the default from Nextflow 26.04 onward, so you will only see this on earlier versions.
    On a version before 26.04 you need to enable the v2 language parser:

    ```bash
    export NXF_SYNTAX_PARSER=v2
    ```

The output will be in the same place as previously, but the contents should be updated with the new text.

??? abstract "File contents"

    ```console title="results/hello/output.txt"
    Hello, Friend!
    ```

Nextflow used the default value of the `name` parameter to create the output.

#### 3.4.3. Override the default value

If you provide the parameter on the command line, the CLI value will override the default value.

Try it out:

```bash
nextflow run hello-world.nf --name 'Celine'
```

??? success "Command output"

    ```console
    N E X T F L O W   ~  version 26.04.4

    Launching `hello-world.nf` [compassionate_yonath] revision: 07609cbd83

    executor >  local (1)
    [7d/69c239] sayHello | 1 of 1 ✔

    Outputs:

      /workspaces/training/hello-nextflow/results

      greetings: hello/output.txt
    ```

Once again, you should find the corresponding updated output in your results directory.

??? abstract "File contents"

    ```console title="results/hello/output.txt"
    Hello, Celine!
    ```

!!! note

    In Nextflow, there are multiple places where you can specify values for parameters.
    If the same parameter is set to different values in multiple places, Nexflow will determine what value to use based on the order of precedence that is described [here](https://www.nextflow.io/docs/latest/config.html).

    The [Configure the pipeline](../nextflow_run/02_configure_pipeline.md) part of Nextflow Run covers this in more detail.

### Takeaway

You know how to insert a variable input into a process command, provide it at runtime via a command-line parameter, and set up, use and override default values.

### What's next?

Take a little break, you've earned it!

When you're ready, move on to [**Part 2: Hello Channels**](./02_hello_channels.md) to learn how to use channels to feed inputs into your workflow, which will allow you to take advantage of Nextflow's built-in dataflow parallelism and other powerful features.

---

## Bonus exercise

This optional exercise lets you practice what you learned in this part.
It can be solved with what this part covered; open the hint if you get stuck.
The solution is also in `solutions/1-hello-world/bonus/`.

### Exercise: Rename the output file

Starting from `hello-world.nf` as you left it at the end of this part, change the `sayHello` process so that it writes the greeting to a file called `greeting.txt` instead of `output.txt`.
Run the workflow and check that `greeting.txt` is published to the `results/hello/` directory.

??? tip "Hint"

    The file name appears in two places in the process: the `output:` declaration and the `script` block.
    If you change only one of them, the run fails with a `Missing output file(s)` error, because the file the process declares is not the file its script wrote.

??? solution

    Change the file name in both places:

    ```groovy title="hello-world.nf" linenums="6" hl_lines="7 11"
    process sayHello {

        input:
        val name

        output:
        path 'greeting.txt'

        script:
        """
        echo 'Hello, ${name}!' > greeting.txt
        """
    }
    ```

    Run the workflow:

    ```bash
    nextflow run hello-world.nf
    ```

    The `Outputs:` section of the console output shows the new file name:

    ```console
    Outputs:

      /workspaces/training/hello-nextflow/results

      greetings: hello/greeting.txt
    ```

    The `output.txt` file in `results/hello/` is left over from earlier runs: Nextflow does not delete files published by previous runs.

    If you changed only the `script` block, the run fails instead:

    ```console
    ERROR ~ Error executing process > 'sayHello'

    Caused by:
      Missing output file(s) `output.txt` expected by process `sayHello`


    Command executed:

      echo 'Hello, Friend!' > greeting.txt
    ```

    The command ran successfully, but Nextflow could not find the file declared in the `output:` block.

---

## Quiz

<quiz>
What are the minimum required components of a Nextflow process?
- [ ] Input and output blocks only
- [x] Output and script blocks
- [ ] Input, output, and script blocks
- [ ] Only a script block

Learn more: [1.1.1. The process definition](#111-the-process-definition)
</quiz>

<quiz>
What is the purpose of the output block in a process?
- [ ] To print results to the console
- [ ] To save files to the work directory
- [x] To declare expected outputs from the process
- [ ] To define environment variables

Learn more: [1.1.1. The process definition](#111-the-process-definition)
</quiz>

<quiz>
Which file in a task directory contains the command that Nextflow ran?
- [ ] `.command.run`
- [x] `.command.sh`
- [ ] `.command.log`
- [ ] `.command.out`

Learn more: [1.2. Run the workflow](#12-run-the-workflow)
</quiz>

<quiz>
How does an entry in the `output` block relate to the workflow?
- [ ] It must have the same name as the process
- [x] It must have the same name as an output listed in the `publish:` section
- [ ] It must have the same name as the output file
- [ ] It can have any name; Nextflow matches outputs by order

Learn more: [2.1.2. Add an `output` block to the script](#212-add-an-output-block-to-the-script)
</quiz>

<quiz>
What is the default mode for publishing workflow outputs?
- [ ] Copy files to the output directory
- [x] Create symbolic links in the output directory
- [ ] Move files to the output directory
- [ ] Compress files in the output directory

Learn more: [2.1.3. Run the workflow](#213-run-the-workflow)
</quiz>

<quiz>
How do you pass a parameter value to a Nextflow workflow from the command line?
- [ ] `-parameter value`
- [ ] `--parameter:value`
- [x] `--parameter value`
- [ ] `-p parameter=value`

Learn more: [3.2. Set up a command-line parameter to capture user input](#32-set-up-a-command-line-parameter-to-capture-user-input)
</quiz>

<quiz>
How do you reference a variable inside a Nextflow script block?
- [ ] Use `%variable%` syntax
- [x] Use `#!groovy ${variable}` syntax
- [ ] Use `{{variable}}` syntax
- [ ] Use `[variable]` syntax

Learn more: [3.1.2. Edit the process command to use the input variable](#312-edit-the-process-command-to-use-the-input-variable)
</quiz>

<quiz>
In `#!groovy params { name: String = 'Friend' }`, what happens when you run the workflow with `--name 'Celine'`?
- [ ] The workflow fails because the parameter already has a value
- [x] The command-line value overrides the default
- [ ] The default value overrides the command-line value
- [ ] Both values are passed to the process

Learn more: [3.4.3. Override the default value](#343-override-the-default-value)
</quiz>
