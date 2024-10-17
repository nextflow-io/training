# Part 1: Hello World

A "Hello, World!" is a minimalist example that is meant to demonstrate the basic syntax and structure of a programming language or software framework. The example typically consists of printing the phrase "Hello, World!" to the output device, such as the console or terminal, or writing it to a file.

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

Learn how to turn that into a step in a Nextflow workflow.

---

## 1. Very first Nextflow run

Now we're going to run a script (named `hello-world.nf`) that does the same thing as before (write 'Hello World!' to a file) but with Nextflow.

!!! note

    We're intentionally not looking at the script yet. Understanding what is the result _before_ we look into the machine will help us understand what each part does.

### 1.1. Run the workflow

```bash
nextflow run hello-world.nf
```

You console should look something like this:

```console title="Output"
N E X T F L O W  ~  version 23.10.1
Launching `hello-world.nf` [mighty_murdock] DSL2 - revision: 80e92a677c
executor >  local (1)
[4e/6ba912] process > sayHello [100%] 1 of 1 ✔
```

Congratulations, you ran your first Nextflow workflow!

The most important output here is the last line (line 4), which reports that the `sayHello` process was successfully executed once.

### 1.2. Explore the work directory

When a Nextflow workflow is run a `work` directory that stores various files is created.

Each task uses a unique directory based on its [hash](https://www.nextflow.io/docs/latest/cache-and-resume.html#task-hash) (e.g., `4e/6ba912`) within the work directory.

When a task is created, Nextflow stages the task input files, script, and other helper files into the task directory. The task writes any output files to this directory during its execution, and Nextflow uses these output files for downstream tasks and/or publishing.

!!! warning

    Your work directory won't necessarily have the same hash as the one shown above.

### 1.3. Explore the log files

Browse the `work` directory in the file explorer to find the log files and any outputs created by the task. You should find the following files:

-   **`.command.begin`**: Metadata related to the beginning of the execution of the process task
-   **`.command.err`**: Error messages (stderr) emitted by the process task
-   **`.command.log`**: Complete log output emitted by the process task
-   **`.command.out`**: Regular output (stdout) by the process task
-   **`.command.sh`**: The command that was run by the process task call
-   **`.exitcode`**: The exit code resulting from the command

In this case, look for your output in the `.command.out` file.

!!! tip

    Some of the specifics will be different in your log output. For example, here `[mighty_murdock]` and `[4e/6ba912]` are randomly generated names, so those will be different every time.

### Takeaway

You know how to run a simple Nextflow script and navigate the work directory.

### What's next?

Learn how to interpret the structure of a Nextflow pipeline by reading the code and identifying the main components.

---

## 2. Interpret the Hello World script

A Nextflow script involves two main types of core components: one or more **processes**, and the **workflow** itself. Each **process** describes what operation(s) the corresponding step in the pipeline should accomplish, while the **workflow** describes the dataflow logic that connects the various steps.

Let's open the `hello-world.nf` script and look at how it's structured.

### 2.1. Open `hello-world.nf` in the editor pane

The file is in the `hello-nextflow` directory, which should be your current working directory. You can either double click on the file in the file explorer, or type `ls` in the terminal and Cmd+Click (MacOS) or Ctrl+Click (PC) on the file to open it.

#### 2.1.1 The `process` definition

The first block of code describes a **process**. The process definition starts with the keyword `process`, followed by process name and finally the process body delimited by curly braces. The process body must contain a script block which represents the command or, more generally, a script that is executed by it.

Here we have a **process** called `sayHello` that writes its **output** to `stdout`:

```groovy title="hello-world.nf"
process sayHello {

    output:
        stdout

    """
    echo 'Hello World!'
    """
}
```

This a very minimal process definition that just contains an output definition and the script itself. In a real-world pipeline, a process usually contains additional blocks such as directives, inputs, and conditional clauses, which we'll introduce in the course of this training series.

!!! note

    The output definition does not _determine_ what output will be created. It simply _declares_ what is the expected output, so that Nextflow can look for it once execution is complete. This is necessary for verifying that the command was executed successfully and for passing the output to downstream processes if needed.

#### 2.1.2 The `workflow` definition

The second block of code describes the **workflow** itself. The workflow definition starts with the keyword `workflow`, followed by an optional name, then the workflow body delimited by curly braces.

Here we have a **workflow** that consists of one call to the `sayHello` process.

```groovy title="hello-world.nf"
workflow {
    sayHello()
}
```

This a very minimal **workflow** definition. In a real-world pipeline, the workflow typically contains multiple calls to **processes** connected by **channels**. You'll learn how to add more processes and connect them by channels in the course of this training.

### 2.2. Add a comment block above the process to document what it does

```groovy title="hello-world.nf" linenums="1"
/*
 * Use echo to print 'Hello World!' to standard out
 */
process sayHello {
```

### 2.3. Add an in-line comment above the process call

```groovy title="hello-world.nf" linenums="14"
workflow {

    // emit a greeting
    sayHello()
}
```

### Takeaway

You know how to interpret the simplest possible Nextflow script and add comments to document it.

### What's next?

Learn how to make the script output a named file.

---

## 3. Send the output to a file

Instead of printing "Hello World!" to the standard output, we'd prefer to save that output to a specific file, just like we did when running in the terminal earlier. This is how most tools that you'll run as part of real-world pipelines typically behave; we'll see examples of that later.

Both the script and the output definition blocks need to be updated.

### 3.1. Change the process command to output a named file

This is the same change we made when we ran the command directly in the terminal earlier.

_Before:_

```groovy title="hello-world.nf" linenums="9"
"""
echo 'Hello World!'
"""
```

_After:_

```groovy title="hello-world.nf" linenums="9"
"""
echo 'Hello World!' > output.txt
"""
```

### 3.2. Change the output declaration in the `sayHello` process

We need to tell Nextflow that it should now look for a specific file to be produced by the process execution.

_Before:_

```groovy title="hello-world.nf" linenums="6"
output:
    stdout
```

_After:_

```groovy title="hello-world.nf" linenums="6"
output:
    path 'output.txt'
```

!!! note

    Inputs and outputs in the process blocks typically require a qualifier and a variable name:

    ```
    <input/output qualifier> <input/output name>
    ```

    The qualifier defines the type of data to be received. This information is used by Nextflow to apply the semantic rules associated with each qualifier, and handle it properly. Common qualifiers include `val` and `path`. In the example above, `stdout` is an exception since it is not associated with a name.

### 3.3. Run the workflow again

```bash
nextflow run hello-world.nf
```

The log output should be very similar to the first time your ran the workflow:

```console title="Output"
N E X T F L O W  ~  version 23.10.1
Launching `scripts/hello-world.nf` [disturbed_cajal] DSL2 - revision: 9512241567
executor >  local (1)
[ab/c61321] process > sayHello [100%] 1 of 1 ✔
```

Like you did before, find the `work` directory in the file explorer. Find the `output.txt` output file and click on it to open it and verify that it contains the greeting as expected.

!!! warning

    This example is brittle because we hardcoded the output filename in two separate places (the script and the output blocks). If we change one but not the other, the script will break. Later, you'll learn how to use variables to avoid this problem.

### 3.4. Add a `publishDir` directive to the process

You'll have noticed that the output is buried in a working directory several layers deep. Nextflow is in control of this directory and we are not supposed to interact with it. To make the output file more accessible, we can utilize the `publishDir` directive. By specifying this directive, we are telling Nextflow to automatically copy the output file to a designated output directory. This allows us to leave the working directory alone, while still having easy access to the desired output file.

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

!!! note

    There is a newer syntax option that makes it possible to declare and publish workflow-level outputs, documented [here](https://www.nextflow.io/docs/latest/workflow.html#publishing-outputs), which makes using `publishDir` at the process level redundant once your pipeline is fully operational. However, `publishDir` is still very useful during pipeline development; that is why we include it in this training series. This will also ensure that you can read and understand the large number of pipelines that have already been written with `publishDir`.

    You'll learn how to use the workflow-level outputs syntax later in this training series.

### 3.5. Run the workflow again

```bash
nextflow run hello-world.nf
```

The log output should start looking very familiar:

```console title="Output"
N E X T F L O W  ~  version 23.10.1
Launching `hello-world.nf` [evil_bose] DSL2 - revision: 6907ac9da2
executor >  local (1)
[46/e4ff05] process > sayHello [100%] 1 of 1 ✔
```

This time, Nextflow will have created a folder called `results/`. In this folder is our `output.txt` file. If you check the contents it should match our existing output. This is how we move results files outside of the working directories.

### Takeaway

You know how to send outputs to a specific named file and use the `publishDir` directive to move files outside of the Nextflow working directory.

### What's next?

Learn how to use resume to re-use the cached results

---

## 4. Use the Nextflow resume feature

Nextflow has an option called `-resume` that allows you to re-run a pipeline you've already launched previously, in such a way that any processes that have already been run with the exact same code, settings and inputs will be skipped. Using this mode means Nextflow will only run processes that are either new, have been modified or are being provided new settings or inputs.

There are two key advantages to doing this:

-   If you're in the middle of developing your pipeline, you can iterate more rapidly since you only effectively have to run the process(es) you're activelyworking on to test your changes.
-   If you're running a pipeline in production and something goes wrong, in many cases you can fix the issue and relaunch the pipeline, and it will resume running from the point of failure, which can save you a lot of time and compute.

### 4.1. Run the workflow again with `-resume`

```bash
nextflow run hello-world.nf -resume
```

The console output should look similar.

```console title="Output"
N E X T F L O W  ~  version 23.10.1
Launching `hello-world.nf` [tiny_elion] DSL2 - revision: 7ad1cd6bfe
executor >  local (1)
[8b/1f9ded] process > sayHello [100%] 1 of 1 ✔, cached: 1 ✔
```

Notice the additional `cached`. Nextflow has cached the process and re-used the result. It will also not replace the output file at `results/output.txt`.

### Takeaway

You know how to to relaunch a pipeline without repeating steps that were already executed in an identical way.

### What's next?

Learn how to add in variable inputs.

---

## 5. Add in variable inputs using a channel

So far, we've been emitting a greeting hardcoded into the process command. Now we're going to add some flexibility by using an input variable, so that we can easily change the greeting. This is going to require us to use a **channel**; more on that in a minute.

### 5.1. Add an input definition to the process block

First we need to adapt the process definition to accept an input.

_Before:_

```groovy title="hello-world.nf" linenums="9"
process sayHello {

    output:
        path "output.txt"
```

_After:_

```groovy title="hello-world.nf" linenums="9"
process sayHello {

    publishDir 'results', mode: 'copy'

    input:
        val greeting

    output:
        path "output.txt"
```

### 5.2. Edit the process command to use the input variable

Now we swap the original hardcoded value for the input variable.

_Before:_

```groovy title="hello-world.nf" linenums="16"
"""
echo 'Hello World!' > output.txt
"""
```

_After:_

```groovy title="hello-world.nf" linenums="16"
"""
echo '$greeting' > output.txt
"""
```

### 5.3. Create an input channel

Now that our process expects an input, we need to set up that input in the workflow body. This is where channels come in: Nextflow uses channels to feed inputs to processes and ferry data between processes that are connected together.

There are multiple ways to do this, but for now, we're just going to use the simplest possible channel, containing a single value.

We're going to create the channel using the `Channel.of()` constructor, which sets up a simple value channel, and give it a hardcoded string to use as greeting by declaring `greeting_ch = Channel.of('Hello world!')`.

_Before:_

```groovy title="hello-world.nf" linenums="21"
workflow {

    // emit a greeting
    sayHello()
}
```

_After:_

```groovy title="hello-world.nf" linenums="21"
workflow {

    // create a channel for inputs
    greeting_ch = Channel.of('Hello world!')

    // emit a greeting
    sayHello()
}
```

### 5.4. Add the channel as input to the process call

Now we need to actually plug our newly created channel into the `sayHello()` process call.

_Before:_

```groovy title="hello-world.nf" linenums="26"
// emit a greeting
sayHello()
```

_After:_

```groovy title="hello-world.nf" linenums="26"
// emit a greeting
sayHello(greeting_ch)
```

### 5.5. Run the workflow command again

```bash
nextflow run hello-world.nf
```

If you made all four edits correctly, you should get another successful execution:

```console title="Output"
N E X T F L O W  ~  version 23.10.1
Launching `hello-world.nf` [maniac_euler] DSL2 - revision: 73bfbe197f
executor >  local (1)
[57/aee130] process > sayHello (1) [100%] 1 of 1 ✔
```

The result is still the same as previously; so far we're just progressively tweaking the internal plumbing to increase the flexibility of our workflow while achieving the same end result.

### Takeaway

You know how to use a simple channel to provide an input to a process.

### What's next?

Learn how to pass inputs from the command line.

---

## 6. Use CLI parameters for inputs

We want to be able to specify the input from the command line, since that is the piece that will almost always be different in subsequent runs of the workflow. Good news: we can use `params`, Nextflow's built-in workflow parameter system, which makes it easy to declare and use CLI parameters.

### 6.1. Edit the input channel declaration to use a parameter

Here we swap out the hardcoded string for `params.greeting` in the channel creation line.

_Before:_

```groovy title="hello-world.nf" linenums="23"
// create a channel for inputs
greeting_ch = Channel.of('Hello world!')
```

_After:_

```groovy title="hello-world.nf" linenums="23"
// create a channel for inputs
greeting_ch = Channel.of(params.greeting)
```

This automatically creates a parameter called `greeting` that you can use to provide a value in the command line.

### 6.2. Run the workflow again with the `--greeting` parameter

To provide a value for this parameter, simply add `--greeting <value>` to your command line.

```bash
nextflow run hello-world.nf --greeting 'Bonjour le monde!'
```

In case you're wondering, yes it's normal to have dreams where the Nextflow log output scrolls endlessly in front of you after running through a training session... Or is that just me?

```console title="Output"
N E X T F L O W  ~  version 23.10.1
Launching `hello-world.nf` [hopeful_laplace] DSL2 - revision: a8ed9a6202
executor >  local (1)
[83/dfbbbc] process > sayHello (1) [100%] 1 of 1 ✔
```

Be sure to open up the output file to check that you now have the new version of the greeting. Voilà!

!!! note

    A double hyphen (`--`) is used to set a `params` item while a single hyphen (`-`) is used to modify a Nextflow setting, _e.g._ the `-resume` feature we used earlier.

### 6.3. Set a default value for a command line parameter

In many cases, it makes sense to supply a default value for a given parameter so that you don't have to specify it for every run.

Let's initialize the `greeting` parameter with a default value by adding the parameter declaration at the top of the script (with a comment block as a free bonus).

```groovy title="hello-world.nf" linenums="1"
/*
 * Pipeline parameters
 */
params.greeting = "Bonjour le monde!"
```

### 6.4. Run the workflow again without specifying the parameter

Now that you have a default value set, you can run the workflow again without having to specify a value in the command line.

```bash
nextflow run hello-world.nf
```

The output should look the same.

```console title="Output"
N E X T F L O W  ~  version 23.10.1
Launching `hello-world.nf` [tiny_elion] DSL2 - revision: 7ad1cd6bfe
executor >  local (1)
[8b/1f9ded] process > sayHello [100%] 1 of 1 ✔
```

Check the output in the results directory, and... Tadaa! It works! Nextflow used the default value to name the output. But wait, what happens now if we provide the parameter in the command line?

### 6.5. Run the workflow again with the `--greeting` parameter on the command line using a different greeting

```bash
nextflow run hello-world.nf --greeting 'Holà!'
```

Nextflow's not complaining, that's a good sign:

```console title="Output"
N E X T F L O W  ~  version 23.10.1
Launching `hello-world.nf` [exotic_lichterman] DSL2 - revision: 7ad1cd6bfe
executor >  local (1)
[36/47354a] process > sayHello [100%] 1 of 1 ✔
```

Check the results directory and look at the contents of `output.txt`. Tadaa again!

The value of the parameter we passed on the command line overrode the value we gave the variable in the script. In fact, parameters can be set in several different ways; if the same parameter is set in multiple places, its value is determined based on the order of precedence that is described [here](https://www.nextflow.io/docs/latest/config.html).

!!! tip

    You can put the parameter declaration inside the workflow block if you prefer. Whatever you choose, try to group similar things in the same place so you don't end up with declarations all over the place.

### Takeaway

You know how to set up an input variable for a process and supply a value in the command line.

### What's next?

Learn how to add in a second process and chain them together.

---

## 7. Add a second step to the workflow

Most real-world workflows involve more than one step. Here we introduce a second process that converts the text to uppercase (all-caps), using the classic UNIX one-liner:

```bash
tr '[a-z]' '[A-Z]'`
```

We're going to run the command by itself in the terminal first to verify that it works as expected without any of the workflow code getting in the way of clarity, just like we did at the start with `echo 'Hello World'`. Then we'll write a process that does the same thing, and finally we'll connect the two processes so the output of the first serves as input to the second.

### 7.1. Run the command in the terminal by itself

```bash
echo 'Hello World' | tr '[a-z]' '[A-Z]'
```

The output is simply the uppercase version of the text string:

```console title="Output"
HELLO WORLD
```

!!! note

    This is a very naive text replacement one-liner that does not account for accented letters, so for example 'Holà' will become 'HOLà'. This is expected.

### 7.2. Make the command take a file as input and write the output to a file

As previously, we want to output results to a dedicated file, which we name by prepending the original filename with `UPPER-`.

```bash
cat output.txt | tr '[a-z]' '[A-Z]' > UPPER-output.txt
```

Now the `HELLO WORLD` output is in the new output file, `UPPER-output.txt`.

### 7.3. Wrap the command in a new Nextflow process definition

We can model our new process on the first one, since we want to use all the same components.

```groovy title="hello-world.nf" linenums="21"
/*
 * Use a text replace utility to convert the greeting to uppercase
 */
process convertToUpper {

    publishDir 'results', mode: 'copy'

    input:
        path input_file

    output:
        path "UPPER-${input_file}"

    """
    cat '$input_file' | tr '[a-z]' '[A-Z]' > UPPER-${input_file}
    """
}
```

!!! note

    As a little bonus, we composed the second output filename based on the first one. Very important to remember: you have to use double quotes around the output filename expression (NOT single quotes) or it will fail.

### 7.4. Add a call to the new process in the workflow body

Don't forget we need to tell Nextflow to actually call the process we just created! To do that, we add it to the `workflow` body.

```groovy title="hello-world.nf" linenums="36"
workflow {

    // create a channel for inputs
    greeting_ch = Channel.of(params.greeting)

    // emit a greeting
    sayHello(greeting_ch)

    // convert the greeting to uppercase
    convertToUpper()
}
```

### 7.5. Pass the output of the first process to the second process

The output of the `sayHello()` process is automatically set up as a channel called `sayHello.out`, so all we need to do is pass that as the input to the `convertToUpper()` process.

```groovy title="hello-world.nf" linenums="44"
// convert the greeting to uppercase
convertToUpper(sayHello.out)
```

For a simple case like this, that's all we need to do to connect two processes!

### 7.6. Run the same workflow command as before

Let's make sure this works:

```bash
nextflow run hello-world.nf --greeting 'Hello World!'
```

Oh, how exciting! There is now an extra line in the log output, which corresponds to the new process we just added:

```console title="Output"
N E X T F L O W  ~  version 23.10.1
Launching `hello-world.nf` [kickass_pasteur] DSL2 - revision: d15b2c482c
executor >  local (2)
[da/8d9221] process > sayHello (1)       [100%] 1 of 1 ✔
[01/2b32ee] process > convertToUpper (1) [100%] 1 of 1 ✔
```

You'll notice that this time the workflow produced two work directories; one per process instance. Check out the work directory of the task from the second process, where you should find two different output files listed. If you look carefully, you'll notice one of them (the output of the first process) has a little arrow icon on the right; that signifies it's a symbolic link. It points to the location where that file lives in the work directory of the first process. By default, Nextflow uses symbolic links to stage input files whenever possible, to avoid making duplicate copies.

!!! note

    Note how all we did was connect the output of `sayHello` to the input of `convertToUpper` and the two processes could be ran in serial. Nextflow did the hard work of handling input and output files and passing them between the two commands for us. This is the power of channels in Nextflow, doing the laborious work of connecting our pipeline steps up together.

    What's more, Nextflow will automatically determine which call needs to be executed first based on how they're connected, so the order in which they're written in the workflow body does not matter. However, we do recommend you be kind to your collaborators and to your future self, and try to write them in a logical order!

### Takeaway

You know how to add a second step that takes the output of the first step as input.

### What's next?

Learn how to make the workflow run on many values for the same input.

---

## 8. Modify the workflow to run on many values for the same input

Workflows typically run on batches of inputs that we want to process in bulk, so we want to upgrade the workflow to accept an input with multiple values.

### 8.1. Modify the channel to contain multiple greetings

For convenience, we go back to hardcoding the greetings instead of using a parameter for the input, but we'll improve on that later on.

_Before:_

```groovy title="hello-world.nf" linenums="38"
// create a channel for inputs
greeting_ch = Channel.of(params.greeting)
```

_After:_

```groovy title="hello-world.nf" linenums="38"
// create a channel for inputs
greeting_ch = Channel.of('Hello','Bonjour','Holà')
```

### 8.2. Modify the first process to generate dynamic filenames so the final filenames will be unique

Since we're going to be publishing all the outputs to the same directory, they can't all be named `output.txt`. We have to give them unique names.

_Before:_

```groovy title="hello-world.nf" linenums="9"
process sayHello {

    publishDir 'results', mode: 'copy'

    input:
        val greeting

    output:
        path "output.txt"

    """
    echo '$greeting' > "output.txt"
    """
}
```

_After:_

```groovy title="hello-world.nf" linenums="9"
process sayHello {

    publishDir 'results', mode: 'copy'

    input:
        val greeting

    output:
        path "${greeting}-output.txt"

    """
    echo '$greeting' > '$greeting-output.txt'
    """
}
```

!!! note

    In practice, naming files based on the data input itself is almost always impractical; the better way to generate dynamic filenames is to use a samplesheet and create a map of metadata (aka metamap) from which we can grab an appropriate identifier to generate the filenames. We'll show how to do that later in this training series.

### 8.3. Run the command and look at the log output

```bash
nextflow run hello-world.nf
```

How many log lines do you expect to see in the terminal? And how many do you actually see?

```console title="Output"
N E X T F L O W  ~  version 23.10.1
Launching `hello-world.nf` [cranky_hypatia] DSL2 - revision: 719dae218c
executor >  local (6)
[6c/91aa50] process > sayHello (3)       [100%] 3 of 3 ✔
[90/80111c] process > convertToUpper (3) [100%] 3 of 3 ✔
```

Something seems wrong! The log lines seem to indicate each process was executed three times (corresponding to the three input elements we provided) but we're only seeing two work directories instead of six.

This is because by default, the ANSI logging system writes the logging from multiple calls to the same process on the same line. Fortunately, we can disable that behavior.

### 8.4. Run the command again with the `-ansi-log false` option

```bash
nextflow run hello-world.nf -ansi-log false
```

This time we see six work directories in the terminal:

```console title="Output"
N E X T F L O W  ~  version 23.10.1
Launching `hello-world.nf` [disturbed_panini] DSL2 - revision: 719dae218c
[8c/77b534] Submitted process > sayHello (1)
[b5/f0bf7e] Submitted process > sayHello (2)
[a8/457f9b] Submitted process > sayHello (3)
[3d/1bb4e6] Submitted process > convertToUpper (2)
[fa/58fbb1] Submitted process > convertToUpper (1)
[90/e88919] Submitted process > convertToUpper (3)
```

That's much better; at least for this number of processes. For a complex workflow, or a large number of inputs, having the full list output to the terminal might get a bit overwhelming.

!!! tip

    Another way to show that all six calls are happening is to delete all the work directories before you run again. Then you'll see the six new ones pop up.

### Takeaway

You know how to feed an input with multiple elements through a channel.

### What's next?

Learn how to make the workflow take a file that contains multiple values for an input.

---

## 9. Modify the workflow to run on a file that contains an input with multiple values

In most cases, when we run on multiple inputs, the input values are contained in a file. Here we're going to use a `.csv` file where the values are separated by commas.

### 9.1. Modify the channel declaration to take an input file through a parameter

Here we introduce a new channel constructor, `Channel.fromPath()`, which has some built-in functionality for handling file paths. We're going to use that instead of the `Channel.of()` constructor we used previously:

`greeting_ch = Channel.fromPath(params.input_file)`

Since our goal is to read in the contents of a `.csv` file, we're going to add the `.splitCsv()` operator to make Nextflow parse the file contents accordingly, as well as the `.flatten()` operator to turn the array element produced by `.splitCsv()` into a channel of individual elements. So the line becomes:

`greeting_ch = Channel.fromPath(params.input_file).splitCsv().flatten()`

And here it is in the context of the workflow body:

_Before:_

```groovy title="hello-world.nf" linenums="38"
// create a channel for inputs
greeting_ch = Channel.of('Hello','Bonjour','Holà')
```

_After:_

```groovy title="hello-world.nf" linenums="38"
// create a channel for inputs from a file
greeting_ch = Channel.fromPath(params.input_file).splitCsv().flatten()
```

!!! tip

    While you're developing your pipeline, you can inspect the contents of any channel by adding the `.view()` operator to the name of the channel. For example, if you add `greeting_ch.view()` anywhere in the workflow body, when you run the script, Nextflow will print the channel contents to standard out.

    You can also use this to inspect the effect of the operators; for example, try comparing the output of `Channel.fromPath(params.input_file).splitCsv().view()` and `Channel.fromPath(params.input_file).splitCsv().flatten().view()`.

### 9.2. Modify the default parameter to point to an input file

_Before:_

```groovy title="hello-world.nf" linenums="38"
/*
 * Pipeline parameters
 */
params.greeting = "Bonjour le monde!"
```

_After:_

```groovy title="hello-world.nf" linenums="38"
/*
 * Pipeline parameters
 */
params.input_file = "data/greetings.txt"
```

### 9.3. Run the workflow with the `-ansi-log false` option and an `--input_file` parameter

```bash
nextflow run hello-world.nf -ansi-log false --input_file data/greetings.txt
```

Once again we see each process get executed three times:

```console title="Output"
N E X T F L O W  ~  version 23.10.1
Launching `hello-world.nf` [small_albattani] DSL2 - revision: 5cea973c3c
[45/18d159] Submitted process > sayHello (1)
[cf/094ea1] Submitted process > sayHello (3)
[27/e3ea5b] Submitted process > sayHello (2)
[7d/63672f] Submitted process > convertToUpper (1)
[62/3184ed] Submitted process > convertToUpper (2)
[02/f0ff38] Submitted process > convertToUpper (3)
```

Looking at the outputs, we see each greeting was correctly extracted and processed through the workflow. We've achieved the same result as the previous step, but now we have a lot more flexibility to add more elements to the channel of greetings we want to process.

!!! tip

    Don't worry if the channel types and operators feel like a lot to grapple with the first time you encounter them. The key learning point here is that we can create a channel from the contents of a file. You'll get more opportunities to practice using these components in various settings in later training modules.

### Takeaway

You know how to provide inputs in a file.

### What's next?

Celebrate your success and take a break!

When you are ready, move on to Part 2 of this training to learn how to apply what you've learned to a more realistic data analysis use case.
