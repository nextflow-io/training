# Part 1: Hello World

A "Hello, World!" is a minimalist example that is meant to demonstrate the basic syntax and structure of a programming language or software framework. The example typically consists of printing the phrase "Hello, World!" to the output device, such as the console or terminal, or writing it to a file.

---

## 0. Warmup: Run Hello World directly

Let's demonstrate this with a simple command that we run directly in the terminal, to show what it does before we wrap it in Nextflow.

#### 1. Make the terminal say hello

```bash
echo 'Hello World!'
```

#### 2. Now make it write the text output to a file

```bash
echo 'Hello World!' > output.txt
```

#### 3. Verify that the output file is there using the `ls` command

```bash
ls
```

#### 4. Show the file contents

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

#### 1. Run the workflow

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

When a Nextflow workflow is run a `work` directory that stores various files is created.

Each task uses a unique directory based on its [hash](https://www.nextflow.io/docs/latest/cache-and-resume.html#task-hash) (e.g., `4e/6ba912`) within the work directory.

When a task is created, Nextflow stages the task input files, script, and other helper files into the task directory. The task writes any output files to this directory during its execution, and Nextflow uses these output files for downstream tasks and/or publishing.

!!! warning

    Your work directory won't necessarily have the same hash as the one shown above.

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

Learn how to interpret the Nextflow code.

---

## 2. Interpret the Hello World script

Nextflow scripts is built up of multiple parts.

A **process** is the basic processing primitive to execute a user script.

The process definition starts with the keyword `process`, followed by process name and finally the process body delimited by curly braces. The process body must contain a script block which represents the command or, more generally, a script that is executed by it.

A process may contain any of the following definition blocks: directives, inputs, outputs, when clauses, and of course, the script.

A **workflow** is a composition of processes and dataflow logic.

The workflow definition starts with the keyword `workflow`, followed by an optional name, and finally the workflow body delimited by curly braces.

Processes are connected through asynchronous first-in, first-out (FIFO) queues, called **channels**. The interaction between processes, and ultimately the workflow execution flow itself, are defined by the process input and output declarations.

Let's open the `hello-world.nf` script and look at how it's structured.

#### 1. Double click on the file in the file explorer to open it in the editor pane

!!! tip

    The file is in the current directory. Optionally, you can type `ls` in the terminal and Ctrl+Click on the file to open it. If you're on macOS, you can use Cmd+Click.

The first block of code describes a **process** called `sayHello` that writes its **output** to `stdout`:

```groovy title="hello-world.nf"
process sayHello {

    output:
        stdout

    """
    echo 'Hello World!'
    """
}
```

The second block of code describes the **workflow** itself, which consists of one call to the `sayHello` process.

```groovy title="hello-world.nf"
workflow {
    sayHello()
}
```

#### 2. Add a comment block above the process to document what it does in plain English

```groovy title="hello-world.nf" linenums="1"
/*
 * Use echo to print 'Hello World!' to standard out
 */
process sayHello {
```

#### 3. Add an in-line comment above the process call

```groovy title="hello-world.nf" linenums="14"
workflow {

    // emit a greeting
    sayHello()
}
```

### Takeaway

You know how to interpret the simplest possible Nextflow script and add comments to document it.

### What's next?

Learn how to make it output a named file.

---

## 3. Send the output to a file

Instead of printing "Hello World!" to the standard output it can be saved to a file (it's the same thing we did when running in the terminal earlier).

In a real-world workflow, this is like having a command that specifies an output file as part of its normal syntax. We'll see examples of that later.

Bother the script and the output definition blocks need to be updated.

!!! note

    Inputs and outputs in the process blocks typically require a qualifier and a variable name:

    ```
    <input/output qualifier> <input/output name>
    ```

    A definition consists of a qualifier and a name. The qualifier defines the type of data to be received. This information is used by Nextflow to apply the semantic rules associated with each qualifier, and handle it properly. Common qualifiers include `val` and `path`.

#### 1. Change the process command to output a named file

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

#### 2. Change the output declaration in the process

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

#### 3. Run the workflow again

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

    This example is brittle because we hardcoded the output filename in two separate places (the script and the output blocks). If we change one but not the other, the script will break.

### Takeaway

You know how to send outputs to a specific named file.

### What's next?

Learn how to pass parameters to the workflow from the command line.

---

## 4. Use a command line parameter for naming the output file

Here we introduce `params` (short for 'parameters') as the construct that can hold command line arguments. This is useful because there will be many parameters, such as filenames and processing options, that you may want to decide at the time you run the workflow. Parameters allow you to do this without editing the script itself every time.

Parameters can be created by prefixing a parameter names with the params scope (e.g., `params.output_file`). When including these in a script block, a `$` must be used to treat it like a variable.

Parameters can be modified when you run your workflow by adding a double hypen (`--`) to the start of the parameter name and including it in the run command (e.g., `nextflow run hello-world --output_file results).

#### 1. Change the output declaration in the process to use a parameter

_Before:_

```groovy title="hello-world.nf" linenums="6"
output:
    path 'output.txt'
```

_After:_

```groovy title="hello-world.nf" linenums="6"
output:
    path params.output_file
```

#### 2. Change the process command to use the parameter too

_Before:_

```groovy title="hello-world.nf" linenums="9"
"""
echo 'Hello World!' > output.txt
"""
```

_After:_

```groovy title="hello-world.nf" linenums="9"
"""
echo 'Hello World!' > $params.output_file
"""
```

#### 3. Run the workflow again with the `--output_file` parameter

```bash
nextflow run hello-world.nf --output_file 'output.txt'
```

The log output should start looking very familiar:

```console title="Output"
N E X T F L O W  ~  version 23.10.1
Launching `hello-world.nf` [evil_bose] DSL2 - revision: 6907ac9da2
executor >  local (1)
[46/e4ff05] process > sayHello [100%] 1 of 1 ✔
```

Follow the same procedure as before to find the `output.txt` output file. If you want to convince yourself that the parameter is working as intended, feel free to repeat this step with a different output filename.

!!! warning

    If you forget to add the output filename parameter, you get a warning and the output file is called `null`. If you add it but don't give it a value, the output file is called `true`.

!!! tip

    Command-line arguments take one dash (-) for Nextflow options, two dashes (--) for workflow parameters.

### Takeaway

You know how to use a command line parameter to set the output filename.

### What's next?

Learn how to set a default value in case we leave out the parameter.

---

## 5. Set a default value for a command line parameter

In many cases, it makes sense to supply a default value for a given parameter so that you don't have to specify it for every run.

Let's initialize the `output_file` parameter with a default value.

#### 1. Add the parameter declaration at the top of the script (with a comment block as a free bonus)

```groovy title="hello-world.nf" linenums="1"
/*
 * Pipeline parameters
 */
params.output_file = 'output.txt'

```

#### 2. Run the workflow again without specifying the parameter

```bash
nextflow run hello-world.nf
```

The console output is expected to look the same...

```console title="Output"
N E X T F L O W  ~  version 23.10.1
Launching `hello-world.nf` [tiny_elion] DSL2 - revision: 7ad1cd6bfe
executor >  local (1)
[8b/1f9ded] process > sayHello [100%] 1 of 1 ✔
```

Check the output in the work directory, and... Tadaa! It works! Nextflow used the default value to name the output. But wait, what happens now if we provide the parameter in the command line?

#### 3. Run the workflow again with the `--output_file` parameter on the command line using a DIFFERENT filename

```bash
nextflow run hello-world.nf --output_file 'output-cli.txt'
```

Nextflow's not complaining, that's a good sign:

```console title="Output"
N E X T F L O W  ~  version 23.10.1
Launching `hello-world.nf` [exotic_lichterman] DSL2 - revision: 7ad1cd6bfe
executor >  local (1)
[36/47354a] process > sayHello [100%] 1 of 1 ✔
```

Check the output directory and look for the output with the new filename. Tadaa again!

The value of the parameter we passed on the command line overrode the value we gave the variable in the script. In fact, parameters can be set in several different ways; if the same parameter is set in multiple places, its value is determined based on the order of precedence that is described [here](https://www.nextflow.io/docs/latest/config.html).

!!! tip

    You can put the parameter declaration inside the workflow block if you prefer. Whatever you choose, try to group similar things in the same place so you don't end up with declarations all over the place.

### Takeaway

You know how to handle command line parameters and set default values.

### What's next?

Learn how to add in variable inputs.

---

## 6. Add in variable inputs

So far, we've been emitting a greeting hardcoded into the process command. Now we're going to add some flexibility by introducing channels.

Nextflow is based on the dataflow programming model in which processes communicate through channels.

Channels are created used channel factories methods. There are several types of channel factories which can be utilized for creating different channel types for different data types.

Importantly, there are two kinds of channels (queue channels and value channels) which behave differently.

**Queue channel**

-   A non-blocking unidirectional first-in first-out queue connecting a producer process (i.e. outputting a value) to a consumer process, or an operators
-   Can be consumed only once

**Value channel**

-   Can be bound (i.e. assigned) with one and only one value
-   Can be consumed any number of times

We're going to start by creating a value channel with the `Channel.of()` channel factory.

!!! tip

    You can build [different kinds of channels](https://www.nextflow.io/docs/latest/channel.html#channel-types) depending on the shape of the input data; we'll cover how to deal with simple inputs later, but more complex input channel types are out of scope for this training.

#### 1. Create an input channel (with a bonus in-line comment)

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

#### 2. Add the channel as input to the process call

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

#### 3. Add an input definition to the process block

_Before:_

```groovy title="hello-world.nf" linenums="9"
process sayHello {

    output:
        path params.output_file
```

_After:_

```groovy title="hello-world.nf" linenums="9"
process sayHello {

    input:
        val greeting

    output:
        path params.output_file
```

#### 4. Edit the process command to use the input variable

_Before:_

```groovy title="hello-world.nf" linenums="16"
"""
echo 'Hello World!' > $params.output_file
"""
```

_After:_

```groovy title="hello-world.nf" linenums="16"
"""
echo '$greeting' > $params.output_file
"""
```

#### 5. Run the workflow command again

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

## 7. Use params for inputs too

We want to be able to specify the input from the command line because that is the piece that will almost always be different in subsequent runs of the workflow. Good news: we can use the same `params` construct we used for the output filename.

#### 1. Edit the input channel declaration to use a parameter

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

#### 2. Run the workflow again with the `--greeting` parameter

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

    The current form of the script doesn't have a variable declaration for `greeting` so that parameter is **REQUIRED** to be included in the command line. If we wanted, we could put in a default value by adding for example `params.greeting = 'Holà el mundo!'` at the top of the script (just like we did for the output filename). But it's less common to want to have a default value set for the input data.

### Takeaway

You know how to set up an input variable for a process and supply a value in the command line.

### What's next?

Learn how to add in a second process and chain them together.

---

## 8. Add a second step to the workflow

Most real-world workflows involve more than one step. Here we introduce a second process that converts the text to uppercase (all-caps), using the classic UNIX one-liner:

```bash
tr '[a-z]' '[A-Z]'`
```

We're going to run the command by itself in the terminal first to verify that it works as expected without any of the workflow code getting in the way of clarity, just like we did at the start with `echo 'Hello World'`. Then we'll write a process that does the same thing, and finally we'll connect the two processes so the output of the first serves as input to the second.

#### 1. Run the command in the terminal by itself

```bash
echo 'Hello World' | tr '[a-z]' '[A-Z]'
```

The output is simply the uppercase version of the text string:

```console title="Output"
HELLO WORLD
```

#### 2. Make the command take a file as input and write the output to a file

```bash
cat output.txt | tr '[a-z]' '[A-Z]' > UPPER-output.txt
```

Now the `HELLO WORLD` output is in the new output file, `UPPER-output.txt`.

#### 3. Turn that into a process definition (documented with a comment block)

```groovy title="hello-world.nf" linenums="21"
/*
 * Use a text replace utility to convert the greeting to uppercase
 */
process convertToUpper {
    input:
        path input_file

    output:
        path "UPPER-${input_file}"

    """
    cat $input_file | tr '[a-z]' '[A-Z]' > UPPER-${input_file}
    """
}
```

#### 4. Add a call to the new process in the workflow block

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

#### 5. Pass the output of the first process to the second process

```groovy title="hello-world.nf" linenums="44"
// convert the greeting to uppercase
convertToUpper(sayHello.out)
```

#### 6. Run the same workflow command as before

```bash
nextflow run hello-world.nf --greeting 'Hello World!'
```

Oh, how exciting! There is now an extra line in the log output, which corresponds to the second process we've added:

```console title="Output"
N E X T F L O W  ~  version 23.10.1
Launching `hello-world.nf` [kickass_pasteur] DSL2 - revision: d15b2c482c
executor >  local (2)
[da/8d9221] process > sayHello (1)       [100%] 1 of 1 ✔
[01/2b32ee] process > convertToUpper (1) [100%] 1 of 1 ✔
```

This time the workflow produced two work directories; one per process instance (task). Check out the work directory of the task from the second process, where you should find two different output files listed. If you look carefully, you'll notice one of them (the output of the first process) has a little arrow icon on the right; that signifies it's a symbolic link. It points to the location where that file lives in the work directory of the first process.

!!! note

    As a little bonus, we composed the second output filename based on the first one. Very important to remember: you have to use double quotes around the filename expression (NOT single quotes) or it will fail.

### Takeaway

You know how to add a second step that takes the output of the first as input.

### What's next?

Learn how to make the workflow run on many values for the same input.

---

## 9. Modify the workflow to run on many values for the same input

Workflows typically run on batches of inputs that we want to process in bulk. Here we upgrade the workflow to accept an input with multiple values. For simplicity, we go back to hardcoding the greetings instead of using a parameter for the input.

#### 1. Modify the channel to contain multiple greetings (hardcoded for now)

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

#### 2. Modify the first process to generate dynamic filenames so the final filenames will be unique

_Before:_

```groovy title="hello-world.nf" linenums="9"
process sayHello {
    input:
        val greeting

    output:
        path params.output_file

    """
    echo '$greeting' > $params.output_file
    """
}
```

_After:_

```groovy title="hello-world.nf" linenums="9"
process sayHello {
    input:
        val greeting

    output:
        path "${greeting}-${params.output_file}"

    """
    echo '$greeting' > '$greeting-$params.output_file'
    """
}
```

!!! note

    In practice, naming files based on the data input itself is almost always impractical; the better way to generate dynamic filenames is to use a samplesheet and create a map of metadata (aka metamap) from which we can grab an appropriate identifier to generate the filenames. We'll show how to do that later in this training.

#### 3. Run the command and look at the log output

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

Something's wrong! The log lines seem to indicate each process was executed three times (corresponding to the three input elements we provided) but we're only seeing two work directories instead of six.

This is because by default, the ANSI logging system writes the logging from multiple calls to the same process on the same line. Fortunately, we can disable that behavior.

#### 4. Run the command again with the `-ansi-log false` option

```bash
nextflow run hello-world.nf -ansi-log false
```

This time it works fine, we see six work directories in the terminal:

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

You know how to feed an input with multiple elements through a queue channel.

### What's next?

Learn how to make the workflow take a file that contains multiple values for an input.

---

## 10. Modify the workflow to run on a file that contains an input with multiple values

In most cases, when we run on multiple inputs, the input values are contained in a file. Here we're going to use a file where each value is on a new line.

Nextflow offers a variety of predefined operators and functions for reading data in from common file formats and applying text transformations. In this example, we used the [`fromPath()`](https://www.nextflow.io/docs/latest/channel.html#frompath) channel factory with the [`splitText()`](https://www.nextflow.io/docs/latest/operator.html#splittext) operator to read each line as a separate value, then we used a [closure](https://www.nextflow.io/docs/latest/script.html) to apply the [`trim()`](https://www.tutorialspoint.com/java/java_string_trim.htm) method to strip the newline (`\n`) character from each element.

#### 1. Modify the channel declaration to take an input file (through a parameter) instead of hardcoded values

_Before:_

```groovy title="hello-world.nf" linenums="38"
// create a channel for inputs
greeting_ch = Channel.of('Hello','Bonjour','Holà')
```

_After:_

```groovy title="hello-world.nf" linenums="38"
// create a channel for inputs from a file
greeting_ch = Channel.fromPath(params.input_file).splitText() { it.trim() }
```

#### 2. Run the workflow with the `-ansi-log false` option and an `--input_file` parameter

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

    But don't worry if this feels like a lot to grapple with all of a sudden! This is just meant to be a little peek at the kind of things you will learn in later training modules.

### Takeaway

You know how to provide inputs in a file.

### What's next?

Celebrate your success and take a break!

When you are ready, move on to Part 2 of this training to learn how to apply what you've learned to a more realistic data analysis use case.
