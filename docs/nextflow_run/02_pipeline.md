# Part 2: Run pipelines

In Part 1 of this course (Run Basic Operations), we started with an example workflow that had only minimal features in order to keep the code complexity low.
However, most real-world pipelines use more sophisticated features in order to enable efficient processing of large amounts of data at scale, and apply multiple processing steps chained together by sometimes complex logic.

In this part of the training, we demonstrate key features of real-world pipelines through a set of example workflows that build on the original Hello World pipeline.

## 1. Processing multiple inputs

Let's start with the question of how to process not a single greeting at a time, but a batch of greetings, to emulate realistic high-throughout data processing.

The `hello-world-plus.nf` workflow we ran in Part 1 used a command-line parameter to provide a single value at a time, which was passed directly to the process call using `sayHello(params.greeting)`.
That was a deliberately simplified approach that won't work for processing multiple values.

In order to process multiple values (experimental data for multiple samples, for example), we have to upgrade the workflow to use Nextflow's powerful system of **channels** and **operators**.

We've prepared a workflow for you that does exactly that, called `channels.nf`, as well as a CSV file called `greetings.csv` containing some input greetings, emulating the kind of columnar data you might want to process in a real data analysis.

```csv title="greetings.csv" linenums="1"
Hello,English,123
Bonjour,French,456
Holà,Spanish,789
```

(The numbers are not significant, they are just there for illustrative purposes.)

Let's run the workflow first, and we'll take a look at what has changed in the code after.

### 1.1. Run the workflow

Run the following command in your terminal.

```bash
nextflow run channels.nf --greeting greetings.csv
```

This should run without error.

```console title="Output" linenums="1"
 N E X T F L O W   ~  version 25.04.3

Launching `channels.nf` [tiny_heisenberg] DSL2 - revision: 845b471427

executor >  local (3)
[1a/1d19ab] sayHello (2) | 3 of 3 ✔
```

Excitingly, this seems to indicate that '3 of 3' calls were made for the process, which is encouraging!
But this only shows us a single run of the process, with one subdirectory path (`1a/1d19ab`).
What's going on?

By default, the ANSI logging system writes the logging from multiple calls to the same process on the same line.
Fortunately, we can disable that behavior to see the full list of process calls.

### 1.2. Run the command again with the `-ansi-log false` option

To expand the logging to display one line per process call, add `-ansi-log false` to the command.

```bash
nextflow run channels.nf -ansi-log false
```

This time we see all three process runs and their associated work subdirectories listed in the output:

```console title="Output" linenums="1"
N E X T F L O W  ~  version 25.04.3
Launching `channels.nf` [pensive_poitras] DSL2 - revision: 778deadaea
[76/f61695] Submitted process > sayHello (1)
[6e/d12e35] Submitted process > sayHello (3)
[c1/097679] Submitted process > sayHello (2)
```

That's much better; at least for a simple workflow.
For a complex workflow, or a large number of inputs, having the full list output to the terminal might get a bit overwhelming, so you might not choose to use `-ansi-log false` in those cases.

!!! note

    The way the status is reported is a bit different between the two logging modes.
    In the condensed mode, Nextflow reports whether calls were completed successfully or not.
    In this expanded mode, it only reports that they were submitted.

### 1.3. Find the outputs

Ok, so this shows us that the process got run three times.
Let's look for the outputs in the individual work directories first, since we've got them listed.

TODO: show example work directory output

There is an output file there but the name has changed, it's no longer just `output.txt`.
File that away in your brain for later.

Now let's look at the 'results' directory to see if our workflow is still writing a copy of our outputs there.

TODO: show results directory contents

Yes! We see all three expected outputs, conveniently with differentiating names.

### 1.4. Examine the code

Now let's take a look at what has changed in the workflow code.

```groovy title="channels.nf" linenums="1" hl_lines="14,18,25,29-32,35"
#!/usr/bin/env nextflow

/*
 * Use echo to print 'Hello World!' to a file
 */
process sayHello {

    publishDir 'results', mode: 'copy'

    input:
        val greeting

    output:
        path "${greeting}-output.txt"

    script:
    """
    echo '$greeting' > '$greeting-output.txt'
    """
}

/*
 * Pipeline parameters
 */
params.greeting = 'greetings.csv'

workflow {

    // create a channel for inputs from a CSV file
    greeting_ch = Channel.fromPath(params.greeting)
                        .splitCsv()
                        .map { line -> line[0] }

    // emit a greeting
    sayHello(greeting_ch)
}
```

#### 1.4.1. Name the outputs dynamically

Let's start with the output naming since that's conceptually the simplest change.

```groovy title="channels.nf" linenums="13"
    output:
        path "${greeting}-output.txt"

    script:
    """
    echo '$greeting' > '$greeting-output.txt'
    """
```

You see that the output declaration and the relevant bit of the command have changed to include the greeting value in the output file name.
This is one way to ensure that the output file names won't collide when they get published to the common `results` directory.

And that's the only change we've had to make inside the process declaration.

#### 1.4.2. Load the inputs from the CSV

This is the really interesting part: how did we switch from taking a single value from the command-line, to taking a CSV file, parsing it and processing the individual greetings it contains?

That is what Nextflow **channels** are for.
Channels are queues designed to handle inputs efficiently and shuttle them from one step to another in multi-step workflows, while providing built-in parallelism and many additional benefits.
They are complemented by **operators** that allow us to transform channel contents as needed.

Confused? Let's break it down.

```groovy title="channels.nf" linenums="25"
params.greeting = 'greetings.csv'

workflow {

    // create a channel for inputs from a CSV file
    greeting_ch = Channel.fromPath(params.greeting)
                        .splitCsv()
                        .map { line -> line[0] }
```

This is where the magic happens, starting at line 30.
Here's what that line means in plain English:

Channel = create a **channel**, i.e. a queue that will hold the data
.fromPath = from the filepath provided in parenthesis
(params.greeting) = the filepath provided with `--greeting` on the command line

Then the next two lines apply **operators** that transform the contents of the newly created channel as follows:

.splitCsv() = parse the CSV file into an array representing rows and columns
.map { line -> line[0] } = for each row (line), take only the element in the first column

So in practice, starting from the following CSV file:

```csv title="greetings.csv" linenums="1"
Hello,English,123
Bonjour,French,456
Holà,Spanish,789
```

We have transformed that into an array that looks like this:

```txt title="Array contents"
[[Hello,English,123],[Bonjour,French,456],[Holà,Spanish,789]]
```

And then we've taken the first element from each of the three rows and loaded them into a Nextflow channel that now contains: `Hello`, `Bonjour`, and `Holà`.

In other words, the result of this very short snippet of code is a channel called `greeting_ch` loaded with the three individual greetings from the CSV file, ready for processing.

#### 1.4.3. Call the process on each greeting

Then in the last line of the workflow block, we call the `sayHello()` process on the loaded `greeting_ch` channel.

```groovy title="channels.nf" linenums="35"
    sayHello(greeting_ch)
}
```

This tells Nextflow to run the process _individually_ on each element in the channel, i.e. on each greeting.

And because Nextflow is smart like that, it will run these process calls in parallel if possible, depending on the available computing infrastructure.

That is how you can achieve efficient and scalable processing of a lot of data (many samples, or data points, whatever is your unit of research) with comparatively very little code.

### 1.5. Optional: Add `view()` to inspect channel contents

If you're interested in getting into the guts of channels and operators, you can use [`view()`](https://www.nextflow.io/docs/latest/reference/operator.html#view) as described below to inspect the contents of the channel.
You can think of `view()` as a debugging tool, like a `print()` statement in Python, or its equivalent in other languages.

In the workflow block, make the following code change:

```groovy title="channels.nf" linenums="29" hl_lines="3,5,7"
    // create a channel for inputs from a CSV file
    greeting_ch = Channel.fromPath(params.greeting)
                        .view { thing -> "Before splitCsv: $thing" }
                        .splitCsv()
                        .view { thing -> "After splitCsv: $thing" }
                        .map { line -> line[0] }
                        .view { thing -> "After map: $thing" }
```

Here we are using an operator **closure**, denoted by the curly brackets, to specify what to do within the scope of the `view()` operator.
This code will be executed for each item in the channel.
We define a temporary variable for the inner value, here called `thing` to be generic (it could be anything), representing each individual item loaded in a channel.
This variable is only used within the scope of that closure.

You can then run the workflow again:

```bash
nextflow run channels.nf --greeting greetings.csv
```

This should once again run without error and produce the following output:

```console title="Output" linenums="1"
 N E X T F L O W   ~  version 25.04.3

Launching `channels.nf` [tiny_heisenberg] DSL2 - revision: 845b471427

executor >  local (3)
[1a/1d19ab] sayHello (2) | 3 of 3 ✔
Before splitCsv: /workspaces/training/nextflow-run/greetings.csv
After splitCsv: [Hello,English,123]
After splitCsv: [Bonjour,French,456]
After splitCsv: [Holà,Spanish,789]
After map: Hello
After map: Bonjour
After map: Holà
```

This time you see the extra lines at the end showing you what are the contents of the channel at each stage.
Feel free to play around with the contents of the CSV and change the number in the `line -> line[0]` bit that controls which column's value the `map()` operator will pull out.
See what happens!

### Takeaway

You understand at a basic level how channels and operators enable us to process multiple inputs efficiently.

### What's next?

Discover how multi-step workflows are constructed and operate.

---

## 2. Multi-step workflows

Most real-world workflows involve more than one step.
Let's build on what we just learned about channels, and look at how Nextflow uses channels and operators to connect processes together in a multi-step workflow.

To that end, we provide you with an example workflow that chains together three separate steps and demonstrates the following:

1. Making data flow from one process to the next
2. Collecting outputs from multiple process calls into a single process call

Specifically, we made a version of the Hello World workflow that takes each input greeting, converts it to uppercase, then collects all the uppercased greetings into a single output file.

As previously, we'll run the workflow first then look at the code to see what's changed.

### 2.1. Run the workflow

Run the following command in your terminal:

```bash
nextflow run flow.nf --greeting greetings.csv
```

Once again this should run successfully.

```console title="Output" linenums="1"
 N E X T F L O W   ~  version 25.04.3

Launching `flow.nf` [soggy_franklin] DSL2 - revision: bc8e1b2726

[d6/cdf466] sayHello (1)       | 3 of 3 ✔
[99/79394f] convertToUpper (2) | 3 of 3 ✔
[1e/83586c] collectGreetings   | 1 of 1 ✔
There were 3 greetings in this batch
```

You see that as promised, multiple steps were run as part of the workflow; the first two (`sayHello` and `convertToUpper`) were presumably run on each individual greeting, and the third (`collectGreetings`) will have been run only once, on the outputs of all three of the `convertToUpper` calls.

### 2.2. Find the outputs

If you'd like to verify that that is in fact what happened (good scientist; have a biscuit), you can take a look in the `results` directory.

```console title="Directory contents"
results
├── Bonjour-output.txt
├── COLLECTED-output.txt
├── COLLECTED-test-batch-output.txt
├── COLLECTED-trio-output.txt
├── Hello-output.txt
├── Holà-output.txt
├── UPPER-Bonjour-output.txt
├── UPPER-Hello-output.txt
└── UPPER-Holà-output.txt
```

Look at the file names and check their contents to confirm that they are what you expect, for example:

```console title="bash"
cat results/COLLECTED-trio-output.txt
```

```console title="Output"
HELLO
BONJOUR
HOLà
```

That is the expected final result of our multi-step pipeline.

### 2.3. Examine the code

Let's look at the code and see what we can tie back to what we just observed.

```groovy title="channels.nf" linenums="1"
#!/usr/bin/env nextflow

/*
 * Use echo to print 'Hello World!' to a file
 */
process sayHello {

    publishDir 'results', mode: 'copy'

    input:
        val greeting

    output:
        path "${greeting}-output.txt"

    script:
    """
    echo '$greeting' > '$greeting-output.txt'
    """
}

/*
 * Use a text replacement tool to convert the greeting to uppercase
 */
process convertToUpper {

    publishDir 'results', mode: 'copy'

    input:
        path input_file

    output:
        path "UPPER-${input_file}"

    script:
    """
    cat '$input_file' | tr '[a-z]' '[A-Z]' > 'UPPER-${input_file}'
    """
}

/*
 * Collect uppercase greetings into a single output file
 */
process collectGreetings {

    publishDir 'results', mode: 'copy'

    input:
        path input_files
        val batch_name

    output:
        path "COLLECTED-${batch_name}-output.txt" , emit: outfile
        val count_greetings , emit: count

    script:
        count_greetings = input_files.size()
    """
    cat ${input_files} > 'COLLECTED-${batch_name}-output.txt'
    """
}

/*
 * Pipeline parameters
 */
params.greeting = 'greetings.csv'
params.batch = 'test-batch'

workflow {

    // create a channel for inputs from a CSV file
    greeting_ch = Channel.fromPath(params.greeting)
                        .splitCsv()
                        .map { line -> line[0] }

    // emit a greeting
    sayHello(greeting_ch)

    // convert the greeting to uppercase
    convertToUpper(sayHello.out)

    // collect all the greetings into one file
    collectGreetings(convertToUpper.out.collect(), params.batch)

    // emit a message about the size of the batch
    collectGreetings.out.count.view { num -> "There were $num greetings in this batch" }
}
```

The most obvious difference compared to the previous version of the workflow is that now there are multiple process definitions, and correspondingly, several process calls in the workflow block.

#### 2.3.1. Multiple process definitions

In addition to the original `sayHello` process, we now also have `convertToUpper` and `collectGreetings`, which match the names of the processes we saw in the console output.

All three are structured in the same way and follow roughly the same logic, though you may notice that the `collectGreetings` process takes two inputs and outputs two outputs.
We won't go into that in detail, but it shows how a process can be given additional parameters and emit multiple outputs.

#### 2.3.2. Processes chained via channels

The really interesting thing to look at here is how the process calls are chained together in the workflow block.

```groovy title="channels.nf" linenums="69"
workflow {

    // create a channel for inputs from a CSV file
    greeting_ch = Channel.fromPath(params.greeting)
                        .splitCsv()
                        .map { line -> line[0] }

    // emit a greeting
    sayHello(greeting_ch)

    // convert the greeting to uppercase
    convertToUpper(sayHello.out)

    // collect all the greetings into one file
    collectGreetings(convertToUpper.out.collect(), params.batch)

    // emit a message about the size of the batch
    collectGreetings.out.count.view { num -> "There were $num greetings in this batch" }
}
```

You can see that the first process call, to `sayHello()`, is unchanged.

Then the next process call, to `convertToUpper`, _refers_ to the output of `sayHello` as `sayHello.out`:

```groovy title="channels.nf" linenums="79"
    // convert the greeting to uppercase
    convertToUpper(sayHello.out)
```

This means 'call `convertToUpper` on the output of `sayHello()`'.

Then the next call is doing the same thing, with a little twist (or two):

```groovy title="channels.nf" linenums="82"
    // collect all the greetings into one file
    collectGreetings(convertToUpper.out.collect(), params.batch)
```

First, you'll note this one has two inputs provided to the `collectGreetings()` call: `convertToUpper.out.collect()` and `params.batch`.
The latter is just a parameter value that the process expects (in second position because it is declared in second position in the process definition).

The other one, `convertToUpper.out.collect()`, is a bit more complicated and deserves its own discussion.

#### 2.3.3. Operators provide plumbing options

What we're seeing in `convertToUpper.out.collect()` is the use of another operator, called `collect()`.
This operator is used to collect the outputs from multiple parallel calls to the same process and package them into a single channel element.

Specifically,
TODO: finish explanation

There are many other operators available to apply transformations to the contents of channels between process calls.

This gives pipeline developers a lot of flexibility for customizing the flow logic of their pipeline.
The downside is that it can sometimes make it harder to decipher what the pipeline is doing.

### 2.4. Use the graph preview

One very helpful tool for understanding what a pipeline does, if it's not adequately documented, is the graph preview functionality available in VSCode. You can see this in the training environment by clicking on the small `DAG preview` link displayed just above the workflow block in any Nextflow script.

TODO: add picture

This does not show operators, but it does give a useful representation of how process calls are connected and what are their inputs.

### Takeaway

You understand at a basic level how multi-step workflows are constructed and operate, using channels and operators, and you can manage their execution.

### What's next?

Learn how Nextflow pipelines are often modularized to promote code reuse and maintainability.

---

## 3. Modular code components

So far, all the workflows we've looked at have consisted of one single workflow file containing all the relevant code.

However, real-world pipelines typically benefit from being _modularized_, meaning that the code is split into different files.
This can make their development and maintenance more efficient and sustainable.

Here we are going to demonstrate the most common form of code modularity in Nextflow, which is the use of **modules**.

In Nextflow, a **module** is a single process definition that is encapsulated by itself in a standalone code file.
To use a module in a workflow, you just add a single-line import statement to your workflow code file; then you can integrate the process into the workflow the same way you normally would.

Putting processes into individual modules makes it possible to reuse process definitions in multiple workflows without producing multiple copies of the code.
This makes the code more shareable, flexible and maintainable.

We have of course once again prepared a suitable workflow for demonstration purposes, called `modular.nf`, along with a set of modules located in the `modules/` directory.

### 3.1. Examine the code

This time we're going to look at the code first.

TODO: show directory contents

Importantly, the processes and workflow logic are exactly the same as in the previous version of the workflow. However the process code is in the modules instead of being in the main workflow file, and there are now import statements in the workflow file telling Nextflow to pull them in at runtime.

```groovy title="hello-modules.nf" linenums="9" hl_lines="4"
// Include modules
include { sayHello } from './modules/sayHello.nf'
include { convertToUpper } from './modules/convertToUpper.nf'
include { collectGreetings } from './modules/collectGreetings.nf'

workflow {
```

You can look inside one of the modules to satisfy yourself that the process definition is unchanged, it's literally just been copy-pasted into a standalone file.

TODO: show module code for `sayHello`

So let's see what it looks like to run this new version.

### 3.2. Run the workflow

Run this command in your terminal, with the `-resume` flag:

```bash
nextflow run modular.nf --greeting greetings.csv -resume
```

Once again this should run successfully.

```console title="Output" linenums="1"
 N E X T F L O W   ~  version 25.04.3

Launching `modular.nf` [soggy_franklin] DSL2 - revision: bc8e1b2726

[j6/cdfa66] sayHello (1)       | 3 of 3, cached: ✔
[95/79484f] convertToUpper (2) | 3 of 3, cached: ✔
[5e/4358gc] collectGreetings   | 1 of 1, cached: ✔
There were 3 greetings in this batch
```

You'll notice that these all cached successfully, meaning that Nextflow recognized that it has already done the requested work, even though the code has been split up and the main workflow file has been renamed.

None of that matters to Nextflow; what matters is the job script that is generated once all the code has been pulled together and evaluated.

!!!note

    It is also possible to encapsulate a section of a workflow as a 'subworkflow' that can be imported into a larger pipeline, but that is outside the scope of this course.

    TODO: add links to learn more about composable workflows

### Takeaway

You know how processes can be stored in standalone modules to promote code reuse and improve maintainability.

### What's next?

Learn to use containers for managing software dependencies.

---

## 4. Using containerized software

So far the workflows we've been using as examples just needed to run very basic text procession operations using UNIX tools available in our environment.

However, real-world pipelines typically require specialized tools and packages that are not included by default in most environments.
Usually, you'd need to install these tools, manage their dependencies, and resolve any conflicts.

That is all very tedious and annoying.
A much better way to address this problem is to use **containers**.

A **container** is a lightweight, standalone, executable unit of software created from a container **image** that includes everything needed to run an application including code, system libraries and settings.

!!! note

    We teach this using the technology [Docker](https://www.docker.com/get-started/), but Nextflow supports [several other container technologies](https://www.nextflow.io/docs/latest/container.html#) as well.

### 4.1. Use a container directly

First, let's try interacting with a container directly.
This will help solidify your understanding of what containers are before we start using them in Nextflow.

TODO: clone the content from hello_containers.md

### 4.2. Use a container in a workflow

TODO: clone the content from hello_containers.md

### Takeaway

You understand what role containers play in managing software tool versions and ensuring reproducibility.

More generally, you have a basic understanding of the most common and most important components of real-world Nexflow pipelines.

### What's next?

Take another break!
TODO: finalize the transition text
