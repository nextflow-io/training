# Part 2: Run real pipelines

In Part 1 of this course (Run Basic Operations), we started with an example workflow that had only minimal features in order to keep the code complexity low.
For example, `1-hello.nf` used a command-line parameter (`--input`) to provide a single value at a time.

However, most real-world pipelines use more sophisticated features in order to enable efficient processing of large amounts of data at scale, and apply multiple processing steps chained together by sometimes complex logic.

In this part of the training, we demonstrate key features of real-world pipelines by trying out expanded versions of the original Hello World pipeline.

## 1. Processing input data from a file

In a real-world pipeline, we typically want to process multiple data points (or data series) contained in one or more input files.
And wherever possible, we want to run the processing of independent data in parallel, to shorten the time spent waiting for analysis.

To demonstrate how Nextflow does this, we've prepared a a CSV file called `greetings.csv` that contains several input greetings, mimicking the kind of columnar data you might want to process in a real data analysis.
Note that the numbers are not meaningful, they are just there for illustrative purposes.

```csv title="data/greetings.csv" linenums="1"
Hello,English,123
Bonjour,French,456
Holà,Spanish,789
```

We've also written an improved version of the original workflow, now called `2a-inputs.nf`, that will read in the CSV file, extract the greetings and write each of them to a separate file.

<figure class="excalidraw">
--8<-- "docs/nextflow_run/img/hello-pipeline-multi-inputs.svg"
</figure>

Let's run the workflow first, and we'll take a look at the relevant Nextflow code afterward.

### 1.1. Run the workflow

Run the following command in your terminal.

```bash
nextflow run 2a-inputs.nf --input data/greetings.csv
```

??? success "Command output"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `2a-inputs.nf` [mighty_sammet] DSL2 - revision: 29fb5352b3

    executor >  local (3)
    [8e/0eb066] sayHello (2) [100%] 3 of 3 ✔
    ```

Excitingly, this seems to indicate that '3 of 3' calls were made for the process, which is encouraging, since there were three rows of data in the CSV we provided as input.
This suggests the `sayHello()` process was called three times, once on each input row.

### 1.2. Find the published outputs in the `results` directory

Let's look at the 'results' directory to see if our workflow is still writing a copy of our outputs there.

??? abstract "Directory contents"

    ```console linenums="1" hl_lines="4-7"
    results
    ├── 1-hello
    |   └── output.txt
    └── 2a-inputs
        ├── Bonjour-output.txt
        ├── Hello-output.txt
        └── Holà-output.txt
    ```

Yes! We see a new directory called `2a-inputs` with three output files with different names, conveniently enough.

You can open each of them to satisfy yourself that they contain the appropriate greeting string.

??? abstract "File contents"

    ```console title="results/2a-inputs/Hello-output.txt"
    Hello
    ```

    ```console title="results/2a-inputs/Bonjour-output.txt"
    Bonjour
    ```

    ```console title="results/2a-inputs/Holà-output.txt"
    Holà
    ```

This confirms each greeting in the input file has been processed appropriately.

### 1.3. Find the original outputs and logs

You may have noticed that the console output above referred to only one task directory.
Does that mean all three calls to `sayHello()` were executed within that one task directory?

#### 1.3.1. Examine the task directory given in the terminal

Let's have a look inside that `8e/0eb066` task directory.

??? abstract "Directory contents"

    ```console title="8e/0eb066"
    work/8e/0eb066071cdb4123906b7b4ea8b047/
    └── Bonjour-output.txt
    ```

We only find the output corresponding to one of the greetings (as well as the accessory files if we enable display of hidden files).

So what's going on here?

By default, the ANSI logging system writes the status information for all calls to the same process on the same line.
As a result, it only showed us one of the three task directory paths (`8e/0eb066`) in the console output.
There are two others that are not listed there.

#### 1.3.2. Make the terminal show more details

We can modify the logging behavior to see the full list of process calls by adding the `-ansi-log false` to the command as follows:

```bash
nextflow run 2a-inputs.nf --input data/greetings.csv -ansi-log false
```

??? success "Command output"

    ```console linenums="1"
    N E X T F L O W  ~  version 25.10.2
    Launching `2a-inputs.nf` [pedantic_hamilton] DSL2 - revision: 6bbc42e49f
    [ab/1a8ece] Submitted process > sayHello (1)
    [0d/2cae24] Submitted process > sayHello (2)
    [b5/0df1d6] Submitted process > sayHello (3)
    ```

This time we see all three process runs and their associated work subdirectories listed in the output.

Notice that the way the status is reported is a bit different between the two logging modes.
In the condensed mode, Nextflow reports whether calls were completed successfully or not.
In this expanded mode, it only reports that they were submitted.

This confirms that the `sayHello()` process gets called three times, and a separate task directory is created for each one.

If we look inside each of the task directories listed there, we can verify that each one corresponds to one of the greetings.

??? abstract "Directory contents"

    ```console title="ab/1a8ece"
    work/ab/1a8ece307e53f03fce689dde904b64/
    └── Hello-output.txt
    ```

    ```console title="0d/2cae24"
    work/0d/2cae2481a53593bc607077c80c9466/
    └── Bonjour-output.txt
    ```

    ```console title="b5/0df1d6"
    work/b5/0df1d642353269909c2ce23fc2a8fa/
    └── Holà-output.txt
    ```

This confirms that each process call is executed in isolation from all the others.
That has many advantages, including avoiding collisions if the process produces any intermediate files with non-unique names.

!!! tip

    For a complex workflow, or a large number of inputs, having the full list output to the terminal might get a bit overwhelming, so people don't normally use `-ansi-log false` in routine usage.

### 1.4. Examine the workflow code

So this version of the workflow is capable of reading in a CSV file of inputs, processing the inputs separately, and naming the outputs uniquely.

Let's take a look at what makes that possible in the workflow code.

??? full-code "Full code file"

    ```groovy title="2a-inputs.nf" linenums="1" hl_lines="31-33 35"
    #!/usr/bin/env nextflow

    /*
    * Use echo to print 'Hello World!' to a file
    */
    process sayHello {

        input:
        val greeting

        output:
        path "${greeting}-output.txt"

        script:
        """
        echo '${greeting}' > '${greeting}-output.txt'
        """
    }

    /*
    * Pipeline parameters
    */
    params {
        input: Path
    }

    workflow {

        main:
        // create a channel for inputs from a CSV file
        greeting_ch = channel.fromPath(params.input)
                            .splitCsv()
                            .map { line -> line[0] }
        // emit a greeting
        sayHello(greeting_ch)

        publish:
        first_output = sayHello.out
    }

    output {
        first_output {
            path '2a-inputs'
            mode 'copy'
        }
    }
    ```

Once again, you don't need to memorize code syntax, but it's good to learn to recognize key components of the workflow that provide important functionality.

#### 1.4.1. Loading the input data from the CSV

This is the most interesting part: how did we switch from taking a single value from the command-line, to taking a CSV file, parsing it and processing the individual greetings it contains?

In Nextflow, we do that with a **channel**: a construct designed to handle inputs efficiently and shuttle them from one step to another in multi-step workflows, while providing built-in parallelism and many additional benefits.

Let's break it down.

```groovy title="2a-inputs.nf" linenums="29" hl_lines="3-5"
    main:
    // create a channel for inputs from a CSV file
    greeting_ch = channel.fromPath(params.input)
                        .splitCsv()
                        .map { line -> line[0] }
    // emit a greeting
    sayHello(greeting_ch)
```

This is where the magic happens, starting at line 31.
Here's what that line means in plain English:

- `channel` creates a **channel**, _i.e._ a queue that will hold the data
- `.fromPath` specifies the data source is a filepath
- `(params.input)` specifies the filepath is provided by `--input` on the command line

In other words, that line tells Nextflow: take the filepath given with `--input` and get ready to treat its contents as input data.

Then the next two lines apply **operators** that do the actual parsing of the file and loading of the data into the appropriate data structure:

- `.splitCsv()` tells Nextflow to parse the CSV file into an array representing rows and columns
- `.map { line -> line[0] }` tells Nextflow to take only the element in the first column from each row

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

The result of this very short snippet of code is a channel called `greeting_ch` loaded with the three individual greetings from the CSV file, ready for processing.

#### 1.4.2. Call the process on each greeting

Next, in the last line of the workflow's `main:` block, we provide the loaded `greeting_ch` channel as input to the `sayHello()` process.

```groovy title="2a-inputs.nf" linenums="29" hl_lines="7"
    main:
    // create a channel for inputs from a CSV file
    greeting_ch = channel.fromPath(params.input)
                        .splitCsv()
                        .map { line -> line[0] }
    // emit a greeting
    sayHello(greeting_ch)
```

This tells Nextflow to run the process individually on each element in the channel, _i.e._ on each greeting.
And because Nextflow is smart like that, it will run these process calls in parallel if possible, depending on the available computing infrastructure.

That is how you can achieve efficient and scalable processing of a lot of data (many samples, or data points, whatever is your unit of research) with comparatively very little code.

#### 1.4.3. How the outputs are named

Finally, it's worth taking a quick look at the process code to see how we get the output files to be named uniquely.

```groovy title="2a-inputs.nf" linenums="6" hl_lines="7 11"
process sayHello {

    input:
    val greeting

    output:
    path "${greeting}-output.txt"

    script:
    """
    echo '${greeting}' > '${greeting}-output.txt'
    """
}
```

You see that, compared to the version of this process in `1-hello.nf`, the output declaration and the relevant bit of the command have changed to include the greeting value in the output file name.

This is one way to ensure that the output file names won't collide when they get published to the common results directory.

And that's the only change we've had to make inside the process declaration!

### Takeaway

You understand at a basic level how channels and operators enable us to process multiple inputs efficiently.

### What's next?

Discover how multi-step workflows are constructed and how they operate.

---

## 2. Running multi-step workflows

Most real-world workflows involve more than one step.
Let's build on what we just learned about channels, and look at how Nextflow uses channels and operators to connect processes together in a multi-step workflow.

To that end, we provide you with an example workflow that chains together three separate steps and demonstrates the following:

1. Making data flow from one process to the next
2. Collecting outputs from multiple process calls into a single process call

Specifically, we made an expanded version of the workflow called `2b-multistep.nf` that takes each input greeting, converts it to uppercase, then collects all the uppercased greetings into a single output file.

<figure class="excalidraw">
--8<-- "docs/nextflow_run/img/hello-pipeline-multi-steps.svg"
</figure>

As previously, we'll run the workflow first then look at the code to see what is new.

### 2.1. Run the workflow

Run the following command in your terminal:

```bash
nextflow run 2b-multistep.nf --input data/greetings.csv
```

??? success "Command output"

    ```console linenums="1"
    N E X T F L O W   ~  version 25.10.2

    Launching `2b-multistep.nf` [soggy_franklin] DSL2 - revision: bc8e1b2726

    [d6/cdf466] sayHello (1)       | 3 of 3 ✔
    [99/79394f] convertToUpper (2) | 3 of 3 ✔
    [1e/83586c] collectGreetings   | 1 of 1 ✔
    ```

You see that as promised, multiple steps were run as part of the workflow; the first two (`sayHello` and `convertToUpper`) were presumably run on each individual greeting, and the third (`collectGreetings`) will have been run only once, on the outputs of all three of the `convertToUpper` calls.

### 2.2. Find the outputs

Let's verify that that is in fact what happened by taking a look in the `results` directory.

??? abstract "Directory contents"

    ```console linenums="1" hl_lines="8-16"
    results
    ├── 1-hello
    |   └── output.txt
    ├── 2a-inputs
    |   ├── Bonjour-output.txt
    |   ├── Hello-output.txt
    |   └── Holà-output.txt
    └── 2b-multistep
        ├── COLLECTED-batch-output.txt
        ├── batch-report.txt
        └── intermediates
            ├── Bonjour-output.txt
            ├── Hello-output.txt
            ├── Holà-output.txt
            ├── UPPER-Bonjour-output.txt
            ├── UPPER-Hello-output.txt
            └── UPPER-Holà-output.txt

    ```

As you can see, we have a new directory called `2b-multistep`, and it contains quite a few more files than before.
Some of the files have been grouped into a subdirectory called `intermediates`, while two files are located at the top level.

Those two are the final results of the multi-step workflow.
Take a minute to look at the file names and check their contents to confirm that they are what you expect.

??? abstract "File contents"

    ```txt title="results/2b-multistep/COLLECTED-batch-output.txt"
    HELLO
    BONJOUR
    HOLà
    ```

    ```txt title="results/2b-multistep/batch-report.txt"
    There were 3 greetings in this batch.
    ```

The first contains our three greetings, uppercased and collected back into a single file as promised.
The second is a report file that summarizes some information about the run.

### 2.3. Examine the code

Let's look at the code and see what we can tie back to what we just observed.

??? full-code "Full code file"

    ```groovy title="2b-multistep.nf" linenums="1" hl_lines="63 75-78 82-84"
    #!/usr/bin/env nextflow

    /*
    * Use echo to print 'Hello World!' to a file
    */
    process sayHello {

        input:
        val greeting

        output:
        path "${greeting}-output.txt"

        script:
        """
        echo '${greeting}' > '${greeting}-output.txt'
        """
    }

    /*
    * Use a text replacement tool to convert the greeting to uppercase
    */
    process convertToUpper {

        input:
        path input_file

        output:
        path "UPPER-${input_file}"

        script:
        """
        cat '${input_file}' | tr '[a-z]' '[A-Z]' > 'UPPER-${input_file}'
        """
    }

    /*
    * Collect uppercase greetings into a single output file
    */
    process collectGreetings {

        input:
        path input_files
        val batch_name

        output:
        path "COLLECTED-${batch_name}-output.txt", emit: outfile
        path "${batch_name}-report.txt", emit: report

        script:
        count_greetings = input_files.size()
        """
        cat ${input_files} > 'COLLECTED-${batch_name}-output.txt'
        echo 'There were ${count_greetings} greetings in this batch.' > '${batch_name}-report.txt'
        """
    }

    /*
    * Pipeline parameters
    */
    params {
        input: Path
        batch: String = 'batch'
    }

    workflow {

        main:
        // create a channel for inputs from a CSV file
        greeting_ch = channel.fromPath(params.input)
                            .splitCsv()
                            .map { line -> line[0] }
        // emit a greeting
        sayHello(greeting_ch)
        // convert the greeting to uppercase
        convertToUpper(sayHello.out)
        // collect all the greetings into one file
        collectGreetings(convertToUpper.out.collect(), params.batch)

        publish:
        first_output = sayHello.out
        uppercased = convertToUpper.out
        collected = collectGreetings.out.outfile
        batch_report = collectGreetings.out.report
    }

    output {
        first_output {
            path '2b-multistep/intermediates'
            mode 'copy'
        }
        uppercased {
            path '2b-multistep/intermediates'
            mode 'copy'
        }
        collected {
            path '2b-multistep'
            mode 'copy'
        }
        batch_report {
            path '2b-multistep'
            mode 'copy'
        }
    }
    ```

There's a lot going on in there, but the most obvious difference compared to the previous version of the workflow is that now there are multiple process definitions, and correspondingly, several process calls in the workflow block.

Let's take a closer look and see if we can identify the most interesting pieces.

#### 2.3.1. Overall wiring of the workflow

If you're using VSCode with the Nextflow extension, you can get a helpful diagram of how the processes are connected by clicking on the small `DAG preview` link displayed just above the workflow block in any Nextflow script.

<figure class="excalidraw">
--8<-- "docs/nextflow_run/img/DAG-multistep.svg"
</figure>

This gives you a nice overview of how the processes are connected and what they produce.

You see that in addition to the original `sayHello` process, we now also have `convertToUpper` and `collectGreetings`, which match the names of the processes we saw in the console output.
The two new process definitions are structured in the same way as the `sayHello` process, except `collectGreetings` takes an additional input parameter called `batch` and produces two outputs.

We won't go into the code for each in detail, but if you're curious, you can look up the details in [Part 2 of Hello Nextflow](../hello_nextflow/03_hello_workflow.md).

For now, let's dig into how the processes are connected to one another.

#### 2.3.2. How the processes are connected

The really interesting thing to look at here is how the process calls are chained together in the workflow's `main:` block.

```groovy title="2b-multistep.nf" linenums="68" hl_lines="9 11"
    main:
    // create a channel for inputs from a CSV file
    greeting_ch = channel.fromPath(params.input)
                        .splitCsv()
                        .map { line -> line[0] }
    // emit a greeting
    sayHello(greeting_ch)
    // convert the greeting to uppercase
    convertToUpper(sayHello.out)
    // collect all the greetings into one file
    collectGreetings(convertToUpper.out.collect(), params.batch)
```

You can see that the first process call, `sayHello(greeting_ch)`, is unchanged.

Then the next process call, to `convertToUpper`, refers to the output of `sayHello` as `sayHello.out`, just like the `publish:` block did.

```groovy title="2b-multistep.nf" linenums="75"
    // convert the greeting to uppercase
    convertToUpper(sayHello.out)
```

This tells Nextflow to provide `sayHello.out`, which represents the channel output by `sayHello()`, as an input to `convertToUpper`.

That is, at its simplest, how we shuttle data from one step to the next in Nextflow.
We take the output channel from the first process, and pass it as an input to the next process.

#### 2.3.3. A process can take multiple inputs

The third process call, to `collectGreetings`, is a little different.

```groovy title="2b-multistep.nf" linenums="77"
    // collect all the greetings into one file
    collectGreetings(convertToUpper.out.collect(), params.batch)
```

You see this call is given two inputs, `convertToUpper.out.collect()` and `params.batch`.
Ignoring the `.collect()` bit for now, we can generalize this as `collectGreetings(input1, input2)`.

That matches the two input declarations in the process module:

```groovy title="2b-multistep.nf" linenums="40"
process collectGreetings {

    input:
    path input_files
    val batch_name
```

When Nextflow parses this, it will assign the first input in the call to `path input_files`, and the second to `val batch_name`.

So now you know a process can take multiple inputs, and what the call looks like in the workflow block.

Now let's take a closer look at that first input, `convertToUpper.out.collect()`.

#### 2.3.4. What `collect()` does in the `collectGreetings` call

To pass the output of `sayHello` to `convertToUpper`, we simply referred to the output channel of `sayHello` as `sayHello.out`. But for the next step, we're seeing a reference to `convertToUpper.out.collect()`.

What is this `collect()` bit and what does it do?

It's an operator, of course. Just like the `splitCsv` and `map` operators we encountered earlier.
This time the operator is called `collect`, and is applied to the output channel produced by `convertToUpper`.

The `collect` operator is used to collect the outputs from multiple calls to the same process and package them into a single channel element.

In the context of this workflow, it's taking the three uppercased greetings in the `convertToUpper.out` channel --which are three separate channel items, and would normally be handled in separate calls by the next process-- and packaging them into a single item.

In more practical terms: if we didn't apply `collect()` to the output of `convertToUpper()` before feeding it to `collectGreetings()`, Nextflow would simply run `collectGreetings()` independently on each greeting, which would not achieve our goal.

<figure class="excalidraw">
--8<-- "docs/nextflow_run/img/without-collect-operator.svg"
</figure>

In contrast, using `collect()` allows us to take all the separate uppercased greetings produced by the second step of the workflow and feed them all together to a single call in the third step of the pipeline.

<figure class="excalidraw">
--8<-- "docs/nextflow_run/img/with-collect-operator.svg"
</figure>

That's how we get all the greetings back into the same file.

There are many other [operators](https://www.nextflow.io/docs/latest/reference/operator.html#operator-page) available to apply transformations to the contents of channels between process calls.

This gives pipeline developers a lot of flexibility for customizing the flow logic of their pipeline.
The downside is that it can sometimes make it harder to decipher what the pipeline is doing.

#### 2.3.5. An input parameter can have a default value

Now that we understand what's going on with the first input to `collectGreetings`, let's look at that `params.batch` input.

```groovy title="2b-multistep.nf" linenums="77"
    // collect all the greetings into one file
    collectGreetings(convertToUpper.out.collect(), params.batch)
```

That automatically causes the workflow to expect a CLI parameter named `--batch` to be included in the command line used to run it.
However, you may recall that when we launched the workflow, we didn't specify a `batch` parameter in the command.

What's going on there?
Have a look at the `params` block.

```groovy title="2b-multistep.nf" linenums="58" hl_lines="6"
    /*
    * Pipeline parameters
    */
    params {
        input: Path
        batch: String = 'batch'
    }
```

Ah! There is a default value configured in the workflow, so we don't have to provide it.
But if we do provide one on the command line, the value we specify will be used instead of the default.

You can try it yourself if you want.

```bash
nextflow run 2b-multistep.nf --input data/greetings.csv --batch test
```

??? success "Command output"

    ```console linenums="1"
    N E X T F L O W   ~  version 25.10.2

    Launching `2b-multistep.nf` [soggy_franklin] DSL2 - revision: bc8e1b2726

    [a5/cdff26] sayHello (1)       | 3 of 3 ✔
    [c5/78794f] convertToUpper (2) | 3 of 3 ✔
    [d3/b4d86c] collectGreetings   | 1 of 1 ✔
    ```

You should see new final outputs named accordingly.

??? abstract "Directory contents"

    ```console linenums="1" hl_lines="10 12"
    results
    ├── 1-hello
    |   └── output.txt
    ├── 2a-inputs
    |   ├── Bonjour-output.txt
    |   ├── Hello-output.txt
    |   └── Holà-output.txt
    └── 2b-multistep
        ├── COLLECTED-batch-output.txt
        ├── COLLECTED-test-output.txt
        ├── batch-report.txt
        ├── test-report.txt
        └── intermediates
            ├── Bonjour-output.txt
            ├── Hello-output.txt
            ├── Holà-output.txt
            ├── UPPER-Bonjour-output.txt
            ├── UPPER-Hello-output.txt
            └── UPPER-Holà-output.txt
    ```

This is an aspect of input configuration, which we'll cover in more detail in Part 3, but for now the important thing is to know that input parameters can be given default values.

#### 2.3.6. A process can produce multiple outputs

In the `collectGreetings` process definition, we see the following output declarations:

```groovy title="2b-multistep.nf" linenums="46"
    output:
    path "COLLECTED-${batch_name}-output.txt", emit: outfile
    path "${batch_name}-report.txt", emit: report
```

Which are then referred to by the name given with `emit:` in the `publish:` block:

```groovy title="2b-multistep.nf" linenums="80" hl_lines="4 5"
    publish:
    first_output = sayHello.out
    uppercased = convertToUpper.out
    collected = collectGreetings.out.outfile
    batch_report = collectGreetings.out.report
```

This makes it easy to then pass specific outputs individually to other processes in the workflow, in combination with various operators.

#### 2.3.7. Published outputs can be organized

In the `output` block, we've used custom paths to group intermediate results in order to make it easier to pick out just the final outputs of the workflow.

```groovy title="2b-multistep.nf" linenums="80" hl_lines="3 7 11 15"
output {
    first_output {
        path '2b-multistep/intermediates'
        mode 'copy'
    }
    uppercased {
        path '2b-multistep/intermediates'
        mode 'copy'
    }
    collected {
        path '2b-multistep'
        mode 'copy'
    }
    batch_report {
        path '2b-multistep'
        mode 'copy'
    }
}
```

There are more sophisticated ways to organize published outputs; we'll touch on a few in the part on configuration.

### Takeaway

You understand at a basic level how multi-step workflows are constructed using channels and operators and how they operate.
You've also seen that processes can take multiple inputs and produce multiple outputs, and that these can be published in a structured way.

### What's next?

Learn how Nextflow pipelines can be modularized to promote code reuse and maintainability.

---

## 3. Running modularized pipelines

So far, all the workflows we've looked at have consisted of one single workflow file containing all the relevant code.

However, real-world pipelines typically benefit from being _modularized_, meaning that the code is split into different files.
This can make their development and maintenance more efficient and sustainable.

Here we are going to demonstrate the most common form of code modularity in Nextflow, which is the use of **modules**.

In Nextflow, a **module** is a single process definition that is encapsulated by itself in a standalone code file.
To use a module in a workflow, you just add a single-line import statement to your workflow code file; then you can integrate the process into the workflow the same way you normally would.
That makes it possible to reuse process definitions in multiple workflows without producing multiple copies of the code.

Until now we've been running workflows that had all their processes included in a monolithic code file.
Now we're going to see what it looks like when the processes are stored in individual modules.

We have of course once again prepared a suitable workflow for demonstration purposes, called `2c-modules.nf`, along with a set of modules located in the `modules/` directory.

<figure class="excalidraw">
--8<-- "docs/nextflow_run/img/modules.svg"
</figure>

??? abstract "Directory contents"

    ```console
    modules/
    ├── collectGreetings.nf
    ├── convertToUpper.nf
    ├── cowpy.nf
    └── sayHello.nf
    ```

You see there are four Nextflow files, each named after one of the processes.
You can ignore the `cowpy.nf` file for now; we'll get to that one later.

### 3.1. Examine the code

This time we're going to look at the code first.
Start by opening the `2c-modules.nf` workflow file.

??? full-code "Full code file"

    ```groovy title="2c-modules.nf" linenums="1"
    #!/usr/bin/env nextflow

    // Include modules
    include { sayHello } from './modules/sayHello.nf'
    include { convertToUpper } from './modules/convertToUpper.nf'
    include { collectGreetings } from './modules/collectGreetings.nf'

    /*
    * Pipeline parameters
    */
    params {
        input: Path
        batch: String = 'batch'
    }

    workflow {

        main:
        // create a channel for inputs from a CSV file
        greeting_ch = channel.fromPath(params.input)
                            .splitCsv()
                            .map { line -> line[0] }
        // emit a greeting
        sayHello(greeting_ch)
        // convert the greeting to uppercase
        convertToUpper(sayHello.out)
        // collect all the greetings into one file
        collectGreetings(convertToUpper.out.collect(), params.batch)

        publish:
        first_output = sayHello.out
        uppercased = convertToUpper.out
        collected = collectGreetings.out.outfile
        batch_report = collectGreetings.out.report
    }

    output {
        first_output {
            path '2c-modules/intermediates'
            mode 'copy'
        }
        uppercased {
            path '2c-modules/intermediates'
            mode 'copy'
        }
        collected {
            path '2c-modules'
            mode 'copy'
        }
        batch_report {
            path '2c-modules'
            mode 'copy'
        }
    }
    ```

You see that the workflow logic is exactly the same as in the previous version of the workflow.
However, the process code is gone from the workflow file, and instead there are `include` statements pointing to separate files under `modules`.

```groovy title="hello-modules.nf" linenums="3"
// Include modules
include { sayHello } from './modules/sayHello.nf'
include { convertToUpper } from './modules/convertToUpper.nf'
include { collectGreetings } from './modules/collectGreetings.nf'
```

Open up one of those files and you'll find the code for the corresponding process.

??? full-code "Full code file"

    ```groovy title="modules/sayHello.nf" linenums="1"
    #!/usr/bin/env nextflow

    /*
    * Use echo to print 'Hello World!' to a file
    */
    process sayHello {

        input:
        val greeting

        output:
        path "${greeting}-output.txt"

        script:
        """
        echo '${greeting}' > '${greeting}-output.txt'
        """
    }
    ```

As you can see, the process code has not changed; it's just been copied into an individual module file instead of being in the main workflow file.
The same applies to the other two processes.

So let's see what it looks like to run this new version.

### 3.2. Run the workflow

Run this command in your terminal, with the `-resume` flag:

```bash
nextflow run 2c-modules.nf --input data/greetings.csv -resume
```

??? success "Command output"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `2c-modules.nf` [soggy_franklin] DSL2 - revision: bc8e1b2726

    [j6/cdfa66] sayHello (1)       | 3 of 3, cached: ✔
    [95/79484f] convertToUpper (2) | 3 of 3, cached: ✔
    [5e/4358gc] collectGreetings   | 1 of 1, cached: ✔
    ```

You'll notice that the process executions all cached successfully, meaning that Nextflow recognized that it has already done the requested work, even though the code has been split up and the main workflow file has been renamed.

None of that matters to Nextflow; what matters is the job script that is generated once all the code has been pulled together and evaluated.

!!! tip

    It is also possible to encapsulate a section of a workflow as a 'subworkflow' that can be imported into a larger pipeline, but that is outside the scope of this course.

    You can learn more about developing composable workflows in the Side Quest on [Workflows of Workflows](https://training.nextflow.io/latest/side_quests/workflows_of_workflows/).

### Takeaway

You know how processes can be stored in standalone modules to promote code reuse and improve maintainability.

### What's next?

Learn to use containers for managing software dependencies.

---

## 4. Using containerized software

So far the workflows we've been using as examples just needed to run very basic text processing operations using UNIX tools available in our environment.

However, real-world pipelines typically require specialized tools and packages that are not included by default in most environments.
Usually, you'd need to install these tools, manage their dependencies, and resolve any conflicts.

That is all very tedious and annoying.
A much better way to address this problem is to use **containers**.

A **container** is a lightweight, standalone, executable unit of software created from a container **image** that includes everything needed to run an application including code, system libraries and settings.

!!! Tip

    We teach this using the technology [Docker](https://www.docker.com/get-started/), but Nextflow supports [several other container technologies](https://www.nextflow.io/docs/latest/container.html#) as well.

### 4.1. Use a container directly

First, let's try interacting with a container directly.
This will help solidify your understanding of what containers are before we start using them in Nextflow.

#### 4.1.1. Pull the container image

To use a container, you usually download or "pull" a container image from a container registry, and then run the container image to create a container instance.

The general syntax is as follows:

```bash title="Syntax"
docker pull '<container>'
```

- `docker pull` is the instruction to the container system to pull a container image from a repository.
- `'<container>'` is the URI address of the container image.

As an example, let's pull a container image that contains [cowpy](https://github.com/jeffbuttars/cowpy), a python implementation of a tool called `cowsay` that generates ASCII art to display arbitrary text inputs in a fun way.

There are various repositories where you can find published containers.
We used the [Seqera Containers](https://seqera.io/containers/) service to generate this Docker container image from the `cowpy` Conda package: `'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'`.

Run the complete pull command:

```bash
docker pull 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'
```

??? success "Command output"

    ```console
    Unable to find image 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273' locally
    131d6a1b707a8e65: Pulling from library/cowpy
    dafa2b0c44d2: Pull complete
    dec6b097362e: Pull complete
    f88da01cff0b: Pull complete
    4f4fb700ef54: Pull complete
    92dc97a3ef36: Pull complete
    403f74b0f85e: Pull complete
    10b8c00c10a5: Pull complete
    17dc7ea432cc: Pull complete
    bb36d6c3110d: Pull complete
    0ea1a16bbe82: Pull complete
    030a47592a0a: Pull complete
    622dd7f15040: Pull complete
    895fb5d0f4df: Pull complete
    Digest: sha256:fa50498b32534d83e0a89bb21fec0c47cc03933ac95c6b6587df82aaa9d68db3
    Status: Downloaded newer image for community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273
    community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273
    ```

This tells the system to download the image specified.
Once the download is complete, you have a local copy of the container image.

#### 4.1.2. Spin up the container

Containers can be run as a one-off command, but you can also use them interactively, which gives you a shell prompt inside the container and allows you to play with the command.

The general syntax is as follows:

```bash title="Syntax"
docker run --rm '<container>' [tool command]
```

- `docker run --rm '<container>'` is the instruction to the container system to spin up a container instance from a container image and execute a command in it.
- `--rm` tells the system to shut down the container instance after the command has completed.

Fully assembled, the container execution command looks like this:

```bash
docker run --rm -it 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'
```

Run that command, and you should see your prompt change to something like `(base) root@b645838b3314:/tmp#`, which indicates that you are now inside the container.

You can verify this by running `ls` to list directory contents:

```bash
ls /
```

??? success "Command output"

    ```console
    bin  boot  dev  etc  home  lib  lib64  media  mnt  opt  proc  root  run  sbin  srv  sys  tmp  usr  var
    ```

You see that the filesystem inside the container is different from the filesystem on your host system.

!!! Tip

    When you run a container, it is isolated from the host system by default.
    This means that the container can't access any files on the host system unless you explicitly allow it to do so by specifying that you want to mount a volume as part of the `docker run` command using the following syntax:

    ```bash title="Syntax"
    -v <outside_path>:<inside_path>
    ```

    This effectively establishes a tunnel through the container wall that you can use to access that part of your filesystem.

    This is covered in more detail in [Part 5 of Hello Nextflow](../hello_nextflow/05_hello_containers.md).

#### 4.1.3. Run the `cowpy` tool

From inside the container, you can run the `cowpy` command directly.

```bash
cowpy "Hello Containers"
```

??? success "Command output"

    ```console
    ______________________________________________________
    < Hello Containers >
    ------------------------------------------------------
        \   ^__^
          \  (oo)\_______
            (__)\       )\/\
              ||----w |
              ||     ||
    ```

This produces ASCII art of the default cow character (or 'cowacter') with a speech bubble containing the text we specified.

Now that you have tested the basic usage, you can try giving it some parameters.
For example, the tool documentation says we can set the character with `-c`.

```bash
cowpy "Hello Containers" -c tux
```

??? success "Command output"

    ```console
    __________________
    < Hello Containers >
    ------------------
      \
        \
            .--.
          |o_o |
          |:_/ |
          //   \ \
        (|     | )
        /'\_   _/`\
        \___)=(___/
    ```

This time the ASCII art output shows the Linux penguin, Tux, because we specified the `-c tux` parameter.

Since you're inside the container, you can run the cowpy command as many times as you like, varying the input parameters, without having to worry about install any libraries on your system itself.

??? tip "Other available characters"

    Use the '-c' flag to pick a different character, including:

    `beavis`, `cheese`, `daemon`, `dragonandcow`, `ghostbusters`, `kitty`, `moose`, `milk`, `stegosaurus`, `turkey`, `turtle`, `tux`

Feel free to play around with this.
When you're done, exit the container using the `exit` command:

```bash
exit
```

You will find yourself back in your normal shell.

### 4.2. Use a container in a workflow

When we run a pipeline, we want to be able to tell Nextflow what container to use at each step, and importantly, we want it to handle all that work we just did: pull the container, spin it up, run the command and tear the container down when it's done.

Good news: that's exactly what Nextflow is going to do for us.
We just need to specify a container for each process.

To demonstrate how this work, we made another version of our workflow that runs `cowpy` on the file of collected greetings produced in the third step.

<figure class="excalidraw">
--8<-- "docs/nextflow_run/img/hello-pipeline-cowpy.svg"
</figure>

This should output a file containing the ASCII art with the three greetings in the speech bubble.

#### 4.2.1. Examine the code

The workflow is very similar to the previous one, plus the extra step to run `cowpy`.

??? full-code "Full code file"

    ```groovy title="2d-container.nf" linenums="1" hl_lines="7 15 32 39 59-62"
    #!/usr/bin/env nextflow

    // Include modules
    include { sayHello } from './modules/sayHello.nf'
    include { convertToUpper } from './modules/convertToUpper.nf'
    include { collectGreetings } from './modules/collectGreetings.nf'
    include { cowpy } from './modules/cowpy.nf'

    /*
    * Pipeline parameters
    */
    params {
        input: Path
        batch: String = 'batch'
        character: String
    }

    workflow {

        main:
        // create a channel for inputs from a CSV file
        greeting_ch = channel.fromPath(params.input)
                            .splitCsv()
                            .map { line -> line[0] }
        // emit a greeting
        sayHello(greeting_ch)
        // convert the greeting to uppercase
        convertToUpper(sayHello.out)
        // collect all the greetings into one file
        collectGreetings(convertToUpper.out.collect(), params.batch)
        // generate ASCII art of the greetings with cowpy
        cowpy(collectGreetings.out.outfile, params.character)

        publish:
        first_output = sayHello.out
        uppercased = convertToUpper.out
        collected = collectGreetings.out.outfile
        batch_report = collectGreetings.out.report
        cowpy_art = cowpy.out
    }

    output {
        first_output {
            path '2d-container/intermediates'
            mode 'copy'
        }
        uppercased {
            path '2d-container/intermediates'
            mode 'copy'
        }
        collected {
            path '2d-container/intermediates'
            mode 'copy'
        }
        batch_report {
            path '2d-container'
            mode 'copy'
        }
        cowpy_art {
            path '2d-container'
            mode 'copy'
        }
    }
    ```

You see that this workflow imports a `cowpy` process from a module file, and calls it on the output of the `collectGreetings()` call, plus an input parameter called `params.character`.

```groovy title="2d-container.nf" linenums="25"
// generate ASCII art with cowpy
cowpy(collectGreetings.out, params.character)
```

The `cowpy` process, which wraps the cowpy command to generate ASCII art, is defined in the `cowpy.nf` module.

??? full-code "Full code file"

    ```groovy title="modules/cowpy.nf" linenums="1"
    #!/usr/bin/env nextflow

    // Generate ASCII art with cowpy (https://github.com/jeffbuttars/cowpy)
    process cowpy {

        container 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'

        input:
        path input_file
        val character

        output:
        path "cowpy-${input_file}"

        script:
        """
        cat ${input_file} | cowpy -c "${character}" > cowpy-${input_file}
        """
    }
    ```

The `cowpy` process requires two inputs: the path to an input file containing the text to put in the speech bubble (`input_file`), and a value for the character variable.

Importantly, it also includes the line `container 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'`, which points to the container URI we used earlier.

#### 4.2.2. Check that Docker is enabled in the configuration

We're going to slightly anticipate Part 3 of this training course by introducing the `nextflow.config` configuration file, which is one of the main ways Nextflow offers for configuring workflow execution.
When a file named `nextflow.config` is present in the current directory, Nextflow will automatically load it in and apply any configuration it contains.

To that end, we included a `nextflow.config` file with a single line of code that enables Docker.

```groovy title="nextflow.config" linenums="1"
docker.enabled = true
```

This configuration tells Nextflow to use Docker for any process that specifies a compatible container.

!!! tip

    It is technically possible to enable Docker execution from the command-line, on a per-run basis, using the `-with-docker <container>` parameter.
    However, that only allows us to specify one container for the entire workflow, whereas the approach we just showed you allows us to specify a different container per process.
    The latter is much better for modularity, code maintenance and reproducibility.

#### 4.2.3. Run the workflow

Just to recap, this is what we are about to run:

<figure class="excalidraw">
--8<-- "docs/hello_nextflow/img/hello_pipeline_complete.svg"
</figure>

Do you think it's going to work?

Let's run the workflow with the `-resume` flag, and specify that we want the character to be the turkey.

```bash
nextflow run 2d-container.nf --input data/greetings.csv --character turkey -resume
```

??? success "Command output"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `2d-container.nf` [elegant_brattain] DSL2 - revision: 028a841db1

    executor >  local (1)
    [95/fa0bac] sayHello (3)       | 3 of 3, cached: 3 ✔
    [92/32533f] convertToUpper (3) | 3 of 3, cached: 3 ✔
    [aa/e697a2] collectGreetings   | 1 of 1, cached: 1 ✔
    [7f/caf718] cowpy              | 1 of 1 ✔
    ```

The first three steps cached since we've already run them before, but the `cowpy` process is new so that actually gets run.

You can find the output of the `cowpy` step in the `results` directory.

??? abstract "File contents"

    ```console title="results/2d-container/cowpy-COLLECTED-batch-output.txt"
    _________
    / HOLà    \
    | HELLO   |
    \ BONJOUR /
    ---------
      \                                  ,+*^^*+___+++_
      \                           ,*^^^^              )
        \                       _+*                     ^**+_
        \                    +^       _ _++*+_+++_,         )
                  _+^^*+_    (     ,+*^ ^          \+_        )
                {       )  (    ,(    ,_+--+--,      ^)      ^\
                { (\@)    } f   ,(  ,+-^ __*_*_  ^^\_   ^\       )
              {:;-/    (_+*-+^^^^^+*+*<_ _++_)_    )    )      /
              ( /  (    (        ,___    ^*+_+* )   <    <      \
              U _/     )    *--<  ) ^\-----++__)   )    )       )
                (      )  _(^)^^))  )  )\^^^^^))^*+/    /       /
              (      /  (_))_^)) )  )  ))^^^^^))^^^)__/     +^^
            (     ,/    (^))^))  )  ) ))^^^^^^^))^^)       _)
              *+__+*       (_))^)  ) ) ))^^^^^^))^^^^^)____*^
              \             \_)^)_)) ))^^^^^^^^^^))^^^^)
              (_             ^\__^^^^^^^^^^^^))^^^^^^^)
                ^\___            ^\__^^^^^^))^^^^^^^^)\\
                      ^^^^^\uuu/^^\uuu/^^^^\^\^\^\^\^\^\^\
                        ___) >____) >___   ^\_\_\_\_\_\_\)
                        ^^^//\\_^^//\\_^       ^(\_\_\_\)
                          ^^^ ^^ ^^^ ^
    ```

You see that the character is saying all the greetings, since it ran on the file of collected uppercased greetings.

More to the point, we were able to run this as part of our pipeline without having to do a proper installation of cowpy and all its dependencies.
And we can now share the pipeline with collaborators and have them run it on their infrastructure without them needing to install anything either, aside from Docker or one of its alternatives (such as Singularity/Apptainer) as mentioned above.

#### 4.2.4. Inspect how Nextflow launched the containerized task

As a final coda to this section, let's take a look at the work subdirectory for one of the `cowpy` process calls to get a bit more insight on how Nextflow works with containers under the hood.

Check the output from your `nextflow run` command to find the path to the work subdirectory for the `cowpy` process.
Looking at what we got for the run shown above, the console log line for the `cowpy` process starts with `[7f/caf718]`.
That corresponds to the following truncated directory path: `work/7f/caf718`.

In that directory, you will find the `.command.run` file that contains all the commands Nextflow ran on your behalf in the course of executing the pipeline.

??? abstract "File contents"

    ```console title="work/7f/caf71890cce1667c094d880f4b6dcc/.command.run"
    #!/bin/bash
    ### ---
    ### name: 'cowpy'
    ### container: 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'
    ### outputs:
    ### - 'cowpy-COLLECTED-batch-output.txt'
    ### ...
    set -e
    set -u
    NXF_DEBUG=${NXF_DEBUG:=0}; [[ $NXF_DEBUG > 1 ]] && set -x
    NXF_ENTRY=${1:-nxf_main}


    nxf_sleep() {
      sleep $1 2>/dev/null || sleep 1;
    }

    nxf_date() {
        local ts=$(date +%s%3N);
        if [[ ${#ts} == 10 ]]; then echo ${ts}000
        elif [[ $ts == *%3N ]]; then echo ${ts/\%3N/000}
        elif [[ $ts == *3N ]]; then echo ${ts/3N/000}
        elif [[ ${#ts} == 13 ]]; then echo $ts
        else echo "Unexpected timestamp value: $ts"; exit 1
        fi
    }

    nxf_env() {
        echo '============= task environment ============='
        env | sort | sed "s/\(.*\)AWS\(.*\)=\(.\{6\}\).*/\1AWS\2=\3xxxxxxxxxxxxx/"
        echo '============= task output =================='
    }

    nxf_kill() {
        declare -a children
        while read P PP;do
            children[$PP]+=" $P"
        done < <(ps -e -o pid= -o ppid=)

        kill_all() {
            [[ $1 != $$ ]] && kill $1 2>/dev/null || true
            for i in ${children[$1]:=}; do kill_all $i; done
        }

        kill_all $1
    }

    nxf_mktemp() {
        local base=${1:-/tmp}
        mkdir -p "$base"
        if [[ $(uname) = Darwin ]]; then mktemp -d $base/nxf.XXXXXXXXXX
        else TMPDIR="$base" mktemp -d -t nxf.XXXXXXXXXX
        fi
    }

    nxf_fs_copy() {
      local source=$1
      local target=$2
      local basedir=$(dirname $1)
      mkdir -p $target/$basedir
      cp -fRL $source $target/$basedir
    }

    nxf_fs_move() {
      local source=$1
      local target=$2
      local basedir=$(dirname $1)
      mkdir -p $target/$basedir
      mv -f $source $target/$basedir
    }

    nxf_fs_rsync() {
      rsync -rRl $1 $2
    }

    nxf_fs_rclone() {
      rclone copyto $1 $2/$1
    }

    nxf_fs_fcp() {
      fcp $1 $2/$1
    }

    on_exit() {
        local last_err=$?
        local exit_status=${nxf_main_ret:=0}
        [[ ${exit_status} -eq 0 && ${nxf_unstage_ret:=0} -ne 0 ]] && exit_status=${nxf_unstage_ret:=0}
        [[ ${exit_status} -eq 0 && ${last_err} -ne 0 ]] && exit_status=${last_err}
        printf -- $exit_status > /workspaces/training/nextflow-run/work/7f/caf71890cce1667c094d880f4b6dcc/.exitcode
        set +u
        docker rm $NXF_BOXID &>/dev/null || true
        exit $exit_status
    }

    on_term() {
        set +e
        docker stop $NXF_BOXID
    }

    nxf_launch() {
        docker run -i --cpu-shares 1024 -e "NXF_TASK_WORKDIR" -v /workspaces/training/nextflow-run/work:/workspaces/training/nextflow-run/work -w "$NXF_TASK_WORKDIR" --name $NXF_BOXID community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273 /bin/bash -ue /workspaces/training/nextflow-run/work/7f/caf71890cce1667c094d880f4b6dcc/.command.sh
    }

    nxf_stage() {
        true
        # stage input files
        rm -f COLLECTED-batch-output.txt
        ln -s /workspaces/training/nextflow-run/work/7f/f435e3f2cf95979b5f3d7647ae6696/COLLECTED-batch-output.txt COLLECTED-batch-output.txt
    }

    nxf_unstage_outputs() {
        true
    }

    nxf_unstage_controls() {
        true
    }

    nxf_unstage() {
        if [[ ${nxf_main_ret:=0} == 0 ]]; then
            (set -e -o pipefail; (nxf_unstage_outputs | tee -a .command.out) 3>&1 1>&2 2>&3 | tee -a .command.err)
            nxf_unstage_ret=$?
        fi
        nxf_unstage_controls
    }

    nxf_main() {
        trap on_exit EXIT
        trap on_term TERM INT USR2
        trap '' USR1

        [[ "${NXF_CHDIR:-}" ]] && cd "$NXF_CHDIR"
        export NXF_BOXID="nxf-$(dd bs=18 count=1 if=/dev/urandom 2>/dev/null | base64 | tr +/ 0A | tr -d '\r\n')"
        NXF_SCRATCH=''
        [[ $NXF_DEBUG > 0 ]] && nxf_env
        touch /workspaces/training/nextflow-run/work/7f/caf71890cce1667c094d880f4b6dcc/.command.begin
        set +u
        set -u
        [[ $NXF_SCRATCH ]] && cd $NXF_SCRATCH
        export NXF_TASK_WORKDIR="$PWD"
        nxf_stage

        set +e
        (set -o pipefail; (nxf_launch | tee .command.out) 3>&1 1>&2 2>&3 | tee .command.err) &
        pid=$!
        wait $pid || nxf_main_ret=$?
        nxf_unstage
    }

    $NXF_ENTRY
    ```

If you search for `nxf_launch` in this file, you should see something like this:

```console
nxf_launch() {
    docker run -i --cpu-shares 1024 -e "NXF_TASK_WORKDIR" -v /workspaces/training/nextflow-run/work:/workspaces/training/nextflow-run/work -w "$NXF_TASK_WORKDIR" --name $NXF_BOXID community.wave.seqera.io/library/pip_cowpy:131d6a1b707a8e65 /bin/bash -ue /workspaces/training/nextflow-run/work/7f/caf7189fca6c56ba627b75749edcb3/.command.sh
}
```

This launch command shows that Nextflow is using a very similar `docker run` command to launch the process call as we did when we ran it manually.
It also mounts the corresponding work subdirectory into the container, sets the working directory inside the container accordingly, and runs our templated bash script in the `.command.sh` file.

This confirms that all the hard work we had to do manually in the previous section is now done for us by Nextflow!

### Takeaway

You understand what role containers play in managing software tool versions and ensuring reproducibility.

More generally, you have a basic understanding of what are the core components of real-world Nextflow pipelines and how they are organized.
You know the fundamentals of how Nextflow can process multiple inputs efficiently, run workflows composed of multiple steps connected together, leverage modular code components, and utilize containers for greater reproducibility and portability.

### What's next?

Take another break! That was a big pile of information about how Nextflow pipelines work.

In the last section of this training, we're going to delve deeper into the topic of configuration.
You will learn how to configure the execution of your pipeline to fit your infrastructure as well as manage configuration of inputs and parameters.

---

## Quiz

<quiz>
Why does Nextflow create a separate task directory for each process call?
- [ ] To improve execution speed
- [ ] To reduce memory usage
- [x] To isolate executions and avoid collisions between outputs
- [ ] To enable parallel file compression
</quiz>

<quiz>
What does the `-ansi-log false` option do when running a workflow?
- [ ] Disables all console output
- [ ] Removes color from the output
- [x] Shows all task directory paths instead of condensing them on one line
- [ ] Enables verbose debugging mode
</quiz>

<quiz>
In the code `channel.fromPath(params.input).splitCsv().map { line -> line[0] }`, what does `.map { line -> line[0] }` do?
- [ ] Filters out empty lines
- [ ] Sorts the lines alphabetically
- [x] Extracts the first column from each CSV row
- [ ] Counts the number of lines
</quiz>

<quiz>
Why is it important to include the input value in output filenames (e.g., `${greeting}-output.txt`)?
- [ ] To improve processing speed
- [ ] To enable resume functionality
- [x] To prevent output files from overwriting each other when processing multiple inputs
- [ ] To make files easier to compress
</quiz>

<quiz>
What is the purpose of the `include` statement in a modularized workflow?
- [ ] To copy process code into the workflow file
- [x] To import a process definition from an external module file
- [ ] To include configuration settings
- [ ] To add documentation comments
</quiz>

<quiz>
When you modularize a workflow and run it with -resume, what happens?
- [ ] Caching is disabled for modular processes
- [ ] All tasks must be re-executed
- [x] Caching works normally based on the generated job scripts
- [ ] Only the main workflow file is cached
</quiz>

<quiz>
What does the `container` directive in a process definition specify?
- [ ] The working directory for the process
- [ ] The maximum memory allocation
- [x] The container image URI to use for running the process
- [ ] The output file format
</quiz>

<quiz>
In the `.command.run` file, what does the `nxf_launch` function contain?
- [ ] The Nextflow version information
- [ ] The workflow parameters
- [x] The docker run command with volume mounts and container settings
- [ ] The process input declarations
</quiz>

<quiz>
What does Nextflow automatically handle when running a containerized process? (Select all that apply)
- [x] Pulling the container image if needed
- [x] Mounting the work directory into the container
- [x] Running the process script inside the container
- [x] Cleaning up the container instance after execution
</quiz>
