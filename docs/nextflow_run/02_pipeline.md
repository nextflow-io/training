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
Disabling ANSI logging also prevented Nextflow from using colours in the terminal output.

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

You don't need to memorize code syntax, but it's good to recognize key components of the workflow.

#### 1.4.1. Channels handle input data

The key to processing multiple inputs is the **channel**: a queue that holds data and shuttles it between workflow steps.

```groovy title="2a-inputs.nf" linenums="29" hl_lines="3-5"
    main:
    // create a channel for inputs from a CSV file
    greeting_ch = channel.fromPath(params.input)
                        .splitCsv()
                        .map { line -> line[0] }
    // emit a greeting
    sayHello(greeting_ch)
```

This code reads the CSV file, parses it, and extracts the first column from each row.
The result is a channel containing `Hello`, `Bonjour`, and `Holà`.

When `sayHello(greeting_ch)` is called, Nextflow automatically runs the process on each element in the channel, in parallel when possible.

!!! tip "Want to learn more about channels?"

    If you want to understand channels and operators in depth, including how to write them yourself, see [Hello Nextflow Part 2: Hello Channels](../hello_nextflow/02_hello_channels.md).

#### 1.4.2. Dynamic output naming

The process uses the input value in the output filename to ensure unique names:

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

Including the input value (`${greeting}`) in the output filename prevents collisions when multiple inputs are processed.

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

Let's look at the code and identify the key patterns for multi-step workflows.

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

#### 2.3.1. Visualizing workflow structure

If you're using VSCode with the Nextflow extension, click the `DAG preview` link above any workflow block to see how processes are connected:

<figure class="excalidraw">
--8<-- "docs/nextflow_run/img/DAG-multistep.svg"
</figure>

#### 2.3.2. Connecting processes

Processes are connected by passing outputs as inputs:

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

The pattern is simple: `processName.out` refers to a process's output channel, which can be passed directly to the next process.

#### 2.3.3. The `collect()` operator

Notice `convertToUpper.out.collect()` in the final process call.
The `collect()` operator gathers all items from a channel into a single element.

Without `collect()`, `collectGreetings` would run separately on each greeting:

<figure class="excalidraw">
--8<-- "docs/nextflow_run/img/without-collect-operator.svg"
</figure>

With `collect()`, all greetings are gathered and processed together:

<figure class="excalidraw">
--8<-- "docs/nextflow_run/img/with-collect-operator.svg"
</figure>

#### 2.3.4. Multiple inputs and outputs

Processes can accept multiple inputs (comma-separated in the call) and produce multiple outputs (using `emit:` to name them):

```groovy title="2b-multistep.nf" linenums="46"
    output:
    path "COLLECTED-${batch_name}-output.txt", emit: outfile
    path "${batch_name}-report.txt", emit: report
```

Named outputs are referenced as `processName.out.outputName` (e.g., `collectGreetings.out.outfile`).

#### 2.3.5. Default parameter values

Parameters can have default values in the `params` block:

```groovy title="2b-multistep.nf" linenums="58" hl_lines="6"
    params {
        input: Path
        batch: String = 'batch'
    }
```

If you don't specify `--batch` on the command line, the default value `'batch'` is used.

!!! tip "Want to learn more about building workflows?"

    For detailed coverage of building multi-step workflows, see [Hello Nextflow Part 3: Hello Workflow](../hello_nextflow/03_hello_workflow.md).

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

Learn more: [1.3. Find the original outputs and logs](#13-find-the-original-outputs-and-logs)
</quiz>

<quiz>
What does the `-ansi-log false` option do when running a workflow?
- [ ] Disables all console output
- [x] Removes color from the output
- [x] Shows all task directory paths instead of condensing them on one line
- [ ] Enables verbose debugging mode

Learn more: [1.3.2. Make the terminal show more details](#132-make-the-terminal-show-more-details)

You can also use either of the following environment variables if you prefer this style:

```bash
export NXF_ANSI_LOG=0
# or
export NO_COLOR=1
```

</quiz>

<quiz>
In the code `#!groovy channel.fromPath(params.input).splitCsv().map { line -> line[0] }`, what does `#!groovy .map { line -> line[0] }` do?
- [ ] Filters out empty lines
- [ ] Sorts the lines alphabetically
- [x] Extracts the first column from each CSV row
- [ ] Counts the number of lines

Learn more: [1.4.1. Loading the input data from the CSV](#141-loading-the-input-data-from-the-csv)
</quiz>

<quiz>
Why is it important to include the input value in output filenames (e.g., `#!groovy "${greeting}-output.txt"`)?
- [ ] To improve processing speed
- [ ] To enable resume functionality
- [x] To prevent output files from overwriting each other when processing multiple inputs
- [ ] To make files easier to compress

Learn more: [1.4.3. How the outputs are named](#143-how-the-outputs-are-named)
</quiz>

<quiz>
What is the purpose of the `include` statement in a modularized workflow?
- [ ] To copy process code into the workflow file
- [x] To import a process definition from an external module file
- [ ] To include configuration settings
- [ ] To add documentation comments

Learn more: [3. Running modularized pipelines](#3-running-modularized-pipelines)
</quiz>

<quiz>
When you modularize a workflow and run it with `-resume`, what happens?
- [ ] Caching is disabled for modular processes
- [ ] All tasks must be re-executed
- [x] Caching works normally based on the generated job scripts
- [ ] Only the main workflow file is cached

Learn more: [3.2. Run the workflow](#32-run-the-workflow)
</quiz>

<quiz>
What does the `container` directive in a process definition specify?
- [ ] The working directory for the process
- [ ] The maximum memory allocation
- [x] The container image URI to use for running the process
- [ ] The output file format

Learn more: [4.2. Use a container in a workflow](#42-use-a-container-in-a-workflow)
</quiz>

<quiz>
In the `.command.run` file, what does the `nxf_launch` function contain?
- [ ] The Nextflow version information
- [ ] The workflow parameters
- [x] The `docker run` command with volume mounts and container settings
- [ ] The process input declarations

Learn more: [4.2.4. Inspect how Nextflow launched the containerized task](#424-inspect-how-nextflow-launched-the-containerized-task)
</quiz>

<quiz>
What does Nextflow automatically handle when running a containerized process? (Select all that apply)
- [x] Pulling the container image if needed
- [x] Mounting the work directory into the container
- [x] Running the process script inside the container
- [x] Cleaning up the container instance after execution

Learn more: [4. Using containerized software](#4-using-containerized-software)
</quiz>
