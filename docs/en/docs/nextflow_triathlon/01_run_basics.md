# Part 1: Run Nextflow

In this first part of the Nextflow Triathlon, we introduce the core concepts of running Nextflow pipelines.
We start with a simple Hello World workflow, then progress to multi-step pipelines that process multiple inputs in parallel using containers.

!!! tip

    Make sure your working directory is set to `nextflow-run/` before starting this part.

    ```bash
    cd /workspaces/training/nextflow-run
    ```

---

## 1. Run a workflow

We provide you with a workflow script named `1-hello.nf` that takes a greeting via a command-line argument named `--input` and writes it to a file.

### 1.1. Launch the workflow

Run the following command in your terminal.

```bash
nextflow run 1-hello.nf --input 'Hello World!'
```

??? success "Command output"

    ```console hl_lines="6"
    N E X T F L O W   ~  version 25.10.4

    Launching `1-hello.nf` [goofy_torvalds] DSL2 - revision: c33d41f479

    executor >  local (1)
    [a3/7be2fa] sayHello | 1 of 1 ✔
    ```

The most important line in the output is the last one:

```console
[a3/7be2fa] sayHello | 1 of 1 ✔
```

This tells us that the `sayHello` process ran successfully once.
The `[a3/7be2fa]` prefix is a truncated path to the task's working directory — more on that below.

### 1.2. Find the output

This workflow is configured to publish its output to a `results` directory.
After running, you should find the output there:

```console title="results/"
results
└── 1-hello
    └── output.txt
```

Open the file to confirm it contains `Hello World!`.

### 1.3. Understand the work/ directory

Behind the scenes, Nextflow creates a unique task directory for every process call inside a directory named `work/`.
The hash shown in the console output (`[a3/7be2fa]`) is the path to that directory.

```bash
ls work/a3/7be2fa*
```

Inside you will find the output file along with several hidden log files:

- **`.command.sh`**: the exact command Nextflow ran
- **`.command.out`** / **`.command.err`**: stdout and stderr from the process
- **`.command.log`**: combined log output
- **`.exitcode`**: the process exit code

The `.command.sh` file is especially useful when debugging — it shows precisely what was executed.

### 1.4. Re-run with -resume

If you re-run the same command, Nextflow will repeat all the work.
Add `-resume` to skip any steps whose inputs and code have not changed:

```bash
nextflow run 1-hello.nf --input 'Hello World!' -resume
```

??? success "Command output"

    ```console hl_lines="5"
    N E X T F L O W   ~  version 25.10.4

    Launching `1-hello.nf` [tiny_noyce] DSL2 - revision: c33d41f479

    [a3/7be2fa] sayHello | 1 of 1, cached: 1 ✔
    ```

The `cached: 1` status tells you Nextflow reused the result from the previous run.
This is particularly valuable when iterating on a long pipeline — you only re-run the steps that changed.

### 1.5. Inspect the execution log

Use `nextflow log` to view the history of all pipeline runs launched from the current directory:

```bash
nextflow log
```

??? success "Command output"

    ```console
    TIMESTAMP               DURATION    RUN NAME              STATUS  REVISION ID     SESSION ID                              COMMAND
    2025-07-04 19:27:09     1.8s        wise_watson           OK      3539118582      a02c9c46-c3c7-4085-9139-d1b9b5b194c8    nextflow run 1-hello.nf --input 'Hello World'
    2025-07-04 19:27:57     2.1s        goofy_wilson          OK      3539118582      5f4b3269-5b53-404a-956c-cac915fbb74e    nextflow run 1-hello.nf --input 'Hello World' -resume
    ```

Runs that used `-resume` share a session ID with the run they resumed from.

### Takeaway

You know how to run a Nextflow pipeline, find its outputs, and use `-resume` to avoid repeating completed work.

### What's next?

Learn what a Nextflow workflow script looks like and how its key components work.

---

## 2. Understand the workflow code

Let's open `1-hello.nf` and look at the three main components.

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

### 2.1. The process block

A **process** defines a single step in the pipeline.
It declares its inputs, outputs and the script to execute:

```groovy title="1-hello.nf" linenums="6"
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

The `val` qualifier means the input is a plain value (string, number, etc.).
The `path` qualifier means the output should be treated as a file path.

### 2.2. The params block

The `params` block declares the command-line parameters the workflow accepts:

```groovy title="1-hello.nf" linenums="22"
params {
    input: String
}
```

Any parameter declared here becomes available on the command line with a double-dash (`--input`).
Supported types include `String`, `Integer`, `Float`, `Boolean`, and `Path`.

!!! tip

    Workflow parameters always use two dashes (`--input`) to distinguish them from Nextflow's own CLI flags, which use one dash (e.g. `-resume`).

### 2.3. The workflow block

The **workflow** block defines the dataflow logic — which processes to run and in what order:

```groovy title="1-hello.nf" linenums="27"
workflow {

    main:
    sayHello(params.input)

    publish:
    first_output = sayHello.out
}
```

The `main:` block calls the `sayHello` process with the `--input` value.
The `publish:` block lists which outputs should be copied to the results directory, where the `output` block at the bottom specifies the destination path and copy mode.

### Takeaway

You understand the three main components of a Nextflow workflow: the `process`, `params`, and `workflow` blocks.

### What's next?

Discover how Nextflow handles multiple inputs efficiently.

---

## 3. Process multiple inputs

Real-world pipelines process many samples, not just one.
Nextflow uses **channels** — queue constructs that shuttle data between steps — to handle multiple inputs efficiently and in parallel.

We provide a workflow called `2a-inputs.nf` that reads from a CSV file and processes each row separately.

### 3.1. Run the workflow

```bash
nextflow run 2a-inputs.nf --input data/greetings.csv
```

??? success "Command output"

    ```console
    N E X T F L O W   ~  version 25.10.4

    Launching `2a-inputs.nf` [mighty_sammet] DSL2 - revision: 29fb5352b3

    executor >  local (3)
    [8e/0eb066] sayHello (2) [100%] 3 of 3 ✔
    ```

The `3 of 3` tells us the `sayHello` process was called three times — once per row in the CSV.

In the `results` directory, you should now see three output files, one per greeting:

```console title="results/2a-inputs/"
2a-inputs
├── Bonjour-output.txt
├── Hello-output.txt
└── Hola-output.txt
```

### 3.2. How the channel works

The key change in `2a-inputs.nf` is in the `main:` block of the workflow:

```groovy title="2a-inputs.nf" linenums="29"
    main:
    // create a channel for inputs from a CSV file
    greeting_ch = channel.fromPath(params.input)
                        .splitCsv()
                        .map { line -> line[0] }
    // emit a greeting
    sayHello(greeting_ch)
```

- `channel.fromPath(params.input)` creates a channel from the file path given with `--input`
- `.splitCsv()` parses the CSV into rows
- `.map { line -> line[0] }` extracts the first column from each row

The result is a channel containing `Hello`, `Bonjour`, and `Hola`.
When passed to `sayHello(greeting_ch)`, Nextflow automatically calls the process once per item, running them in parallel when resources allow.

### Takeaway

You understand how channels load data from files and enable automatic parallelism.

### What's next?

Learn how multiple processes are connected in a multi-step workflow.

---

## 4. Run a multi-step workflow

Most real-world pipelines chain several processes together.
The workflow `2b-multistep.nf` demonstrates this: it takes each greeting, converts it to uppercase, then collects all the results into a single output file.

### 4.1. Run the workflow

```bash
nextflow run 2b-multistep.nf --input data/greetings.csv
```

??? success "Command output"

    ```console
    N E X T F L O W   ~  version 25.10.4

    Launching `2b-multistep.nf` [soggy_franklin] DSL2 - revision: bc8e1b2726

    [d6/cdf466] sayHello (1)       | 3 of 3 ✔
    [99/79394f] convertToUpper (2) | 3 of 3 ✔
    [1e/83586c] collectGreetings   | 1 of 1 ✔
    ```

Three processes ran: `sayHello` and `convertToUpper` each ran once per input (3 of 3), while `collectGreetings` ran once on all results combined (1 of 1).

### 4.2. How the processes are connected

In the `main:` block, each process call feeds its output to the next:

```groovy title="2b-multistep.nf" linenums="68"
    main:
    greeting_ch = channel.fromPath(params.input)
                        .splitCsv()
                        .map { line -> line[0] }
    sayHello(greeting_ch)
    convertToUpper(sayHello.out)
    collectGreetings(convertToUpper.out.collect(), params.batch)
```

The pattern `processName.out` refers to a process's output channel.
Passing it as input to the next process is how data flows from step to step.

The `.collect()` operator in `collectGreetings(convertToUpper.out.collect(), ...)` gathers all individual outputs from `convertToUpper` into a single channel item, which is how `collectGreetings` receives all three greetings at once rather than running three times separately.

### Takeaway

You know how multi-step workflows chain processes together using output channels and operators.

### What's next?

Learn how containerized software is used in Nextflow pipelines.

---

## 5. Use containerized software

Real pipelines depend on specialized tools that would normally require complex installation.
Containers solve this by bundling a tool and all its dependencies into a portable unit.

The workflow `2d-container.nf` adds a fourth step that uses the `cowpy` tool (packaged in a Docker container) to generate ASCII art from the collected greetings.

### 5.1. Examine the container directive

Open `modules/cowpy.nf` and notice the `container` directive:

```groovy title="modules/cowpy.nf" linenums="1"
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

The `container` directive specifies the URI of a Docker image.
When Nextflow runs this process, it automatically pulls the image, spins up a container, runs the script inside it, and cleans up afterward.
Docker is enabled for this project in `nextflow.config`:

```groovy title="nextflow.config"
docker.enabled = true
```

### 5.2. Run the workflow

```bash
nextflow run 2d-container.nf --input data/greetings.csv --character turkey -resume
```

??? success "Command output"

    ```console
    N E X T F L O W   ~  version 25.10.4

    Launching `2d-container.nf` [elegant_brattain] DSL2 - revision: 028a841db1

    executor >  local (1)
    [95/fa0bac] sayHello (3)       | 3 of 3, cached: 3 ✔
    [92/32533f] convertToUpper (3) | 3 of 3, cached: 3 ✔
    [aa/e697a2] collectGreetings   | 1 of 1, cached: 1 ✔
    [7f/caf718] cowpy              | 1 of 1 ✔
    ```

The first three steps were cached from the previous run.
Only the new `cowpy` step actually ran.

Check the output in `results/2d-container/` for a text file containing the ASCII art.

### Takeaway

You understand how the `container` directive enables Nextflow to manage software dependencies automatically, ensuring reproducible execution without manual installation.

### What's next?

You've covered the fundamentals of Nextflow.
Move on to Part 2 to discover the nf-core community pipeline ecosystem.

---

## Summary

In this part you learned to:

- Run a Nextflow workflow and find its outputs
- Use `-resume` to skip completed steps
- Process multiple inputs from a CSV file using channels
- Connect processes in a multi-step pipeline
- Use containerized software with the `container` directive
