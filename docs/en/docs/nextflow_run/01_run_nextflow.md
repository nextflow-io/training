# Part 1: Run Nextflow

In this part, we introduce the core concepts of running Nextflow pipelines.
We start with a simple Hello World workflow, then progress to a complete multi-step pipeline that processes multiple inputs in parallel using containers.

---

## 1. Hello World

The workflow `1-hello.nf` takes a greeting via a command-line argument and writes it to a file.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello_world.svg"
</figure>

### 1.1. Launch the workflow

Run the following command in your terminal.

```bash
nextflow run 1-hello.nf --input 'Hello World!'
```

??? success "Command output"

    ```console hl_lines="6"
    N E X T F L O W   ~  version 25.10.4

    Launching `1-hello.nf` [infallible_volhard] DSL2 - revision: 82d40dbf94

    executor >  local (1)
    [6d/740edd] sayHello | 1 of 1 ✔
    ```

The most important line in the output is the last one:

```console
[6d/740edd] sayHello | 1 of 1 ✔
```

This tells us that the `sayHello` process ran successfully once.
The `[6d/740edd]` prefix is a truncated path to the task's working directory — more on that below.

### 1.2. Find the output

This workflow is configured to publish its output to a `results` directory.
After running, you should find the output there:

```console title="results/"
results
└── 1-hello
    └── Hello World!-output.txt
```

Open the file to confirm it contains `Hello World!`.

### 1.3. Explore the `work/` directory

Behind the scenes, Nextflow creates a unique task directory for every process call inside a directory named `work/`.
The hash shown in the console output (`[6d/740edd]`) is the path to that directory.

```bash
ls work/6d/740edd*
```

Inside you will find the output file along with several hidden log files:

- **`.command.sh`**: the exact command Nextflow ran
- **`.command.out`** / **`.command.err`**: stdout and stderr from the process
- **`.command.log`**: combined log output
- **`.exitcode`**: the process exit code

The `.command.sh` file is especially useful when debugging — it shows precisely what was executed.

### 1.4. Optional: Code walkthrough

??? optional "Click to explore the code associated with this exercise"

    Let's open `1-hello.nf` and look at its main components.

    ```groovy title="1-hello.nf" linenums="1"
    #!/usr/bin/env nextflow

    include { sayHello } from './modules/sayHello.nf'

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

    We see the following:

    - an `include` statement pointing to a `process` module
    - a `params` block defining pipeline parameters
    - a `workflow` block describing the work to be done
    - an `output` block describing what to do with the outputs

    Let's take a look at each in turn.

    ### The `process` module

    The `include` statement tells Nextflow to load something called `sayHello` from a separate code file.

    ```groovy title="1-hello.nf" linenums="3"
    include { sayHello } from './modules/sayHello.nf'
    ```

    In that file, we find the definition for a process called `sayHello`:

    ```groovy title="modules/sayHello.nf" linenums="4"
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

    A **process** defines a single step in the pipeline.
    It declares its inputs, outputs, and the script to execute.
    The `val` qualifier means the input is a plain value (string, number, etc.).
    The `path` qualifier means the output is a file path.

    It is possible to write the process definition in the main workflow file, but keeping them in separate module files makes them reusable: the same module can be imported by multiple workflow scripts.

    ### The `params` block

    The `params` block declares the command-line parameters the workflow accepts:

    ```groovy title="1-hello.nf" linenums="8"
    params {
        input: String
    }
    ```

    Any parameter declared here becomes available on the command line with a double-dash (`--input`).
    Supported types include `String`, `Integer`, `Float`, `Boolean`, and `Path`.

    !!! tip

        Workflow parameters always use two dashes (`--input`) to distinguish them from Nextflow's own CLI flags, which use one dash (e.g. `-resume`).

    ### The `workflow` block

    The **workflow** block defines the dataflow logic: which processes to run and in what order.

    ```groovy title="1-hello.nf" linenums="12"
    workflow {

        main:
        // emit a greeting
        sayHello(params.input)

        publish:
        first_output = sayHello.out
    }
    ```

    Here there is only one process being called so it's very simple; we will cover more realistic examples later.

    The `main:` section calls the `sayHello` process with the `--input` value.
    The `publish:` section lists which outputs should be copied to the results directory.

    ### The `output` block

    The `output` block at the bottom of the file specifies the destination path and copy mode.

    ```groovy title="1-hello.nf" linenums="22"
    output {
        first_output {
            path '1-hello'
            mode 'copy'
        }
    }
    ```

    Each named entry corresponds to a `publish:` label in the workflow and maps it to a subdirectory under `results/`.

### Takeaway

You know how to run a Nextflow pipeline and find its outputs, and you know that the work is executed in task directories under `work/`.

### What's next?

Find out how Nextflow handles multiple inputs efficiently.

---

## 2. Process multiple inputs

Real-world pipelines typically process many pieces of data, not just one.
The workflow `2-inputs.nf` reads from a CSV file and runs `sayHello` once per row, in parallel.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-pipeline-multi-inputs-csv.svg"
</figure>

Let's run the workflow first, then we'll take a look at what mechanism Nextflow uses to handle these multiple inputs.

### 2.1. Run the workflow

Run the following command in your terminal.

```bash
nextflow run 2-inputs.nf --input data/greetings.csv
```

??? success "Command output"

    ```console hl_lines="6"
    N E X T F L O W   ~  version 25.10.4

    Launching `2-inputs.nf` [nauseous_babbage] DSL2 - revision: b90778224d

    executor >  local (3)
    [66/de7844] sayHello (3) | 3 of 3 ✔
    ```

The `3 of 3` tells us the `sayHello` process was called three times, once per row in the CSV.

In the `results` directory, you should now see three output files, one per greeting:

```console title="results/2-inputs/"
2-inputs
├── Bonjour-output.txt
├── Hello-output.txt
└── Hola-output.txt
```

Open any of the output files to confirm each one contains a greeting.

### 2.2. Run the workflow again with `-ansi-log false`

By default, Nextflow condenses the output to a single summary line per process.
To see each process call listed individually, add `-ansi-log false`:

```bash
nextflow run 2-inputs.nf --input data/greetings.csv -ansi-log false
```

??? success "Command output"

    ```console
    N E X T F L O W   ~  version 25.10.4

    Launching `2-inputs.nf` [extravagant_bardeen] DSL2 - revision: b90778224d
    [43/0bac1c] Submitted process > sayHello (1)
    [2d/99f604] Submitted process > sayHello (2)
    [6d/7578d7] Submitted process > sayHello (3)
    ```

This shows all three process calls and the unique work subdirectory created for each one.

### 2.3. Use `-resume` to skip completed work

Now switch to the extended input file, which adds two more greetings, and add `-resume` to the command line:

```bash
nextflow run 2-inputs.nf --input data/greetings-extended.csv -resume
```

??? success "Command output"

    ```console hl_lines="6"
    N E X T F L O W   ~  version 25.10.4

    Launching `2-inputs.nf` [adoring_mayer] DSL2 - revision: b90778224d

    executor >  local (2)
    [84/2f3067] sayHello (5) | 5 of 5, cached: 3 ✔
    ```

Nextflow ran only the two new inputs.
The three greetings processed in the previous run were cached and reused automatically.

This also works to skip executing processes for steps that have already been run successfully in a multi-step pipeline.
For example, if a pipeline run was interrupted by a system error, or if you added new steps to a pipeline in development.

The `-resume` capability is especially valuable in long pipelines where recovering from failure can save critical time and resources.

### 2.4. Optional: Code walkthrough

??? optional "Click to explore the code associated with this exercise"

    The key change in `2-inputs.nf` is in the `main:` section of the workflow:

    ```groovy title="2-inputs.nf" linenums="14" hl_lines="3 4 5"
        main:
        // create a channel for inputs from a CSV file
        greeting_ch = channel.fromPath(params.input)
                            .splitCsv()
                            .map { line -> line[0] }
        // emit a greeting
        sayHello(greeting_ch)
    ```

    What you see here is called a **channel**: a queue construct that handles input data in a way that makes it easy to parallelize operations.

    - `channel.fromPath(params.input)` creates a channel from the file path given with `--input`
    - `.splitCsv()` parses the CSV into rows
    - `#!groovy .map { line -> line[0] }` extracts the first column from each row

    The result is a channel containing `Hello`, `Bonjour`, and `Hola`.
    When passed to `sayHello(greeting_ch)`, Nextflow automatically calls the process once per item, running them in parallel when resources allow.

### Takeaway

You know how to process multiple inputs from a CSV file in parallel, and how to use `-resume` to avoid repeating completed work.

### What's next?

Learn how a complete multi-step pipeline chains processes together using channels, and how to use containers to manage analysis tools and their dependencies.

---

## 3. Run a multi-step pipeline

The workflow `main.nf` chains four processes into a complete pipeline.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello_pipeline_complete.svg"
</figure>

Each input greeting flows through all four steps: `sayHello` writes it to a file, `convertToUpper` converts the text to uppercase, `collectGreetings` merges all results into one file, and `cowpy` generates ASCII art from the merged output using a containerized tool.

Note that this workflow uses modules: each process is defined in its own file under `modules/`, and `main.nf` imports them with `include` statements instead of defining them inline.
This makes each process reusable across multiple workflows without duplicating code. To learn more, see the code exploration section further below.

### 3.1. Run the workflow

Run the following command in your terminal.

```bash
nextflow run main.nf --input data/greetings.csv --character turkey
```

??? success "Command output"

    ```console hl_lines="6 7 8 9"
    N E X T F L O W   ~  version 25.10.4

    Launching `main.nf` [reverent_jepsen] DSL2 - revision: c3c85dec78

    executor >  local (8)
    [87/35824b] sayHello (3)       | 3 of 3 ✔
    [0f/31c37e] convertToUpper (1) | 3 of 3 ✔
    [49/c25ca1] collectGreetings   | 1 of 1 ✔
    [a5/29e13f] cowpy              | 1 of 1 ✔
    ```

Four processes ran: `sayHello` and `convertToUpper` each ran once per input (3 of 3), `collectGreetings` gathered all outputs into one file (1 of 1), and `cowpy` generated the ASCII art (1 of 1).

Check `results/full_pipeline/` for the ASCII art file.

??? abstract "File contents"

    ```console title="results/full_pipeline/cowpy-COLLECTED-batch-output.txt"
     _________
    / HOLA    \
    | BONJOUR |
    \ HELLO   /
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

!!! note

    The `cowpy` step runs inside a Docker container rather than relying on software installed locally.
    A container packages an application together with everything it needs to run, so you don't have to install and manage dependencies yourself, and the pipeline behaves the same way on any machine that can run the container.
    Nextflow also supports Conda as an alternative to containers; see [Part 2](./02_configure_pipeline.md) for how to switch between them.

### 3.2. Optional: Code walkthrough

??? optional "Click to explore the code associated with this exercise"

    ## How data flows from one step to the next

    Each process passes its output channel to the next:

    ```groovy title="main.nf" linenums="19" hl_lines="7 8 9"
        main:
        // create a channel for inputs from a CSV file
        greeting_ch = channel.fromPath(params.input)
                            .splitCsv()
                            .map { line -> line[0] }
        sayHello(greeting_ch)
        convertToUpper(sayHello.out)
        collectGreetings(convertToUpper.out.collect(), params.batch)
        cowpy(collectGreetings.out.outfile, params.character)
    ```

    The pattern `processName.out` refers to a process's output channel.

    The `.collect()` operator gathers all individual outputs from `convertToUpper` into a single channel item before passing them to `collectGreetings`.

    ## Using process modules

    `main.nf` doesn't define any process code directly.
    Instead, it imports each process from its own file under `modules/`:

    ```groovy title="main.nf" linenums="3"
    include { sayHello } from './modules/sayHello.nf'
    include { convertToUpper } from './modules/convertToUpper.nf'
    include { collectGreetings } from './modules/collectGreetings.nf'
    include { cowpy } from './modules/cowpy.nf'
    ```

    Each module file contains a single process definition, structured the same way as the `sayHello` module in [1.4](#14-optional-code-walkthrough).
    Keeping processes in separate files makes them reusable across multiple workflows without duplicating code.

    <figure class="excalidraw">
    --8<-- "docs/en/docs/hello_nextflow/img/modules.svg"
    </figure>

    ## Using containerized software

    The `cowpy` process runs inside a Docker container specified in its module file:

    ```groovy title="modules/cowpy.nf" linenums="2" hl_lines="3"
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

    Nextflow automatically pulls the image, runs the script inside the container, and cleans up afterward.
    Docker is enabled for this project in `nextflow.config`:

    ```groovy title="nextflow.config"
    docker.enabled = true
    ```

    This single line enables Docker for any process in the pipeline that has a container specified.

### Takeaway

You've run a complete multi-step pipeline that processes multiple inputs in parallel using a containerized tool.

### What's next?

Head on to [Part 2](./02_configure_pipeline.md), where you'll learn how to configure pipeline behavior using `nextflow.config`.

---

## Summary

In this part you learned to:

- Run a Nextflow workflow and find its outputs
- Explore the `work/` directory and its log files
- Process multiple inputs from a CSV file in parallel
- Use `-resume` to skip completed work when adding new inputs
- Run a multi-step pipeline that uses a containerized tool
