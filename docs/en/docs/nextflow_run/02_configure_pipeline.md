# Part 2: Configure the pipeline

In [Part 1](./01_run_nextflow.md), you ran a complete multi-step pipeline that processes multiple inputs in parallel using containers.
Now we look at how to configure pipeline behavior using `nextflow.config`.

---

## 1. Configure the pipeline

Nextflow automatically picks up `nextflow.config` from the working directory and applies its settings to every run.

We provide you with a configuration file that covers four areas: software packaging, process settings, pipeline parameters and execution profiles.
We'll review each briefly then zoom in on two points of special interest: how to generate an execution report, and how to use a test profile.

### 1.1. The configuration file

The full file contains four configuration sections.

??? full-code "nextflow.config"

    ```groovy title="nextflow.config" linenums="1"
    /*
     * Software packaging
     */
    docker.enabled = true

    /*
     * Process settings
     */
    process {
        memory = 1.GB
        // withName: 'cowpy' {
        //     memory = 2.GB
        //     cpus = 2
        // }
    }

    /*
     * Pipeline parameters
     */
    params {
        input = 'data/greetings.csv'
        batch = 'batch'
        character = 'turkey'
    }

    /*
     * Profiles
     */
    profiles {
        test {
            params.input = 'data/greetings.csv'
            params.batch = 'test'
            params.character = 'tux'
        }
        conda {
            docker.enabled = false
            conda.enabled = true
        }
    }
    ```

**Software packaging** (`docker.enabled = true`): enables Docker for all processes.
Any process that declares a `container` directive runs inside the specified image.
To use Conda instead, activate the `conda` profile.

**Process settings**: sets a 1 GB memory limit for all processes.
The commented `withName: 'cowpy'` block is included as a syntax example showing how to apply separate limits to a specific process.
It is commented out because `cowpy` only appears in `main.nf`, and an unmatched selector produces a warning when running the other scripts.
Feel free to uncomment it, run `main.nf` again, and experiment with different resource values.

**Parameter defaults**: provide fallback values for parameters not supplied on the command line.
Running `nextflow run main.nf` with no flags uses these values.

**Profiles**: group settings that activate together when you pass `-profile <name>`.
The `test` profile overrides three parameters to run the pipeline with a small, well-defined input set.
The `conda` profile switches software packaging from Docker to Conda.

!!! note

    This config covers local execution on a single machine.
    Nextflow also supports HPC schedulers (SLURM, PBS, LSF) and cloud executors (AWS Batch, Google Cloud Batch, Azure Batch), all configured through the same `nextflow.config` mechanism.
    See the archived [Nextflow Run: Configuration](../archive/nextflow_run/03_config.md) lesson for a full walkthrough of these options.

### 1.2. Generate an execution report

Add `-with-report` to any `nextflow run` command to generate an HTML report after the pipeline completes:

```bash
nextflow run main.nf --input data/greetings.csv --character turkey -with-report
```

??? success "Command output"

    ```console hl_lines="6 7 8 9"
    N E X T F L O W   ~  version 26.04.4

    Launching `main.nf` [lonely_aryabhata] DSL2 - revision: c3c85dec78

    executor >  local (8)
    [9e/f3fa5e] sayHello (3)       | 3 of 3 ✔
    [08/4a93d7] convertToUpper (2) | 3 of 3 ✔
    [c5/c8595a] collectGreetings   | 1 of 1 ✔
    [2a/1e5abe] cowpy              | 1 of 1 ✔
    ```

Nextflow writes the report to a file named `report-<timestamp>.html` in the working directory.
Open it in a browser to see an execution summary, a table of every task with its status and runtime, and resource usage charts broken down by process.

The report is especially useful when a pipeline takes longer than expected or a task fails — the task table shows exactly where time was spent and which tasks succeeded or failed.

### 1.3. Run with a profile

The `test` profile is a standard nf-core convention. Every nf-core pipeline ships with one for quick validation.

Run the pipeline using the `test` profile:

```bash
nextflow run main.nf -profile test
```

??? success "Command output"

    ```console hl_lines="6 7 8 9"
    N E X T F L O W   ~  version 26.04.4

    Launching `main.nf` [silly_goodall] DSL2 - revision: c3c85dec78

    executor >  local (8)
    [06/9614da] sayHello (2)       | 3 of 3 ✔
    [8d/5f1b8a] convertToUpper (3) | 3 of 3 ✔
    [1e/cd52db] collectGreetings   | 1 of 1 ✔
    [2c/5f41d8] cowpy              | 1 of 1 ✔
    ```

The pipeline runs with `batch = 'test'` and `character = 'tux'`.
Check `results/full_pipeline/`. The batch name appears in the output file names, and the ASCII art features the tux penguin instead of a turkey.

To inspect the fully resolved settings for any profile combination, run `nextflow config` with the relevant profiles:

```bash
nextflow config -profile test
```

??? success "Command output"

    ```console hl_lines="3 4"
    params {
       input = 'data/greetings.csv'
       batch = 'test'
       character = 'tux'
    }

    docker {
       enabled = true
    }

    process {
       memory = '1 GB'
    }
    ```

This shows the merged result: the `test` profile's overrides applied on top of the base config.
Use it to confirm what settings will be active before running a pipeline.

### Takeaway

You know how to configure pipeline behavior using `nextflow.config`, how to generate an HTML execution report with `-with-report`, and how profiles bundle related settings that activate together with a single flag.

### What's next?

Head on to [Part 3](./03_remote_repositories.md), where you'll learn how to run pipelines directly from remote repositories such as GitHub.

---

## Summary

In this part you learned to:

- Generate an HTML execution report with `-with-report`
- Configure pipeline behavior using `nextflow.config` and profiles
