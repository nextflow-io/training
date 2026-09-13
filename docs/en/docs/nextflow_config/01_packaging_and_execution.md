# Part 1: Adapt to your compute environment

In [Nextflow Run](../nextflow_run/index.md), you configured a pipeline's inputs, parameters, and outputs.
This course covers the other half of the picture: adapting a pipeline's execution to whatever compute environment it happens to run on, without changing the workflow code.

!!! example "Scenario"

    You developed and tested your pipeline on your laptop using Docker.
    Now you need to hand it off: a collaborator only has Conda set up, and your institution's HPC cluster expects jobs to go through its own scheduler with its own resource limits.
    None of that should require rewriting the pipeline itself.

The same pipeline code can run in all of these places, because none of that is baked into the workflow.
Software packaging, execution platform, and resource allocation are all controlled through configuration, layered on top of the code, and that's what this course covers: how to adapt the same pipeline to a new environment by changing config, not code.

---

## 1. Select a software packaging technology

In [Nextflow Run](../nextflow_run/index.md), you saw a `conda` profile already set up in `nextflow.config` as an alternative to Docker.
Here you'll build that same switch yourself, and see what it takes to make a process actually usable with Conda.

### 1.1. Disable Docker and enable Conda

Switch `docker.enabled` to `false` and add a directive enabling Conda.

=== "After"

    ```groovy title="nextflow.config" linenums="1" hl_lines="1-2"
    docker.enabled = false
    conda.enabled = true
    ```

=== "Before"

    ```groovy title="nextflow.config" linenums="1"
    docker.enabled = true
    ```

This lets Nextflow create and use Conda environments for any process that has a Conda package specified.
The `cowpy` process doesn't have one yet, so let's add one, entirely from config.

### 1.2. Add a Conda package via config

A `conda` directive can be set in the process definition itself, the same way `container` already is in `modules/cowpy.nf`, but it doesn't have to be: `withName` lets you set it from config instead, scoped to just the `cowpy` process.

=== "After"

    ```groovy title="nextflow.config" linenums="6" hl_lines="3-5"
    process {
        memory = 1.GB
        withName: 'cowpy' {
            conda = 'conda-forge::cowpy==1.1.5'
        }
    }
    ```

=== "Before"

    ```groovy title="nextflow.config" linenums="6"
    process {
        memory = 1.GB
    }
    ```

This doesn't replace the `container` directive already in the pipeline code, it adds an alternative alongside it, without touching that code at all.

!!! tip

    The [Seqera Containers](https://seqera.io/containers/) search is a convenient way to look up the Conda package URI for a given tool, even if you're not planning to build a container from it.

### 1.3. Run the workflow to verify that it can use Conda

```bash
nextflow run main.nf --batch conda
```

??? success "Command output"

    ```console
    N E X T F L O W   ~  version 26.04.4

    Launching `main.nf` [extravagant_mccarthy] DSL2 - revision: c3c85dec78

    executor >  local (8)
    [c3/1a514c] sayHello (2)       | 3 of 3 ✔
    [14/615655] convertToUpper (3) | 3 of 3 ✔
    [78/4d519c] collectGreetings   | 1 of 1 ✔
    Creating env using conda: conda-forge::cowpy==1.1.5 [cache /workspaces/training/nextflow-config/work/conda/env-898314d566668b6587ad714ae06b8520]
    [0a/2e4f16] cowpy              | 1 of 1 ✔
    ```

This produces the same output as running with Docker, even though the mechanics are different behind the scenes: Nextflow retrieves the Conda package and builds an environment from it, instead of pulling a container image.

!!! info

    Building a new Conda environment can take a bit longer than pulling a container the first time around, but the package used here is small so it should be quick.

Now switch back to Docker for the rest of this course.

```groovy title="nextflow.config" linenums="1"
docker.enabled = true
```

??? tip "Mixing and matching Docker and Conda"

    Because these settings are scoped per process, you can mix and match: some processes use Docker, others use Conda, depending on what's available for each tool.
    If both a `container` directive (in the pipeline code) and a `conda` directive (here, from config) are set for the same process and both packaging systems are enabled, Nextflow prioritizes containers.

### Takeaway

You know how to configure which software packaging technology a process should use, and how to switch between Docker and Conda.

### What's next?

Learn how to change the execution platform Nextflow uses to actually run your tasks.

---

## 2. Select an execution platform

Every pipeline you've run so far has used the local executor: each task runs on the same machine as Nextflow itself.
Nextflow checks the available CPUs and memory, and holds tasks back until enough resources free up.

The local executor is convenient, but it doesn't scale past a single machine.
Nextflow supports [many other execution backends](https://nextflow.io/docs/latest/executor.html), including HPC schedulers (Slurm, LSF, SGE, PBS, and others) and cloud platforms (AWS Batch, Google Cloud Batch, Azure Batch, Kubernetes, and more).

### 2.1. Target a different backend

The executor is set by a process directive called `executor`.
By default it's `local`, so the following is implied:

```groovy title="Built-in configuration"
process {
    executor = 'local'
}
```

To target a different backend, set the directive to the executor you want.

```groovy title="nextflow.config"
process {
    executor = 'slurm'
}
```

!!! warning

    The training environment isn't connected to an HPC cluster, so this isn't something you can run here.

### 2.2. Backend-specific syntax is abstracted away

Most HPC platforms require job submissions to specify resource requests, such as CPUs, memory, and a queue name, using their own syntax.
The same request for 8 CPUs and 4 GB of RAM on a queue called `my-science-work` looks completely different depending on the scheduler.

??? abstract "Examples"

    ```bash title="Config for SLURM / submit using sbatch"
    #SBATCH -o /path/to/my/task/directory/my-task-1.log
    #SBATCH --no-requeue
    #SBATCH -c 8
    #SBATCH --mem 4096M
    #SBATCH -p my-science-work
    ```

    ```bash title="Config for PBS / submit using qsub"
    #PBS -o /path/to/my/task/directory/my-task-1.log
    #PBS -j oe
    #PBS -q my-science-work
    #PBS -l nodes=1:ppn=8
    #PBS -l mem=4gb
    ```

    ```bash title="Config for SGE / submit using qsub"
    #$ -o /path/to/my/task/directory/my-task-1.log
    #$ -j y
    #$ -terse
    #$ -notify
    #$ -q my-science-work
    #$ -l slots=8
    #$ -l h_rss=4096M,mem_free=4096M
    ```

Nextflow abstracts all of this away: you specify standardized properties such as `cpus`, `memory`, and `queue` once (see [process directives](https://nextflow.io/docs/latest/reference/process.html#process-directives) for the full list), and Nextflow translates them into the appropriate backend-specific scripts at runtime.

### 2.3. See what Nextflow actually runs

That translation isn't just a config-file convenience: it's backed by something concrete you can inspect right now, even with the local executor.
In [Nextflow Run, section 1.3](../nextflow_run/01_run_nextflow.md#13-explore-the-work-directory), you looked inside a task directory under `work/` and found `.command.sh`, the exact command Nextflow ran.
That same directory also contains a file you didn't look at yet: `.command.run`.

```bash
cat work/0a/0df4a1*/.command.run
```

??? success "Command output (excerpt)"

    ```console
    #!/bin/bash
    ### ---
    ### name: 'convertToUpper (3)'
    ### container: 'null'
    ### outputs:
    ### - 'UPPER-Bonjour-output.txt'
    ### ...
    set -e
    set -u
    ...
    nxf_launch() {
        /bin/bash -ue /workspaces/training/nextflow-run/work/0a/0df4a1028c2001758b1841cff92fc7/.command.sh
    }
    ...
    ```

`.command.run` is the real script Nextflow hands off for execution.
It wraps `.command.sh` with everything needed to actually run it: environment setup, input/output staging, and reporting the result back to Nextflow.
With the `local` executor, Nextflow simply runs this script on the same machine.

This is exactly what changes when you set a different `executor`.
For an HPC scheduler such as Slurm or PBS, Nextflow generates that same kind of wrapper script, adds the scheduler-specific header you saw in [2.2](#22-backend-specific-syntax-is-abstracted-away) (translated from your `cpus`, `memory`, and `queue` settings), and hands the result to that scheduler's own submission command, for example `sbatch` for Slurm.
From there, Nextflow polls the scheduler for job status instead of watching a local process directly.
Cloud batch backends work a little differently, since they're driven by API calls rather than a submission command, but the same underlying idea applies: the same task script runs, only how it gets launched and tracked changes.

### Takeaway

You know how to change the executor to target different compute infrastructure, that Nextflow abstracts away backend-specific submission syntax, and what actually happens behind the scenes when a task runs on a different backend.

### What's next?

Head on to [Part 2](./02_resources_and_retries.md), where you'll learn how to profile and allocate compute resources, and handle task failures with retries.

---

## Summary

In this part you learned to:

- Switch software packaging technology between Docker and Conda
- Add a `conda` directive to a process definition
- Change the execution platform with the `executor` directive
- Inspect what Nextflow actually generates and runs for a task, and how that changes across executors
