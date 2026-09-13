# Part 1: Adapt to your compute environment

In [Nextflow Run](../nextflow_run/index.md), you configured a pipeline's inputs, parameters, and outputs.
This course covers the other half of the picture: adapting a pipeline's execution to whatever compute environment it happens to run on, without changing the workflow code.

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
The `cowpy` process doesn't have one yet, so let's add it.

### 1.2. Specify a Conda package in the process definition

Add the `conda` directive to the `cowpy` process, alongside its existing `container` directive.

=== "After"

    ```groovy title="modules/cowpy.nf" linenums="1" hl_lines="4"
    // Generate ASCII art with cowpy (https://github.com/jeffbuttars/cowpy)
    process cowpy {

        container 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'
        conda 'conda-forge::cowpy==1.1.5'

        input:
    ```

=== "Before"

    ```groovy title="modules/cowpy.nf" linenums="1"
    // Generate ASCII art with cowpy (https://github.com/jeffbuttars/cowpy)
    process cowpy {

        container 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'

        input:
    ```

This doesn't replace the `container` directive, it adds an alternative alongside it.

!!! tip

    The [Seqera Containers](https://seqera.io/containers/) search is a convenient way to look up the Conda package URI for a given tool, even if you're not planning to build a container from it.

### 1.3. Run the workflow to verify that it can use Conda

```bash
nextflow run main.nf --batch conda
```

??? success "Command output"

    ```console
    N E X T F L O W   ~  version 26.04.4

    Launching `main.nf` [suspicious_hopper] DSL2 - revision: c3c85dec78

    executor >  local (8)
    [6b/a3cf85] sayHello (2)       | 3 of 3 ✔
    [2c/f5eba3] convertToUpper (3) | 3 of 3 ✔
    [e9/1285e8] collectGreetings   | 1 of 1 ✔
    [5e/7cc4f9] cowpy              | 1 of 1 ✔
    ```

This produces the same output as running with Docker, even though the mechanics are different behind the scenes: Nextflow retrieves the Conda package and builds an environment from it, instead of pulling a container image.

!!! info

    Building a new Conda environment can take a bit longer than pulling a container the first time around, but the package used here is small so it should be quick.

??? info "Mixing and matching Docker and Conda"

    Because these directives are set per process, you can mix and match: some processes use Docker, others use Conda, depending on what's available for each tool.
    If both a `container` and a `conda` directive are set for a process and both packaging systems are enabled, Nextflow prioritizes containers.

Now switch back to Docker for the rest of this course.

```groovy title="nextflow.config" linenums="1"
docker.enabled = true
```

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

### Takeaway

You know how to change the executor to target different compute infrastructure, and that Nextflow abstracts away backend-specific submission syntax.

### What's next?

Learn how to evaluate and set compute resource allocations.

---

## 3. Control compute resource allocations

By default, Nextflow allocates a single CPU to each process via the `cpus` directive, and does not impose a memory limit unless you set one:

```groovy title="Built-in configuration"
process {
    cpus = 1
}
```

You already know from [Nextflow Run](../nextflow_run/index.md) that this pipeline's configuration sets `memory` to 1 GB for all processes.
But how do you know what values to actually use for your own pipelines?

### 3.1. Generate a resource utilization report

You already generated an execution report with `-with-report` in [Nextflow Run](../nextflow_run/02_configure_pipeline.md).
That same report is how you find out how much CPU and memory your processes actually need: run the workflow with some default allocations, record actual usage, then adjust from there.

```bash
nextflow run main.nf -with-report report-config-1.html
```

The report is an HTML file you can open in a browser.
It breaks down runtime and resource utilization per process, including what percentage of the allocated resources was actually used.
Here's what it shows for `cowpy` with the current defaults (1 CPU, 1 GB memory):

| Metric           | Value  |
| ---------------- | ------ |
| CPU usage        | 116%   |
| Peak memory used | 6.4 MB |
| Allocated memory | 1 GB   |

`cowpy` uses well under 1% of its 1 GB allocation; the `%cpu` above 100% just means it briefly uses more than one CPU's worth of processing inside the container, in short bursts.

See [Reports](https://nextflow.io/docs/latest/reports.html) for the full list of available features.

### 3.2. Set resource allocations for a specific process

The report above shows `cowpy` comfortably within its current allocation, but say you wanted to give it more headroom anyway, for example because you expect larger inputs in production.
You can override the defaults for a single process with `withName`.

=== "After"

    ```groovy title="nextflow.config" linenums="6" hl_lines="3-6"
    process {
        memory = 1.GB
        withName: 'cowpy' {
            memory = 2.GB
            cpus = 2
        }
    }
    ```

=== "Before"

    ```groovy title="nextflow.config" linenums="6"
    process {
        memory = 1.GB
    }
    ```

With this in place, every process requests 1 GB of memory and a single CPU, except `cowpy`, which requests 2 GB and 2 CPUs.

!!! info

    If your machine has few CPUs and you allocate a high number per process, task calls may queue up behind each other, since Nextflow won't request more CPUs than are available.

Run it again with a different report filename, so you can compare before and after.

```bash
nextflow run main.nf -with-report report-config-2.html
```

??? success "Command output"

    ```console
    N E X T F L O W   ~  version 26.04.4

    Launching `main.nf` [furious_roentgen] DSL2 - revision: c3c85dec78

    executor >  local (8)
    [e8/1584eb] sayHello (2)       | 3 of 3 ✔
    [0f/b88085] convertToUpper (3) | 3 of 3 ✔
    [11/d81bfa] collectGreetings   | 1 of 1 ✔
    [a7/ca287f] cowpy              | 1 of 1 ✔
    ```

Comparing the two reports for `cowpy`:

| Metric           | Before (1 CPU, 1 GB) | After (2 CPUs, 2 GB) |
| ---------------- | -------------------- | -------------------- |
| Peak memory used | 6.4 MB               | 6.4 MB               |
| CPU usage        | 116%                 | 118%                 |

Doubling the allocation didn't change actual usage at all, which tells you the original 1 GB / 1 CPU was already generous for this toy workload.
On a real pipeline processing non-trivial data, you'd expect the numbers themselves to differ meaningfully between processes, which is exactly why you profile before deciding what to allocate, rather than guessing.

!!! tip

    Nextflow also has built-in [dynamic retry logic](https://nextflow.io/docs/latest/process.html#dynamic-task-resources) to automatically retry tasks that fail due to resource limits, with increased allocations.

### 3.3. Add resource limits

Depending on your compute infrastructure, there may be hard constraints on what you can request, for example a cluster-wide cap.
The `resourceLimits` directive lets you set those limits:

```groovy title="Syntax example"
process {
    resourceLimits = [
        memory: 750.GB,
        cpus: 200,
        time: 30.d
    ]
}
```

Nextflow translates these into whatever the target executor expects.
If a process requests more than the limit, the request gets capped rather than rejected.

!!! warning

    This isn't something you can run in the training environment, since it requires HPC infrastructure to have an effect.

??? info "Institutional reference configurations"

    The nf-core project maintains a [collection of configuration files](https://nf-co.re/configs/) shared by institutions worldwide, covering a wide range of HPC and cloud executors.
    They're a useful starting point whether or not your own institution is among them.

### Takeaway

You know how to generate a profiling report to assess resource utilization, override resource allocations for a specific process, and cap allocations with `resourceLimits`.

### What's next?

Head on to [Part 2](./02_profiles.md), where you'll learn how to bundle configuration like this into switchable profiles.

---

## Summary

In this part you learned to:

- Switch software packaging technology between Docker and Conda
- Add a `conda` directive to a process definition
- Change the execution platform with the `executor` directive
- Generate a resource profiling report and set per-process resource allocations
- Cap resource requests with `resourceLimits`
