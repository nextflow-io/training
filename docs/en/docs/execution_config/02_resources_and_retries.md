# Part 2: Manage compute resources and failures

In [Part 1](./01_packaging_and_execution.md), you adapted where and how a pipeline's tasks run.
Here you'll adapt how much compute each task gets, and what happens when a task fails despite your best guess at an allocation.

---

## 1. Control compute resource allocations

By default, Nextflow allocates a single CPU to each process via the `cpus` directive, and does not impose a memory limit unless you set one:

```groovy title="Built-in configuration"
process {
    cpus = 1
}
```

You already know from [Nextflow Run](../nextflow_run/index.md) that this pipeline's configuration sets `memory` to 1 GB for all processes.
But how do you know what values to actually use for your own pipelines?

### 1.1. Generate a resource utilization report

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

### 1.2. Set resource allocations for a specific process

The report above shows `cowpy` comfortably within its current allocation, but say you wanted to give it more headroom anyway, for example because you expect larger inputs in production.
You can override the defaults for a single process with `withName`.

=== "After"

    ```groovy title="nextflow.config" linenums="6" hl_lines="5-6"
    process {
        memory = 1.GB
        withName: 'cowpy' {
            conda = 'conda-forge::cowpy==1.1.5'
            memory = 2.GB
            cpus = 2
        }
    }
    ```

=== "Before"

    ```groovy title="nextflow.config" linenums="6"
    process {
        memory = 1.GB
        withName: 'cowpy' {
            conda = 'conda-forge::cowpy==1.1.5'
        }
    }
    ```

With this in place, every process requests 1 GB of memory and a single CPU, except `cowpy`, which requests 2 GB and 2 CPUs (on top of the `conda` setting from [Part 1](./01_packaging_and_execution.md)).

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

### 1.3. Add resource limits

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

Learn how to make a pipeline recover automatically when a task fails, whether or not your resource allocation guess was right.

---

## 2. Handle task failures with retries

Profiling tells you what a process needs most of the time, but real workloads vary: an allocation that's comfortable for most inputs can still be too tight for an unusually large one, and guesses can simply be wrong.
Rather than letting a single failed task bring down the whole run, Nextflow can retry a failed task automatically, optionally giving it more resources on each attempt.

### 2.1. Retry a failed task automatically

To see this in action, deliberately set `cowpy`'s memory allocation below what it actually needs: recall from [1.1](#11-generate-a-resource-utilization-report) that it peaks at around 6.4 MB, so 6 MB should be just short of enough.

=== "After"

    ```groovy title="nextflow.config" linenums="6" hl_lines="5-7"
    process {
        memory = 1.GB
        withName: 'cowpy' {
            conda = 'conda-forge::cowpy==1.1.5'
            memory = 6.MB
            errorStrategy = 'retry'
            maxRetries = 2
        }
    }
    ```

=== "Before"

    ```groovy title="nextflow.config" linenums="6"
    process {
        memory = 1.GB
        withName: 'cowpy' {
            conda = 'conda-forge::cowpy==1.1.5'
            memory = 2.GB
            cpus = 2
        }
    }
    ```

`errorStrategy` tells Nextflow what to do when a task fails: `'retry'` resubmits the task instead of stopping the whole pipeline.
`maxRetries` caps how many extra attempts it gets before Nextflow gives up.

```bash
nextflow run main.nf
```

??? failure "Command output (abridged)"

    ```console
    [PROCESS 20/24eec4] cowpy
    [ERROR] cowpy
    exit: 137
    cmd: cat COLLECTED-batch-output.txt | cowpy -c "turkey" > cowpy-COLLECTED-batch-output.txt
    workdir: .../work/20/24eec4...
    [PROCESS 7f/5814b3] cowpy
    [ERROR] cowpy
    exit: 1
    cmd: cat COLLECTED-batch-output.txt | cowpy -c "turkey" > cowpy-COLLECTED-batch-output.txt
    workdir: .../work/7f/5814b3...
    [PROCESS 12/4c9601] cowpy
    [ERROR] ERROR ~ Error executing process > 'cowpy'

    Caused by:
      Process `cowpy` terminated with an error exit status (137)

    Command exit status:
      137

    Work dir:
      .../work/12/4c9601...

    Tip: when you have fixed the problem you can continue the execution adding the option `-resume` to the run command line

    [FAILED] completed=10 failed=3 cached=0
    ```

Exit code 137 is the standard signal for an out-of-memory kill: the container didn't have enough memory to run `cowpy` at all.
Nextflow retried the task twice, three attempts in total, matching `maxRetries = 2`.
Since the memory allocation never changed between attempts, every attempt hit the same wall; once retries are exhausted, Nextflow reports the failure in full and stops the pipeline, exiting with a non-zero status.

Retrying on its own doesn't fix anything if the underlying cause doesn't change between attempts.

### 2.2. Increase resources on each retry

Inside a process directive, `task.attempt` holds the current attempt number, starting at 1.
You can use it in a closure to scale a resource allocation up with each retry.

=== "After"

    ```groovy title="nextflow.config" linenums="6" hl_lines="5 7"
    process {
        memory = 1.GB
        withName: 'cowpy' {
            conda = 'conda-forge::cowpy==1.1.5'
            memory = { 6.MB * task.attempt }
            errorStrategy = 'retry'
            maxRetries = 3
        }
    }
    ```

=== "Before"

    ```groovy title="nextflow.config" linenums="6"
    process {
        memory = 1.GB
        withName: 'cowpy' {
            conda = 'conda-forge::cowpy==1.1.5'
            memory = 6.MB
            errorStrategy = 'retry'
            maxRetries = 2
        }
    }
    ```

Run the workflow again:

```bash
nextflow run main.nf
```

??? success "Command output (abridged)"

    ```console
    [PROCESS 99/b9c7f4] cowpy
    [ERROR] cowpy
    exit: 137
    cmd: cat COLLECTED-batch-output.txt | cowpy -c "turkey" > cowpy-COLLECTED-batch-output.txt
    workdir: .../work/99/b9c7f4...
    [PROCESS d6/d3627d] cowpy

    Outputs:

      ...
      cowpy_art: full_pipeline/cowpy-COLLECTED-batch-output.txt

    [FAILED] completed=9 failed=1 cached=0
    ```

The first attempt still fails at 6 MB, but the retry runs with 12 MB (`6.MB * 2`) and succeeds, and the pipeline completes with all outputs published.

!!! warning

    The console summary tag above still reads `[FAILED]`, even though the pipeline as a whole succeeded: that tag reflects individual task attempts, not overall outcome, and one attempt did fail along the way.
    Check for the `Outputs:` listing, or the command's exit status, to see whether the run actually succeeded.

See [Dynamic computing resources](https://nextflow.io/docs/latest/process.html#dynamic-task-resources) in the Nextflow documentation for more advanced retry patterns, including scaling based on which specific error occurred.

### Takeaway

You know how to make a pipeline automatically retry failed tasks, and how to scale resource allocations with each retry using `task.attempt`.

### What's next?

Head on to [Part 3](./03_profiles.md), where you'll learn how to bundle configuration like this into switchable profiles.

---

## Summary

In this part you learned to:

- Generate a resource profiling report and set per-process resource allocations
- Cap resource requests with `resourceLimits`
- Automatically retry a failed task with `errorStrategy` and `maxRetries`
- Scale a resource allocation up with each retry using `task.attempt`
