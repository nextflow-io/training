# Part 7: Tune compute resources

In [Part 3](./03_run_production_pipeline.md), the `nf-core/rnaseq` test run failed because a process requested more memory than the machine had.
You fixed it with `laptop.config`, which sets every resource label to 6 GB of memory.

That fix got the pipeline running, but it was based on what the machine could provide, not on what each process needs.
nf-core pipelines record how much CPU and memory every task actually used, so you can do better.

In this part, you read the execution report from your Part 3 run, use it to write a configuration file that matches resource requests to actual usage, and see how that affects the run.

---

## 1. Read the execution report

### 1.1. Find the reports

nf-core pipelines generate execution reports automatically; you don't need to add `-with-report`.
They are written to the `pipeline_info/` directory inside the output directory:

```bash
ls rnaseq-results/pipeline_info/
```

```console title="Directory contents"
execution_report_<timestamp>.html
execution_timeline_<timestamp>.html
execution_trace_<timestamp>.txt
nf_core_rnaseq_software_mqc_versions.yml
params_<timestamp>.json
pipeline_dag_<timestamp>.html
```

- **The execution report** summarizes resource usage per process.
- **The execution timeline** shows when each task started and ended.
- **The execution trace** is a tab-separated table with one row per task, convenient for further analysis.

If you ran the pipeline more than once, there is one set of files per run, identified by its timestamp.

### 1.2. Compare usage with allocation

Open the execution report in your browser.
In the **Resource Usage** section, the **Memory** plot has two tabs:

- **Raw usage** shows the peak memory each process used, in absolute terms.
- **% Allocated** shows that peak as a percentage of the memory requested.

<!-- TODO: screenshot of the % Allocated memory plot from a Codespaces run -->

With `laptop.config`, most processes use a small fraction of their allocation.
Here are a few representative processes:

<!-- TODO: replace with values from a Codespaces run -->

| Process                 | Label            | Requested | Peak usage |
| ----------------------- | ---------------- | --------- | ---------- |
| `FASTQC`                | `process_low`    | 6 GB      | ~0.3 GB    |
| `STAR_ALIGN`            | `process_high`   | 6 GB      | ~0.9 GB    |
| `PICARD_MARKDUPLICATES` | `process_medium` | 6 GB      | ~1.8 GB    |

Your numbers will differ somewhat, but the pattern is the same: the requests are far higher than the usage.

### 1.3. Understand the cost of over-allocation

Over-allocating memory does not make a task run any faster.
It does limit how many tasks can run at the same time.

When Nextflow runs tasks locally, it only starts a task if enough CPUs and memory remain unreserved on the machine.
On a machine with about 8 GB of memory, a task that requests 6 GB leaves no room for a second one, so the pipeline runs one task at a time even when the CPUs are idle.

Open the execution timeline to see this.
Each bar is a task; with `laptop.config`, the bars form a mostly sequential staircase rather than overlapping stacks.

<!-- TODO: screenshot of the laptop.config timeline -->

The same problem appears on HPC clusters and in the cloud in a different form: oversized requests wait longer in the scheduler queue and cost more.

### Takeaway

The execution report shows how much of its allocation each process actually used.
Over-allocation limits how many tasks can run in parallel.

### What's next?

Write a configuration that matches resource requests to actual usage.

---

## 2. Right-size the resource allocations

### 2.1. Check what the machine provides

Find out how many CPUs and how much memory your machine has:

```bash
nproc
free -g
```

<!-- TODO: confirm Codespaces values -->

The training environment provides 4 CPUs and about 8 GB of memory.

### 2.2. Set a ceiling with `resourceLimits`

`resourceLimits` caps every resource request at the values you give it, whatever the pipeline asks for.
A request above the cap is reduced to the cap instead of failing.

The `test` profile of `nf-core/rnaseq` already sets `resourceLimits`, but to 15 GB of memory, which is more than the training environment has.
That is why the test run in Part 3 failed.

Create a file named `tuned.config` with the following content:

```groovy title="tuned.config" linenums="1"
process {
    // Never request more than the machine has
    resourceLimits = [cpus: 4, memory: 7.GB]
}
```

This alone would have prevented the failure in Part 3, but every process would still request 7 GB, so tasks would still run one at a time.

### 2.3. Right-size each label

The pipeline assigns resources through labels defined in `conf/base.config`.
Each label definition multiplies its values by `task.attempt`:

```groovy title="conf/base.config (excerpt)"
withLabel:process_low {
    cpus   = { 2     * task.attempt }
    memory = { 12.GB * task.attempt }
    time   = { 4.h   * task.attempt }
}
```

`task.attempt` is 1 on the first attempt and increases each time a task is retried.
nf-core pipelines retry tasks that fail with an out-of-memory exit code, so a task that runs out of memory is resubmitted automatically with double the memory.

`laptop.config` replaces these with fixed values such as `memory = 6.GB`, which disables that safety net: a retried task gets the same 6 GB as the first attempt.
When you override resource values, keep the `task.attempt` multiplier.

Based on the execution report, update `tuned.config`:

=== "After"

    ```groovy title="tuned.config" linenums="1" hl_lines="5-21"
    --8<-- "nfcore-use/solutions/tuned.config"
    ```

=== "Before"

    ```groovy title="tuned.config" linenums="1"
    process {
        // Never request more than the machine has
        resourceLimits = [cpus: 4, memory: 7.GB]
    }
    ```

Each value leaves headroom above the peak usage in the report, and the `task.attempt` multiplier covers inputs that turn out larger than expected.

`PICARD_MARKDUPLICATES` uses more memory than the other `process_medium` processes, so it gets its own `withName` selector.
A `withName` selector takes precedence over a `withLabel` selector, so `PICARD_MARKDUPLICATES` gets 3 GB while the other `process_medium` processes get 2 GB.

### 2.4. Run the pipeline with the tuned configuration

Run the test profile again, using `tuned.config` instead of `laptop.config`, with a new output directory:

```bash
nextflow run nf-core/rnaseq -profile test,docker -c tuned.config --outdir rnaseq-tuned
```

Compare the run time with the Part 3 run, then open the new execution timeline in `rnaseq-tuned/pipeline_info/`.
Tasks now overlap, because several of them fit in memory at once.

In our tests on a machine with 4 CPUs and 8 GB of memory, the run with `tuned.config` took about 30% less time than the run with `laptop.config`, with no task running out of memory.
The gain is limited by the 4 CPUs, since the processes with the `process_high` label each take all of them.

??? tip "Applying this to a real analysis"

    The test dataset is tiny, so its resource usage does not represent a real analysis.
    For real data, run the pipeline on one or two representative samples first, read the execution report, then tune the configuration before launching the full dataset.
    Memory needs often scale with input size and reference genome size, so keep generous headroom.

??? exercise "Apply it to the pipeline you chose in Part 4"

    Return to the published results of the pipeline you chose in [Part 4](./04_choose_pipeline.md) and open the execution report in `pipeline_info/`.
    Identify the processes with the highest memory usage, and compare that usage with the memory requested by their labels in the pipeline's `conf/base.config` file.

    Which processes would you right-size for your dataset, and in which direction?

### Takeaway

`resourceLimits` keeps requests within what the machine can provide, and label-level overrides based on the execution report let more tasks run in parallel.
Keeping the `task.attempt` multiplier preserves the pipeline's automatic retry with more resources.

### What's next?

You've completed the Use nf-core course.
See [Next steps](next_steps.md) for where to go from here.

---

## Summary

In this part you learned to:

- Find and read the execution report and timeline generated by nf-core pipelines
- Recognize over-allocation and its effect on parallel execution
- Cap resource requests with `resourceLimits` and right-size them by label and by process name, while preserving automatic retries
