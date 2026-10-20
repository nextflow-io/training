# Part 6: Scale up and tune compute resources

In [Part 3](./03_run_production_pipeline.md), you got `nf-core/rnaseq` running on the default machine by constraining its resource requests with `laptop.config`.
That is fine for testing, but real analyses need more resources, not fewer.

In this part, you move to a bigger machine and run the pipeline with its default resource requests.
You then read the execution report from that run to see how much of those resources each process actually used, write a configuration that matches requests to usage, and see how that affects the run.

This part uses a bigger codespace machine, which uses more of your GitHub Codespaces quota, as explained below.

---

## 1. Run on a bigger machine

### 1.1. Change the machine type

You can change the machine type of your codespace from the editor:

1. Open the Command Palette with ++ctrl+shift+p++ (++cmd+shift+p++ on macOS).
2. Type `Codespaces: Change Machine Type` and select it.
3. Select the **8-core** machine type.
4. Confirm that you want to restart the codespace when prompted.

Alternatively, go to [github.com/codespaces](https://github.com/codespaces), click the **...** menu next to your codespace, and select **Change machine type**.

The change takes effect when the codespace restarts, which can take a few minutes.
Your files in `/workspaces` are preserved, but anything running in the terminal is stopped.

The machine types available for this repository include the following (a 16-core machine with 64 GB is also available):

| Machine type | CPUs | Memory | Quota use relative to the default |
| ------------ | ---- | ------ | --------------------------------- |
| 2-core       | 2    | 8 GB   | 1x (default)                      |
| 4-core       | 4    | 16 GB  | 2x                                |
| 8-core       | 8    | 32 GB  | 4x                                |

Your quota is counted in core-hours, so a bigger machine uses it up faster: the free allowance of 120 core-hours per month gives you about 60 hours on the default machine, but only 15 hours on an 8-core machine.
The rest of this part uses the 8-core machine.

### 1.2. Check the new resources

Once the codespace has restarted, open a terminal and return to the working directory:

```bash
cd /workspaces/training/nfcore-use
```

Check the number of CPUs and the amount of memory available:

```bash
nproc
free -g
```

```console title="Command output"
8
               total        used        free      shared  buff/cache   available
Mem:              31           2           0           0          28          28
```

### 1.3. Run the pipeline with its default resource requests

Run the pipeline without `laptop.config`, so that each process gets the resources the pipeline requests by default:

```bash
nextflow run nf-core/rnaseq -r 3.27.0 -profile test,docker --outdir rnaseq-default
```

??? success "Command output (excerpt)"

    ```console
    ...
    -[nf-core/rnaseq] Pipeline completed successfully-
    ```

With 8 CPUs and 32 GB, the default requests fit, so the run proceeds without any overrides.
In our tests on a newly created 8-core codespace, the run took about 14 minutes.
Much of that time goes to downloading the container image for each tool the first time it is used; once the images are downloaded, a full run takes about 6.5 minutes.

Four times as many CPUs as the default machine does not make the run four times faster, for two reasons.
The `test` profile limits each task to 4 CPUs and 15 GB of memory, and the test dataset is so small that many tasks finish in seconds regardless of the resources they get.
On a real dataset, the long-running tasks such as read alignment benefit much more from additional resources.

Real RNA-seq experiments typically involve dozens of samples and run for hours or days.
Nextflow supports HPC schedulers (SLURM, PBS, LSF) and cloud platforms (AWS, Google Cloud, Azure), which can dramatically reduce wall-clock time by distributing work across many machines.
Setting up those environments, however, adds significant complexity.

The Seqera platform (developed by the creators of Nextflow) provides a web-based interface for launching Nextflow pipelines on HPC or cloud infrastructure (either your own or one managed for you), with compute and data management capabilities that streamline the process of running pipelines at scale.

??? tip "Seqera academic program"

    Academic researchers can access Seqera Platform free of charge through the [Seqera academic program](https://seqera.io/academic-program/).

### Takeaway

You can give a pipeline more resources by running it on a bigger machine, and local execution is a starting point rather than a destination for real-scale analysis.

### What's next?

Check how much of those resources the pipeline actually used.

---

## 2. Read the execution report

### 2.1. Find the reports

nf-core pipelines generate execution reports automatically; you don't need to add `-with-report`.
They are written to the `pipeline_info/` directory inside the output directory:

```bash
ls rnaseq-default/pipeline_info/
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

There is one set of files per run, identified by its timestamp.

### 2.2. Compare usage with allocation

Open the execution report in your browser.
In the **Resource Usage** section, the **Memory** plot has two tabs:

- **Raw usage** shows the peak memory each process used, in absolute terms.
- **% Allocated** shows that peak as a percentage of the memory requested.

<!-- TODO: screenshot of the % Allocated memory plot from a Codespaces run -->

With the pipeline's default requests, most processes use a small fraction of their allocation.
Here are a few representative processes from our own test run on the 8-core machine:

| Process                 | Label            | Requested     | Peak usage |
| ----------------------- | ---------------- | ------------- | ---------- |
| `FASTQC`                | `process_low`    | 2 CPUs, 12 GB | 0.5 GB     |
| `STAR_ALIGN`            | `process_high`   | 4 CPUs, 15 GB | 0.9 GB     |
| `BBMAP_BBSPLIT`         | `process_high`   | 4 CPUs, 15 GB | 1.3 GB     |
| `PICARD_MARKDUPLICATES` | `process_medium` | 4 CPUs, 15 GB | 4.2 GB     |
| `MULTIQC`               | `process_single` | 1 CPU, 6 GB   | 1.0 GB     |

The requests of 15 GB and 4 CPUs are lower than what `conf/base.config` specifies for these labels because the `test` profile caps every request at those values.
Your numbers will differ somewhat, but the pattern is the same: the requests are far higher than the usage.

### 2.3. Understand the cost of over-allocation

Over-allocating memory does not make a task run any faster.
It does limit how many tasks can run at the same time.

When Nextflow runs tasks locally, it only starts a task if enough CPUs and memory remain unreserved on the machine.
In our test run, the tasks running at the same time reserved up to 30 GB of the 32 GB available, so no more than 4 tasks ever ran at once, and 2.2 on average, even though the machine has 8 CPUs.
Memory, not CPU, was the bottleneck.

Open the execution timeline to see this.
Each bar is a task; many tasks wait for memory to be released before they can start, even when CPUs are idle.

<!-- TODO: screenshot of the timeline from the default 8-core run -->

The same problem appears on HPC clusters and in the cloud in a different form: oversized requests wait longer in the scheduler queue and cost more.

### Takeaway

The execution report shows how much of its allocation each process actually used.
Over-allocation limits how many tasks can run in parallel.

### What's next?

Write a configuration that matches resource requests to actual usage.

---

## 3. Right-size the resource allocations

### 3.1. Set a ceiling with `resourceLimits`

`resourceLimits` caps every resource request at the values you give it, whatever the pipeline asks for.
A request above the cap is reduced to the cap instead of failing.

The `test` profile of `nf-core/rnaseq` sets `resourceLimits` to 4 CPUs and 15 GB of memory.
That is why the run in Part 3 failed on the default machine, which has less than 8 GB: a 12 GB request is within the profile's limits, but not within the machine's.

The 8-core machine provides 8 CPUs and 32 GB of memory.
Create a file named `tuned.config` that sets the limits to what it can provide:

```groovy title="tuned.config" linenums="1"
process {
    // Never request more than the machine has
    resourceLimits = [cpus: 8, memory: 30.GB]
}
```

On its own, this does not help parallelism: every process still requests far more memory than it uses.

### 3.2. Right-size each label

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

    ```groovy title="tuned.config" linenums="1" hl_lines="5-24"
    --8<-- "nfcore-use/solutions/tuned.config"
    ```

=== "Before"

    ```groovy title="tuned.config" linenums="1"
    process {
        // Never request more than the machine has
        resourceLimits = [cpus: 8, memory: 30.GB]
    }
    ```

Each value leaves headroom above the peak usage in the report, and the `task.attempt` multiplier covers inputs that turn out larger than expected.
`process_single` gets 2 GB rather than 1 GB because `MULTIQC`, which has that label, peaked at about 1 GB.

`PICARD_MARKDUPLICATES` uses more memory than the other `process_medium` processes, so it gets its own `withName` selector.
A `withName` selector takes precedence over a `withLabel` selector, so `PICARD_MARKDUPLICATES` gets 6 GB while the other `process_medium` processes get 2 GB.
The simple name is enough here because, unlike `TRIMGALORE` in Part 3, the pipeline does not set resources for `PICARD_MARKDUPLICATES` by name.

### 3.3. Run the pipeline with the tuned configuration

Run the test profile again with `tuned.config` and a new output directory:

```bash
nextflow run nf-core/rnaseq -r 3.27.0 -profile test,docker -c tuned.config --outdir rnaseq-tuned
```

Note the run time, then open the new execution timeline in `rnaseq-tuned/pipeline_info/`.
The run in section 1.3 also had to download the container images, so its run time is not directly comparable.
Tasks now overlap more, because several of them fit in memory at once.

In our tests on the 8-core machine, with the container images already downloaded, the run took about 5.5 minutes with `tuned.config` compared to about 6.5 minutes with the default requests, and no task ran out of memory.
At most 6 tasks ran at once instead of 4, and the memory reserved at any one time peaked at 16 GB instead of 30 GB.

The time saved is modest here because the test dataset is tiny: most tasks finish within seconds, so starting them sooner gains little.
On real data, where individual tasks run for minutes or hours, running more of them in parallel makes a much larger difference.

??? tip "Applying this to a real analysis"

    The test dataset is tiny, so its resource usage does not represent a real analysis.
    For real data, run the pipeline on one or two representative samples first, read the execution report, then tune the configuration before launching the full dataset.
    Memory needs often scale with input size and reference genome size, so keep generous headroom.

??? exercise "Apply it to the pipeline you chose in Part 4"

    Return to the published results of the pipeline you chose in [Part 4](./04_choose_pipeline.md) and open the execution report in `pipeline_info/`.
    Identify the processes with the highest memory usage, and compare that usage with the memory requested by their labels in the pipeline's `conf/base.config` file.

    Which processes would you right-size for your dataset, and in which direction?

### 3.4. Switch back to the default machine

When you're done, switch your codespace back to the 2-core machine type, following the same steps as in section 1.1, to avoid using up your quota.

### Takeaway

`resourceLimits` keeps requests within what the machine can provide, and label-level overrides based on the execution report let more tasks run in parallel.
Keeping the `task.attempt` multiplier preserves the pipeline's automatic retry with more resources.

### What's next?

You've completed the Use nf-core course.
See [Next steps](next_steps.md) for where to go from here.

---

## Summary

In this part you learned to:

- Run a pipeline with more compute resources by changing the machine type
- Find and read the execution report and timeline generated by nf-core pipelines
- Recognize over-allocation and its effect on parallel execution
- Cap resource requests with `resourceLimits` and right-size them by label and by process name, while preserving automatic retries
