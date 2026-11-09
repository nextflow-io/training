# Part 3: Run a production pipeline

In [Part 2](./02_configure_execution.md), you learned how to set parameters and customize configuration for nf-core/demo.
Now we apply what you've learned to a real production pipeline, nf-core/rnaseq.

You first get it running on the default machine by constraining its resource requests.
You then move to a bigger machine, run it with its default requests, and use its execution report to tune them.

---

## 1. Pull and run nf-core/rnaseq

So far we have used `nf-core/demo`, which is a minimal pipeline designed for training.
Now we pull a real production pipeline and run it with its test profile.

The `nf-core/rnaseq` pipeline performs the core steps of bulk RNA sequencing analysis: quality control, adapter trimming, read alignment, and gene-level quantification.
It is probably the most widely used nf-core pipeline to date.

### 1.1. Pull the pipeline

Run the following command to download it:

```bash
nextflow pull nf-core/rnaseq -r 3.27.0
```

??? success "Command output"

    ```console
    Checking nf-core/rnaseq ...
     downloaded from https://github.com/nf-core/rnaseq.git - revision: a1fcdddd3b [3.27.0]
    ```

The pipeline is now cached locally and ready to run.

As mentioned in [Part 1](./01_run_demo.md#22-run-the-pipeline), the `-r` option pins the pipeline to a specific version, here release 3.27.0.
When you ran `nf-core/demo` without `-r`, Nextflow used whatever version of the pipeline's default branch it had downloaded, which changes as the pipeline is updated.

Pipelines change from one release to the next: parameters, default resource requests, error messages and even results can differ.
For any real analysis, always pin the version with `-r`, so that your results are reproducible and you can report exactly which version you used.
In this course, we pin `nf-core/rnaseq` to the release the materials were tested with, so that what you see in your terminal matches this page.

### 1.2. Run the test profile

Run it with the test profile and Docker:

```bash
nextflow run nf-core/rnaseq -r 3.27.0 -profile test,docker --outdir rnaseq-results
```

??? failure "Command output"

    ```console
     N E X T F L O W   ~  version 26.04.4

    Launching `https://github.com/nf-core/rnaseq` [suspicious_dijkstra] revision: a1fcdddd3b [3.27.0]

    ------------------------------------------------------
                                            ,--./,-.
            ___     __   __   __   ___     /,-._.--~'
      |\ | |__  __ /  ` /  \ |__) |__         }  {
      | \| |       \__, \__/ |  \ |___     \`-._,-`-,
                                            `._,._,'
      nf-core/rnaseq 3.27.0
    ------------------------------------------------------
    ...
    [-        ] NFCORE_RNASEQ:RNASEQ:FASTQ_QC_TRIM_FILTER_SETSTRANDEDNESS:FQ_LINT  -
    Plus 47 more processes waiting for tasks…

    Execution cancelled -- Finishing pending tasks before exit
    -[nf-core/rnaseq] Pipeline completed with errors-
    ERROR ~ Error executing process > 'NFCORE_RNASEQ:RNASEQ:FASTQ_QC_TRIM_FILTER_SETSTRANDEDNESS:FQ_LINT (WT_REP2)'

    Caused by:
      Process requirement exceeds available memory -- req: 12 GB; avail: 7.8 GB


    Command executed:

      fq lint \
          --disable-validator P001 \
          SRR6357072_1.fastq.gz SRR6357072_2.fastq.gz > WT_REP2.fq_lint.txt

    Command exit status:
      -

    Command output:
      (empty)

    Work dir:
      /workspaces/training/nfcore-use/work/xx/xxxxxxxxxxxxxxxxxxxxxx

    Container:
      quay.io/biocontainers/fq:0.12.0--h9ee0642_0

    Tip: you can replicate the issue by changing to the process work dir and entering the command `bash .command.run`

     -- Check '.nextflow.log' file for details
    ERROR ~ Pipeline failed. Please refer to troubleshooting docs for common issues: https://nf-co.re/docs/running/troubleshooting
    ```

The key line in that error is:

```console
Process requirement exceeds available memory -- req: 12 GB; avail: 7.8 GB
```

The default Codespaces machine has 8 GB of RAM.
The pipeline is requesting 12 GB for the `FQ_LINT` process — more than the machine can provide.

That 12 GB comes from the `process_low` resource label defined in `conf/base.config`:

```groovy title="conf/base.config"
withLabel:process_low {
    cpus   = { 2     * task.attempt }
    memory = { 12.GB * task.attempt }
    time   = { 4.h   * task.attempt }
}
```

There are two ways to deal with this.
You can make the pipeline's requests fit the machine, which is fine for testing on a small dataset, or you can give the pipeline a bigger machine, which is what you need for real analyses.
In this part, you do both: first you make the requests fit the default machine, then in [section 2](#2-run-on-a-bigger-machine) you move to a bigger machine.

### 1.3. Re-run with a custom configuration

We provide you with a custom config file that overrides the resource defaults to fit within the default Codespaces machine (2 CPUs, about 8 GB of memory).

??? full-code "laptop.config"

    ```groovy title="laptop.config" hl_lines="21-23"
    --8<-- "nfcore-use/laptop.config"
    ```

[Part 2](./02_configure_execution.md) introduced `withName:` to target a single process by name.
Here we mostly use `withLabel:` to target all processes that share a label at once.

The highlighted `withName` block is needed because one process, `TRIMGALORE`, does not get its CPUs from a label.
The pipeline sets them by name in its `conf/modules/trimgalore.config` file:

```groovy title="conf/modules/trimgalore.config (excerpt)"
withName: '.*:FASTQ_FASTQC_UMITOOLS_TRIMGALORE:TRIMGALORE' {
    cpus       = { 8 * task.attempt }
```

A `withName` setting takes precedence over a `withLabel` setting, so changing the labels does not affect this process.
Among `withName` selectors, one that matches only the simple name of the process, such as `withName: 'TRIMGALORE'`, has lower priority than the pipeline's selector, which matches the full name.
`withName: '.*:TRIMGALORE'` also matches the full name, and because your config file is applied last, it takes effect.

This file is already present in your working directory.
Pass it with `-c` to apply the overrides:

```bash
nextflow run nf-core/rnaseq -r 3.27.0 -profile test,docker -c laptop.config --outdir rnaseq-results
```

??? success "Command output (pipeline launching)"

    ```console
     N E X T F L O W   ~  version 26.04.4

    Launching `https://github.com/nf-core/rnaseq` [romantic_faraday] revision: a1fcdddd3b [3.27.0]

    ------------------------------------------------------
                                            ,--./,-.
            ___     __   __   __   ___     /,-._.--~'
      |\ | |__  __ /  ` /  \ |__) |__         }  {
      | \| |       \__, \__/ |  \ |___     \`-._,-`-,
                                            `._,._,'
      nf-core/rnaseq 3.27.0
    ------------------------------------------------------
    ...
    executor >  local (7)
    [xx/xxxxxx] NFCORE_RNASEQ:RNASEQ:FASTQ_QC_TRIM_FILTER_SETSTRANDEDNESS:FQ_LINT (RAP1_IAA_30M_REP1)                                 | 3 of 5
    [xx/xxxxxx] NFCORE_RNASEQ:RNASEQ:FASTQ_QC_TRIM_FILTER_SETSTRANDEDNESS:FASTQ_FASTQC_UMITOOLS_TRIMGALORE:FASTQC (RAP1_IAA_30M_REP1) | 2 of 5
    ...
    ```

The pipeline is now running, and you can watch tasks completing one by one.
On the default Codespaces machine, this minimal test dataset takes about 20 minutes to process, executing over 200 tasks in total.

Watch it for a minute or two, then move on to the next step.

??? warning "Disk space on the default machine"

    If you let this run go to completion on the default 2-core machine, it may fail near the end with a `no space left on device` error.
    The pipeline uses a separate container image for each tool, and together they take up most of the 32 GB disk of the default machine.
    The 8-core machine used in [section 2](#2-run-on-a-bigger-machine) has a 64 GB disk, so it is not affected.

    To free up space, remove the container images that are not in use:

    ```bash
    docker image prune -a
    ```

    Answer `y` when asked to confirm.
    Images that are needed again are downloaded on the next run.

### 1.4. Cancel the run

You don't need to wait for this run to finish.
Stopping a run is something you will do regularly, for example when you notice a mistake in the parameters, so it is worth seeing what happens.

Press ++ctrl+c++ once in the terminal where the pipeline is running.

??? failure "Command output"

    ```console
    ^C-[nf-core/rnaseq] Pipeline completed with errors-
    ERROR ~ Pipeline failed. Please refer to troubleshooting docs for common issues: https://nf-co.re/docs/running/troubleshooting

     -- Check '.nextflow.log' file for details
    ```

Nextflow stops submitting new tasks, terminates the tasks that are running, and exits within a few seconds.
The `Pipeline failed` message is printed by nf-core pipelines whenever a run does not complete, including when you stop it yourself, so it is expected here.

Press ++ctrl+c++ only once, and give Nextflow a few seconds to shut down.
Pressing it again forces Nextflow to exit immediately, which on a cluster or in the cloud can leave tasks running without Nextflow to manage them.

The tasks that completed before you cancelled the run are cached, so you can pick up where you left off with `-resume`.

??? tip "Stopping a run launched in the background"

    For long runs, you can launch Nextflow in the background with the `-bg` option, so that it keeps running after you close the terminal:

    ```bash
    nextflow run nf-core/rnaseq -r 3.27.0 -profile test,docker -c laptop.config --outdir rnaseq-results -bg
    ```

    Nextflow records the ID of its process in a file named `.nextflow.pid` in the launch directory.
    To stop the run, send that process a termination signal:

    ```bash
    kill $(cat .nextflow.pid)
    ```

    Nextflow shuts down the same way as with ++ctrl+c++, terminating the running tasks before it exits.
    The output goes to `.nextflow.log`, where you can see `Killing running tasks` followed by the same `Pipeline failed` message.

??? tip "Managing runs with Seqera Platform"

    Keeping track of runs from the terminal gets harder as you run more pipelines, especially in the background or on remote infrastructure.
    [Seqera Platform](https://seqera.io/platform/) provides a web interface to launch Nextflow pipelines, follow their progress task by task, and cancel or resume runs, on your own HPC or cloud infrastructure.
    See the [Scale with Seqera](../seqera_scale/index.md) course to get started.

### Takeaway

You have pulled `nf-core/rnaseq`, seen how nf-core resource labels and process-specific settings work, and overridden them with a custom config file to fit the pipeline on a small machine.
You also know how to stop a run cleanly.

### What's next?

Give the pipeline a bigger machine, so that it can run with its default resource requests.

---

## 2. Run on a bigger machine

Running with constrained resources is fine for testing, but real analyses need more resources, not fewer.
In this section, you move to a bigger machine and run the pipeline with its default resource requests.
The bigger machine uses more of your GitHub Codespaces quota, as explained below.

### 2.1. Change the machine type

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
The rest of this part, and Part 4, use the 8-core machine.

### 2.2. Check the new resources

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

### 2.3. Run the pipeline with its default resource requests

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

## 3. Read the execution report

### 3.1. Find the reports

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

### 3.2. Compare usage with allocation

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

### 3.3. Understand the cost of over-allocation

Over-allocating memory does not make a task run any faster.
It does limit how many tasks can run at the same time.

When Nextflow runs tasks locally, it only starts a task if enough CPUs and memory remain unreserved on the machine.
In our test run, the tasks running at the same time reserved up to 30 GB of the 32 GB available, so no more than 4 tasks ever ran at once, and 2.2 on average, even though the machine has 8 CPUs.
Memory, not CPU, was the bottleneck.

Open the execution timeline to see this.
Each bar is a task; many tasks wait for memory to be released before they can start, even when CPUs are idle.

<!-- TODO: screenshot of the timeline from the default 8-core run -->

The same problem appears on HPC clusters and in the cloud in a different form: oversized requests wait longer in the scheduler queue and cost more.

Generic labels are a common cause of over-allocation.
For example, `nf-core/fetchngs` requests 12 GB of memory for each download by default, because its download process uses the `process_low` label from the nf-core pipeline template, while each download actually uses about 15 MB.
You will run it in Part 4.

### Takeaway

The execution report shows how much of its allocation each process actually used.
Over-allocation limits how many tasks can run in parallel.

### What's next?

Write a configuration that matches resource requests to actual usage.

---

## 4. Right-size the resource allocations

### 4.1. Set a ceiling with `resourceLimits`

`resourceLimits` caps every resource request at the values you give it, whatever the pipeline asks for.
A request above the cap is reduced to the cap instead of failing.

The `test` profile of `nf-core/rnaseq` sets `resourceLimits` to 4 CPUs and 15 GB of memory.
That is why the run in section 1.2 failed on the default machine, which has 8 GB, of which about 7.8 GB is available: a 12 GB request is within the profile's limits, but not within the machine's.

The 8-core machine provides 8 CPUs and 32 GB of memory.
Create a file named `tuned.config` that sets the limits to what it can provide:

```groovy title="tuned.config" linenums="1"
process {
    // Never request more than the machine has
    resourceLimits = [cpus: 8, memory: 30.GB]
}
```

On its own, this does not help parallelism: every process still requests far more memory than it uses.

### 4.2. Right-size each label

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
The simple name is enough here because, unlike `TRIMGALORE` in section 1.3, the pipeline does not set resources for `PICARD_MARKDUPLICATES` by name.

### 4.3. Run the pipeline with the tuned configuration

Run the test profile again with `tuned.config` and a new output directory:

```bash
nextflow run nf-core/rnaseq -r 3.27.0 -profile test,docker -c tuned.config --outdir rnaseq-tuned
```

Note the run time, then open the new execution timeline in `rnaseq-tuned/pipeline_info/`.
The run in section 2.3 also had to download the container images, so its run time is not directly comparable.
Tasks now overlap more, because several of them fit in memory at once.

In our tests on the 8-core machine, with the container images already downloaded, the run took about 5.5 minutes with `tuned.config` compared to about 6.5 minutes with the default requests, and no task ran out of memory.
At most 6 tasks ran at once instead of 4, and the memory reserved at any one time peaked at 16 GB instead of 30 GB.

The time saved is modest here because the test dataset is tiny: most tasks finish within seconds, so starting them sooner gains little.
On real data, where individual tasks run for minutes or hours, running more of them in parallel makes a much larger difference.

??? tip "Applying this to a real analysis"

    The test dataset is tiny, so its resource usage does not represent a real analysis.
    For real data, run the pipeline on one or two representative samples first, read the execution report, then tune the configuration before launching the full dataset.
    Memory needs often scale with input size and reference genome size, so keep generous headroom.

### Takeaway

`resourceLimits` keeps requests within what the machine can provide, and label-level overrides based on the execution report let more tasks run in parallel.
Keeping the `task.attempt` multiplier preserves the pipeline's automatic retry with more resources.

### What's next?

Run the pipeline on your own inputs.
Head on to [Part 4](./04_custom_inputs.md), and stay on the 8-core machine: Part 4 uses it too.

---

## Summary

In this part you learned to:

- Pull and run a production-scale pipeline (nf-core/rnaseq), and override its default resource allocations by label and by process name
- Pin a pipeline to a specific release with `-r`
- Cancel a run cleanly
- Run a pipeline with more compute resources by changing the machine type
- Find and read the execution report and timeline generated by nf-core pipelines
- Recognize over-allocation and its effect on parallel execution
- Cap resource requests with `resourceLimits` and right-size them by label and by process name, while preserving automatic retries
