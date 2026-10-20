# Part 3: Run a production pipeline

In [Part 2](./02_configure_execution.md), you learned how to set parameters and customize configuration for nf-core/demo.
Now we apply what you've learned to a real production pipeline, nf-core/rnaseq.

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
We start with the first approach; section 2 covers the second.

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

!!! warning "Disk space on the default machine"

    If you let this run go to completion on the default 2-core machine, it may fail near the end with a `no space left on device` error.
    The pipeline uses a separate container image for each tool, and together they take up most of the 32 GB disk of the default machine.
    The 8-core machine used in section 2 has a 64 GB disk, so it is not affected.

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

### Takeaway

You have pulled `nf-core/rnaseq`, seen how nf-core resource labels and process-specific settings work, and overridden them with a custom config file to fit the pipeline on a small machine.
You also know how to stop a run cleanly.

### What's next?

Give the pipeline more resources instead of constraining it.

---

## 2. Optional: Run with more compute resources

Constraining resource requests, as in section 1, gets a pipeline running on a small machine, which is useful for testing.
For real analyses, you need the opposite: more compute resources.

This section shows how to switch your codespace to a bigger machine and run the pipeline with its default resource requests.
It is optional, and it uses more of your GitHub Codespaces quota, as explained below.

### 2.1. Change the machine type

You can change the machine type of your codespace from the editor:

1. Open the Command Palette with ++ctrl+shift+p++ (++cmd+shift+p++ on macOS).
2. Type `Codespaces: Change Machine Type` and select it.
3. Select the **8-core** machine type.
4. Confirm that you want to restart the codespace when prompted.

Alternatively, go to [github.com/codespaces](https://github.com/codespaces), click the **...** menu next to your codespace, and select **Change machine type**.

The change takes effect when the codespace restarts, which can take a few minutes.
Your files in `/workspaces` are preserved, including the work directory of the run you cancelled, but anything running in the terminal is stopped.

The machine types available for this repository include the following (a 16-core machine with 64 GB is also available):

| Machine type | CPUs | Memory | Quota use relative to the default |
| ------------ | ---- | ------ | --------------------------------- |
| 2-core       | 2    | 8 GB   | 1x (default)                      |
| 4-core       | 4    | 16 GB  | 2x                                |
| 8-core       | 8    | 32 GB  | 4x                                |

Your quota is counted in core-hours, so a bigger machine uses it up faster: the free allowance of 120 core-hours per month gives you about 60 hours on the default machine, but only 15 hours on an 8-core machine.
The rest of this section uses the 8-core machine.

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

Run the pipeline again, this time without `laptop.config`, and add `-resume` to reuse the tasks that completed before you cancelled the previous run:

```bash
nextflow run nf-core/rnaseq -r 3.27.0 -profile test,docker --outdir rnaseq-results -resume
```

??? success "Command output (excerpt)"

    ```console
    [52/0bb6df] Cached process > NFCORE_RNASEQ:PREPARE_GENOME_REFERENCES:GUNZIP_GTF (genes_with_empty_tid.gtf.gz)
    [73/fbf95a] Cached process > NFCORE_RNASEQ:RNASEQ:FASTQ_QC_TRIM_FILTER_SETSTRANDEDNESS:FQ_LINT (RAP1_UNINDUCED_REP1)
    [b8/0bda15] Cached process > NFCORE_RNASEQ:PREPARE_GENOME_INDICES:UNTAR_SALMON_INDEX (salmon_2.7.0.tar.gz)
    ...
    -[nf-core/rnaseq] Pipeline completed successfully-
    ```

The tasks marked `Cached` completed before you cancelled the previous run, so Nextflow reuses their results instead of running them again.
How many tasks are cached depends on how long you let the previous run go.

With 8 CPUs and 32 GB, the default requests fit, so the run proceeds without any overrides.
In our tests, the resumed run took about 8 minutes on the 8-core machine, compared to about 20 minutes for a full run on the default machine with `laptop.config`.
Part of that time goes to downloading the container image for each tool the first time it is used; once the images are downloaded, a full run on the 8-core machine takes about 6.5 minutes.

Four times as many CPUs does not make the run four times faster, for two reasons.
The `test` profile limits each task to 4 CPUs and 15 GB of memory, and the test dataset is so small that many tasks finish in seconds regardless of the resources they get.
On a real dataset, the long-running tasks such as read alignment benefit much more from additional resources.

Real RNA-seq experiments typically involve dozens of samples and run for hours or days.
Nextflow supports HPC schedulers (SLURM, PBS, LSF) and cloud platforms (AWS, Google Cloud, Azure), which can dramatically reduce wall-clock time by distributing work across many machines.
Setting up those environments, however, adds significant complexity.

The Seqera platform (developed by the creators of Nextflow) provides a web-based interface for launching Nextflow pipelines on HPC or cloud infrastructure (either your own or one managed for you), with compute and data management capabilities that streamline the process of running pipelines at scale.

??? tip "Seqera academic program"

    Academic researchers can access Seqera Platform free of charge through the [Seqera academic program](https://seqera.io/academic-program/).

### 2.4. Switch back to the default machine

When you're done, switch your codespace back to the 2-core machine type, following the same steps as in 2.1, to avoid using up your quota.
Parts 4 and 5 run on the default machine; section 2.4 of Part 6 asks you to switch to the 8-core machine again.

### Takeaway

You can give a pipeline more resources by running it on a bigger machine, and see why local execution is a starting point rather than a destination for real-scale analysis.

### What's next?

Find and evaluate a pipeline for a scientific question of your own.
Head on to [Part 4](./04_choose_pipeline.md).

---

## Summary

In this part you learned to:

- Pull and run a production-scale pipeline (nf-core/rnaseq), and override its default resource allocations by label and by process name
- Pin a pipeline to a specific release with `-r`
- Cancel a run cleanly, and resume it later
- Run a pipeline with more compute resources by changing the machine type
