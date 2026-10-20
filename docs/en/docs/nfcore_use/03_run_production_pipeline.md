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
In this part, you use the first approach, which is fine for testing.
Later in the course, in Part 6, you do a full run on a bigger machine.

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
    The 8-core machine used in Part 6 has a 64 GB disk, so it is not affected.

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

Find and evaluate a pipeline for a scientific question of your own.
Head on to [Part 4](./04_choose_pipeline.md).

---

## Summary

In this part you learned to:

- Pull and run a production-scale pipeline (nf-core/rnaseq), and override its default resource allocations by label and by process name
- Pin a pipeline to a specific release with `-r`
- Cancel a run cleanly
