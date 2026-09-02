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

Run the following command to download it.

```bash
nextflow pull nf-core/rnaseq
```

??? success "Command output"

    ```console
    Checking nf-core/rnaseq ...
    downloaded from https://github.com/nf-core/rnaseq.git - revision: e7ca46272c [master]
    ```

The pipeline is now cached locally and ready to run.

### 1.2. Run the test profile

Run it with the test profile and Docker:

```bash
nextflow run nf-core/rnaseq -profile test,docker --outdir rnaseq-results
```

??? failure "Command output"

    ```console
     N E X T F L O W   ~  version 26.04.4

    Launching `https://github.com/nf-core/rnaseq` [suspicious_dijkstra] DSL2 - revision: e7ca46272c [master]

    ------------------------------------------------------
                                            ,--./,-.
            ___     __   __   __   ___     /,-._.--~'
      |\ | |__  __ /  ` /  \ |__) |__         }  {
      | \| |       \__, \__/ |  \ |___     \`-._,-`-,
                                            `._,._,'
      nf-core/rnaseq 3.26.0
    ------------------------------------------------------
    ...
    [-        ] NFCORE_RNASEQ:RNASEQ:FASTQ_QC_TRIM_FILTER_SETSTRANDEDNESS:FQ_LINT  -
    Plus 47 more processes waiting for tasks…

    Execution cancelled -- Finishing pending tasks before exit
    -[nf-core/rnaseq] Pipeline completed with errors-
    ERROR ~ Error executing process > 'NFCORE_RNASEQ:RNASEQ:FASTQ_QC_TRIM_FILTER_SETSTRANDEDNESS:FQ_LINT (WT_REP2)'

    Caused by:
      Process requirement exceeds available memory -- req: 12 GB; avail: 7.7 GB

    Command executed:
      fq lint \
          --disable-validator P001 \
          SRR6357072_1.fastq.gz SRR6357072_2.fastq.gz > WT_REP2.fq_lint.txt

    Command exit status:
      -

    Command output:
      (empty)

    Work dir:
      /workspaces/training/nfcore-run/work/xx/xxxxxxxxxxxxxxxxxxxxxx

    Container:
      quay.io/biocontainers/fq:0.12.0--h9ee0642_0

    Tip: you can replicate the issue by changing to the process work dir and entering the command `bash .command.run`
     -- Check '.nextflow.log' file for details
    ERROR ~ Pipeline failed. Please refer to troubleshooting docs for common issues: https://nf-co.re/docs/running/troubleshooting
     -- Check '.nextflow.log' file for details
    ```

The key line in that error is:

```console
Process requirement exceeds available memory -- req: 12 GB; avail: 7.7 GB
```

The default Codespaces machine has 8 GB of RAM, which is also the typical default for Docker Desktop.
The pipeline is requesting 12 GB for the `FQ_LINT` process — more than the machine can provide.

That 12 GB comes from the `process_low` resource label defined in `conf/base.config`:

```groovy title="conf/base.config"
withLabel:process_low {
    cpus   = { 2     * task.attempt }
    memory = { 12.GB * task.attempt }
    time   = { 4.h   * task.attempt }
}
```

One option would be to use a larger machine type, but for testing purposes we want to be able to run on whatever hardware is available.
The better approach is to override the resource defaults in a custom config file.

### 1.3. Re-run with a custom configuration

We provide you with a custom config file that overrides the label-based resource defaults.

??? full-code "laptop.config"

    ```groovy title="laptop.config"
    process {
        withLabel: 'process_low' {
            cpus   = 2
            memory = 6.GB
        }
        withLabel: 'process_medium' {
            cpus   = 4
            memory = 6.GB
        }
        withLabel: 'process_high' {
            cpus   = 6
            memory = 6.GB
        }
        withLabel: 'process_high_memory' {
            memory = 6.GB
        }
    }
    ```

[Part 2](./02_configure_execution.md) introduced `withName:` to target a single process by name.
Here we use `withLabel:` to target all processes that share a label at once.

This file is already present in your working directory.
Pass it with `-c` to apply the overrides:

```bash
nextflow run nf-core/rnaseq -profile test,docker -c laptop.config --outdir rnaseq-results
```

??? success "Command output (pipeline launching)"

    ```console
     N E X T F L O W   ~  version 26.04.4

    Launching `https://github.com/nf-core/rnaseq` [romantic_faraday] DSL2 - revision: e7ca46272c [master]

    ------------------------------------------------------
                                            ,--./,-.
            ___     __   __   __   ___     /,-._.--~'
      |\ | |__  __ /  ` /  \ |__) |__         }  {
      | \| |       \__, \__/ |  \ |___     \`-._,-`-,
                                            `._,._,'
      nf-core/rnaseq 3.26.0
    ------------------------------------------------------
    ...
    executor >  local
    [xx/xxxxxx] NFCORE_RNASEQ:RNASEQ:FQ_LINT (RAP1_IAA_30M_REP1)   | 3 of 5, running
    [xx/xxxxxx] NFCORE_RNASEQ:RNASEQ:FASTQC (RAP1_IAA_30M_REP1)    | 2 of 5, running
    ...
    ```

The pipeline is now running, and you can watch tasks completing one by one.
On this minimal test dataset it will complete in 15–20 minutes, executing over 200 tasks in total.

Real RNA-seq experiments typically involve dozens of samples and run for hours or days.
Nextflow supports HPC schedulers (SLURM, PBS, LSF) and cloud platforms (AWS, Google Cloud, Azure), which can dramatically reduce wall-clock time by distributing work across many nodes.
Setting up those environments, however, adds significant complexity.

The Seqera platform (developed by the creators of Nextflow) provides a web-based interface for launching Nextflow pipelines on HPC or cloud infrastructure (either your own or one managed for you), with compute and data management capabilities that streamline the process of running pipelines at scale.

!!! tip

    Academic researchers can access Seqera Platform free of charge through the [Seqera academic program](https://seqera.io/academic-program/).

### Takeaway

You have pulled `nf-core/rnaseq`, seen how nf-core resource labels work, and learned to override them with a custom config file.
More importantly, you have seen why local execution is a starting point rather than a destination for real-scale analysis.

### What's next?

You've covered the fundamentals of running nf-core pipelines.
See [Next steps](next_steps.md) for where to go from here.

---

## Summary

In this part you learned to:

- Pull and run a production-scale pipeline (nf-core/rnaseq), and override its default resource labels
