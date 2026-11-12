# Part 2: Run nf-core

In this second part of the Nextflow Triathlon, we show you how to find and try out an nf-core pipeline, configure and customize its execution, and understand how input validation protects against common errors.

We are going to use a pipeline called `nf-core/demo` that is maintained by the nf-core project as part of its inventory of pipelines for demonstration and training purposes.

!!! tip

    Switch to the `triathlon/nf-core/` directory for this part:

    ```bash
    cd /workspaces/training/triathlon/nf-core
    ```

    If you are coming directly from Part 1, you can use the relative path instead:

    ```bash
    cd ../nf-core
    ```

    Then enable the v1 syntax parser (required for nf-core pipelines on Nextflow 25.10 and later):

    ```bash
    export NXF_SYNTAX_PARSER=v1
    ```

---

## 1. Find and retrieve the nf-core/demo pipeline

Let's start by locating the `nf-core/demo` pipeline on the project website at [nf-co.re](https://nf-co.re), which centralizes all information such as: general documentation and help articles, documentation for each of the pipelines, blog posts, event announcements and so forth.

### 1.1. Find the pipeline on the website

In your web browser, go to [https://nf-co.re/pipelines/](https://nf-co.re/pipelines/) and type `demo` in the search bar.

Click on the pipeline name, `demo`, to access the pipeline documentation page.

Each released pipeline has a dedicated page that includes the following documentation sections:

- **Introduction:** An introduction and overview of the pipeline
- **Usage:** Descriptions of how to execute the pipeline
- **Parameters:** Grouped pipeline parameters with descriptions
- **Output:** Descriptions and examples of the expected output files
- **Results:** Example output files generated from the full test dataset
- **Releases & Statistics:** Pipeline version history and statistics

Whenever you are considering adopting a new pipeline, read the documentation carefully to understand what it does and how it should be configured before attempting to run it.

Have a look now and see if you can find out:

- Which tools the pipeline will run (Check the tab: `Introduction`)
- Which inputs and parameters the pipeline accepts or requires (Check the tab: `Parameters`)
- What outputs the pipeline produces (Check the tab: `Output`)

#### 1.1.1. Pipeline overview

The `Introduction` tab provides an overview of the pipeline, including a visual representation (called a subway map) and a list of tools that are run as part of the pipeline.

1. Read QC (`FASTQC`)
2. Adapter and quality trimming (`SEQTK_TRIM`)
3. Present QC for raw reads (`MULTIQC`)

#### 1.1.2. Example command line

The documentation also provides an example command line:

```bash
nextflow run nf-core/demo \
  -profile <docker/singularity/.../institute> \
  --input samplesheet.csv \
  --outdir <OUTDIR>
```

Notice that the command references the pipeline repository (`nf-core/demo`) rather than a local file.
When invoked this way, Nextflow assumes a standard code organization.
Let's retrieve the code to examine it.

### 1.2. Retrieve the pipeline code

#### 1.2.1. Use `nextflow pull`

Run the following command to download a local copy of the pipeline code.

```bash
nextflow pull nf-core/demo
```

??? success "Command output"

    ```console
    Checking nf-core/demo ...
    downloaded from https://github.com/nf-core/demo.git - revision: 45904cb9d1 [master]
    ```

Nextflow downloads the full repository to your local drive.
This works with any Nextflow pipeline that is appropriately set up in GitHub, not just nf-core pipelines.
However, nf-core is the largest open-source collection of Nextflow pipelines.

#### 1.2.2. Find your pipelines in `$NXF_HOME/assets/`

Nextflow saves pulled pipelines to `$NXF_HOME/assets`, keeping the source code out of your working directory on the principle that these pipelines should be used more like libraries.

```bash
tree -L 2 $NXF_HOME/assets/
```

```console title="Directory contents"
/workspaces/.nextflow/assets/
└── nf-core
    └── demo
```

Each pipeline is stored in its own subdirectory, named after the GitHub organisation and repository.

#### 1.2.3. Overview of the code organization

Run the following command to see the top-level layout of the pipeline code.

```bash
tree -L 1 $NXF_HOME/assets/nf-core/demo
```

??? abstract "Directory contents"

    ```console
    /workspaces/.nextflow/assets/nf-core/demo
    ├── assets
    ├── CHANGELOG.md
    ├── CITATIONS.md
    ├── CODE_OF_CONDUCT.md
    ├── conf
    ├── docs
    ├── LICENSE
    ├── main.nf
    ├── modules
    ├── modules.json
    ├── nextflow.config
    ├── nextflow_schema.json
    ├── nf-test.config
    ├── README.md
    ├── ro-crate-metadata.json
    ├── subworkflows
    ├── tests
    ├── tower.yml
    └── workflows
    ```

At the top level you'll find the main pipeline file (`main.nf`), configuration (`nextflow.config`, `conf/`), and the modular code components (`modules/`, `subworkflows/`, `workflows/`).

!!! tip

    You can browse any nf-core pipeline's source code on GitHub, e.g. [github.com/nf-core/demo](https://github.com/nf-core/demo).
    Every nf-core pipeline follows the same directory layout, so once you know the structure, you can navigate any pipeline the same way.

### Takeaway

You now know how to find an nf-core pipeline and retrieve a local copy of the source code.

### What's next?

Learn how to try out an nf-core pipeline with minimal effort.

---

## 2. Try out the pipeline with its test profile

Every nf-core pipeline comes with a test profile: a minimal configuration that runs the pipeline on a small hosted test dataset.
It's a quick way to verify a pipeline works in your environment.

### 2.1. Examine the test profile

The test profile for `nf-core/demo` lives in `conf/test.config`.
Check it before running:

```groovy title="conf/test.config" linenums="1" hl_lines="8 26"
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Nextflow config file for running minimal tests
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Defines input files and everything required to run a fast and simple pipeline test.

    Use as follows:
        nextflow run nf-core/demo -profile test,<docker/singularity> --outdir <OUTDIR>

----------------------------------------------------------------------------------------
*/

process {
    resourceLimits = [
        cpus: 2,
        memory: '4.GB',
        time: '1.h'
    ]
}

params {
    config_profile_name        = 'Test profile'
    config_profile_description = 'Minimal test dataset to check pipeline function'

    // Input data
    input  = 'https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/samplesheet/samplesheet_test_illumina_amplicon.csv'

}
```

The comment block shows you how to run with this profile.
The `input` parameter is already set to a hosted test samplesheet — no data needed on your end.

### 2.2. Run the pipeline

We'll use Docker for the container system and `demo-results` as the output directory:

```bash
nextflow run nf-core/demo -profile docker,test --outdir demo-results
```

??? success "Command output"

    ```console
     N E X T F L O W   ~  version 25.10.4

    Launching `https://github.com/nf-core/demo` [magical_pauling] DSL2 - revision: 45904cb9d1 [master]

    ------------------------------------------------------
                                            ,--./,-.
            ___     __   __   __   ___     /,-._.--~'
      |\ | |__  __ /  ` /  \ |__) |__         }  {
      | \| |       \__, \__/ |  \ |___     \`-._,-`-,
                                            `._,._,'
      nf-core/demo 1.1.0
    ------------------------------------------------------
    Input/output options
      input                     : https://...samplesheet_test_illumina_amplicon.csv
      outdir                    : demo-results
    ...

    executor >  local (7)
    [ff/a6976b] NFCORE_DEMO:DEMO:FASTQC (SAMPLE3_SE)     | 3 of 3 ✔
    [39/731ab7] NFCORE_DEMO:DEMO:SEQTK_TRIM (SAMPLE3_SE) | 3 of 3 ✔
    [7c/78d96e] NFCORE_DEMO:DEMO:MULTIQC                 | 1 of 1 ✔
    -[nf-core/demo] Pipeline completed successfully-
    ```

Notice the fuller process names such as `NFCORE_DEMO:DEMO:FASTQC`, which include the parent workflow names and reflect the pipeline's modular structure.

For reproducible production runs, pin a specific version with `-r`:

```bash
nextflow run nf-core/demo -r 1.1.0 -profile docker,test --outdir demo-results
```

### 2.3. Examine the pipeline's outputs

Run the following command to see the output structure.

```bash
tree -L 2 demo-results
```

??? abstract "Directory contents"

    ```console
    demo-results
    ├── fastqc
    │   ├── SAMPLE1_PE
    │   ├── SAMPLE2_PE
    │   └── SAMPLE3_SE
    ├── fq
    │   ├── SAMPLE1_PE
    │   ├── SAMPLE2_PE
    │   └── SAMPLE3_SE
    ├── multiqc
    │   ├── multiqc_data
    │   ├── multiqc_plots
    │   └── multiqc_report.html
    └── pipeline_info
        ├── execution_report_*.html
        ├── execution_timeline_*.html
        ├── execution_trace_*.txt
        └── pipeline_dag_*.html
    ```

Results are organized by module.
The `pipeline_info/` directory contains automatically generated execution reports, including a timeline showing what ran, in what order, and how long it took.

### Takeaway

You know how to run an nf-core pipeline using its built-in test profile and where to find its outputs.

### What's next?

Learn how to configure the pipeline to customize its execution.

---

## 3. Configure pipeline execution

nf-core pipelines distinguish two kinds of configuration:

- **Pipeline parameters** (set via `--param_name` or `-params-file`): input files, tool behavior flags, analysis parameters.
- **Configuration** (set via `-c`): execution logistics such as resource allocation and tool arguments.

<figure class="excalidraw">
    --8<-- "docs/en/docs/hello_nf-core/img/params_vs_config.excalidraw.svg"
</figure>

### 3.1. Pipeline parameters

#### 3.1.1. Get the list of parameters with `--help`

All nf-core pipelines support `--help` because every pipeline defines its parameters in a JSON schema file (`nextflow_schema.json`):

```bash
nextflow run nf-core/demo --help
```

??? success "Command output"

    ```console
    Input/output options
      --input                       [string]           Path to a metadata file containing information about the samples in the experiment.
      --outdir                      [string]           The output directory where the results will be saved.
      --email                       [string]           Email address for completion summary.
      --multiqc_title               [string]           MultiQC report title.

    Process skipping options
      --skip_trim                   [boolean]          Skip trimming fastq files with seqtk
    ...
    ```

Use `--help` whenever you want a quick reference to a pipeline's available parameters and their types.

#### 3.1.2. Set parameter values

Pass parameters on the command line with `--param_name`, or collect them in a YAML file and pass it with `-params-file`.

For example, to skip the trimming step:

```bash
nextflow run nf-core/demo -profile docker,test --outdir demo-results-notrim --skip_trim
```

??? success "Command output"

    ```console
    executor >  local (4)
    [3f/a82c91] NFCORE_DEMO:DEMO:FASTQC (SAMPLE3_SE) | 3 of 3 ✔
    [7d/c5e014] NFCORE_DEMO:DEMO:MULTIQC             | 1 of 1 ✔
    -[nf-core/demo] Pipeline completed successfully-
    ```

The `SEQTK_TRIM` step no longer appears in the output.

#### 3.1.3. Validation

Because parameters are defined in a schema, nf-core pipelines are able to validate parameters and input files at launch time.
If you pass an unknown parameter, a value of the wrong type, or a malformed samplesheet, the pipeline reports all errors upfront and stops before any work is done, saving you from wasted compute time.

For more details, see the more detailed sections on [Parameter validation](../hello_nf-core/01_run_demo.md#313-parameter-validation) and [Input validation](../hello_nf-core/01_run_demo.md#314-input-validation) in the [Hello nf-core](../hello_nf-core/index.md) course.

### 3.2. Configuration

Configuration in the strict sense controls how the pipeline runs: resource allocation, tool-specific arguments, executor settings.

List the available config files:

```bash
tree $NXF_HOME/assets/nf-core/demo/conf/
```

```console
.nextflow/assets/nf-core/demo/conf/
├── base.config
├── igenomes.config
├── igenomes_ignored.config
├── modules.config
├── test.config
└── test_full.config

1 directory, 6 files
```

The key files are:

- **`conf/base.config`**: Default resource labels (`process_low`, `process_medium`, `process_high`) used to assign CPUs, memory, and time.
- **`conf/modules.config`**: Per-process tool arguments (`ext.args`) and output publishing settings.

Never modify these files directly.
Instead, create your own config file and pass it with `-c` to override only what you need.

#### 3.2.1. Change resource allocation for a process

Create an empty config file:

```bash
touch custom.config
```

Open it and add the following:

```groovy title="custom.config" linenums="1"
process {
    withName: 'FASTQC' {
        cpus = 2
        memory = 4.GB
    }
}
```

Run the pipeline with the custom config:

```bash
nextflow run nf-core/demo -profile docker,test --outdir demo-results-custom -c custom.config -resume
```

??? success "Command output"

    ```console
    N E X T F L O W   ~  version 25.10.4

    Launching `https://github.com/nf-core/demo` [jolly_curie] DSL2 - revision: 45904cb9d1 [master]

    ...

    executor >  local (7)
    [ff/a6976b] NFCORE_DEMO:DEMO:FASTQC (SAMPLE3_SE)     | 3 of 3 ✔
    [39/731ab7] NFCORE_DEMO:DEMO:SEQTK_TRIM (SAMPLE3_SE) | 3 of 3 ✔
    [7c/78d96e] NFCORE_DEMO:DEMO:MULTIQC                 | 1 of 1 ✔
    -[nf-core/demo] Pipeline completed successfully-
    ```

To verify the resource override took effect, open the execution report that nf-core generates automatically for every run.
It is saved to `demo-results-custom/pipeline_info/execution_report_<timestamp>.html`.
Open it in a browser and check the `FASTQC` process — it will show the CPU and memory limits you set.

#### 3.2.2. Set tool arguments with `ext.args`

Many nf-core modules accept additional tool arguments through an `ext.args` configuration convention.
For example, to trim 5 bases from the start of each read during the `SEQTK_TRIM` step, update `custom.config` to add a second `withName` block:

```groovy title="custom.config" linenums="1" hl_lines="6 7 8"
process {
    withName: 'FASTQC' {
        cpus = 2
        memory = 4.GB
    }
    withName: 'SEQTK_TRIM' {
        ext.args = '-b 5'
    }
}
```

Run the pipeline again with the updated config:

```bash
nextflow run nf-core/demo -profile docker,test --outdir demo-results-extargs -c custom.config -resume
```

??? success "Command output"

    ```console
    N E X T F L O W   ~  version 25.10.4

    Launching `https://github.com/nf-core/demo` [sharp_volta] DSL2 - revision: 45904cb9d1 [master]

    ...

    executor >  local (7)
    [ff/a6976b] NFCORE_DEMO:DEMO:FASTQC (SAMPLE3_SE)     | 3 of 3 ✔
    [39/731ab7] NFCORE_DEMO:DEMO:SEQTK_TRIM (SAMPLE3_SE) | 3 of 3 ✔
    [7c/78d96e] NFCORE_DEMO:DEMO:MULTIQC                 | 1 of 1 ✔
    -[nf-core/demo] Pipeline completed successfully-
    ```

To verify the argument was applied, find the `SEQTK_TRIM` work directory from the run output and inspect the `.command.sh` file.

??? abstract ".command.sh"

    ```bash
    printf "%s\n" SAMPLE1_PE_R1.fastq.gz SAMPLE1_PE_R2.fastq.gz | while read f;
    do
        seqtk \
            trimfq \
            -b 5 \
            $f \
            | gzip --no-name > SAMPLE1_PE_$(basename $f)
    done
    ```

    The `-b 5` argument from `ext.args` is injected directly into the `seqtk trimfq` call.

!!! note

    If a module already has a default `ext.args` set in `conf/modules.config`, your value will **replace** it entirely.
    To keep the default and add to it, include the original value in your override.

### Takeaway

You know how to get help from an nf-core pipeline, set and validate parameters, and customize configuration through config files.

### What's next?

Try a real analysis pipeline!

---

## 4. Pull and run nf-core/rnaseq

So far we have used `nf-core/demo`, which is a minimal pipeline designed for training.
Now we pull a real production pipeline and run it with its test profile.

The `nf-core/rnaseq` pipeline performs the core steps of bulk RNA sequencing analysis: quality control, adapter trimming, read alignment, and gene-level quantification.
It is probably the most widely used nf-core pipeline to date.

### 4.1. Pull the pipeline

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

### 4.2. Run the test profile

Run it with the test profile and Docker:

```bash
nextflow run nf-core/rnaseq -profile docker,test --outdir rnaseq-results
```

??? failure "Command output"

    ```console
     N E X T F L O W   ~  version 25.10.4

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
      /workspaces/training/triathlon/nf-core/work/xx/xxxxxxxxxxxxxxxxxxxxxx

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

### 4.3. Re-run with a custom configuration

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

Section 3.2 introduced `withName:` to target a single process by name.
Here we use `withLabel:` to target all processes that share a label at once.

This file is already present in your working directory.
Pass it with `-c` to apply the overrides:

```bash
nextflow run nf-core/rnaseq -profile docker,test -c laptop.config --outdir rnaseq-results
```

??? success "Command output (pipeline launching)"

    ```console
     N E X T F L O W   ~  version 25.10.4

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

Move on to Part 3, where you will run `nf-core/rnaseq` at scale with Seqera.
