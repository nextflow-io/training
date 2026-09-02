# Part 2: Configure pipeline execution

In [Part 1](./01_run_demo.md), you found and ran the nf-core/demo pipeline using its test profile.
Now we look at how to configure pipeline execution: setting parameters, understanding validation, and customizing resource allocation and tool arguments.

---

## 1. Configure pipeline execution

As explained in [Hello Config](../hello_nextflow/06_hello_config.md), we want to be able to change what data our pipeline will run on and how it will run without changing the pipeline code itself.
To that end, Nextflow supports multiple ways of controlling pipeline configuration, which can be a bit overwhelming.

The nf-core project specifies conventions for organizing configuration elements, distinguishing two kinds of configuration at the top level: **pipeline parameters** and **configuration** in the strict sense.

- **Pipeline parameters** (set through the `params` system) typically include things like input files, tool behavior flags and analysis parameters.
- **Configuration** in the strict sense refers to the logistics of how the pipeline gets run, i.e. the executor, compute resource allocations and so on.

<figure class="excalidraw">
    --8<-- "docs/en/docs/nfcore_run/img/params_vs_config.excalidraw.svg"
</figure>

Let's start by tackling pipeline parameters, then we'll look at configuration in the strict sense.

### 1.1. Pipeline parameters

For all nf-core pipelines, you can obtain a full list of pipeline parameters directly from the command line by using the `--help` flag, which is itself a pipeline parameter.

#### 1.1.1. Get the list of parameters with `--help`

Run the help command for the demo pipeline:

```bash
nextflow run nf-core/demo --help
```

??? success "Command output"

    ```console
     N E X T F L O W   ~  version 26.04.4

    Launching `https://github.com/nf-core/demo` [adoring_meucci] revision: 32893afef8 [master]


    ------------------------------------------------------
                                            ,--./,-.
            ___     __   __   __   ___     /,-._.--~'
      |\ | |__  __ /  ` /  \ |__) |__         }  {
      | \| |       \__, \__/ |  \ |___     \`-._,-`-,
                                            `._,._,'
      nf-core/demo 1.2.0
    ------------------------------------------------------
    Typical pipeline command:

      nextflow run nf-core/demo -profile <docker/singularity/.../institute> --input samplesheet.csv --outdir <OUTDIR>

    Input/output options
      --input                       [string] Path to a metadata file containing information about the samples in the experiment.
      --outdir                      [string] The output directory where the results will be saved. You have to use absolute paths to storage on Cloud infrastructure.

      --email                       [string] Email address for completion summary.
      --multiqc_title               [string] MultiQC report title. Printed as page header, used for filename if not otherwise specified.

    Reference genome options
      --genome                      [string] Name of iGenomes reference.
      --fasta                       [string] Path to FASTA genome file.

    Process skipping options
      --skip_trim                   [boolean] Skip trimming fastq files with seqtk

    Generic options
      --multiqc_methods_description [string]          Custom MultiQC yaml file containing HTML including a methods description.
      --help                        [boolean, string] Display the help message.
      --help_full                   [boolean]         Display the full detailed help message.
      --show_hidden                 [boolean]         Display hidden parameters in the help message (only works when --help or --help_full are provided).
    !! Hiding 19 param(s), use the `--showHidden` parameter to show them !!
    ------------------------------------------------------

    * The pipeline
        https://doi.org/10.5281/zenodo.12192442

    * The nf-core framework
        https://doi.org/10.1038/s41587-020-0439-x

    * Software dependencies
        https://github.com/nf-core/demo/blob/master/CITATIONS.md
    ```

As you can see, the output groups parameters into categories (Input/output options, Reference genome options, etc.) with types and descriptions for each one.

This categorization is determined by a schema file, which is covered further below.
In plain Nextflow pipelines, `--help` only works if the developer implemented it manually.

!!! tip

    Use `--help --show_hidden` to see additional parameters that are hidden by default, such as `--publish_dir_mode` or `--monochrome_logs`.

#### 1.1.2. Set parameter values

As covered in [Hello Config](../hello_nextflow/06_hello_config.md), you can set parameter values on the command line with `--param_name` or collect a set of parameters in a YAML file and pass it with `-params-file`.
Both approaches work the same way with nf-core pipelines.

For example, to skip the trimming step, we want to set the boolean parameter `skip_trim` to `true`.
A params file called `my_params.yml` is provided in your working directory with that value already set:

```yaml title="my_params.yml"
skip_trim: true
```

Pass it with `-params-file`:

```bash
nextflow run nf-core/demo -profile test,docker --outdir demo-results-notrim -params-file my_params.yml
```

??? success "Command output"

    ```console
     N E X T F L O W   ~  version 26.04.4

    Launching `https://github.com/nf-core/demo` [focused_heisenberg] revision: 32893afef8 [master]


    ------------------------------------------------------
                                            ,--./,-.
            ___     __   __   __   ___     /,-._.--~'
      |\ | |__  __ /  ` /  \ |__) |__         }  {
      | \| |       \__, \__/ |  \ |___     \`-._,-`-,
                                            `._,._,'
      nf-core/demo 1.2.0
    ------------------------------------------------------

    Input/output options
      input                     : https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/samplesheet/samplesheet_test_illumina_amplicon.csv
      outdir                    : demo-results-notrim

    Process skipping options
      skip_trim                 : true

    Institutional config options
      config_profile_name       : Test profile
      config_profile_description: Minimal test dataset to check pipeline function

    Generic options
      trace_report_suffix       : 2026-07-03_22-08-47

    Core Nextflow options
      revision                  : master
      runName                   : focused_heisenberg
      containerEngine           : docker
      launchDir                 : /workspaces/training/hello-nf-core
      workDir                   : /workspaces/training/hello-nf-core/work
      projectDir                : /workspaces/.nextflow/assets/.repos/nf-core/demo/clones/32893afef8076a03a2767a020b3f0cab2e0b40b2
      userName                  : root
      profile                   : test,docker
      configFiles               : /workspaces/.nextflow/assets/.repos/nf-core/demo/clones/32893afef8076a03a2767a020b3f0cab2e0b40b2/nextflow.config

    !! Only displaying parameters that differ from the pipeline defaults !!
    ------------------------------------------------------

    * The pipeline
        https://doi.org/10.5281/zenodo.12192442

    * The nf-core framework
        https://doi.org/10.1038/s41587-020-0439-x

    * Software dependencies
        https://github.com/nf-core/demo/blob/master/CITATIONS.md

    executor >  local (5)
    [7a/f3599e] NFCORE_DEMO:DEMO:FASTQC (SAMPLE3_SE) [100%] 3 of 3 ✔
    [b0/2f0bdc] NFCORE_DEMO:DEMO:COWPY               [100%] 1 of 1 ✔
    [c3/3c2278] NFCORE_DEMO:DEMO:MULTIQC (demo)      [100%] 1 of 1 ✔
    -[nf-core/demo] Pipeline completed successfully-
    ```

The `SEQTK_TRIM` process no longer appears in the output.

!!! warning "Important limitations about parameter inputs"

    **Setting boolean parameters on the command line**

    Starting with Nextflow version 26.04, all values supplied on the command line are typed as strings.
    For a boolean parameter like `skip_trim`, passing it as a bare flag (`--skip_trim`) or as `--skip_trim true` is evaluated as the **string** `"true"`, which fails schema validation:

    ```console
    * --skip_trim (true): Value is [string] but should be [boolean]
    ```

    To set a boolean parameter to a genuine `true`/`false` value, use a `-params-file` as shown above, or set it in a config file.
    String, integer and file-path parameters are unaffected and can still be set directly on the command line.
    This course uses this pattern throughout for boolean parameters.

    **Using custom configuration files**

    Although it is technically possible to set pipeline parameters in a custom configuration file passed with `-c`, this may not override defaults already set in the pipeline's own `nextflow.config`, depending on Nextflow's configuration precedence rules.
    Using `--param_name` on the command line or `-params-file` is more reliable, as these always take precedence.

    As a rule of thumb: If it appears in the `--help` output, set it via the command line or a params file rather than a config file.

#### 1.1.3. Parameter validation

Fun fact: the `--help` command works for all nf-core pipelines because the nf-core project requires developers to define all pipeline parameters formally in a JSON schema file (`nextflow_schema.json`).
This schema records each parameter's type, description, default value, and grouping.

In addition to powering the `--help` output, the schema file also enables automated validation at launch time.
This means that Nextflow can check that every parameter you pass exists and has been given an appropriate value (of appropriate type, within the allowed range of values etc).

We cover this in more detail in [input validation section](../hello_nf-core/04_input_validation.md), but you can already see it in action by giving the demo pipeline some invalid parameter input.

##### 1.1.3.1. Unrecognized parameters

Try passing a parameter that does not exist:

```bash
nextflow run nf-core/demo -profile test,docker --outdir demo-results --foobar "invalid"
```

The console output includes a warning:

```console
WARN: The following invalid input values have been detected:

* --foobar: invalid
```

The pipeline still runs, but the warning alerts you right away that `--foobar` is not a recognized parameter.
This is meant to draw your attention to non-breaking typos, like `--outDir` being used instead of `--outdir`, which can help you avoid wasting time and compute.

##### 1.1.3.2. Invalid parameter values

Validation also checks parameter **values**.
The `--skip_trim` parameter is a boolean flag, so passing a string value causes the pipeline to fail immediately:

```bash
nextflow run nf-core/demo -profile test,docker --outdir demo-results --skip_trim yes
```

```console
ERROR ~ Validation of pipeline parameters failed!

 -- Check '.nextflow.log' file for details
The following invalid input values have been detected:

* --skip_trim (yes): Value is [string] but should be [boolean]
```

The pipeline stops before any processes run, saving you from a failed or incorrect execution.
As noted in section 1.1.2, boolean parameters should be set to a genuine `true`/`false` value in a params file rather than passed on the command line, since command-line values are typed as strings.

#### 1.1.4. Input validation

The same validation logic can also be used to check the validity of input files.
For example, if a pipeline expects a samplesheet as its main data input (which is the case of many if not most nf-core pipelines), the developer can provide an input schema (distinct from the parameters schema) describing how the input file should be structured.

Then, at runtime, Nextflow can check that the input file provided is valid.

We also cover this in more detail in [input validation section](../hello_nf-core/04_input_validation.md), but you can already see it in action by giving the demo pipeline an invalid input samplesheet.

The `nf-core/demo` pipeline expects a CSV file with columns `sample`, `fastq_1`, and `fastq_2`.
This is defined in a schema file (`assets/schema_input.json`) that specifies the expected structure, column types, and constraints.

??? abstract "Schema file for inputs"

    ```json title="assets/schema_input.json"
    {
        "$schema": "https://json-schema.org/draft/2020-12/schema",
        "$id": "https://raw.githubusercontent.com/nf-core/demo/master/assets/schema_input.json",
        "title": "nf-core/demo pipeline - params.input schema",
        "description": "Schema for the file provided with params.input",
        "type": "array",
        "items": {
            "type": "object",
            "properties": {
                "sample": {
                    "type": "string",
                    "pattern": "^\\S+$",
                    "errorMessage": "Sample name must be provided and cannot contain spaces",
                    "meta": ["id"]
                },
                "fastq_1": {
                    "type": "string",
                    "format": "file-path",
                    "exists": true,
                    "pattern": "^([\\S\\s]*\\/)?[^\\s\\/]+\\.f(ast)?q\\.gz$",
                    "errorMessage": "FastQ file for reads 1 must be provided, cannot contain spaces and must have extension '.fq.gz' or '.fastq.gz'"
                },
                "fastq_2": {
                    "type": "string",
                    "format": "file-path",
                    "exists": true,
                    "pattern": "^([\\S\\s]*\\/)?[^\\s\\/]+\\.f(ast)?q\\.gz$",
                    "errorMessage": "FastQ file for reads 2 cannot contain spaces and must have extension '.fq.gz' or '.fastq.gz'"
                }
            },
            "required": ["sample", "fastq_1"]
        }
    }
    ```

The schema specifies that `sample` and `fastq_1` are required, while `fastq_2` is optional (supporting both paired-end and single-end data).
File paths are validated for existence and extension pattern.

To demonstrate this, we provide a malformed samplesheet called `malformed_samplesheet.csv` in your working directory:

```csv title="malformed_samplesheet.csv"
sample,fastq_2
SAMPLE1,/not/a/real/file.fastq.gz
```

This samplesheet is missing the required `fastq_1` column and has a non-existent file path in `fastq_2`.

Run the demo pipeline using `malformed_samplesheet.csv` as the input:

```bash
nextflow run nf-core/demo -profile test,docker --outdir demo-results --input malformed_samplesheet.csv
```

```console
ERROR ~ Validation of pipeline parameters failed!

 -- Check '.nextflow.log' file for details
The following invalid input values have been detected:

* --input (malformed_samplesheet.csv): Validation of file failed:
    -> Entry 1: Error for field 'fastq_2' (/not/a/real/file.fastq.gz): the file or directory
       '/not/a/real/file.fastq.gz' does not exist (FastQ file for reads 2 cannot contain spaces
       and must have extension '.fq.gz' or '.fastq.gz')
    -> Entry 1: Missing required field(s): fastq_1
```

As you can see, the pipeline fails immediately and reports **all** validation errors at once.
nf-schema does not stop at the first error — it collects every problem and lists them together, so you can fix everything in one go rather than discovering issues one by one.

Each error identifies the exact entry and field that caused the problem, so you can fix your samplesheet then re-launch the pipeline with confidence that it's not going to fail at some later point when Nextflow actually goes to access the file path.

For developers, all of this is covered in more detail in [Part 4 of Build with nf-core](../hello_nf-core/04_input_validation.md).

### 1.2. Configuration

Configuration in the strict sense controls **how** the pipeline runs: resource allocation, tool-specific arguments, where jobs execute, and which software packaging system to use.

nf-core pipelines include default configuration in `nextflow.config` and the `conf/` directory.
Before overriding anything, it helps to know where the defaults live.

You already saw in [Part 1](./01_run_demo.md) that the pipeline source code lives under `$NXF_HOME/assets`.
Using the `pipelines` symlink you created in [Part 1](./01_run_demo.md), list the config files to see what's available:

```bash
ls pipelines/nf-core/demo/conf/
```

```console
base.config
containers_conda_lock_files_amd64.config
containers_conda_lock_files_arm64.config
containers_docker_amd64.config
containers_docker_arm64.config
containers_singularity_https_amd64.config
containers_singularity_https_arm64.config
containers_singularity_oras_amd64.config
containers_singularity_oras_arm64.config
igenomes.config
igenomes_ignored.config
modules.config
test.config
test_full.config
```

<figure class="excalidraw">
--8<-- "docs/en/docs/nfcore_run/img/nfcore_config_files.excalidraw.svg"
</figure>

The most important configuration files are:

- **`conf/base.config`**: Defines resource labels (`process_low`, `process_medium`, `process_high`) that assign CPUs, memory, and time to processes. When you see a process using more resources than expected, this is where those defaults come from.
- **`conf/modules.config`**: Sets per-process tool arguments (`ext.args`) and output publishing settings (`publishDir`). Open this file to see what arguments each tool receives by default.
- **`conf/test.config`**: The test profile you used in [Part 1](./01_run_demo.md), which caps resources via `resourceLimits` and sets a test samplesheet. Activated with `-profile test`.
  There is also a `conf/test_full.config` for running with a full-sized test dataset, useful for benchmarking.

The central `nextflow.config` loads all of the above and sets the appropriate default values for everything.

If you wish to modify any of the settings specified in these files, do not modify any of them files directly.
Instead, create your own config file and pass it with `-c`.
The values you specify will override the default values set in those other files.

Let's try this in practice.

#### 1.2.1. Customize process resources and tool arguments

nf-core modules support two common types of configuration override: **resource allocation** (CPUs, memory, time) and **tool arguments** via `ext.args`.

Many command-line tools have arguments that are not commonly enough used to be exposed as pipeline parameters.
The `ext.args` convention lets you pass these arguments to the underlying tool through a config file instead.

The `custom.config` file provided in your working directory demonstrates both overrides:

```groovy title="custom.config" linenums="1"
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

The first block overrides `FASTQC` resource allocation.
By default, `FASTQC` uses the `process_medium` label from `base.config`, which allocates 6 CPUs and 36 GB of memory; here we cap it at 2 CPUs and 4 GB.

The second block passes an extra argument to `SEQTK_TRIM` via `ext.args`.
The `-b 5` flag tells `seqtk trimfq` to trim 5 bases from the beginning of each read in addition to quality trimming.

Run the pipeline with this config:

```bash
nextflow run nf-core/demo -profile test,docker --outdir demo-results-custom -c custom.config
```

??? success "Command output"

    ```console
    executor >  local (8)
    [95/b32876] NFCORE_DEMO:DEMO:FASTQC (SAMPLE1_PE)     | 3 of 3 ✔
    [17/428668] NFCORE_DEMO:DEMO:SEQTK_TRIM (SAMPLE1_PE) | 3 of 3 ✔
    [cf/85991a] NFCORE_DEMO:DEMO:COWPY                   | 1 of 1 ✔
    [3c/94a7a0] NFCORE_DEMO:DEMO:MULTIQC (demo)          | 1 of 1 ✔
    -[nf-core/demo] Pipeline completed successfully-
    ```

The `-c` flag adds your config on top of the pipeline's built-in configuration.

To verify the `ext.args` override took effect, find the `SEQTK_TRIM` work directory hash from the run output (e.g. `work/17/428668...`) and check the `.command.sh` file inside it:

```bash
cat work/17/428668/.command.sh
```

??? success "Command output"

    ```console
    #!/usr/bin/env bash -e -u -o pipefail
    printf "%s\n" sample1_R1.fastq.gz sample1_R2.fastq.gz | while read f;
    do
        seqtk \
            trimfq \
            -b 5 \
            $f \
            | gzip --no-name > SAMPLE1_PE_$(basename $f)
    done
    ...
    ```

You should see `-b 5` in the `seqtk trimfq` command.

One important thing to know about `ext.args`: if a module already has a default value set, your value will **completely replace** it rather than append to it.
For example, `FASTQC` has `ext.args = '--quiet'` set by default in `conf/modules.config`:

```groovy title="conf/modules.config" linenums="21" hl_lines="2"
    withName: FASTQC {
        ext.args   = '--quiet'
        publishDir = [
            path: { "${params.outdir}/fastqc/${meta.id}" },
            mode: params.publish_dir_mode,
            pattern: "*.{html,json}",
        ]
    }
```

If you set `ext.args = '--kmers 8'` for `FASTQC`, the `--quiet` flag will no longer be applied.
To keep both, set `ext.args = '--quiet --kmers 8'`.

You should always check a module's default configuration before overriding `ext.args`.

### Takeaway

You know how to get help from an nf-core pipeline, set parameters and understand how they are validated, and customize configuration through config files.

### What's next?

Head on to [Part 3](./03_run_production_pipeline.md), where you'll apply what you've learned to a real production pipeline.

---

## Summary

In this part you learned to:

- Get help, set parameters, and understand parameter and input validation
- Customize resource allocation and tool arguments through configuration files
