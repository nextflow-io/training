# Part 11: Hello nf-core pipeline

nf-core is a community effort to develop and maintain a curated set of analysis pipelines built using Nextflow.

![nf-core logo](../nf_customize/img/nf-core-logo.png)

The community provides tooling to create pipeline templates and use ready-made components that they have developed. 

## Create a basic pipeline from template

The nf-core pipeline template is a standardized framework designed to streamline the development of Nextflow-based bioinformatics pipelines.

Creating a pipeline using the nf-core template is greatly simplified by the nf-core tooling. It will help you create a pipeline using the set framework that can be modified to suit your own purposes.

Here, you will use the nf-core template to kickstart your pipeline development using the latest version of Nextflow and the nf-core tooling.

### Creating your pipeline

nf-core tooling has commands for pipeline users and developers.

View all of the tooling using the `nf-core --help` argument.

```bash
nf-core --help
```

Here we will focus on the tooling to assist pipeline developers, starting with the `nf-core pipelines create` command.

The `nf-core pipelines create` command creates a new pipeline using the nf-core base template with a pipeline name, description, and author. It is the first and most important step for creating a pipeline that will integrate with the wider Nextflow ecosystem.

```bash
nf-core pipelines create
```

Running this command will open a Text User Interface (TUI) for pipeline creation.

<!-- TODO: Change this clip to what we'll do -->
<div style="text-align: center;">
    <iframe width="560" height="315" src="https://www.youtube.com/embed/VwjXNXONHlY?si=d0HkFSISnKn76TeI" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" allowfullscreen="" data-ruffle-polyfilled=""></iframe>
</div>

Template features can be flexibly included or excluded at the time of creation:

Follow these steps create your first pipeline using the `nf-core pipelines create` TUI:

    1. Run the `nf-core pipelines create` command
    2. Select **Let's go!** on the welcome screen
    3. Select **Custom** on the Choose pipeline type screen
    4. Enter your pipeline details, replacing < YOUR NAME > with your own name, then select **Next**
        - **GitHub organisation:** myorg
        - **Workflow name:** myfirstpipeline
        - **A short description of your pipeline:** My first pipeline
        - **Name of the main author / authors:** < YOUR NAME >
    5. On the Template features screen, turn **off**:
        - `Use a GitHub repository`
        - `Add Github CI tests`
        - `Use reference genomes`
        - `Add Github badges`
        - `Include citations`
        - `Include a gitpod environment`
        - `Include GitHub Codespaces`
        - `Use fastqc`
        - `Add a changelog`
        - `Support Microsoft Teams notifications`
        - `Support Slack notifications`
    6. Select **Finish** on the Final details screen
    7. Wait for the pipeline to be created, then select **Continue**
    8. Select **Finish without creating a repo** on the Create GitHub repository screen
    9. Select **Close** on the HowTo create a GitHub repository page

If run successfully, you will see a new folder in your current directory named `myorg-myfirstpipeline`.

###  Testing your pipeline

Let's try it:

```bash
cd /workspace/gitpod/hello-nextflow/hello-nf-core-template/nf-core-training-firstpipeline
nextflow run . -profile docker,test --outdir results
```

The pipeline should run successfully!

Here's the console output from the pipeline:

```console title="Output"
Launching `./main.nf` [evil_mestorf] DSL2 - revision: 11a3012ba7

Input/output options
input                     : https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/samplesheet/samplesheet_test_illumina_amplicon.csv
outdir                    : results

Institutional config options
config_profile_name       : Test profile
config_profile_description: Minimal test dataset to check pipeline function

Core Nextflow options
runName                   : evil_mestorf
containerEngine           : docker
launchDir                 : /workspace/gitpod/hello-nextflow/hello-nf-core-template/nf-core-training-firstpipeline
workDir                   : /workspace/gitpod/hello-nextflow/hello-nf-core-template/nf-core-training-firstpipeline/work
projectDir                : /workspace/gitpod/hello-nextflow/hello-nf-core-template/nf-core-training-firstpipeline
userName                  : gitpod
profile                   : docker,test
configFiles               : 

!! Only displaying parameters that differ from the pipeline defaults !!
------------------------------------------------------
executor >  local (1)

[89/bce167] process > NFCORETRAINING_FIRSTPIPELINE:FIRSTPIPELINE:MULTIQC [100%] 1 of 1 ✔

-[nf-core-training/firstpipeline] Pipeline completed successfully-
Completed at: 12-Nov-2024 15:53:37
Duration    : 1m 6s
CPU hours   : (a few seconds)
Succeeded   : 1

```

Let's dissect what we are seeing:

The nf-core pipeline template is a working pipeline and comes pre-configured with some modules, here we only chose MultiQC:

-   [MultiQC](https://multiqc.info/): A modular tool to aggregate results from bioinformatics analyses across many samples into a single report.

!!! note "The template can be granularly configured"

    From nf-core tools 3.0 onwards many features can be removed during template creation. This is what we did when deselecting 
    features earlier.

At the top you see, all parameters displayed that differ from the pipeline defaults. These were configured with the `test` profile.

You can use the `test` profile to check if your pipeline is still working during your development cycle. 

The default template `test` profile leverages small test files that are stored in the nf-core [test data GitHub repository](https://github.com/nf-core/test-datasets) as inputs for the pipeline.

Additionally, the template comes with profiles for the management of software dependencies (e.g., `docker`, `singularity`, and `conda`). nf-core modules come with containers/images/recipes and profiles can be used to change the way dependencies are handled when you execute your pipeline.

!!! warning

    If `-profile` for managing software dependencies is not specified, the pipeline will run locally and expect all software to be installed and available on `PATH`. **This is not recommended.**

Additional test profiles can be created to test different parts of your pipeline.

### Template tour

The nf-core pipeline template comes packed with a lot of files and folders.

Here, we selected a subset (what we did in step 5 during the template creation).

While the template may feel overwhelming, a complete understanding isn't required to start developing your pipeline. Let's look at the important places that we need to get touch during pipeline development.

#### Workflows, subworkflows, and modules

The nf-core pipeline template has a `main.nf` script that calls `myfirstpipeline.nf` from the `workflows` folder. The `myfirstpipeline.nf` file inside the workflows folder is the central pipeline file that is used to bring everything else together.

Instead of having one large monolithic pipeline script, it's broken up into smaller script components, namely, modules and subworkflows:

-   **Modules:** Wrappers around a single process
-   **Subworkflows:** Two or more modules that are packaged together as a mini workflow

<figure class="excalidraw">
--8<-- "docs/nf_develop/img/nested.excalidraw.svg"
</figure>

Within your pipeline repository, `modules` and `subworkflows` are stored within `local` and `nf-core` folders. The `nf-core` folder is for components that have come from the nf-core GitHub repository while the `local` folder is for components that have been developed independently:

```console
modules/
├── local
│   └── <toolname>.nf
│   .
│
└── nf-core
    ├── <tool name>
    │   ├── environment.yml
    │   ├── main.nf
    │   ├── meta.yml
    │   └── tests
    │       ├── main.nf.test
    │       ├── main.nf.test.snap
    │       └── tags.yml
    .
```

Modules from nf-core follow a similar same structure and contain a small number of additional files that are used for testing using [nf-test](https://www.nf-test.com/) and documentation about the module.

!!!note

    Some nf-core modules are also split into command specific directories:

    ```console
    │
    └── <tool name>
        └── <command>
            ├── environment.yml
            ├── main.nf
            ├── meta.yml
            └── tests
                ├── main.nf.test
                ├── main.nf.test.snap
                └── tags.yml
    ```

!!!note

    The nf-core template does not come with a local modules folder by default.

#### Configuration files

The nf-core pipeline template utilizes Nextflows flexible customization options and has a series of configuration files throughout the template.

In the template, the `nextflow.config` file is a central configuration file and is used to set default values for parameters and other configuration options. The majority of these configuration options are applied by default while others (e.g., software dependency profiles) are included as optional profiles.

There are several configuration files that are stored in the `conf` folder and are added to the configuration by default or optionally as profiles:

-   `base.config`: A 'blank slate' config file, appropriate for general use on most high performance compute environments.
-   `modules.config`: Additional module directives and arguments.
-   `test.config`: A profile to run the pipeline with minimal test data.
-   `test_full.config`: A profile to run the pipeline with a full-sized test dataset.

#### `nextflow_schema.json`

The `nextflow_schema.json` is a file used to store parameter related information including type, description and help text in a machine readable format. The schema is used for various purposes, including automated parameter validation, help text generation, and interactive parameter form rendering in UI interfaces.

### Takeaway

You have now created a template pipeline, and learned about important template files

### What's next?

Congratulations! You have now created a template pipeline. In the next step, we will start adding new tools to it. 

---

## Check the input data

Above, we said that the `test` profile comes with small test files that are stored in the nf-core. Let's check what type of files we are dealing with to plan our expansion. We can inspect any channel content by using the `dump` operator: 

```groovy title="firstpipeline.nf" linenums="27"
ch_samplesheet.dump(pretty: true)
```

We can then run the pipeline with the additional argument `-dump-channels`. This will run `dump()`. It is great for development, because a normal user will not add this flag and thus not see the channel content:
```bash
nextflow run . -profile docker,test --outdir results -dump-channels
```

The output should look like this: We see that we have FastQ files as input and each pair of files is accompanied by some metadata: the `id` and whether or not it is single end:

```console title="Output"
Launching `./main.nf` [evil_mestorf] DSL2 - revision: 11a3012ba7

...

!! Only displaying parameters that differ from the pipeline defaults !!
------------------------------------------------------
executor >  local (1)

[50/de1541] NFCORETRAINING_FIRSTPIPELINE:FIRSTPIPELINE:MULTIQC [100%] 1 of 1 ✔

[DUMP] [
    {
        "id": "SAMPLE1_PE",
        "single_end": "false"
    },
    [
        "https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/illumina/amplicon/sample1_R1.fastq.gz",
        "https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/illumina/amplicon/sample1_R2.fastq.gz"
    ]
]
[DUMP] [
    {
        "id": "SAMPLE2_PE",
        "single_end": "false"
    },
    [
        "https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/illumina/amplicon/sample2_R1.fastq.gz",
        "https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/illumina/amplicon/sample2_R2.fastq.gz"
    ]
]
[DUMP] [
    {
        "id": "SAMPLE3_SE",
        "single_end": "true"
    },
    [
        "https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/illumina/amplicon/sample1_R1.fastq.gz",
        "https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/illumina/amplicon/sample2_R1.fastq.gz"
    ]
]
```

### Takeaway

You have now created a template pipeline, and learned about important template files

### What's next?

In the next step we will start changing the code and add new tools to the pipeline.

---

## Add an nf-core module

nf-core provides a large library of modules and subworkflows: pre-made nextflow wrappers around tools that can be installed into nextflow pipelines. They are designed to be flexible but may require additional configuration to suit different use cases. Currently, there are more than [1300 nf-core modules](https://nf-co.re/modules) and [60 nf-core subworkflows](https://nf-co.re/subworkflows) (November 2024) available. Modules and subworkflows can be listed, installed, updated, removed, and patched using nf-core tooling. 

While you could develop a module for this tool independently, you can save a lot of time and effort by leveraging nf-core modules and subworkflows.

Let's see which modules are there:

```console
nf-core modules list remote
```

This command lists all currently available modules, > 1300. An easier way to find them, is to go to the nf-core website and visit the modules subpage [https://nf-co.re/modules](https://nf-co.re/modules). Here you can search for modules by name or tags, find documentation for each module and which nf-core pipeline uses it:

<!-- TODO add screen grab -->

### Install an nf-core module

Now let's add another tool to the pipeline.

`Seqtk` is a fast and lightweight tool for processing sequences in the FASTA or FASTQ format. Here, you will use the [`seqtk trim`](https://github.com/lh3/seqtk) command to trim FASTQ files.

In your pipeline, you will add a new step that will take FASTQ files from the sample sheet as inputs and will produce trimmed fastq files that can be used as an input for other tools and version information about the seqtk tools to mix into the inputs for the MultiQC process.

<!-- <figure class="excalidraw">
--8<-- "docs/nf_template/img/pipeline.excalidraw.svg"
</figure> -->

The `nf-core modules install` command can be used to install the `seqtk/trim` module directly from the nf-core repository:

```
nf-core modules install
```

!!!warning

    You need to be in the myorg-myfirstpipeline directory when executing `nf-core modules install`

You can follow the prompts to find and install the module you are interested in:

```console
? Tool name: seqtk/trim
```

Once selected, the tooling will install the module in the `modules/nf-core/` folder and suggest code that you can add to your main workflow file (`workflows/mypipeline.nf`).

```console
INFO     Installing 'seqtk/trim'
INFO     Use the following statement to include this module:

include { SEQTK_TRIM } from '../modules/nf-core/seqtk/trim/main'
```

To enable reporting and reproducibility, modules and subworkflows from the nf-core repository are tracked using hashes in the `modules.json` file. When modules are installed or removed using the nf-core tooling the `modules.json` file will be automatically updated.

When you open the `modules.json`, you will see an entry for each module that is currently installed from the nf-core modules repository:

```console
"nf-core": {
    "multiqc": {
        "branch": "master",
        "git_sha": "cf17ca47590cc578dfb47db1c2a44ef86f89976d",
        "installed_by": ["modules"]
    },
    "seqtk/trim": {
        "branch": "master",
        "git_sha": "666652151335353eef2fcd58880bcef5bc2928e1",
        "installed_by": ["modules"]
    }
}
```

### Add the module to your pipeline

Although the module has been installed in your local pipeline repository, it is not yet added to your pipeline.

The suggested `include` statement needs to be added to your `workflows/mypipeline.nf` file and the process call (with inputs) needs to be added to the workflow block.

```groovy title="workflows/firstpipeline.nf" linenums="6"
include { FASTQC                 } from '../modules/nf-core/fastqc/main'
include { SEQTK_TRIM             } from '../modules/nf-core/seqtk/trim/main' 
include { MULTIQC                } from '../modules/nf-core/multiqc/main'
```

To add the `SEQTK_TRIM` module to your workflow you will need to check what inputs are required.

You can view the input channels for the module by opening the `./modules/nf-core/seqtk/trim/main.nf` file.

```groovy title="/modules/nf-core/seqtk/trim/main.nf" linenums="11"
input:
tuple val(meta), path(reads)
```

Each nf-core module also has a `meta.yml` file which describes the inputs and outputs. This meta file is rendered on the [nf-core website](https://nf-co.re/modules/seqtk_trim), or can be viewed using the `nf-core modules info` command:

```console
nf-core modules info seqtk/trim
```

It outputs a table with all defined inputs and outputs of the module: 

```console title="Output"

╭─ Module: seqtk/trim  ─────────────────────────────────────────────────────────────────────────────╮
│ Location: modules/nf-core/seqtk/trim                                                              │
│ 🔧 Tools: seqtk                                                                                   │
│ 📖 Description: Trim low quality bases from FastQ files                                           │
╰───────────────────────────────────────────────────────────────────────────────────────────────────╯
               ╷                                                                       ╷
 📥 Inputs     │Description                                                            │     Pattern
╺━━━━━━━━━━━━━━┿━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┿━━━━━━━━━━━━╸
 input[0]      │                                                                       │
╶──────────────┼───────────────────────────────────────────────────────────────────────┼────────────╴
  meta  (map)  │Groovy Map containing sample information e.g. [ id:'test',             │
               │single_end:false ]                                                     │
╶──────────────┼───────────────────────────────────────────────────────────────────────┼────────────╴
  reads  (file)│List of input FastQ files                                              │*.{fastq.gz}
               ╵                                                                       ╵
                      ╷                                                                ╷
 📥 Outputs           │Description                                                     │     Pattern
╺━━━━━━━━━━━━━━━━━━━━━┿━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┿━━━━━━━━━━━━╸
 reads                │                                                                │
╶─────────────────────┼────────────────────────────────────────────────────────────────┼────────────╴
  meta  (map)         │Groovy Map containing sample information e.g. [ id:'test',      │
                      │single_end:false ]                                              │
╶─────────────────────┼────────────────────────────────────────────────────────────────┼────────────╴
  *.fastq.gz  (file)  │Filtered FastQ files                                            │*.{fastq.gz}
╶─────────────────────┼────────────────────────────────────────────────────────────────┼────────────╴
 versions             │                                                                │
╶─────────────────────┼────────────────────────────────────────────────────────────────┼────────────╴
  versions.yml  (file)│File containing software versions                               │versions.yml
                      ╵                                                                ╵

 Use the following statement to include this module:

 include { SEQTK_TRIM } from '../modules/nf-core/seqtk/trim/main'
```

Using this module information you can work out what inputs are required for the `SEQTK_TRIM` process:

1.  `tuple val(meta), path(reads)`

    -   A tuple with a meta _map_ and a list of FASTQ _files_
    -   The channel `ch_samplesheet` used by the `FASTQC` process can be used as the reads input.

As only one input channel required, and it already exists, it can be added to your `firstpipeline.nf` file without any additional channel creation or modifications.

```groovy title="workflows/myfirstpipeline.nf" linenums="30"
//
// MODULE: Run SEQTK_TRIM
//
SEQTK_TRIM (
    ch_samplesheet
)
```

Let's test, that it works:

```bash
nextflow run . -profile docker,test --outdir results
```

```console title="Output"
Launching `./main.nf` [grave_lagrange] DSL2 - revision: 11a3012ba7

Input/output options
  input                     : https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/samplesheet/samplesheet_test_illumina_amplicon.csv
  outdir                    : results

Institutional config options
  config_profile_name       : Test profile
  config_profile_description: Minimal test dataset to check pipeline function

Core Nextflow options
  runName                   : grave_lagrange
  containerEngine           : docker
  launchDir                 : /workspace/gitpod/hello-nextflow/hello-nf-core-template/nf-core-training-firstpipeline
  workDir                   : /workspace/gitpod/hello-nextflow/hello-nf-core-template/nf-core-training-firstpipeline/work
  projectDir                : /workspace/gitpod/hello-nextflow/hello-nf-core-template/nf-core-training-firstpipeline
  userName                  : gitpod
  profile                   : docker,test
  configFiles               : 

!! Only displaying parameters that differ from the pipeline defaults !!
------------------------------------------------------
executor >  local (4)
[0d/725202] process > NFCORETRAINING_FIRSTPIPELINE:FIRSTPIPELINE:SEQTK_TRIM (SAMPLE3_SE) [100%] 3 of 3 ✔
[0a/977475] process > NFCORETRAINING_FIRSTPIPELINE:FIRSTPIPELINE:MULTIQC                 [100%] 1 of 1 ✔
-[nf-core-training/firstpipeline] Pipeline completed successfully-
```

### Inspect results folder

nf-core by default, publishes the output of each process into the `<outdir>/<TOOL>`. After running the previous command, you 
should have a `results` folder that looks like this:

```console
results
├── multiqc
│   ├── multiqc_data
│   └── multiqc_report.html
├── pipeline_info
│   ├── execution_report_2024-11-14_10-31-26.html
│   ├── execution_timeline_2024-11-14_10-31-26.html
│   ├── execution_trace_2024-11-14_10-31-26.txt
│   ├── params_2024-11-14_10-31-27.json
│   ├── pipeline_dag_2024-11-14_10-31-26.html
│   └── pipeline_software_mqc_versions.yml
└── seqtk
    ├── SAMPLE1_PE_sample1_R1.fastq.gz
    ├── SAMPLE1_PE_sample1_R2.fastq.gz
    ├── SAMPLE2_PE_sample2_R1.fastq.gz
    ├── SAMPLE2_PE_sample2_R2.fastq.gz
    ├── SAMPLE3_SE_sample1_R1.fastq.gz
    └── SAMPLE3_SE_sample2_R1.fastq.gz
```

The resulting files of `multiqc` and `seqtk` are published respective subdirectories. In addition, `nf-core` pipelines by default have all sorts of reporting switched on. These files are stored in the `pipeline_info` subdirectory and time-stamped so that multiple runs don't overwrite them.

### Handle modules output

As with the inputs, you can view the outputs for the module by opening the `/modules/nf-core/seqtk/trim/main.nf` file and viewing the module metadata.

```groovy title="/modules/nf-core/seqtk/trim/main.nf" linenums="13"
output:
tuple val(meta), path("*.fastq.gz"), emit: reads
path "versions.yml"                , emit: versions
```

To help with organization and readability it is beneficial to create named output channels.

For `SEQTK_TRIM`, the `reads` output could be put into a channel named `ch_trimmed`.

```groovy title="workflows/mypipeline.nf"
ch_trimmed  = SEQTK_TRIM.out.reads
```

Similarly, it is beneficial immediately mix the versions of tools into the `ch_versions` channel so they can be used as an input for the `MULTIQC` process.

```groovy title="workflows/mypipeline.nf"
ch_versions = ch_versions.mix(SEQTK_TRIM.out.versions.first())
```

!!! note

    The `first` operator is used to emit the first item from `SEQTK_TRIM.out.versions` to avoid duplication.

### Add a parameter to the `seqtk/trim` tool

To prevent changing the nf-core modules, additional configuration options can be applied to a module using scopes within configuration files.

The configuration of modules is commonly added to the `modules.conf` file in the `conf` folder. Process selectors (e.g., `withName`) are used to apply configuration to modules selectively. Process selectors must be used within the `process` scope.

Extra configuration may also be applied as directives by using `args`. You can find many examples of how arguments are added to modules in nf-core pipelines, for example, the nf-core/rnaseq [modules.config](https://github.com/nf-core/rnaseq/blob/master/conf/modules.config) file.

Add this snippet to your `conf/modules.config` file to call the tool with an additional argument: `-b 5` trims 5bp from the left end of each read:

```console title="conf/modules.config" linenums="21"
withName: 'SEQTK_TRIM' {
    ext.args = "-b 5"
}
```

Run the pipeline again and check if the new parameter is applied:

```bash
nextflow run . -profile docker,test --outdir results

[fb/8a18ce] process > NFCORETRAINING_FIRSTPIPELINE:FIRSTPIPELINE:SEQTK_TRIM (SAMPLE3_SE) [100%] 3 of 3 ✔
[e4/039743] process > NFCORETRAINING_FIRSTPIPELINE:FIRSTPIPELINE:MULTIQC                 [100%] 1 of 1 ✔
```

Copy the hash, that you see in your console output (here `fb/8a18ce`, it is different for each run). Use tab-completion to expand the complete hash.
In this folder you will find various log files. The `.command.sh` file contains the resolved command:

```bash
less work/fb/8a18cedc5127f9a2c26eb6579c6887/.command.sh
```

We can see, that the parameter `-b 5`, that we set in the `modules.config` is applied to the task:

```console title="Output"
#!/usr/bin/env bash

set -e # Exit if a tool returns a non-zero status/exit code
set -u # Treat unset variables and parameters as an error
set -o pipefail # Returns the status of the last command to exit with a non-zero status or zero if all successfully execute
set -C # No clobber - prevent output redirection from overwriting files.

printf "%s\n" sample1_R1.fastq.gz sample2_R1.fastq.gz | while read f;
do
    seqtk \
        trimfq \
        -b 5 \
        $f \
        | gzip --no-name > SAMPLE3_SE_$(basename $f)
done

cat <<-END_VERSIONS > versions.yml
"NFCORETRAINING_FIRSTPIPELINE:FIRSTPIPELINE:SEQTK_TRIM":
    seqtk: $(echo $(seqtk 2>&1) | sed 's/^.*Version: //; s/ .*$//')
END_VERSIONS
```

### Takeaway

You have now added a nf-core/module to your pipeline, configured it with a particular parameter, and made the output available in the workflow.

### What's next?

In the next step we will add a pipeline parameter to allow users to skip the trimming step.

---

## Adding parameters to your pipeline

Parameters that can be overridden, either using the command line or the Nextflow configuration file, and should be used for anything that a pipeline user may want to configure regularly.

Here, as a simple example, you will add a new parameter to your pipeline that will skip the `SEQTK_TRIM` process.

Parameters are accessible in the pipeline script.

### Default values

In the nf-core template the default values for parameters are set in the `nextflow.config` in the base repository.

Any new parameters should be added to the `nextflow.config` with a default value within the `params` scope.

Parameter names should be unique and easily identifiable.

We can a new parameter `skip_trim` to your `nextflow.config` file and set it to `false`.

```groovy title="nextflow.config" linenums="21"
// Trimming
skip_trim                   = false
```

### Adding parameters to your pipeline 

Here, an `if` statement that is depended on the `skip_trim` parameter can be used to control the execution of the `SEQTK_TRIM` process. An `!` can be used to imply the logical "not".

Thus, if the `skip_trim` parameter is **not** `true`, the `SEQTK_TRIM` will be be executed.

```groovy title="workflows/mypipeline.nf" linenums="37"
//
// MODULE: Run SEQTK_TRIM
//
if (!params.skip_trim) {
    SEQTK_TRIM (
        ch_samplesheet
    )
    ch_trimmed  = SEQTK_TRIM.out.reads
    ch_versions = ch_versions.mix(SEQTK_TRIM.out.versions.first())
}
```

Now your if statement has been added to your main workflow file and has a default setting in your `nextflow.config` file, you will be able to flexibly skip the new trimming step using the `skip_trim` parameter.

We can now run the pipeline with the new `skip_trim` parameter to check it is working:

```console
cd /workspace/gitpod/nf-develop/
nextflow run myorg-myfirstpipeline -profile test,docker --outdir results --skip_trim
```

```console title="Output"
!! Only displaying parameters that differ from the pipeline defaults !!
------------------------------------------------------
WARN: The following invalid input values have been detected:

* --skip_trim: true

executor >  local (1)
[fe/cd373d] process > NFCORETRAINING_FIRSTPIPELINE:FIRSTPIPELINE:MULTIQC [100%] 1 of 1 ✔
-[nf-core-training/firstpipeline] Pipeline completed successfully-
```

You should see that the `SEQTK_TRIM` process has been skipped in your execution.

### Validate input parameters

When we ran the pipeline, we saw a warning message:

```console
WARN: The following invalid input values have been detected:

* --skip_trim: true
```

Parameters are validated through the `nextflow_schema.json` file. This file is also used by the nf-core website (for example in [nf-core/mag](https://nf-co.re/mag/3.2.1/parameters/)) to render the parameter documentation, and to print the pipeline help message (`nextflow run . --help`). If you have added parameters and they have not been documented in the `nextflow_schema.json` file then the input validation does not recognize the parameter.

The `nextflow_schema.json` file can get very big and very complicated very quickly.

The `nf-core pipelines schema build` command is designed to support developers write, check, validate, and propose additions to your `nextflow_schema.json` file.

```console
nf-core pipelines schema build
```

It will enable you to launch a web builder to edit this file in your web browser rather than trying to edit this file manually.

```console
INFO     [✓] Default parameters match schema validation
INFO     [✓] Pipeline schema looks valid (found 20 params)
✨ Found 'params.skip_trim' in the pipeline config, but not in the schema. Add to pipeline schema? [y/n]: y
INFO     Writing schema with 21 params: 'nextflow_schema.json'
🚀  Launch web builder for customization and editing? [y/n]: y
```

Using the web builder you can add add details about your new parameters.

The parameters that you have added to your pipeline will be added to the bottom of the `nf-core schema build` file. Some information about these parameters will be automatically filled based on the default value from your `nextflow.config`. You will be able to categorize your new parameters into a group, add icons, and add descriptions for each.

![Pipeline parameters](img/schemabuild.png)

!!!note

    Ungrouped parameters in schema will cause a warning.

Once you have made your edits you can click `Finished` and all changes will be automatically added to your `nextflow_schema.json` file.

If you rerun the previous command, the warning should disappear:

```console
cd /workspace/gitpod/nf-develop/
nextflow run myorg-myfirstpipeline -profile test,docker --outdir results --skip_trim
```

### Takeaway

You have added a new parameter to the pipeline.

### What's next?

In the next step we will take a look at how we track additional information to an input file.

---

## Meta maps 

Datasets often have additional information that is relevant for the analysis, such as a sample name, information about sequencing protocols, or other conditions that are needed in the pipeline to process certain samples together, determine their output name, or adjust parameters.

nf-core tracks this type of information in `meta` maps. These are `key`-`value` pairs that are passed into modules together with the files. We already saw this briefly, when inspecting the `input` for `seqtk`:

```groovy title="/modules/nf-core/seqtk/trim/main.nf" linenums="11"
input:
tuple val(meta), path(reads)
```

If we run the pipeline again with `-dump-channels`, we can take a look at the current content of the `meta` maps:

```console
    {
        "id": "SAMPLE1_PE",
        "single_end": "false"
    },
```

You can add any field, that you like to the `meta` map. By default, nf-core modules expect an `id` field. 

### Takeaway

You know that a `meta` map is used to pass along additional information for a sample. 

### What's next?

In the next step we will take a look how we can add a new key to the `meta` map using the samplesheet.

---

## Simple Samplesheet adaptations

nf-core pipelines typically use samplesheets as inputs to the pipelines. This allows us to:

    - validate each entry and print specific error messages
    - attach information to each input file
    - track which datasets are processed 

Samplesheets are comma-separated text files with a header row specifying the column names, followed by one entry per row. For example, the samplesheet that we have been using during this session looks like this:

```console "samplesheet_test_illumina_amplicon.csv"
sample,fastq_1,fastq_2
SAMPLE1_PE,https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/illumina/amplicon/sample1_R1.fastq.gz,https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/illumina/amplicon/sample1_R2.fastq.gz
SAMPLE2_PE,https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/illumina/amplicon/sample2_R1.fastq.gz,https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/illumina/amplicon/sample2_R2.fastq.gz
SAMPLE3_SE,https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/illumina/amplicon/sample1_R1.fastq.gz,
SAMPLE3_SE,https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/illumina/amplicon/sample2_R1.fastq.gz,
```

The structure of the samplesheet is specified in its own schema file in `assets/schema_input.json`. Each column has its own entry together with information about the column:

```console title="schema_input.json"
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
        "pattern": "^\\S+\\.f(ast)?q\\.gz$",
        "errorMessage": "FastQ file for reads 1 must be provided, cannot contain spaces and must have extension '.fq.gz' or '.fastq.gz'"
    },
    "fastq_2": {
        "type": "string",
        "format": "file-path",
        "exists": true,
        "pattern": "^\\S+\\.f(ast)?q\\.gz$",
        "errorMessage": "FastQ file for reads 2 cannot contain spaces and must have extension '.fq.gz' or '.fastq.gz'"
    }
},
"required": ["sample", "fastq_1"]
```

This validates that the samplesheet has at least two columns: `sample` and `fastq1` (`"required": ["sample", "fastq_1"]`). It also checks that `fastq1` and `fastq2` are files, and that the file endings match a particular pattern. 
Lastly, `sample` is information about the files that we want to attach and pass along the pipeline. nf-core uses `meta` maps for this: objects that have a key and a value. We can indicate this in the schema file directly by using the meta field:

```console 
    "sample": {
        "type": "string",
        "pattern": "^\\S+$",
        "errorMessage": "Sample name must be provided and cannot contain spaces",
        "meta": ["id"]
    },
```

This sets the key name as `id` and the value that is in the `sample` column, for example `SAMPLE1_PE`:

```console title=meta
[id: SAMPLE1_PE] 
```

By adding a new entry into the json schema, we can attach additional meta information that we want to track. This will automtically validate it for us and add it to the meta map.

Let's add some new meta information, like the `machineid` as an optional column:

```console title="schema_input.json"
"properties": {
    "sample": {
        "type": "string",
        "pattern": "^\\S+$",
        "errorMessage": "Sample name must be provided and cannot contain spaces",
        "meta": ["id"]
    },
    "machineid": {
        "type": "string",
        "pattern": "^\\S+$",
        "meta": ["machineid"]
    },
    "fastq_1": {
        "type": "string",
        "format": "file-path",
        "exists": true,
        "pattern": "^\\S+\\.f(ast)?q\\.gz$",
        "errorMessage": "FastQ file for reads 1 must be provided, cannot contain spaces and must have extension '.fq.gz' or '.fastq.gz'"
    },
    "fastq_2": {
        "type": "string",
        "format": "file-path",
        "exists": true,
        "pattern": "^\\S+\\.f(ast)?q\\.gz$",
        "errorMessage": "FastQ file for reads 2 cannot contain spaces and must have extension '.fq.gz' or '.fastq.gz'"
    }
},
"required": ["sample", "fastq_1"]
```

We can now run our normal tests with the old samplesheet. Let's add `-dump-channels` again to inspect the channel content:

```console
nextflow run . -profile docker,test --outdir results -dump-channels
```

The meta map now has a new key `machineid`, that is empty because we did not specify a value yet:

```console title="Output"
[DUMP] [
    {
        "id": "SAMPLE1_PE",
        "machineid": [
            
        ],
        "single_end": "false"
    },
    [
        "https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/illumina/amplicon/sample1_R1.fastq.gz",
        "https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/illumina/amplicon/sample1_R2.fastq.gz"
    ]
```

We have also prepared a new samplesheet, that has the `machineid` column. You can overwrite the existing input with this command:

```console
nextflow run . -profile docker,test --outdir results --input -dump-channels
```

This populates the `machineid` and we could access it in the pipeline:

```console
[DUMP] [
    {
        "id": "SAMPLE3_SE",
        "machineid": "myfavorite_machine",
        "single_end": "true"
    },
    [
        "https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/illumina/amplicon/sample1_R1.fastq.gz",
        "https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/illumina/amplicon/sample2_R1.fastq.gz"
    ]
```

### Use the new meta key in the pipeline

<!-- TODO No good idea yet on what to do here -->

### Takeaway

You know how to adapt the samplesheet to add new meta information to your files. 

### What's next?

In the next step we will add a module that is not yet in nf-core.

---

## Create a custom module for your pipeline

nf-core offers a comprehensive set of modules that have been created and curated by the community. However, as a developer, you may be interested in bespoke pieces of software that are not apart of the nf-core repository or customizing a module that already exists.

In this instance, we will write a local module for the QC Tool [FastQE](https://fastqe.com/), which computes stats for FASTQ files and print those stats as emoji.

This section should feel familiar to the `hello_modules` section.

### Create the module

Start by using the nf-core tooling to create a sceleton local module. It will prompt you to type in the tool name `fastqe`, for the remaining fields press `enter` to accpet the default: 

```console
nf-core modules create
```

```console title="Output"
INFO     Repository type: pipeline                                                                                                               
INFO     Press enter to use default values (shown in brackets) or type your own responses. ctrl+click underlined text to open links.             
Name of tool/subtool: fastqe
INFO     Using Bioconda package: 'bioconda::fastqe=0.3.3'                                                                                        
INFO     Using Docker container: 'biocontainers/fastqe:0.3.3--pyhdfd78af_0'                                                                      
INFO     Using Singularity container: 'https://depot.galaxyproject.org/singularity/fastqe:0.3.3--pyhdfd78af_0'                                   
GitHub Username: (@<your-name>): 
INFO     Provide an appropriate resource label for the process, taken from the nf-core pipeline template.                                        
         For example: process_single, process_low, process_medium, process_high, process_long, process_high_memory                               
? Process resource label: process_single
INFO     Where applicable all sample-specific information e.g. 'id', 'single_end', 'read_group' MUST be provided as an input via a Groovy Map    
         called 'meta'. This information may not be required in some instances, for example indexing reference genome files.                     
Will the module require a meta map of sample information? [y/n] (y): 
INFO     Created component template: 'fastqe'                                                                                                    
INFO     Created following files:                                                                                                                
           modules/local/fastqe.nf     
```

This will create a new file in `modules/local/fastqe.nf` that already contains the container and conda definitions, the general structure of the process, and a number of TODO statements to guide you through the adaptation. 

You will notice, that it still calls `samtools` and the input are `bam`.

From our sample sheet, we know we have fastq files instead, so let's change the input definition accordingly:

```groovy title="fastqe.nf" linenums="38"
tuple val(meta), path(reads)
```

The output of this tool is a tsv file with the emoji annotation, let's adapt the output as well:

```groovy title="fastqe.nf" linenums="42"
tuple val(meta), path("*.tsv"), emit: tsv
```

The script section still calls `samtools`. Let's change this to the proper call of the tool:

```groovy title="fastqe.nf" linenums="62"
    fastqe \\
        $args \\
        $reads \\
        --output ${prefix}.tsv
```

And at last, we need to adapt the version retrieval. This tool does not have a version command, so we will add the release number manualy:

```groovy title="fastqe.nf" linenums="52"
    def VERSION = '0.3.3'
```

and write it to a file in the script section:

```groovy title="fastqe.nf" linenums="68"
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastqe: $VERSION
    END_VERSIONS
```

We will not cover `stubs` in this training, but look at them at a later point. They are not necessary to run a module, so let's remove them for now and delete:

```groovy title="fastqe.nf" linenums="74"
stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // TODO nf-core: A stub section should mimic the execution of the original module as best as possible
    //               Have a look at the following examples:
    //               Simple example: https://github.com/nf-core/modules/blob/818474a292b4860ae8ff88e149fbcda68814114d/modules/nf-core/bcftools/annotate/main.nf#L47-L63
    //               Complex example: https://github.com/nf-core/modules/blob/818474a292b4860ae8ff88e149fbcda68814114d/modules/nf-core/bedtools/split/main.nf#L38-L54
    """
    touch ${prefix}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastqe: \$(samtools --version |& sed '1!d ; s/samtools //')
    END_VERSIONS
    """
```

### Include the module into the pipeline

The module is now ready in your `modules/local` folder, but not yet included into the pipeline. Similar to `seqtk/trim` we need to add it to `workflows/myfirstpipeline.nf`:

```groovy title="myfirstpipeline.nf" linenums="6"
    include { FASTQE                 } from '../modules/local/fastqe'
```

and call it on our input data:

```groovy title="myfirstpipeline.nf" linenums="42"
    FASTQE(ch_samplesheet)
    ch_versions = ch_versions.mix(FASTQE.out.versions.first())
```

Let's run the pipeline again:

```console
nextflow run . -profile docker,test --outdir results
```

In the results folder, you should see a new subdirectory `fastqe/`, with the mean read qualities:

```console title="SAMPLE1_PE.tsv"
Filename	Statistic	Qualities
sample1_R1.fastq.gz	mean	😝 😝 😝 😝 😝 😁 😁 😁 😁 😁 😁 😁 😁 😁 😁 😁 😁 😁 😁 😁 😁 😁 😁 😁 😁 😁 😁 😁 😁 😁 😁 😁 😁 😁 😁 😁 😁 😁 😁 😁 😁 😁 😁 😁 😁 😁 😁 😁 😁 😁 😁 😁 😁 😁 😁 😁 😁 😁 😁 😁 😁 😁 😁 😁 😁 😁 😁 😁 😁 😁 😁 😁 😁 😁 😁 😁 😁 😁 😁 😁 😁 😁 😁 😁 😁 😁 😁 😁 😁 😁 😁 😁 😁 😁 😁 😁 😁 😁 😁 😁 😁 😁 😁 😁 😁 😁 😁 😁 😉 😉 😜 😜 😜 😉 😉 😁 😁 😁 😁 😁 😁 😁 😁 😁 😁 😁 😁 😁 😁 😁 😁 😁 😁 😁 😁 😁 😁 😁 😁 😁 😁 😁 😁 😁 😁 😁 😁 😁 😁 😁 😁 😁 😁 😁 😁 😁 😁 😁 😁 😁 😁 😁 😁 😁 😁 😁 😁 😁 😁 😁 😁 😉 😛 😜 😉 😉 😉 😉 😜 😜 😉 😉 😉 😉 😉 😁 😁 😁 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😜 😉 😉 😉 😉 😉 😜 😜 😜 😜 😜 😜 😜 😜 😜 😜 😜 😜 😜 😜 😛 😜 😜 😛 😛 😛 😚
sample1_R2.fastq.gz	mean	😌 😌 😌 😝 😝 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😜 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😉 😜 😉 😉 😜 😜 😉 😜 😜 😜 😜 😜 😜 😜 😜 😜 😜 😜 😜 😛 😜 😜 😜 😛 😜 😜 😜 😜 😛 😜 😛 😛 😛 😛 😛 😛 😛 😛 😛 😛 😛 😛 😝 😛 😝 😝 😝 😝 😝 😝 😝 😝 😝 😝 😝 😝 😝 😝 😌 😌 😌 😌 😌 😌 😌 😌 😌 😌 😌 😌 😌 😌 😌 😌 😌 😌 😌 😋 😋 😋 😋 😋 😋 😋 😋 😀
```

### Takeaway

You know how to add a new module that is not yet available in nf-core.

!!! note "New module contributions are always welcome!"

    If you have a module that you would like to contribute back to the commmunity, reach out on the nf-core slack or just open a pull request to the modules repository.

---

## Takeaway

You know how to use the nf-core tooling to create a new pipeline, add modulea to it, apply tool and pipeline parameters, and adapt the samplesheet. 

## What's next?

Celebrate and take another break!