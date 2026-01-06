# Part 1: Run a demo pipeline

In this first part of the Hello nf-core training course, we show you how to find and try out an nf-core pipeline, understand how the code is organized, and recognize how it differs from plain Nextflow code as shown in [Hello Nextflow](../hello_nextflow/index.md).

We are going to use a pipeline called nf-core/demo that is maintained by the nf-core project as part of its inventory of pipelines for demonstrating code structure and tool operations.

Make sure your working directory is set to `hello-nf-core/` as instructed on the [Getting started](./00_orientation.md) page.

---

## 1. Find and retrieve the nf-core/demo pipeline

Let's start by locating the nf-core/demo pipeline on the project website at [nf-co.re](https://nf-co.re), which centralizes all information such as: general documentation and help articles, documentation for each of the pipelines, blog posts, event announcements and so forth.

### 1.1. Find the pipeline on the website

In your web browser, go to [https://nf-co.re/pipelines/](https://nf-co.re/pipelines/) and type `demo` in the search bar.

![search results](./img/search-results.png)

Click on the pipeline name, `demo`, to access the pipeline documentation page.

Each released pipeline has a dedicated page that includes the following documentation sections:

- **Introduction:** An introduction and overview of the pipeline
- **Usage:** Descriptions of how to execute the pipeline
- **Parameters:** Grouped pipeline parameters with descriptions
- **Output:** Descriptions and examples of the expected output files
- **Results:** Example output files generated from the full test dataset
- **Releases & Statistics:** Pipeline version history and statistics

Whenever you are considering adopting a new pipeline, you should read the pipeline documentation carefully first to understand what it does and how it should be configured before attempting to run it.

Have a look now and see if you can find out:

- Which tools the pipeline will run (Check the tab: `Introduction`)
- Which inputs and parameters the pipeline accepts or requires (Check the tab: `Parameters`)
- What are the outputs produced by the pipeline (Check the tab: `Output`)

#### 1.1.1. Pipeline overview

The `Introduction` tab provides an overview of the pipeline, including a visual representation (called a subway map) and a list of tools that are run as part of the pipeline.

![pipeline subway map](./img/nf-core-demo-subway-cropped.png)

1. Read QC (FASTQC)
2. Adapter and quality trimming (SEQTK_TRIM)
3. Present QC for raw reads (MULTIQC)

#### 1.1.2. Example command line

The documentation also provides an example input file (discussed further below) and an example command line.

```bash
nextflow run nf-core/demo \
  -profile <docker/singularity/.../institute> \
  --input samplesheet.csv \
  --outdir <OUTDIR>
```

You'll notice that the example command does NOT specify a workflow file, just the reference to the pipeline repository, `nf-core/demo`.

When invoked this way, Nextflow will assume that the code is organized in a certain way.
Let's retrieve the code so we can examine this structure.

### 1.2. Retrieve the pipeline code

Once we've determined that the pipeline appears to be suitable for our purposes, let's try it out.
Fortunately Nextflow makes it easy to retrieve pipelines from correctly-formatted repositories without having to download anything manually.

Let's return to the terminal and run the following:

```bash
nextflow pull nf-core/demo
```

Nextflow will `pull` the pipeline code, meaning it will download the full repository to your local drive.

```console title="Output"
Checking nf-core/demo ...
 downloaded from https://github.com/nf-core/demo.git - revision: 04060b4644 [master]
```

To be clear, you can do this with any Nextflow pipeline that is appropriately set up in GitHub, not just nf-core pipelines.
However nf-core is the largest open-source collection of Nextflow pipelines.

You can get Nextflow to give you a list of what pipelines you have retrieved in this way:

```bash
nextflow list
```

```console title="Output"
nf-core/demo
```

You'll notice that the files are not in your current work directory.
By default, Nextflow saves them to `$NXF_HOME/assets`.

```bash
tree -L 2 $NXF_HOME/assets/
```

```console title="Directory contents"
/workspaces/.nextflow/assets/
└── nf-core
    └── demo

2 directories, 0 files
```

!!! note

    The full path may differ on your system if you're not using our training environment.

Nextflow keeps the downloaded source code intentionally 'out of the way' on the principle that these pipelines should be used more like libraries than code that you would directly interact with.

However, for the purposes of this training, we want to be able to poke around and see what's in there.
So to make that easier, let's create a symbolic link to that location from our current working directory.

```bash
ln -s $NXF_HOME/assets pipelines
```

This creates a shortcut that makes it easier to explore the code we just downloaded.

```bash
tree -L 2 pipelines
```

```console title="Directory contents"
pipelines
└── nf-core
    └── demo

2 directories, 0 files
```

Now we can more easily peek into the source code as needed.

But first, let's try running our first nf-core pipeline!

### Takeaway

You now know how to find a pipeline via the nf-core website and retrieve a local copy of the source code.

### What's next?

Learn how to try out an nf-core pipeline with minimal effort.

---

## 2. Try out the pipeline with its test profile

Conveniently, every nf-core pipeline comes with a test profile.
This is a minimal set of configuration settings for the pipeline to run using a small test dataset hosted in the [nf-core/test-datasets](https://github.com/nf-core/test-datasets) repository.
It's a great way to quickly try out a pipeline at small scale.

!!! note

    Nextflow's configuration profile system allows you to easily switch between different container engines or execution environments.
    For more details, see [Hello Nextflow Part 6: Configuration](../hello_nextflow/06_hello_configuration.md).

### 2.1. Examine the test profile

It's good practice to check what a pipeline's test profile specifies before running it.
The `test` profile for `nf-core/demo` lives in the configuration file `conf/test.config` and is shown below.

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
        cpus: 4,
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

You'll notice right away that the comment block at the top includes a usage example showing how to run the pipeline with this test profile.

```groovy title="conf/test.config" linenums="7"
Use as follows:
        nextflow run nf-core/demo -profile test,<docker/singularity> --outdir <OUTDIR>
```

The only things we need to supply are what's shown between carets in the example command: `<docker/singularity>` and `<OUTDIR>`.

As a reminder, `<docker/singularity>` refers to the choice of container system. All nf-core pipelines are designed to be usable with containers (Docker, Singularity, etc.) to ensure reproducibility and eliminate software installation issues.
So we'll need to specify whether we want to use Docker or Singularity to test the pipeline.

The `--outdir <OUTDIR>` part refers to the directory where Nextflow will write the pipeline's outputs.
We need to provide a name for it, which we can just make up.
If it does not exist already, Nextflow will create it for us at runtime.

Moving on to the section after the comment block, the test profile shows us what has been pre-configured for testing: most notably, the `input` parameter is already set to point to a test dataset, so we don't need to provide our own data.
If you follow the link to the pre-configured input, you'll see it is a csv file containing sample identifiers and file paths for several experimental samples.

```csv title="samplesheet_test_illumina_amplicon.csv"
sample,fastq_1,fastq_2
SAMPLE1_PE,https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/illumina/amplicon/sample1_R1.fastq.gz,https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/illumina/amplicon/sample1_R2.fastq.gz
SAMPLE2_PE,https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/illumina/amplicon/sample2_R1.fastq.gz,https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/illumina/amplicon/sample2_R2.fastq.gz
SAMPLE3_SE,https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/illumina/amplicon/sample1_R1.fastq.gz,
SAMPLE3_SE,https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/illumina/amplicon/sample2_R1.fastq.gz,
```

This is called a samplesheet, and is the most common form of input to nf-core pipelines.

!!! note

    Don't worry if you're not familiar with the data formats and types, it's not important for what follows.

So this confirms that we have everything we need to try out the pipeline.

### 2.2. Run the pipeline

Let's decide to use Docker for the container system and `demo-results` as the output directory, and we're ready to run the test command:

```bash
nextflow run nf-core/demo -profile docker,test --outdir demo-results
```

??? success "Command output"

    ```console
     N E X T F L O W   ~  version 25.04.3

    Launching `https://github.com/nf-core/demo` [magical_pauling] DSL2 - revision: db7f526ce1 [master]


    ------------------------------------------------------
                                            ,--./,-.
            ___     __   __   __   ___     /,-._.--~'
      |\ | |__  __ /  ` /  \ |__) |__         }  {
      | \| |       \__, \__/ |  \ |___     \`-._,-`-,
                                            `._,._,'
      nf-core/demo 1.0.2
    ------------------------------------------------------
    Input/output options
      input                     : https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/samplesheet/samplesheet_test_illumina_amplicon.csv
      outdir                    : demo-results

    Institutional config options
      config_profile_name       : Test profile
      config_profile_description: Minimal test dataset to check pipeline function

    Generic options
      trace_report_suffix       : 2025-11-21_04-57-41

    Core Nextflow options
      revision                  : master
      runName                   : magical_pauling
      containerEngine           : docker
      launchDir                 : /workspaces/training/hello-nf-core
      workDir                   : /workspaces/training/hello-nf-core/work
      projectDir                : /workspaces/.nextflow/assets/nf-core/demo
      userName                  : root
      profile                   : docker,test
      configFiles               : /workspaces/.nextflow/assets/nf-core/demo/nextflow.config

    !! Only displaying parameters that differ from the pipeline defaults !!
    ------------------------------------------------------
    * The pipeline
        https://doi.org/10.5281/zenodo.12192442

    * The nf-core framework
        https://doi.org/10.1038/s41587-020-0439-x

    * Software dependencies
        https://github.com/nf-core/demo/blob/master/CITATIONS.md


    executor >  local (7)
    [ff/a6976b] NFCORE_DEMO:DEMO:FASTQC (SAMPLE3_SE)     | 3 of 3 ✔
    [39/731ab7] NFCORE_DEMO:DEMO:SEQTK_TRIM (SAMPLE3_SE) | 3 of 3 ✔
    [7c/78d96e] NFCORE_DEMO:DEMO:MULTIQC                 | 1 of 1 ✔
    -[nf-core/demo] Pipeline completed successfully-
    ```

If your output matches that, congratulations! You've just run your first nf-core pipeline.

You'll notice that there is a lot more console output than when you run a basic Nextflow pipeline.
There's a header that includes a summary of the pipeline's version, inputs and outputs, and a few elements of configuration.

!!! note

    Your output will show different timestamps, execution names, and file paths, but the overall structure and process execution should be similar.

Moving on to the execution output, let's have a look at the lines that tell us what processes were run:

```console title="Output (subset)"
    [ff/a6976b] NFCORE_DEMO:DEMO:FASTQC (SAMPLE3_SE)     | 3 of 3 ✔
    [39/731ab7] NFCORE_DEMO:DEMO:SEQTK_TRIM (SAMPLE3_SE) | 3 of 3 ✔
    [7c/78d96e] NFCORE_DEMO:DEMO:MULTIQC                 | 1 of 1 ✔
```

This tells us that three processes were run, corresponding to the three tools shown in the pipeline documentation page on the nf-core website: FASTQC, SEQTK_TRIM and MULTIQC.

The full process names as shown here, such as `NFCORE_DEMO:DEMO:MULTIQC`, are longer than what you may have seen in the introductory Hello Nextflow material.
These include the names of their parent workflows and reflect the modularity of the pipeline code.
We'll go into more detail about that in a little bit.

### 2.3. Examine the pipeline's outputs

Finally, let's have a look at the `demo-results` directory produced by the pipeline.

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
        ├── execution_report_2025-11-21_04-57-41.html
        ├── execution_timeline_2025-11-21_04-57-41.html
        ├── execution_trace_2025-11-21_04-57-41.txt
        ├── nf_core_demo_software_mqc_versions.yml
        ├── params_2025-11-21_04-57-46.json
        └── pipeline_dag_2025-11-21_04-57-41.html
    ```

That might seem like a lot.
To learn more about the `nf-core/demo` pipeline's outputs, check out its [documentation page](https://nf-co.re/demo/1.0.2/docs/output/).

At this stage, what's important to observe is that the results are organized by module, and there is additionally a directory called `pipeline_info` containing various timestamped reports about the pipeline execution.

For example, the `execution_timeline_*` file shows you what processes were run, in what order and how long they took to run:

![execution timeline report](./img/execution_timeline.png)

!!! note

    Here the tasks were not run in parallel because we are running on a minimalist machine in Github Codespaces.
    To see these run in parallel, try increasing the CPU allocation of your codespace and the resource limits in the test configuration.

These reports are generated automatically for all nf-core pipelines.

### Takeaway

You know how to run an nf-core pipeline using its built-in test profile and where to find its outputs.

### What's next?

Learn how the pipeline code is organized.

---

## 3. Examine the pipeline code structure

Now that we've successfully run the pipeline as users, let's shift our perspective to look at how nf-core pipelines are structured internally.

The nf-core project enforces strong guidelines for how pipelines are structured, and for how the code is organized, configured and documented.
Understanding how this is all organized is the first step toward developing your own nf-core-compatible pipelines, which we will tackle in Part 2 of this course.

Let's have a look at how the pipeline code is organized in the `nf-core/demo` repository, using the `pipelines` symlink we created earlier.

You can either use `tree` or use the file explorer to find and open the `nf-core/demo` directory.

```bash
tree -L 1 pipelines/nf-core/demo
```

??? abstract "Directory contents"

    ```console
    pipelines/nf-core/demo
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

There's a lot going on in there, so we'll tackle this step by step.

First, let's note that at the top level, you can find a README file with summary information, as well as accessory files that summarize project information such as licensing, contribution guidelines, citation and code of conduct.
Detailed pipeline documentation is located in the `docs` directory.
All of this content is used to generate the web pages on the nf-core website programmatically, so they're always up to date with the code.

Now, for the rest, we're going to divide our exploration in three stages:

1. Pipeline code components (`main.nf`, `workflows`, `subworkflows`, `modules`)
2. Pipeline configuration
3. Inputs and validation

Let's start with the pipeline code components.
We're going to focus on the file hierarchy and structural organization, rather than diving into the code within individual files.

### 3.1. Pipeline code components

The standard nf-core pipeline code organization follows a modular structure that is designed to maximize code reuse, as introduced in [Hello Modules](../hello_nextflow/04_hello_modules.md), Part 4 of the [Hello Nextflow](../hello_nextflow/index.md) course, although in true nf-core fashion, this is implemented with a bit of additional complexity.
Specifically, nf-core pipelines make abundant use of subworkflows, i.e. workflow scripts that are imported by a parent workflow.

That may sound a bit abstract, so let's take a look how this is used in practice in the `nf-core/demo` pipeline.

!!! note

    We won't go over the actual code for _how_ these modular components are connected, because there is some additional complexity associated with the use of subworkflows that can be confusing, and understanding that is not necessary at this stage of the training.
    For now, we're going to focus on the overall organization and logic.

#### 3.1.1. General overview

Here is what the relationships between the relevant code components look like for the `nf-core/demo` pipeline:

<figure class="excalidraw">
    --8<-- "docs/hello_nf-core/img/nf-core_demo_code_organization.svg"
</figure>

There is a so-called _entrypoint_ script called `main.nf`, which acts as a wrapper for two kinds of nested workflows: the workflow containing the actual analysis logic, located under `workflows/` and called `demo.nf`, and a set of housekeeping workflows located under `subworkflows/`.
The `demo.nf` workflow calls on **modules** located under `modules/`; these contain the **processes** that will perform the actual analysis steps.

!!! note

    Subworkflows are not limited to housekeeping functions, and they can make use of process modules.

    The `nf-core/demo` pipeline shown here happens to be on the simpler side on the spectrum, but other nf-core pipelines (such as `nf-core/rnaseq`) utilize subworkflows that are involved in the actual analysis.

Now, let's review these components in turn.

#### 3.1.2. The entrypoint script: `main.nf`

The `main.nf` script is the entrypoint that Nextflow starts from when we execute `nextflow run nf-core/demo`.
That means when you run `nextflow run nf-core/demo` to run the pipeline, Nextflow automatically finds and executes the `main.nf` script.
This works for any Nextflow pipeline that follows this conventional naming and structure, not just nf-core pipelines.

Using an entrypoint script makes it easy to run standardized 'housekeeping' subworkflows before and after the actual analysis script gets run.
We'll go over those after we've reviewed the actual analysis workflow and its modules.

#### 3.1.3. The analysis script: `workflows/demo.nf`

The `workflows/demo.nf` workflow is where the central logic of the pipeline is stored.
It is structured much like a normal Nextflow workflow, except it is designed to be called from a parent workflow, which requires a few extra features.
We'll cover the relevant differences in the next part of this course, when we tackle the conversion of the simple Hello pipeline from Hello Nextflow into an nf-core-compatible form.

The `demo.nf` workflow calls on **modules** located under `modules/`, which we'll review next.

!!! note

    Some nf-core analysis workflows display additional levels of nesting by calling on lower-level subworkflows.
    This is mostly used for wrapping two or more modules that are commonly used together into easily reusable pipeline segments.
    You can see some examples by browsing available [nf-core subworkflows](https://nf-co.re/subworkflows/) on the nf-core website.

    When the analysis script uses subworkflows, those are stored under the `subworkflows/` directory.

#### 3.1.4. The modules

The modules are where the process code lives, as described in [Part 4 of the Hello Nextflow training course](../hello_nextflow/04_hello_modules.md).

In the nf-core project, modules are organized using a multi-level nested structure that reflect both their origin and their contents.
At the top level, modules are differentiated as either `nf-core` or `local` (not part of the nf-core project), and then further placed into a directory named after the tool(s) they wrap.
If the tool belongs to a toolkit (i.e. a package containing multiple tools) then there is an intermediate directory level named after the toolkit.

You can see this applied in practice to the `nf-core/demo` pipeline modules:

```bash
tree -L 3 pipelines/nf-core/demo/modules
```

??? abstract "Directory contents"

    ```console
    pipelines/nf-core/demo/modules
    └── nf-core
        ├── fastqc
        │   ├── environment.yml
        │   ├── main.nf
        │   ├── meta.yml
        │   └── tests
        ├── multiqc
        │   ├── environment.yml
        │   ├── main.nf
        │   ├── meta.yml
        │   └── tests
        └── seqtk
            └── trim

    7 directories, 6 files
    ```

Here you see that the `fastqc` and `multiqc` modules sit at the top level within the `nf-core` modules, whereas the `trim` module sits under the toolkit that it belongs to, `seqtk`.
In this case there are no `local` modules.

The module code file describing the process is always called `main.nf`, and is accompanied by tests and `.yml` files which we'll ignore for now.

Taken together, the entrypoint workflow, analysis workflow and modules are sufficient for running the 'interesting' parts of the pipeline.
However, we know there are also housekeeping subworkflows in there, so let's look at those now.

#### 3.1.5. The housekeeping subworkflows

Like modules, subworkflows are differentiated into `local` and `nf-core` directories, and each subworkflow has its own nested directory structure with its own `main.nf` script, tests and `.yml` file.

```bash
tree -L 3 pipelines/nf-core/demo/subworkflows
```

??? abstract "Directory contents"

    ```console
    pipelines/nf-core/demo/subworkflows
    ├── local
    │   └── utils_nfcore_demo_pipeline
    │       └── main.nf
    └── nf-core
        ├── utils_nextflow_pipeline
        │   ├── main.nf
        │   ├── meta.yml
        │   └── tests
        ├── utils_nfcore_pipeline
        │   ├── main.nf
        │   ├── meta.yml
        │   └── tests
        └── utils_nfschema_plugin
            ├── main.nf
            ├── meta.yml
            └── tests

    9 directories, 7 files
    ```

As noted above, the `nf-core/demo` pipeline does not include any analysis-specific subworkflows, so all the subworkflows we see here are so-called 'housekeeping' or 'utility' workflows, as denoted by the `utils_` prefix in their names.
These subworkflows are what produces the fancy nf-core header in the console output, among other accessory functions.

!!! tip

    Aside from their naming pattern, another indication that these subworkflows do not perform any truly analysis-related function is that they do not call any processes at all.

This completes the round-up of core code components that constitute the `nf-core/demo` pipeline.
Now let's take a look at the remaining elements that you should know a little bit about before diving into development: pipeline configuration and input validation.

### 3.2. Pipeline configuration

You've learned previously that Nextflow offers many options for configuring pipeline execution, be it in terms of inputs and parameters, computing resources, and other aspects of orchestration.
The nf-core project applies highly standardized guidelines for pipeline configuration that aim to build on Nextflow's flexible customization options in a way that provides greater consistency and maintainability across pipelines.

The central configuration file `nextflow.config` is used to set default values for parameters and other configuration options.
The majority of these configuration options are applied by default while others (e.g., software dependency profiles) are included as optional profiles.

There are several additional configuration files that are stored in the `conf` folder and which can be added to the configuration by default or optionally as profiles:

- `base.config`: A 'blank slate' config file, appropriate for general use on most high-performance computing environments. This defines broad bins of resource usage, for example, which are convenient to apply to modules.
- `modules.config`: Additional module directives and arguments.
- `test.config`: A profile to run the pipeline with minimal test data, which we used when we ran the demo pipeline.
- `test_full.config`: A profile to run the pipeline with a full-sized test dataset.

We will touch a few of those files later in the course.

### 3.3. Inputs and validation

As we noted earlier, when we examined the `nf-core/demo` pipeline's test profile, it is designed to take as input a samplesheet containing file paths and sample identifiers.
The file paths linked to real data located in the `nf-core/test-datasets` repository.

An example samplesheet is also provided under the `assets` directory, although the paths in this one are not real.

```csv title="assets/samplesheet.csv" linenums="1"
sample,fastq_1,fastq_2
SAMPLE_PAIRED_END,/path/to/fastq/files/AEG588A1_S1_L002_R1_001.fastq.gz,/path/to/fastq/files/AEG588A1_S1_L002_R2_001.fastq.gz
SAMPLE_SINGLE_END,/path/to/fastq/files/AEG588A4_S4_L003_R1_001.fastq.gz,

```

This particular samplesheet is fairly simple, but some pipelines run on samplesheets that are more complex, with a lot more metadata associated with the primary inputs.

Unfortunately, because these files can be difficult to check by eye, improper formatting of input data is a very common source of pipeline failures.
A related problem is when parameters are provided incorrectly.

The solution to these problems is to run automated validation checks on all input files to ensure they contain the expected types of information, formatted correctly, and on parameters to ensure they are of the expected type.
This is called input validation, and should ideally be done _before_ trying to run a pipeline, rather than waiting for the pipeline to fail to find out there was a problem with the inputs.

Just like for configuration, the nf-core project is very opinionated about input validation, and recommends the use of the [nf-schema plugin](https://nextflow-io.github.io/nf-schema/latest/), a Nextflow plugin that provides comprehensive validation capabilities for Nextflow pipelines.

We'll cover this topic in more detail in Part 5 of this course.
For now, just be aware that there are two JSON files provided for that purpose, `nextflow_schema.json` and `assets/schema_input.json`.

The `nextflow_schema.json` is a file used to store information about the pipeline parameters including type, description and help text in a machine readable format.
This is used for various purposes, including automated parameter validation, help text generation, and interactive parameter form rendering in UI interfaces.

The `schema_input.json` is a file used to define the input samplesheet structure.
Each column can have a type, pattern, description and help text in a machine readable format.
The schema is used for various purposes, including automated validation, and providing helpful error messages.

### Takeaway

You know what are the main components of an nf-core pipeline and how the code is organized; where the main elements of configuration are located; and you're aware of what input validation is for.

### What's next?

Take a break! That was a lot. When you're feeling refreshed and ready, move on to the next section to apply what you've learned to write an nf-core compatible pipeline.

!!! tip

    If you would like to learn how to compose workflows with subworkflows before moving on to the next part, check out the [Workflows of Workflows](../side_quests/workflows_of_workflows.md) Side Quest.
