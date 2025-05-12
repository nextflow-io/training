# Part 2: Rewrite Hello for nf-core

In this second part of the Hello nf-core training course, we show you how to create an nf-core-compliant pipeline version of the pipeline produced by the [Hello Nextflow](../hello_nextflow/index.md) course.

You'll have noticed in the first section of the training that nf-core pipelines follow a fairly elaborate structure with a lot of accessory files.
Creating all that from scratch would be very tedious, so the nf-core community has developed tooling to do it from a template instead, to bootstrap the process.

We are going to show you how to use this tooling to create a pipeline scaffold, then adapt existing 'regular' pipeline code onto the nf-core scaffold.

!!! note
The nf-core-tools package is pre-installed for you in our training environment.
If you are using a different environment, you need to check whether the package is installed (run `nf-core --help` in your terminal) and if not, install it as described here: https://nf-co.re/docs/nf-core-tools/installation.

---

## 1. Create a new pipeline project

First, we create the scaffold for the new pipeline.

!!! note
Make sure you are in the `hello_nf-core` directory in your terminal.

### 1.1. Run the template-based pipeline creation tool

Let's start by creating a new pipeline with the `nf-core pipelines create` command.
This will create a new pipeline scaffold using the nf-core base template, customized with a pipeline name, description, and author.

```bash
nf-core pipelines create
```

Running this command will open a Text User Interface (TUI) for pipeline creation:

<div style="text-align: center;">
    <iframe width="560" height="315" src="https://www.youtube.com/embed/VwjXNXONHlY?si=d0HkFSISnKn76TeI" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" allowfullscreen="" data-ruffle-polyfilled=""></iframe>
</div>

This TUI will ask you to provide basic information about your pipeline and will provide you with a choice of features to include or exclude in your pipeline scaffold.

1. Select **Let's go!** on the welcome screen
2. Select **Custom** on the `Choose pipeline type` screen
3. Enter your pipeline details, replacing < YOUR NAME > with your own name, then select **Next**

- **GitHub organisation:** myorg
- **Workflow name:** hello
- **A short description of your pipeline:** nf-core compliant version of Hello Nextflow
- **Name of the main author / authors:** < YOUR NAME >

4. On the Template features screen, set "Toggle all features" to **off**, then selectively **enable** the following:

- `Add configuration files`
- `Use nf-core components`
- `Use nf-schema`
- `Add documentation`
- `Add testing profiles`

5. Select **Finish** on the Final details screen
6. Wait for the pipeline to be created, then select **Continue**
7. Select **Finish without creating a repo** on the Create GitHub repository screen
8. Select **Close** on the HowTo create a GitHub repository page

Once the TUI closes, you should see the following console output.

```console title="Output"
                                          ,--./,-.
          ___     __   __   __   ___     /,-._.--~\
    |\ | |__  __ /  ` /  \ |__) |__         }  {
    | \| |       \__, \__/ |  \ |___     \`-._,-`-,
                                          `._,._,'

    nf-core/tools version 3.2.1 - https://nf-co.re


INFO     Launching interactive nf-core pipeline creation tool.
```

There is no explicit confirmation in the console output that the pipeline creation worked, but you should see a new directory called `myorg-hello`.

View the contents of the new directory to see how much work you saved yourself by using the template.

```bash
tree myorg-hello
```

```console title="Output"
myorg-hello
├── assets
│   ├── samplesheet.csv
│   └── schema_input.json
├── conf
│   ├── base.config
│   ├── modules.config
│   ├── test.config
│   └── test_full.config
├── docs
│   ├── output.md
│   ├── README.md
│   └── usage.md
├── main.nf
├── modules.json
├── nextflow.config
├── nextflow_schema.json
├── README.md
├── subworkflows
│   ├── local
│   │   └── utils_nfcore_hello_pipeline
│   │       └── main.nf
│   └── nf-core
│       ├── utils_nextflow_pipeline
│       │   ├── main.nf
│       │   ├── meta.yml
│       │   └── tests
│       │       ├── main.function.nf.test
│       │       ├── main.function.nf.test.snap
│       │       ├── main.workflow.nf.test
│       │       ├── nextflow.config
│       │       └── tags.yml
│       ├── utils_nfcore_pipeline
│       │   ├── main.nf
│       │   ├── meta.yml
│       │   └── tests
│       │       ├── main.function.nf.test
│       │       ├── main.function.nf.test.snap
│       │       ├── main.workflow.nf.test
│       │       ├── main.workflow.nf.test.snap
│       │       ├── nextflow.config
│       │       └── tags.yml
│       └── utils_nfschema_plugin
│           ├── main.nf
│           ├── meta.yml
│           └── tests
│               ├── main.nf.test
│               ├── nextflow.config
│               └── nextflow_schema.json
└── workflows
    └── hello.nf

14 directories, 36 files
```

That's a lot of files!

<!-- TODO: add some commentary tying this back to what we covered in Part 1 -->

!!! note
One important difference compared to the `nf-core/demo` pipeline we examined in the first part of this training is that there is no `modules` directory.
This is because we didn't include any of the default nf-core modules.

### 1.2. Test that the scaffold is functional

Believe it or not, even though you haven't yet added any modules to make it do real work, the pipeline scaffold can actually be run using the test profile, the same way we ran the `nf-core/demo` pipeline.

```bash
nextflow run myorg-hello -profile docker,test --outdir results
```

```console title="Output"
 N E X T F L O W   ~  version 24.10.4

Launching `myorg-hello/main.nf` [naughty_babbage] DSL2 - revision: c0376c97f3

Input/output options
  input                     : https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/samplesheet/samplesheet_test_illumina_amplicon.csv
  outdir                    : test_hello

Institutional config options
  config_profile_name       : Test profile
  config_profile_description: Minimal test dataset to check pipeline function

Generic options
  trace_report_suffix       : 2025-05-12_14-46-59

Core Nextflow options
  runName                   : naughty_babbage
  launchDir                 : /workspaces/training/hello-nf-core
  workDir                   : /workspaces/training/hello-nf-core/work
  projectDir                : /workspaces/training/hello-nf-core/myorg-hello
  userName                  : root
  profile                   : docker,test
  configFiles               :

!! Only displaying parameters that differ from the pipeline defaults !!
------------------------------------------------------
-[myorg/hello] Pipeline completed successfully-
```

This shows you that all the basic wiring is in place.

### 1.3. Examine the placeholder workflow

If you look inside the `main.nf` file, you'll see it imports a workflow called `HELLO` from `workflows/hello`.
This is a placeholder workflow for our workflow of interest, with some nf-core functionality already in place.

```groovy title="conf/test.config" linenums="1" hl_lines="15,17,19,35"
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow HELLO {

    take:
    ch_samplesheet // channel: samplesheet read in from --input
    main:

    ch_versions = Channel.empty()

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name:  'hello_software_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }


    emit:
    versions       = ch_versions                 // channel: [ path(versions.yml) ]

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
```

Compared to a basic Nextflow workflow like the one developed in Hello Nextflow, you'll notice a few things that are new here:

- The workflow block has a name
- Workflow inputs are declared using the `take:` keyword
- Workflow content is placed inside a `main:` block
- Outputs are declared using the `emit:` keyword

These are optional features of Nextflow that make the workflow **composable**, meaning that it can be called from within another workflow.

We are going to need to plug the relevant logic from our workflow of interest into that structure.

### Takeaway

You now know how to create a pipeline scaffold using nf-core tools.

### What's next?

Learn how to make a simple workflow composable.

---

## 2. Make the original Hello Nextflow workflow composable

We provide you with a clean, fully functional copy of the completed Hello Nextflow workflow in the directory `original-hello` along with its modules and the default CSV file it expects to use as input.

```bash
tree original-hello/
```

```console title="Output"
original-hello/
├── hello.nf
├── modules
│   ├── collectGreetings.nf
│   ├── convertToUpper.nf
│   ├── cowpy.nf
│   └── sayHello.nf
└── nextflow.config
```

Feel free to run it to satisfy yourself that it works:

```bash
nextflow run original-hello/hello.nf
```

```console title="Output"
 N E X T F L O W   ~  version 24.10.4

Launching `hello-original/hello.nf` [goofy_babbage] DSL2 - revision: e9e72441e9

executor >  local (8)
[a4/081cec] sayHello (1)       | 3 of 3 ✔
[e7/7e9058] convertToUpper (3) | 3 of 3 ✔
[0c/17263b] collectGreetings   | 1 of 1 ✔
[94/542280] cowpy              | 1 of 1 ✔
There were 3 greetings in this batch
```

### 2.1. Add `take`, `main` and `emit` statements

TODO: instructions

### 2.2. Test the workflow

TODO: instructions

### Takeaway

You know how to [...].

!!!note
If you're interested in digging deeper into options for composing workflows of workflows, check out the [Workflow of Workflows](https://training.nextflow.io/latest/side_quests/workflows_of_workflows) (a.k.a. WoW) side quest.

### What's next?

Learn how to [...].

---

## 3. Fit the updated workflow logic into the placeholder workflow

TODO: overview

### 3.1. Plug in the main workflow logic

TODO: instructions

### 3.2. Copy over the modules

TODO: instructions

### 3.3. Set up the module imports

TODO: instructions

### Takeaway

You know how to [...].

### What's next?

Learn how to [...].

---

## 4. Set up inputs, parameters and validation

TODO: instructions

### Takeaway

You know how to [...].

### What's next?

[...].

---

## 5. Wire up tool version output <!-- is this a separate point? -->

TODO: instructions

### Takeaway

You know how to [...].

### What's next?

[...].

---

## 6. Test the pipeline

TODO: instructions

### Takeaway

You know how to [...].

### What's next?

[...].
