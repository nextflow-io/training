# Part 2: Rewrite Hello for nf-core

In this second part of the Hello nf-core training course, we show you how to create an nf-core compatible version of the pipeline produced by the [Hello Nextflow](../hello_nextflow/index.md) course.

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

- On the welcome screen, click **Let's go!**.
- On the `Choose pipeline type` screen, click **Custom**.
- Enter your pipeline details as follows (replacing `< YOUR NAME >` with your own name), then click **Next**.

```
[ ] GitHub organisation: core
[ ] Workflow name: hello
[ ] A short description of your pipeline: A basic nf-core style version of Hello Nextflow
[ ] Name of the main author(s): < YOUR NAME >
```

- On the Template features screen, set `Toggle all features` to **off**, then selectively **enable** the following. Check your selections and click **Continue**.

```
[ ] Add configuration files
[ ] Use nf-core components
[ ] Use nf-schema
[ ] Add documentation
[ ] Add testing profiles
```

- On the `Final details` screen, click **Finish**. Wait for the pipeline to be created, then click **Continue**.
- On the Create GitHub repository screen, click **Finish without creating a repo**. This will display instructions for creating a GitHub repository later. Ignore these and click **Close**.

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

There is no explicit confirmation in the console output that the pipeline creation worked, but you should see a new directory called `core-hello`.

View the contents of the new directory to see how much work you saved yourself by using the template.

```bash
tree core-hello
```

```console title="Output"
core-hello/
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

Don't worry too much right now about what they all are; we are going to walk through the important parts together in the course of this training.

!!! note

    One important difference compared to the `nf-core/demo` pipeline we examined in the first part of this training is that there is no `modules` directory.
    This is because we didn't include any of the default nf-core modules.

### 1.2. Test that the scaffold is functional

Believe it or not, even though you haven't yet added any modules to make it do real work, the pipeline scaffold can actually be run using the test profile, the same way we ran the `nf-core/demo` pipeline.

```bash
nextflow run ./core-hello -profile docker,test --outdir core-hello-results
```

```console title="Output"
 N E X T F L O W   ~  version 24.10.4

Launching `core-hello/main.nf` [special_ride] DSL2 - revision: c31b966b36

Downloading plugin nf-schema@2.2.0
Input/output options
  input                     : https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/samplesheet/samplesheet_test_illumina_amplicon.csv
  outdir                    : core-hello-results

Institutional config options
  config_profile_name       : Test profile
  config_profile_description: Minimal test dataset to check pipeline function

Generic options
  trace_report_suffix       : 2025-05-14_10-01-18

Core Nextflow options
  runName                   : special_ride
  containerEngine           : docker
  launchDir                 : /workspaces/training/hello-nf-core
  workDir                   : /workspaces/training/hello-nf-core/work
  projectDir                : /workspaces/training/hello-nf-core/core-hello
  userName                  : root
  profile                   : docker,test
  configFiles               :

!! Only displaying parameters that differ from the pipeline defaults !!
------------------------------------------------------
-[core/hello] Pipeline completed successfully
```

This shows you that all the basic wiring is in place.
You can take a look at the reports in the `pipeline_info` directory to see what was run; not much at all!

!!! note

    The nf-core pipeline template includes an example samplesheet, but at time of writing it is very domain-specific.
    Future work will aim to produce something more generic.

### 1.3. Examine the placeholder workflow

If you look inside the `main.nf` file, you'll see it imports a workflow called `HELLO` from `workflows/hello`.
This is a placeholder workflow for our workflow of interest, with some nf-core functionality already in place.

```groovy title="core-hello/workflows/hello.nf" linenums="1" hl_lines="15 17 19 35"
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

Compared to a basic Nextflow workflow like the one developed in Hello Nextflow, you'll notice a few things that are new here (highlighted lines above):

- The workflow block has a name
- Workflow inputs are declared using the `take:` keyword and the channel construction is moved up to the parent workflow
- Workflow content is placed inside a `main:` block
- Outputs are declared using the `emit:` keyword

These are optional features of Nextflow that make the workflow **composable**, meaning that it can be called from within another workflow.

We are going to need to plug the relevant logic from our workflow of interest into that structure.
The first step for that is to make our original workflow composable.

### Takeaway

You now know how to create a pipeline scaffold using nf-core tools.

### What's next?

Learn how to make a simple workflow composable as a prelude to making it nf-core compatible.

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

Launching `original-hello/hello.nf` [goofy_babbage] DSL2 - revision: e9e72441e9

executor >  local (8)
[a4/081cec] sayHello (1)       | 3 of 3 ✔
[e7/7e9058] convertToUpper (3) | 3 of 3 ✔
[0c/17263b] collectGreetings   | 1 of 1 ✔
[94/542280] cowpy              | 1 of 1 ✔
There were 3 greetings in this batch
```

Open the `hello.nf` workflow file to inspect the code, which is shown in full below (not counting the processes, which are in modules):

```groovy title="original-hello/hello.nf" linenums="1"
#!/usr/bin/env nextflow

/*
* Pipeline parameters
*/
params.greeting = 'greetings.csv'
params.batch = 'test-batch'
params.character = 'turkey'

// Include modules
include { sayHello } from './modules/sayHello.nf'
include { convertToUpper } from './modules/convertToUpper.nf'
include { collectGreetings } from './modules/collectGreetings.nf'
include { cowpy } from './modules/cowpy.nf'

workflow {

  // create a channel for inputs from a CSV file
  greeting_ch = Channel.fromPath(params.greeting)
                      .splitCsv()
                      .map { line -> line[0] }

  // emit a greeting
  sayHello(greeting_ch)

  // convert the greeting to uppercase
  convertToUpper(sayHello.out)

  // collect all the greetings into one file
  collectGreetings(convertToUpper.out.collect(), params.batch)

  // emit a message about the size of the batch
  collectGreetings.out.count.view { "There were $it greetings in this batch" }

  // generate ASCII art of the greetings with cowpy
  cowpy(collectGreetings.out.outfile, params.character)
}
```

Let's walk through the necessary changes one by one.

### 2.1. Name the workflow

First, let's give the workflow a name so we can refer to it from a parent workflow.

=== "After"

    ```groovy title="original-hello/hello.nf" linenums="16"
    workflow HELLO {
    ```

=== "Before"

    ```groovy title="original-hello/hello.nf" linenums="16"
    workflow {
    ```

The same conventions apply to workflow names as to module names.

### 2.2. Replace channel construction with `take`

Now, replace the channel construction with a simple `take` statement declaring expected inputs.

=== "After"

    ```groovy title="original-hello/hello.nf" linenums="18"
        take:
        // channel of greetings
        greeting_ch
    ```

=== "Before"

    ```groovy title="original-hello/hello.nf" linenums="18"
        // create a channel for inputs from a CSV file
        greeting_ch = Channel.fromPath(params.greeting)
                            .splitCsv()
                            .map { line -> line[0] }
    ```

This leaves the details of how the inputs are provided up to the parent workflow.

As part of this change, you should also delete the line `params.greeting = 'greetings.csv'` from the block of parameter definitions (line 6).
That will also be left to the parent workflow to declare.

### 2.3. Preface workflow operations with `main` statement

Next, add a `main` statement before the rest of the operations called in the body of the workflow.

=== "After"

    ```groovy title="original-hello/hello.nf" linenums="21" hl_lines="1"
        main:

        // emit a greeting
        sayHello(greeting_ch)

        // convert the greeting to uppercase
        convertToUpper(sayHello.out)

        // collect all the greetings into one file
        collectGreetings(convertToUpper.out.collect(), params.batch)

        // emit a message about the size of the batch
        collectGreetings.out.count.view { "There were $it greetings in this batch" }

        // generate ASCII art of the greetings with cowpy
        cowpy(collectGreetings.out.outfile, params.character)
    ```

=== "Before"

    ```groovy title="original-hello/hello.nf" linenums="21"
        // emit a greeting
        sayHello(greeting_ch)

        // convert the greeting to uppercase
        convertToUpper(sayHello.out)

        // collect all the greetings into one file
        collectGreetings(convertToUpper.out.collect(), params.batch)

        // emit a message about the size of the batch
        collectGreetings.out.count.view { "There were $it greetings in this batch" }

        // generate ASCII art of the greetings with cowpy
        cowpy(collectGreetings.out.outfile, params.character)
    ```

This basically says 'this is what this workflow _does_'.

### 2.4. Add `emit` statement

Finally, add an `emit` statement declaring what are the final outputs of the workflow.

```groovy title="original-hello/hello.nf" linenums="37"
    emit:
    cowpy_hellos = cowpy.out
```

This is a net new addition to the code compared to the original workflow.

### 2.5. Recap of the completed changes

If you've done all the changes as described, your workflow should now look like this:

```groovy title="original-hello/hello.nf" linenums="1" hl_lines="15 17-19 21 37-38"
#!/usr/bin/env nextflow

/*
* Pipeline parameters
*/
params.batch = 'test-batch'
params.character = 'turkey'

// Include modules
include { sayHello } from './modules/sayHello.nf'
include { convertToUpper } from './modules/convertToUpper.nf'
include { collectGreetings } from './modules/collectGreetings.nf'
include { cowpy } from './modules/cowpy.nf'

workflow HELLO {

    take:
    // channel of greetings
    greeting_ch

    main:

    // emit a greeting
    sayHello(greeting_ch)

    // convert the greeting to uppercase
    convertToUpper(sayHello.out)

    // collect all the greetings into one file
    collectGreetings(convertToUpper.out.collect(), params.batch)

    // emit a message about the size of the batch
    collectGreetings.out.count.view { "There were $it greetings in this batch" }

    // generate ASCII art of the greetings with cowpy
    cowpy(collectGreetings.out.outfile, params.character)

    emit:
    cowpy_hellos = cowpy.out
}
```

This describes everything Nextflow needs EXCEPT what to feed into the input channel.
That is going to be defined in the parent workflow, also called the **entrypoint** workflow.

### 2.6. Make a dummy entrypoint workflow

We can make a dummy entrypoint workflow to test the composable workflow without yet having to deal with the rest of the complexity of the nf-core pipeline scaffold.

Create a blank file named `main.nf` in the same`original-hello` directory.

```bash
touch original-hello/main.nf
```

Copy the following code into the `main.nf` file.

```groovy title="original-hello/main.nf" linenums="1"
#!/usr/bin/env nextflow

// import the workflow code from the hello.nf file
include { HELLO } from './hello.nf'

// declare input parameter
params.greeting = 'greetings.csv'

workflow {
  // create a channel for inputs from a CSV file
  greeting_ch = Channel.fromPath(params.greeting)
                      .splitCsv()
                      .map { line -> line[0] }

  // call the imported workflow on the channel of greetings
  HELLO(greeting_ch)

  // view the outputs emitted by the workflow
  HELLO.out.view { "Output: $it" }
}
```

There are two important observations to make here:

- The syntax for calling the imported workflow (line 16) is essentially the same as the syntax for calling modules.
- Everything that is related to pulling the inputs into the workflow (input parameter and channel construction) is now declared in this parent workflow.

!!! note

    Naming the entrypoint workflow file `main.nf` is a convention, not a requirement.

    If you follow this convention, you can omit specifying the workflow file name in your `nextflow run` command.
    Nextflow will automatically look for a file named `main.nf` in the execution directory.

    However, you can name the entrypoint workflow file something else if you prefer.
    In that case, be sure to specify the workflow file name in your `nextflow run` command.

### 2.7. Test that the workflow runs

We finally have all the pieces we need to verify that the composable workflow works.

```bash
nextflow run ./original-hello
```

!!! note

    Here you see the advantage of using the `main.nf` naming convention.
    If we had named the entrypoint workflow `something_else.nf`, we would have had to do `nextflow run original-hello/something_else.nf`.

If you made all the changes correctly, this should run to completion.

```console title="Output"
 N E X T F L O W   ~  version 24.10.4

Launching `original-hello/main.nf` [friendly_wright] DSL2 - revision: 1ecd2d9c0a

executor >  local (8)
[24/c6c0d8] HELLO:sayHello (3)       | 3 of 3 ✔
[dc/721042] HELLO:convertToUpper (3) | 3 of 3 ✔
[48/5ab2df] HELLO:collectGreetings   | 1 of 1 ✔
[e3/693b7e] HELLO:cowpy              | 1 of 1 ✔
There were 3 greetings in this batch
Output: /workspaces/training/hello-nf-core/work/e3/693b7e48dc119d0c54543e0634c2e7/cowpy-COLLECTED-test-batch-output.txt
```

This means we've successfully upgraded our HELLO workflow to be composable.

### Takeaway

You know how to make a workflow composable by giving it a name and adding `take`, `main` and `emit` statements, and how to call it from an entrypoint workflow.

!!! note

    If you're interested in digging deeper into options for composing workflows of workflows, check out the [Workflow of Workflows](https://training.nextflow.io/latest/side_quests/workflows_of_workflows) (a.k.a. WoW) side quest.

### What's next?

Learn how to graft a basic composable workflow onto the nf-core scaffold.

---

## 3. Fit the updated workflow logic into the placeholder workflow

This is the current content of the `HELLO` workflow in `core-hello/workflows/hello.nf`.
Overall this code does very little aside from some housekeeping that has to do with capturing the version of any software tools that get run in the pipeline.

We need to add the relevant code from the version of the original workflow that we made composable.

```groovy title="core-hello/workflows/hello.nf" linenums="1"
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

We're going to tackle this in the following stages:

1. Copy over the modules and set up module imports
2. Leave the `take` declaration as is
3. Add the workflow logic to the `main` block
4. Update the `emit` block

!!! note

    We're going to ignore the version capture for this first pass and will look at how to wire that up in a later section.

### 3.1. Copy the modules and set up module imports

In the original workflow, the four processes are stored in modules, so we need to copy those over to this new project (into a new `local` directory) and add import statements to the workflow file.

First let's copy the module files over:

```bash
mkdir -p core-hello/modules/local/
cp original-hello/modules/* core-hello/modules/local/.
```

You should now see the directory of modules listed under `core-hello/`.

```bash
tree core-hello/modules
```

```console title="Output"
core-hello/modules
└── local
    ├── collectGreetings.nf
    ├── convertToUpper.nf
    ├── cowpy.nf
    └── sayHello.nf
```

Now let's set up the module import statements.

These were the import statements in the `original-hello/hello.nf` workflow:

```groovy title="original-hello/hello.nf" linenums="9"
// Include modules
include { sayHello } from './modules/sayHello.nf'
include { convertToUpper } from './modules/convertToUpper.nf'
include { collectGreetings } from './modules/collectGreetings.nf'
include { cowpy } from './modules/cowpy.nf'
```

Open the `core-hello/workflows/hello.nf` file and transpose those import statements into it as shown below.

=== "After"

    ```groovy title="core-hello/workflows/hello.nf" linenums="1" hl_lines="8-11"
    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */
    include { paramsSummaryMap       } from 'plugin/nf-schema'
    include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
    include { sayHello               } from '../modules/local/sayHello.nf'
    include { convertToUpper         } from '../modules/local/convertToUpper.nf'
    include { collectGreetings       } from '../modules/local/collectGreetings.nf'
    include { cowpy                  } from '../modules/local/cowpy.nf'
    ```

=== "Before"

    ```groovy title="core-hello/workflows/hello.nf" linenums="1"
    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */
    include { paramsSummaryMap       } from 'plugin/nf-schema'
    include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
    ```

Two more interesting observations here:

- We've adapted the formatting of the import statements to follow the nf-core style convention.
- We've updated the relative paths to the modules to reflect that they're now stored at a different level of nesting.

### 3.2. Leave the `take` declaration as is

The nf-core project has a lot of prebuilt functionality around the concept of the samplesheet, which is typically a CSV file containing columnar data.
Since that is essentially what our `greetings.csv` file is, we'll keep the current `take` declaration as is, and simply update the name of the input channel in the next step.

```groovy title="core-hello/workflows/hello.nf" linenums="21"
    take:
    ch_samplesheet // channel: samplesheet read in from --input
```

The input handling will be done upstream of this workflow (not in this code file).

### 3.3. Add the workflow logic to the `main` block

Now that our modules are available to the workflow, we can plug the workflow logic into the `main` block.

As a reminder, this is the relevant code in the original workflow, which didn't change much when we made it composable (we just added the `main:` line):

```groovy title="original-hello/hello.nf" linenums="21"
    main:

    // emit a greeting
    sayHello(greeting_ch)

    // convert the greeting to uppercase
    convertToUpper(sayHello.out)

    // collect all the greetings into one file
    collectGreetings(convertToUpper.out.collect(), params.batch)

    // emit a message about the size of the batch
    collectGreetings.out.count.view { "There were $it greetings in this batch" }

    // generate ASCII art of the greetings with cowpy
    cowpy(collectGreetings.out.outfile, params.character)
```

We need to copy this code into the new version of the workflow (minus the `main:` keyword which is already there).

There is already some code in there that has to do with capturing the versions of the tools that get run by the workflow. We're going to leave that alone for now (we'll deal with the tool versions later) and simply insert our code right after the `main:` line.

=== "After"

    ```groovy title="core-hello/workflows/hello.nf" linenums="23"

        main:

        // emit a greeting
        sayHello(greeting_ch)

        // convert the greeting to uppercase
        convertToUpper(sayHello.out)

        // collect all the greetings into one file
        collectGreetings(convertToUpper.out.collect(), params.batch)

        // emit a message about the size of the batch
        collectGreetings.out.count.view { "There were $it greetings in this batch" }

        // generate ASCII art of the greetings with cowpy
        cowpy(collectGreetings.out.outfile, params.character)

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

    ```

=== "Before"

    ```groovy title="core-hello/workflows/hello.nf" linenums="23"
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

    ```

This looks great, but we still need to update the name of the channel we're passing to the `sayHello()` process from `greeting_ch` to `ch_samplesheet` (see highlighted lines), to match what is written under the `take:` keyword.

=== "After"

    ```groovy title="core-hello/workflows/hello.nf" linenums="26"
        // emit a greeting (updated to use the nf-core convention for samplesheets)
        sayHello(ch_samplesheet)
    ```

=== "Before"

    ```groovy title="core-hello/workflows/hello.nf" linenums="26"
        // emit a greeting
        sayHello(greeting_ch)
    ```

Now the workflow logic is correctly wired up.

### 3.4. Update the `emit` block

Finally, we need to update the `emit` block to include the declaration of the workflow's final outputs.

=== "After"

    ```groovy title="core-hello/workflows/hello.nf" linenums="55" hl_lines="2"
        emit:
        cowpy_hellos   = cowpy.out
        versions       = ch_versions                 // channel: [ path(versions.yml) ]
    ```

=== "Before"

    ```groovy title="core-hello/workflows/hello.nf" linenums="55"
        emit:
        versions       = ch_versions                 // channel: [ path(versions.yml) ]
    ```

This concludes the modifications we need to make to the HELLO workflow itself.

### Takeaway

You know how to fit the core pieces of a composable workflow into an nf-core placeholder workflow.

### What's next?

Learn how to adapt how the inputs are handle in the nf-core pipeline scaffold.

---

## 4. Adapt the input handling

Now that the HELLO workflow is ready to go, we need to adapt how the inputs are handled to make sure our `greetings.csv` will be handled appropriately.

### 4.1. Identify where inputs are handled

The first step is to figure out where the input handling is done.

You may recall that when we rewrote the Hello Nextflow workflow to be composable, we moved the input parameter declaration up one level, in the `main.nf` entrypoint workflow.
So let's have a look at the top level `main.nf` entrypoint workflow that was created as part of the pipeline scaffold:

```groovy title="core-hello/main.nf" linenums="1"
#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    core/hello
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/core/hello
----------------------------------------------------------------------------------------
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { HELLO  } from './workflows/hello'
include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_nfcore_hello_pipeline'
include { PIPELINE_COMPLETION     } from './subworkflows/local/utils_nfcore_hello_pipeline'
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOWS FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// WORKFLOW: Run main analysis pipeline depending on type of input
//
workflow CORE_HELLO {

    take:
    samplesheet // channel: samplesheet read in from --input

    main:

    //
    // WORKFLOW: Run pipeline
    //
    HELLO (
        samplesheet
    )
}
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {

    main:
    //
    // SUBWORKFLOW: Run initialisation tasks
    //
    PIPELINE_INITIALISATION (
        params.version,
        params.validate_params,
        params.monochrome_logs,
        args,
        params.outdir,
        params.input
    )

    //
    // WORKFLOW: Run main workflow
    //
    CORE_HELLO (
        PIPELINE_INITIALISATION.out.samplesheet
    )
    //
    // SUBWORKFLOW: Run completion tasks
    //
    PIPELINE_COMPLETION (
        params.outdir,
        params.monochrome_logs,
    )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
```

The nf-core project makes heavy use of nested subworkflows, so this bit can be a little confusing on first approach.

What matters here is that there are two workflows defined:

- `CORE_HELLO` is a thin wrapper for running the HELLO workflow we just finished adapting in `core-hello/workflows/hello.nf`.
- An unnamed workflow that calls `CORE_HELLO` as well as two other subworkflows, `PIPELINE_INITIALISATION` and `PIPELINE_COMPLETION`.

Importantly, we cannot find any code constructing an input channel at this level, only references to a samplesheet provided via the `--input` parameter.

A bit of poking around reveals that the input handling is done by the `PIPELINE_INITIALISATION` subworkflow, appropriately enough.

If we open up `core-hello/subworkflows/local/utils_nfcore_hello_pipeline/main.nf` and scroll down, we come to this chunk of code:

```groovy title="core-hello/subworkflows/local/utils_nfcore_hello_pipeline/main.nf" linenums="64"
    //
    // Create channel from input file provided through params.input
    //

    Channel
        .fromList(samplesheetToList(params.input, "${projectDir}/assets/schema_input.json"))
        .map {
            meta, fastq_1, fastq_2 ->
                if (!fastq_2) {
                    return [ meta.id, meta + [ single_end:true ], [ fastq_1 ] ]
                } else {
                    return [ meta.id, meta + [ single_end:false ], [ fastq_1, fastq_2 ] ]
                }
        }
        .groupTuple()
        .map { samplesheet ->
            validateInputSamplesheet(samplesheet)
        }
        .map {
            meta, fastqs ->
                return [ meta, fastqs.flatten() ]
        }
        .set { ch_samplesheet }

    emit:
    samplesheet = ch_samplesheet
    versions    = ch_versions
```

This is the channel factory that parses the samplesheet and passes it on in a form that is ready to be consumed by the HELLO workflow.
It is quite complex because it does a lot of parsing and validation work.

!!! note

    The syntax above is a little different from what we've used previously, but basically this:

    ```groovy
    Channel.<...>.set { ch_samplesheet }
    ```

    is equivalent to this:

    ```groovy
    ch_samplesheet = Channel.<...>
    ```

### 4.2. Replace the templated input channel code

The good news is that our pipeline's needs are much simpler, so we can replace all of that by the channel construction code we developed in the original Hello Nextflow workflow.

As a reminder, this is what the channel construction looked like (as seen in the solutions directory):

```groovy title="solutions/composable-hello/main.nf" linenums="10" hl_lines="4"
    // create a channel for inputs from a CSV file
    greeting_ch = Channel.fromPath(params.greeting)
                        .splitCsv()
                        .map { line -> line[0] }
```

So we just need to plug that into the initialisation workflow, with minor changes: we update the channel name from `greeting_ch` to `ch_samplesheet`, and the parameter name from `params.greeting` to `params.input` (see highlighted line).

=== "After"

    ```groovy title="core-hello/subworkflows/local/utils_nfcore_hello_pipeline/main.nf" linenums="64" hl_lines="4"
        //
        // Create channel from input file provided through params.input
        //
        ch_samplesheet = Channel.fromPath(params.input)
                            .splitCsv()
                            .map { line -> line[0] }

        emit:
        samplesheet = ch_samplesheet
        versions    = ch_versions
    ```

=== "Before"

    ```groovy title="core-hello/subworkflows/local/utils_nfcore_hello_pipeline/main.nf" linenums="64"
        //
        // Create channel from input file provided through params.input
        //

        Channel
            .fromList(samplesheetToList(params.input, "${projectDir}/assets/schema_input.json"))
            .map {
                meta, fastq_1, fastq_2 ->
                    if (!fastq_2) {
                        return [ meta.id, meta + [ single_end:true ], [ fastq_1 ] ]
                    } else {
                        return [ meta.id, meta + [ single_end:false ], [ fastq_1, fastq_2 ] ]
                    }
            }
            .groupTuple()
            .map { samplesheet ->
                validateInputSamplesheet(samplesheet)
            }
            .map {
                meta, fastqs ->
                    return [ meta, fastqs.flatten() ]
            }
            .set { ch_samplesheet }

        emit:
        samplesheet = ch_samplesheet
        versions    = ch_versions
    ```

That completes the changes we need to make the input processing work.

In its current form, this won't let us take advantage of nf-core's built-in capabilities for schema validation, but we can add that in later.
For now, we're focused on keeping it as simple as possible to get to something we can run successfully on test data.

### 4.3. Update the test profile

Speaking of test data and parameters, let's update the test profile for this pipeline to use the `greetings.csv` mini-samplesheet instead of the example samplesheet provided in the template.

Under `core-hello/config`, we find two templated test profiles: `test.config` and `test_full.config`, which are meant to test a small data sample and a full-size one.
Given the purpose of our pipeline, there's not really a point to setting up a full-size test profile, so feel free to ignore or delete `test_full.config`.
We're going to focus on setting up `test.config` to run on our `greetings.csv` file with a few default parameters.

First we need to copy the `greetings.csv` file to an appropriate place in our pipeline project.
Typically small test files are stored in the `assets` directory, so let's copy the file over from our working directory.

```bash
cp greetings.csv core-hello/assets/.
```

Now we can update the `test.config` file as follows:

=== "After"

    ```groovy title="core-hello/conf/test.config" linenums="21" hl_lines="5-10"
        params {
            config_profile_name        = 'Test profile'
            config_profile_description = 'Minimal test dataset to check pipeline function'

            // Input data
            input  = 'core-hello/assets/greetings.csv'

            // Other parameters
            batch     = 'test'
            character = 'tux'
        }
    ```

=== "Before"

    ```groovy title="core-hello/config/test.config" linenums="21"
        params {
            config_profile_name        = 'Test profile'
            config_profile_description = 'Minimal test dataset to check pipeline function'

            // Input data
            // TODO nf-core: Specify the paths to your test data on nf-core/test-datasets
            // TODO nf-core: Give any required params for the test so that command line flags are not needed
            input  = params.pipelines_testdata_base_path + 'viralrecon/samplesheet/samplesheet_test_illumina_amplicon.csv'
        }
    ```

And while we're at it, let's lower the default resource limitations:

=== "After"

    ```groovy title="core-hello/config/test.config" linenums="13"
    process {
        resourceLimits = [
            cpus: 2,
            memory: '4.GB',
            time: '1.h'
        ]
    }
    ```

=== "Before"

    ```groovy title="core-hello/config/test.config" linenums="13"
    process {
        resourceLimits = [
            cpus: 4,
            memory: '15.GB',
            time: '1.h'
        ]
    }
    ```

This completes the code modifications we need to do.

### 4.4. Run the pipeline with the test profile

That was a lot, but we can finally try running the pipeline!
Note that we have to add `--validate_params false` to the command line because we didn't set up the validation yet (that will come later).

```bash
nextflow run core-hello --outdir core-hello-results -profile test,docker --validate_params false
```

If you've done all of the modifications correctly, it should run to completion.

```console title="Output"
 N E X T F L O W   ~  version 24.10.4

Launching `core-hello/main.nf` [agitated_noyce] DSL2 - revision: c31b966b36

Input/output options
  input                     : core-hello/assets/greetings.csv
  outdir                    : core-hello-results

Institutional config options
  config_profile_name       : Test profile
  config_profile_description: Minimal test dataset to check pipeline function

Generic options
  validate_params           : false
  trace_report_suffix       : 2025-05-14_11-10-22

Core Nextflow options
  runName                   : agitated_noyce
  containerEngine           : docker
  launchDir                 : /workspaces/training/hello-nf-core
  workDir                   : /workspaces/training/hello-nf-core/work
  projectDir                : /workspaces/training/hello-nf-core/core-hello
  userName                  : root
  profile                   : test,docker
  configFiles               :

!! Only displaying parameters that differ from the pipeline defaults !!
------------------------------------------------------
executor >  local (8)
[d6/b59dca] CORE_HELLO:HELLO:sayHello (1)       | 3 of 3 ✔
[0b/42f9a1] CORE_HELLO:HELLO:convertToUpper (2) | 3 of 3 ✔
[73/bec621] CORE_HELLO:HELLO:collectGreetings   | 1 of 1 ✔
[3f/e0a67a] CORE_HELLO:HELLO:cowpy              | 1 of 1 ✔
-[core/hello] Pipeline completed successfully-
```

As you can see, this produced the typical nf-core summary at the start thanks to the initialisation subworkflow, and the lines for each module now show the full PIPELINE:WORKFLOW:module names.

### 4.5. Find the pipeline outputs

The question now is: where are the outputs of the pipeline?
And the answer is quite interesting: there are now two different places to look for the results.

We didn't change anything to the modules themselves, so the outputs handled by module-level `publishDir` directives are still going to a `results` directory as specified in the original pipeline.

```bash
tree results
```

```console title="Output"
results
├── Bonjour-output.txt
├── COLLECTED-test-batch-output.txt
├── COLLECTED-test-output.txt
├── cowpy-COLLECTED-test-batch-output.txt
├── cowpy-COLLECTED-test-output.txt
├── Hello-output.txt
├── Holà-output.txt
├── UPPER-Bonjour-output.txt
├── UPPER-Hello-output.txt
└── UPPER-Holà-output.txt
```

Anything that is hooked up to the nf-core template code gets put into a directory generated automatically, called `core-hello-results/`.
This includes the various reports produced by the nf-core utility subworkflows, which you can find under `core-hello-results/pipeline_info`.

```bash
tree core-hello-results
```

```console title="Output"
core-hello-results
└── pipeline_info
    ├── execution_report_2025-06-03_18-22-28.html
    ├── execution_report_2025-06-03_20-11-39.html
    ├── execution_timeline_2025-06-03_18-22-28.html
    ├── execution_timeline_2025-06-03_20-11-39.html
    ├── execution_trace_2025-06-03_18-22-28.txt
    ├── execution_trace_2025-06-03_20-10-11.txt
    ├── execution_trace_2025-06-03_20-11-39.txt
    ├── hello_software_versions.yml
    ├── params_2025-06-03_18-22-32.json
    ├── params_2025-06-03_20-10-15.json
    ├── params_2025-06-03_20-11-43.json
    ├── pipeline_dag_2025-06-03_18-22-28.html
    └── pipeline_dag_2025-06-03_20-11-39.html
```

In our case, we didn't explicitly mark anything else as an output, so there's nothing else there.

And there it is! It may seem like a lot of work to accomplish the same result as the original pipeline, but you do get all those lovely reports generated automatically, and you now have a solid foundation for taking advantage of additional features of nf-core, including input validation and some neat metadata handling capabilities that we'll cover in a later section.

---

### Takeaway

You know how to convert a regular Nextflow pipeline into an nf-core style pipeline using the nf-core template. As part of that, you learned how to make a workflow composable, and identify the most common elements of the nf-core template that need to be adapted when developing a custom nf-core style pipeline.

### What's next?

Take a big break, that was hard work! Your brain deserves to chill out and you could probably use some hydration and a bit of stretching. When you're ready, move on to the next section to learn how to add an nf-core module to an existing nf-core style pipeline. (COMING SOON)
