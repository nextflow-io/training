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

```groovy title="conf/test.config" linenums="1" hl_lines="15 17 19 35"
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
The first step for that is to make our original workflow composable.

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

Launching `original-hello/hello.nf` [goofy_babbage] DSL2 - revision: e9e72441e9

executor >  local (8)
[a4/081cec] sayHello (1)       | 3 of 3 ✔
[e7/7e9058] convertToUpper (3) | 3 of 3 ✔
[0c/17263b] collectGreetings   | 1 of 1 ✔
[94/542280] cowpy              | 1 of 1 ✔
There were 3 greetings in this batch
```

For reference, here is the complete workflow code (not counting the processes, which are in modules):

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

### 2.2. Replace channel construction with `take` statement

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

    ```groovy title="original-hello/hello.nf" linenums="21" hl_lines="22"
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
  final_result = cowpy.out
```

This is a net new addition to the code compared to the original workflow.

### 2.5. Recap of the completed changes

If you've done all the changes as described, your workflow should now look like this:

```groovy title="original-hello/hello.nf" linenums="1" hl_lines="15 18-20 22 38-39"
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
    final_result = cowpy.out
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
  HELLO.out.view { "Outputs: $it" }
}
```

You can see that the syntax for calling the imported workflow is essentially the same as the syntax for calling modules.

You should also note that everything that has to do with pulling the inputs into the workflow (input parameter and channel construction) is now declared in this parent workflow.

!!!note
You can name the entrypoint workflow file whatever you want, it does not have to be named `main.nf`.
The advantage of naming it `main.nf` is that if you don't specify a workflow file, Nextflow will automatically look for a file named `main.nf` in the specified directory.

### 2.7. Test that the workflow runs

We finally have all the pieces we need to verify that the composable workflow works.

```bash
nextflow run original-hello
```

!!! note
Here you see the advantage of using the `main.nf` naming convention, which allows us to omit including the name of the workflow file in the command.
If we had named it `something_else.nf`, we would have had to do `nextflow run original-hello/something_else.nf`.

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
Outputs: /workspaces/training/hello-nf-core/work/e3/693b7e48dc119d0c54543e0634c2e7/cowpy-COLLECTED-test-batch-output.txt
```

This means we've successfully upgraded our HELLO workflow to be composable.

### Takeaway

You know how to make a workflow composable by giving it a name and adding `take`, `main` and `emit` statements, and how to call it from an entrypoint workflow.

!!!note
If you're interested in digging deeper into options for composing workflows of workflows, check out the [Workflow of Workflows](https://training.nextflow.io/latest/side_quests/workflows_of_workflows) (a.k.a. WoW) side quest.

### What's next?

Learn how to graft a basic composable workflow onto the nf-core scaffold.

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
