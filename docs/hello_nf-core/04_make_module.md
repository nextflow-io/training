# Part 4: Make an nf-core module

In this fourth part of the Hello nf-core training course, we show you how to create an nf-core module by applying the key conventions that make modules portable and maintainable.

The nf-core project provides a command (`nf-core modules create`) that generates properly structured module templates automatically, similar to what we used for the workflow in Part 2.
However, for teaching purposes, we're going to start by doing it manually: transforming the local `cowpy` module in your `core-hello` pipeline into an nf-core-style module step-by-step.
After that, we'll show you how to use the template-based module creation to work more efficiently in the future.

??? info "How to begin from this section"

    This section assumes you have completed [Part 3: Use an nf-core module](./03_use_module.md) and have integrated the `CAT_CAT` module into your pipeline.

    If you did not complete Part 3 or want to start fresh for this part, you can use the `core-hello-part3` solution as your starting point.
    Run these commands from inside the `hello-nf-core/` directory:

    ```bash
    cp -r solutions/core-hello-part3 core-hello
    cd core-hello
    ```

    This gives you a pipeline with the `CAT_CAT` module already integrated.
    You can test that it runs successfully by running the following command:

    ```bash
    nextflow run . --outdir core-hello-results -profile test,docker --validate_params false
    ```

---

## 1. Transform `cowpy` into an nf-core module

In this section, we'll apply nf-core conventions to the local `cowpy` module in your `core-hello` pipeline, transforming it into a module that follows nf-core community standards.

This is the current code for the `cowpy` process module:

```groovy title="core-hello/modules/local/cowpy.nf" linenums="1"
#!/usr/bin/env nextflow

// Generate ASCII art with cowpy (https://github.com/jeffbuttars/cowpy)
process cowpy {

    publishDir 'results', mode: 'copy'

    container 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'
    conda 'conda-forge::cowpy==1.1.5'

    input:
    path input_file
    val character

    output:
    path "cowpy-${input_file}"

    script:
    """
    cat $input_file | cowpy -c "$character" > cowpy-${input_file}
    """
}
```

We'll apply the following nf-core conventions incrementally:

1. **Uppercase the process name to `COWPY`** to follow convention.
2. **Update `COWPY` to use metadata tuples** to propagate sample metadata through the workflow.
3. **Centralize tool argument configuration with `ext.args`** to increase module versatility while keeping the interface minimal.
4. **Standardize output naming with `ext.prefix`** to promote consistency.
5. **Centralize the publishing configuration** to promote consistency.

After each step, we'll run the pipeline to test that everything works as expected.

!!! warning "Working directory"

    Make sure you're in the `core-hello` directory (your pipeline root) for all file edits and command executions in this section.

    ```bash
    cd core-hello
    ```

### 1.1. Uppercase the process name

This is purely a stylistic convention (there is no technical justification) but since it is the norm for nf-core modules, let's comply.

We need to make three sets of changes:

1. Update the process name in the module
2. Update the module import statement in the workflow header
3. Update the process call and emit declaration in the workflow body

Let's get started!

#### 1.1.1. Update the process name in the module

Open the `cowpy.nf` module file (under `core-hello/modules/local/`) and modify the process name to uppercase:

=== "After"

    ```groovy title="core-hello/modules/local/cowpy.nf" linenums="3" hl_lines="2"
    // Generate ASCII art with cowpy (https://github.com/jeffbuttars/cowpy)
    process COWPY {
    ```

=== "Before"

    ```groovy title="core-hello/modules/local/cowpy.nf" linenums="3" hl_lines="2"
    // Generate ASCII art with cowpy (https://github.com/jeffbuttars/cowpy)
    process cowpy {
    ```

In this case the uppercasing is completely straightforward.

If the process name was composed of several words, for example if we had a process called MyCowpyTool originally in camel case, the nf-core convention would be to use underscores to separate them, yielding MY_COWPY_TOOL.

#### 1.1.2. Update the module import statement

Process names are case-sensitive, so now that we've changed the process name, we need to update the module import statement accordingly in the workflow header of `hello.nf`:

=== "After"

    ```groovy title="core-hello/workflows/hello.nf" linenums="1" hl_lines="10"
    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */
    include { paramsSummaryMap       } from 'plugin/nf-schema'
    include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
    include { sayHello               } from '../modules/local/sayHello.nf'
    include { convertToUpper         } from '../modules/local/convertToUpper.nf'
    include { COWPY                  } from '../modules/local/cowpy/main.nf'
    include { CAT_CAT                } from '../modules/nf-core/cat/cat/main'
    ```

=== "Before"

    ```groovy title="core-hello/workflows/hello.nf" linenums="1" hl_lines="10"
    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */
    include { paramsSummaryMap       } from 'plugin/nf-schema'
    include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
    include { sayHello               } from '../modules/local/sayHello.nf'
    include { convertToUpper         } from '../modules/local/convertToUpper.nf'
    include { cowpy                  } from '../modules/local/cowpy/main.nf'
    include { CAT_CAT                } from '../modules/nf-core/cat/cat/main'
    ```

We could use an alias in the import statement to avoid having to update calls to the process, but that would somewhat defeat the point of adopting the uppercasing convention.

#### 1.1.3. Update the process call and emit declaration

So now let's update the two references to the process in the workflow block of `hello.nf`:

=== "After"

    ```groovy title="core-hello/workflows/hello.nf" linenums="43" hl_lines="2 17"
    // generate ASCII art of the greetings with cowpy
    COWPY(CAT_CAT.out.file_out)

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
    cowpy_hellos   = COWPY.out.cowpy_output
    versions       = ch_versions
    ```

=== "Before"

    ```groovy title="core-hello/workflows/hello.nf" linenums="43" hl_lines="2 17"
    // generate ASCII art of the greetings with cowpy
    cowpy(CAT_CAT.out.file_out)

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
    cowpy_hellos   = cowpy.out.cowpy_output
    versions       = ch_versions
    ```

Be sure to make **both** changes, otherwise you will get an error when you run this.

#### 1.1.4. Run the pipeline to test it

Let's run the workflow to test that everything is working correctly after these changes.

```bash
nextflow run . --outdir core-hello-results -profile test,docker --validate_params false
```

??? success "Command output"

    ```console
     N E X T F L O W   ~  version 25.04.3

    Launching `./main.nf` [elegant_plateau] DSL2 - revision: b9e9b3b8de

    Input/output options
      input                     : /workspaces/training/hello-nf-core/core-hello/assets/greetings.csv
      outdir                    : core-hello-results

    Institutional config options
      config_profile_name       : Test profile
      config_profile_description: Minimal test dataset to check pipeline function

    Generic options
      validate_params           : false
      trace_report_suffix       : 2026-01-06_04-51-29

    Core Nextflow options
      runName                   : elegant_plateau
      containerEngine           : docker
      launchDir                 : /workspaces/training/hello-nf-core/core-hello
      workDir                   : /workspaces/training/hello-nf-core/core-hello/work
      projectDir                : /workspaces/training/hello-nf-core/core-hello
      userName                  : root
      profile                   : test,docker
      configFiles               : /workspaces/training/hello-nf-core/core-hello/nextflow.config

    !! Only displaying parameters that differ from the pipeline defaults !!
    ------------------------------------------------------
    executor >  local (8)
    [7b/66ceb5] CORE_HELLO:HELLO:sayHello (3)       | 3 of 3 ✔
    [8e/1bafb9] CORE_HELLO:HELLO:convertToUpper (3) | 3 of 3 ✔
    [bb/203575] CORE_HELLO:HELLO:CAT_CAT (test)     | 1 of 1 ✔
    [39/715489] CORE_HELLO:HELLO:COWPY              | 1 of 1 ✔
    -[core/hello] Pipeline completed successfully-
    ```

Alright, this works! Now let's move on to making more substantial changes.

### 1.2. Update `COWPY` to use metadata tuples

In the current version of the `core-hello` pipeline, we're extracting the file from `CAT_CAT`'s output tuple to pass to `COWPY`, as shown in the top half of the diagram below.

<figure class="excalidraw">
    --8<-- "docs/hello_nf-core/img/cowpy-inputs.svg"
</figure>

It would be better to have `COWPY` accept metadata tuples directly, allowing metadata to flow on through the workflow, as shown in the bottom half of the diagram.

To the end, we'll need to make the following changes:

1. Update the input and output definitions
2. Update the process call in the workflow
3. Update the emit block in the workflow

Once we've done all that, we'll run the pipeline to test that everything still works as before.

#### 1.2.1. Update the input and output definitions

Return to the `cowpy.nf` module file and modify it to accept metadata tuples as shown below.

=== "After"

    ```groovy title="core-hello/modules/local/cowpy.nf" linenums="11" hl_lines="2 6"
        input:
            tuple val(meta), path(input_file)
            val character

        output:
            tuple val(meta), path("cowpy-${input_file}"), emit: cowpy_output
    ```

=== "Before"

    ```groovy title="core-hello/modules/local/cowpy.nf" linenums="11" hl_lines="2 6"
        input:
            path input_file
            val character

        output:
            path "cowpy-${input_file}"
    ```

As you can see, we changed both the **main input** and the **output** to a tuple that follows the `tuple val(meta), path(input_file)` pattern introduced in Part 3 of this training.
For the output, we also took this opportunity to add `emit: cowpy_output` in order to give the output channel a descriptive name.

Now that we've changed what the process expects, we need to update what we provide to it in the process call.

#### 1.2.2. Update the process call in the workflow

The good news is that this change will simplify the process call.
Now that the output of `CAT_CAT` and the input of `COWPY` are the same 'shape', i.e. they both consist of a `tuple val(meta), path(input_file)` structure, we can simply connect them directly instead of having to extract the file explicitly from the output of the `CAT_CAT` process.

Open the `hello.nf` workflow file (under `core-hello/workflows/`) and update the call to `COWPY` as shown below.

=== "After"

    ```groovy title="core-hello/workflows/hello.nf" linenums="43" hl_lines="2"
        // generate ASCII art of the greetings with cowpy
        COWPY(CAT_CAT.out.file_out, params.character)
    ```

=== "Before"

    ```groovy title="core-hello/workflows/hello.nf" linenums="43" hl_lines="1-2 5"
        // extract the file from the tuple since cowpy doesn't use metadata yet
        ch_for_cowpy = CAT_CAT.out.file_out.map{ meta, file -> file }

        // generate ASCII art of the greetings with cowpy
        COWPY(ch_for_cowpy, params.character)
    ```

We now call `COWPY` on `CAT_CAT.out.file_out` directly.

As a result, we no longer need to construct the `ch_for_cowpy` channel, so that line (and its comment line) can be deleted entirely.

#### 1.2.3. Update the emit block in the workflow

Since `COWPY` now emits a named output, `cowpy_output`, we can update the `hello.nf` workflow's `emit:` block to use that.

=== "After"

    ```groovy title="core-hello/workflows/hello.nf" linenums="60" hl_lines="2"
        emit:
        cowpy_hellos   = COWPY.out.cowpy_output
        versions       = ch_versions
    ```

=== "Before"

    ```groovy title="core-hello/workflows/hello.nf" linenums="60" hl_lines="2"
        emit:
        cowpy_hellos   = COWPY.out
        versions       = ch_versions
    ```

This is technically not required, but it's good practice to refer to named outputs whenever possible.

#### 1.2.4. Run the pipeline to test it

Let's run the workflow to test that everything is working correctly after these changes.

```bash
nextflow run . --outdir core-hello-results -profile test,docker --validate_params false
```

??? success "Command output"

    ```console
     N E X T F L O W   ~  version 25.04.3

    Launching `./main.nf` [modest_saha] DSL2 - revision: b9e9b3b8de

    Downloading plugin nf-schema@2.5.1
    Input/output options
      input                     : /workspaces/training/hello-nf-core/core-hello/assets/greetings.csv
      outdir                    : core-hello-results

    Institutional config options
      config_profile_name       : Test profile
      config_profile_description: Minimal test dataset to check pipeline function

    Generic options
      validate_params           : false
      trace_report_suffix       : 2025-12-27_06-16-55

    Core Nextflow options
      runName                   : modest_saha
      containerEngine           : docker
      launchDir                 : /workspaces/training/hello-nf-core/core-hello
      workDir                   : /workspaces/training/hello-nf-core/core-hello/work
      projectDir                : /workspaces/training/hello-nf-core/core-hello
      userName                  : root
      profile                   : test,docker
      configFiles               : /workspaces/training/hello-nf-core/core-hello/nextflow.config

    !! Only displaying parameters that differ from the pipeline defaults !!
    ------------------------------------------------------
    executor >  local (8)
    [a8/447993] CORE_HELLO:HELLO:sayHello (3)       [100%] 3 of 3 ✔
    [00/1fc59c] CORE_HELLO:HELLO:convertToUpper (3) [100%] 3 of 3 ✔
    [57/ac800d] CORE_HELLO:HELLO:CAT_CAT (test)     [100%] 1 of 1 ✔
    [b7/092f2b] CORE_HELLO:HELLO:COWPY              [100%] 1 of 1 ✔
    -[core/hello] Pipeline completed successfully-
    ```

The pipeline should run successfully, with metadata now flowing from `CAT_CAT` through `COWPY`.

That completes what we needed to do to make `COWPY` handle metadata tuples.
Now, let's look at what else we can do to take advantage of nf-core module patterns.

### 1.3. Centralize tool argument configuration with `ext.args`

In its current state, the `COWPY` process expects to receive a value for the `character` parameter.
As a result, we have to provide a value every time we call the process, even if we'd be happy with the defaults set by the tool.
For `COWPY` this is admittedly not a big problem, but for tools with many optional parameters, it can get quite cumbersome.

The nf-core project recommends using a Nextflow feature called [`ext.args`](https://www.nextflow.io/docs/latest/reference/process.html#ext) to manage tool arguments more conveniently via configuration files.

Instead of declaring process inputs for every tool option, you write the module to reference `ext.args` in the construction of its command line.
Then it's just a matter of setting up the `ext.args` variable to hold the arguments and values you want to use in the `modules.config` file, which consolidates configuration details for all modules.
Nextflow will add those arguments with their values into the tool command line at runtime.

Let's apply this approach to the `COWPY` module.
We're going to need to make the following changes:

1. Update the `COWPY` module
2. Configure `ext.args` in the `modules.config` file
3. Update the `hello.nf` workflow

Once we've done all that, we'll run the pipeline to test that everything still works as before.

#### 1.3.1. Update the `COWPY` module

Let's do it.
Open the `cowpy.nf` module file (under `core-hello/modules/local/`) and modify it to reference `ext.args` as shown below.

=== "After"

    ```groovy title="modules/local/cowpy.nf" linenums="1" hl_lines="18 20"
    #!/usr/bin/env nextflow

    // Generate ASCII art with cowpy (https://github.com/jeffbuttars/cowpy)
    process COWPY {

        publishDir 'results', mode: 'copy'

        container 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'
        conda 'conda-forge::cowpy==1.1.5'

        input:
            tuple val(meta), path(input_file)

        output:
            tuple val(meta), path("cowpy-${input_file}"), emit: cowpy_output

        script:
        def args = task.ext.args ?: ''
        """
        cat $input_file | cowpy $args > cowpy-${input_file}
        """
    }
    ```

=== "Before"

    ```groovy title="core-hello/modules/local/cowpy.nf" linenums="1" hl_lines="13 20"
    #!/usr/bin/env nextflow

    // Generate ASCII art with cowpy (https://github.com/jeffbuttars/cowpy)
    process COWPY {

        publishDir 'results', mode: 'copy'

        container 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'
        conda 'conda-forge::cowpy==1.1.5'

        input:
            tuple val(meta), path(input_file)
            val character

        output:
            tuple val(meta), path("cowpy-${input_file}"), emit: cowpy_output

        script:
        """
        cat $input_file | cowpy -c "$character" > cowpy-${input_file}
        """
    }
    ```

You can see we made three changes.

1. **In the `input:` block, we removed the `val character` input.**
   Going forward, we'll supply that argument via the `ext.args` configuration as described further below.

2. **In the `script:` block, we added the line `def args = task.ext.args ?: ''`.**
   That line uses the `?:` operator to determine the value of the `args` variable: the content of `task.ext.args` if it is not empty, or an empty string if it is.
   Note that while we generally refer to `ext.args`, this code must reference `task.ext.args` to pull out the module-level `ext.args` configuration.

3. **In the command line, we replaced `-c "$character"` with `$args`.**
   This is where Nextflow will inject any tool arguments set in `ext.args` in the `modules.config` file.

As a result, the module interface is now simpler: it only expects the essential metadata and file inputs.

!!! note

    The `?:` operator is often called the 'Elvis operator' because it looks like a sideways Elvis Presley face, with the `?` character symbolizing the wave in his hair.

#### 1.3.2. Configure `ext.args` in the `modules.config` file

Now that we've taken the `character` declaration out of the module, we've got to add it to `ext.args` in the `modules.config` configuration file.

Specifically, we're going to add this little chunk of code to the `process {}` block:

```groovy title="Code to add"
withName: 'COWPY' {
    ext.args = { "-c ${params.character}" }
}
```

The `withName:` syntax assigns this configuration to the `COWPY` process only, and `ext.args = { "-c ${params.character}" }` simply composes a string that will include the value of the `character` parameter.
Note the use of curly braces, which tell Nextflow to evaluate the value of the parameter at runtime.

Makes sense? Let's add it in.

Open `conf/modules.config` and add the configuration code inside the `process {}` block as shown below.

=== "After"

    ```groovy title="core-hello/conf/modules.config" linenums="13" hl_lines="8-10"
    process {
        publishDir = [
            path: { "${params.outdir}/${task.process.tokenize(':')[-1].tokenize('_')[0].toLowerCase()}" },
            mode: params.publish_dir_mode,
            saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
        ]

        withName: 'COWPY' {
            ext.args = { "-c ${params.character}" }
        }
    }
    ```

=== "Before"

    ```groovy title="core-hello/conf/modules.config" linenums="13"
    process {
        publishDir = [
            path: { "${params.outdir}/${task.process.tokenize(':')[-1].tokenize('_')[0].toLowerCase()}" },
            mode: params.publish_dir_mode,
            saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
        ]
    }
    ```

Hopefully you can imagine having all the modules in a pipeline have their `ext.args` specified in this file, with the following benefits:

- The **module interface stays simple** - It only accepts the essential metadata and file inputs
- The **pipeline still exposes `params.character`** - End-users can still configure it as before
- The **module is now portable** - It can be reused in other pipelines without expecting a specific parameter name
- The configuration is **centralized** in `modules.config`, keeping workflow logic clean

By using the `modules.config` file as the place where all pipelines centralize per-module configuration, we make our modules more reusable across different pipelines.

#### 1.3.3. Update the `hello.nf` workflow

Since the `COWPY` module no longer requires the `character` parameter as an input, we need to update the workflow call accordingly.

Open the `hello.nf` workflow file (under `core-hello/workflows/`) and update the call to `COWPY` as shown below.

=== "After"

    ```groovy title="core-hello/workflows/hello.nf" linenums="39" hl_lines="2"
        // generate ASCII art of the greetings with cowpy
        COWPY(CAT_CAT.out.file_out)
    ```

=== "Before"

    ```groovy title="core-hello/workflows/hello.nf" linenums="39" hl_lines="2"
        // generate ASCII art of the greetings with cowpy
        COWPY(CAT_CAT.out.file_out, params.character)
    ```

The workflow code is now cleaner: we don't need to pass `params.character` directly to the process.
The module interface is kept minimal, making it more portable, while the pipeline still provides the explicit option through configuration.

#### 1.3.4. Run the pipeline to test it

Let's test that the workflow still works as expected, specifying a different character to verify that the `ext.args` configuration is working.

Run this command using `kosh`, one of the more... enigmatic options:

```bash
nextflow run . --outdir core-hello-results -profile test,docker --validate_params false --character kosh
```

??? success "Command output"

    ```console
     N E X T F L O W   ~  version 25.04.3

    Launching `./main.nf` [exotic_planck] DSL2 - revision: b9e9b3b8de

    Input/output options
      input                     : /workspaces/training/hello-nf-core/core-hello/assets/greetings.csv
      outdir                    : core-hello-results

    Institutional config options
      config_profile_name       : Test profile
      config_profile_description: Minimal test dataset to check pipeline function

    Generic options
      validate_params           : false
      trace_report_suffix       : 2025-12-27_06-23-13

    Core Nextflow options
      runName                   : exotic_planck
      containerEngine           : docker
      launchDir                 : /workspaces/training/hello-nf-core/core-hello
      workDir                   : /workspaces/training/hello-nf-core/core-hello/work
      projectDir                : /workspaces/training/hello-nf-core/core-hello
      userName                  : root
      profile                   : test,docker
      configFiles               : /workspaces/training/hello-nf-core/core-hello/nextflow.config

    !! Only displaying parameters that differ from the pipeline defaults !!
    ------------------------------------------------------
    executor >  local (8)
    [13/9e3c0e] CORE_HELLO:HELLO:sayHello (3)       [100%] 3 of 3 ✔
    [e2/5b0ee5] CORE_HELLO:HELLO:convertToUpper (3) [100%] 3 of 3 ✔
    [b6/4fb569] CORE_HELLO:HELLO:CAT_CAT (test)     [100%] 1 of 1 ✔
    [38/eb29ea] CORE_HELLO:HELLO:COWPY              [100%] 1 of 1 ✔
    -[core/hello] Pipeline completed successfully-
    ```

This should run successfully as previously.

Let's verify that the `ext.args` configuration worked by checking the output.
Find the output in the file browser or use the task hash (the `38/eb29ea` part in the example above) to look at the output file:

```bash
cat work/38/eb29ea*/cowpy-test.txt
```

??? success "Command output"

    ```console
    _________
    / HELLO   \
    | HOLà    |
    \ BONJOUR /
    ---------
        \
        \
          \
      ___       _____     ___
    /   \     /    /|   /   \
    |     |   /    / |  |     |
    |     |  /____/  |  |     |
    |     |  |    |  |  |     |
    |     |  | {} | /   |     |
    |     |  |____|/    |     |
    |     |    |==|     |     |
    |      \___________/      |
    |                         |
    |                         |
    ```

You should see the ASCII art displayed with the `kosh` character, confirming that the `ext.args` configuration worked!

??? info "(Optional) Inspect the command file"

    If you want to see exactly how the configuration was applied, you can inspect the `.command.sh` file:

    ```bash
    cat work/38/eb29ea*/.command.sh
    ```

    You'll see the `cowpy` command with the `-c kosh` argument:

    ```console
    #!/usr/bin/env bash
    ...
    cat test.txt | cowpy -c kosh > cowpy-test.txt
    ```

    This shows that the `.command.sh` file was generated correctly based on the `ext.args` configuration.

Take a moment to think about what we achieved here.
This approach keeps the module interface focused on essential data (files, metadata, and any mandatory per-sample parameters), while options that control the behavior of the tool are handled separately through configuration.

This may seem unnecessary for a simple tool like `cowpy`, but it can make a big difference for data analysis tools that have a lot of optional arguments.

To summarize the benefits of this approach:

- **Clean interface**: The module focuses on essential data inputs (metadata and files)
- **Flexibility**: Users can specify tool arguments via configuration, including sample-specific values
- **Consistency**: All nf-core modules follow this pattern
- **Portability**: Modules can be reused without hardcoded tool options
- **No workflow changes**: Adding or changing tool options doesn't require updating workflow code

!!! note

    The `ext.args` system has powerful additional capabilities not covered here, including switching argument values dynamically based on metadata. See the [nf-core module specifications](https://nf-co.re/docs/guidelines/components/modules) for more details.

### 1.4. Standardize output naming with `ext.prefix`

Now that we've given the `COWPY` process access to the metamap, we can start taking advantage of another useful nf-core pattern: naming output files based on metadata.

Here we're going to use a Nextflow feature called `ext.prefix` that will allow us to standardize output file naming across modules using `meta.id` (the identifier included in the metamap), while still being able to configure modules individually if desired.

This will be similar to what we did with `ext.args`, with a few differences that we'll detail as we go.

Let's apply this approach to the `COWPY` module.
We're going to need to make the following changes:

1. Update the `COWPY` module
2. Configure `ext.prefix` in the `modules.config` file

(No changes needed to the workflow.)

Once we've done that, we'll run the pipeline to test that everything still works as before.

#### 1.4.1. Update the `COWPY` module

Open the `cowpy.nf` module file (under `core-hello/modules/local/`) and modify it to reference `ext.prefix` as shown below.

=== "After"

    ```groovy title="core-hello/modules/local/cowpy.nf" linenums="1" hl_lines="2 6 8"
        output:
            tuple val(meta), path("${prefix}.txt"), emit: cowpy_output

        script:
        def args = task.ext.args ?: ''
        prefix = task.ext.prefix ?: "${meta.id}"
        """
        cat $input_file | cowpy $args > ${prefix}.txt
        """
    }
    ```

=== "Before"

    ```groovy title="core-hello/modules/local/cowpy.nf" linenums="1" hl_lines="2 7"
        output:
            tuple val(meta), path("cowpy-${input_file}"), emit: cowpy_output

        script:
        def args = task.ext.args ?: ''
        """
        cat $input_file | cowpy $args > cowpy-${input_file}
        """
    }
    ```

You can see we made three changes.

1. **In the `script:` block, we added the line `prefix = task.ext.prefix ?: "${meta.id}"`.**
   That line uses the `?:` operator to determine the value of the `prefix` variable: the content of `task.ext.prefix` if it is not empty, or the identifier from the metamap (`meta.id`) if it is.
   Note that while we generally refer to `ext.prefix`, this code must reference `task.ext.prefix` to pull out the module-level `ext.prefix` configuration.

2. **In the command line, we replaced `cowpy-${input_file}` with `${prefix}.txt`.**
   This is where Nextflow will inject the value of `prefix` determined by the line above.

3. **In the `output:` block, we replaced `path("cowpy-${input_file}")` with `path("${prefix}.txt")`.\*\***
   This simply reiterates what the file path will be according to what is written in the command line.

As a result, the output file name is now constructed using a sensible default (the identifier from the metamap) combined with the appropriate file format extension.

#### 1.4.2. Configure `ext.prefix` in the `modules.config` file

In this case the sensible default is not sufficiently expressive for our taste; we want to use a custom naming pattern that includes the tool name, `cowpy-<id>.txt`, like we had before.

We'll do that by configuring `ext.prefix` in `modules.config`, just like we did for the `character` parameter with `ext.args`, except this time the `withName: 'COWPY' {}` block already exists, and we just need to add the following line:

```groovy title="Code to add"
ext.prefix = { "cowpy-${meta.id}" }
```

This will compose the string we want.
Note that once again we use curly braces, this time to tell Nextflow to evaluate the value of `meta.id` at runtime.

Let's add it in.

Open `conf/modules.config` and add the configuration code inside the `process {}` block as shown below.

=== "After"

    ```groovy title="core-hello/conf/modules.config" linenums="21" hl_lines="3"
        withName: 'COWPY' {
            ext.args = { "-c ${params.character}" }
            ext.prefix = { "cowpy-${meta.id}" }
        }
    ```

=== "Before"

    ```groovy title="core-hello/conf/modules.config" linenums="21"
        withName: 'COWPY' {
            ext.args = { "-c ${params.character}" }
        }
    ```

In case you're wondering, the `ext.prefix` closure has access to the correct piece of metadata because the configuration is evaluated in the context of the process execution, where metadata is available.

#### 1.4.3. Run the pipeline to test it

Let's test that the workflow still works as expected.

```bash
nextflow run . --outdir core-hello-results -profile test,docker --validate_params false
```

??? success "Command output"

    ```console
     N E X T F L O W   ~  version 25.04.3

    Launching `./main.nf` [admiring_turing] DSL2 - revision: b9e9b3b8de

    Input/output options
      input                     : /workspaces/training/hello-nf-core/core-hello/assets/greetings.csv
      outdir                    : core-hello-results

    Institutional config options
      config_profile_name       : Test profile
      config_profile_description: Minimal test dataset to check pipeline function

    Generic options
      validate_params           : false
      trace_report_suffix       : 2025-12-27_06-29-02

    Core Nextflow options
      runName                   : admiring_turing
      containerEngine           : docker
      launchDir                 : /workspaces/training/hello-nf-core/core-hello
      workDir                   : /workspaces/training/hello-nf-core/core-hello/work
      projectDir                : /workspaces/training/hello-nf-core/core-hello
      userName                  : root
      profile                   : test,docker
      configFiles               : /workspaces/training/hello-nf-core/core-hello/nextflow.config

    !! Only displaying parameters that differ from the pipeline defaults !!
    ------------------------------------------------------
    executor >  local (8)
    [b2/e08524] CORE_HELLO:HELLO:sayHello (3)       [100%] 3 of 3 ✔
    [13/88939f] CORE_HELLO:HELLO:convertToUpper (3) [100%] 3 of 3 ✔
    [23/4554e1] CORE_HELLO:HELLO:CAT_CAT (test)     [100%] 1 of 1 ✔
    [a3/c6cbe9] CORE_HELLO:HELLO:COWPY              [100%] 1 of 1 ✔
    -[core/hello] Pipeline completed successfully-
    ```

Take a look at the output in the results directory.
You should see the cowpy output file with the same naming as before: `cowpy-test.txt`, based on the default batch name.

??? abstract "Directory contents"

    ```console hl_lines="3"
    results
    ├── Bonjour-output.txt
    ├── cowpy-test.txt
    ├── Hello-output.txt
    ├── Holà-output.txt
    ├── UPPER-Bonjour-output.txt
    ├── UPPER-Hello-output.txt
    └── UPPER-Holà-output.txt
    ```

Feel free to change the `ext.prefix` configuration in `conf/modules.config` to satisfy yourself that you can change the naming pattern without having to make any changes to the module or workflow code.

Alternatively, you can also try running this again with a different `--batch` parameter specified on the command line to satisfy yourself that that part is still customizable on the fly.

This demonstrates how `ext.prefix` allows you to maintain your preferred naming convention while keeping the module interface flexible.

To summarize the benefits of this approach:

- **Standardized naming**: Output files are typically named using sample IDs from metadata
- **Configurable**: Users can override the default naming if needed
- **Consistent**: All nf-core modules follow this pattern
- **Predictable**: Easy to know what output files will be called

Pretty good, right?
Well, there's one more important change we need to make to improve our module to fit the nf-core guidelines.

### 1.5. Centralize the publishing configuration

You may have noticed that we've been publishing outputs to two different directories:

- **`results`** — The original output directory we've been using from the beginning for our local modules, set individually using per-module `publishDir` directives;
- **`core-hello-results`** — The output directory set with `--outdir` on the command line, which has been receiving the nf-core logs and the results published by `CAT_CAT`.

This is messy and suboptimal; it would be better to have one location for everything.
Of course, we could go into each of our local modules and update the `publishDir` directive manually to use the `core-hello-results` directory, but what about next time we decide to change the output directory?

Having individual modules make publishing decisions is clearly not the way to go, especially in a world where the same module might be used in a lot of different pipelines, by people who have different needs or preferences.
We want to be able to control where outputs get published at the level of the workflow configuration.

"Hey," you might say, "`CAT_CAT` is sending its outputs to the `--outdir`. Maybe we should copy its `publishDir` directive?"

Yes, that's a great idea.

Except it doesn't have a `publishDir` directive. (Go ahead, look at the module code.)

That's because nf-core pipelines centralize control at the workflow level by configuring `publishDir` in `conf/modules.config` rather than in individual modules.
Specifically, the nf-core template declares a default `publishDir` directive (with a predefined directory structure) that applies to all modules unless an overriding directive is provide.

Doesn't that sound awesome? Could it be that to take advantage of this default directive, all we need to do is remove the current `publishDir` directive from our local modules?

Let's try that out on `COWPY` to see what happens, then we'll look at the code for the default configuration to understand how it works.

Finally, we'll demonstrate how to override the default behavior if desired.

#### 1.5.1. Remove the `publishDir` directive from `COWPY`

Let's do this.
Open the `cowpy.nf` module file (under `core-hello/modules/local/`) and remove the `publishDir` directive as shown below.

=== "After"

    ```groovy title="core-hello/modules/local/cowpy.nf (excerpt)" linenums="1"
    #!/usr/bin/env nextflow

    // Generate ASCII art with cowpy (https://github.com/jeffbuttars/cowpy)
    process COWPY {

        container 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'
        conda 'conda-forge::cowpy==1.1.5'
    ```

=== "Before"

    ```groovy title="core-hello/modules/local/cowpy.nf (excerpt)" linenums="1" hl_lines="6"
    #!/usr/bin/env nextflow

    // Generate ASCII art with cowpy (https://github.com/jeffbuttars/cowpy)
    process COWPY {

        publishDir 'results', mode: 'copy'

        container 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'
        conda 'conda-forge::cowpy==1.1.5'

    ```

That's it!

#### 1.5.2. Run the pipeline to test it

Let's have a look at what happens if we run the pipeline now.

```bash
nextflow run . --outdir core-hello-results -profile test,docker --validate_params false
```

??? success "Command output"

    ```console
     N E X T F L O W   ~  version 25.04.3

    Launching `./main.nf` [silly_caravaggio] DSL2 - revision: b9e9b3b8de

    Input/output options
      input                     : /workspaces/training/hello-nf-core/core-hello/assets/greetings.csv
      outdir                    : core-hello-results

    Institutional config options
      config_profile_name       : Test profile
      config_profile_description: Minimal test dataset to check pipeline function

    Generic options
      validate_params           : false
      trace_report_suffix       : 2025-12-27_06-35-56

    Core Nextflow options
      runName                   : silly_caravaggio
      containerEngine           : docker
      launchDir                 : /workspaces/training/hello-nf-core/core-hello
      workDir                   : /workspaces/training/hello-nf-core/core-hello/work
      projectDir                : /workspaces/training/hello-nf-core/core-hello
      userName                  : root
      profile                   : test,docker
      configFiles               : /workspaces/training/hello-nf-core/core-hello/nextflow.config

    !! Only displaying parameters that differ from the pipeline defaults !!
    ------------------------------------------------------
    executor >  local (8)
    [db/39978e] CORE_HELLO:HELLO:sayHello (3)       [100%] 3 of 3 ✔
    [b5/bf6a8d] CORE_HELLO:HELLO:convertToUpper (3) [100%] 3 of 3 ✔
    [b7/c61842] CORE_HELLO:HELLO:CAT_CAT (test)     [100%] 1 of 1 ✔
    [46/5839d6] CORE_HELLO:HELLO:COWPY              [100%] 1 of 1 ✔
    -[core/hello] Pipeline completed successfully-
    ```

Have a look at your current working directory.
Now the `core-hello-results` also contains the outputs of the `COWPY` module.

??? abstract "Directory contents"

    ```console hl_lines="4-5"
    core-hello-results/
    ├── cat
    │   └── test.txt
    ├── cowpy
    │   └── cowpy-test.txt
    └── pipeline_info
        ├── execution_report_2025-12-27_06-16-55.html
        ├── execution_report_2025-12-27_06-23-13.html
        ├── execution_report_2025-12-27_06-29-02.html
        ├── execution_report_2025-12-27_06-35-56.html
        ├── execution_timeline_2025-12-27_06-16-55.html
        ├── execution_timeline_2025-12-27_06-23-13.html
        ├── execution_timeline_2025-12-27_06-29-02.html
        ├── execution_timeline_2025-12-27_06-35-56.html
        ├── execution_trace_2025-12-27_06-16-55.txt
        ├── execution_trace_2025-12-27_06-23-13.txt
        ├── execution_trace_2025-12-27_06-29-02.txt
        ├── execution_trace_2025-12-27_06-35-56.txt
        ├── hello_software_versions.yml
        ├── params_2025-12-27_06-17-00.json
        ├── params_2025-12-27_06-23-17.json
        ├── params_2025-12-27_06-29-07.json
        ├── params_2025-12-27_06-36-01.json
        ├── pipeline_dag_2025-12-27_06-16-55.html
        ├── pipeline_dag_2025-12-27_06-23-13.html
        ├── pipeline_dag_2025-12-27_06-29-02.html
        └── pipeline_dag_2025-12-27_06-35-56.html
    ```

You can see that Nextflow created this hierarchy of directories based on the names of the workflow and of the module.

The code responsible lives in the `conf/modules.config` file.
This is the default `publishDir` configuration that is part of the nf-core template and applies to all processes:

```groovy
process {
    publishDir = [
        path: { "${params.outdir}/${task.process.tokenize(':')[-1].tokenize('_')[0].toLowerCase()}" },
        mode: params.publish_dir_mode,
        saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
    ]
}
```

This may look complicated, so let's look at each of the three components:

- **`path:`** Determines the output directory based on the process name.
  The full name of a process contained in `task.process` includes the hierarchy of workflow and module imports (such as `CORE_HELLO:HELLO:CAT_CAT`).
  The `tokenize` operations strip away that hierarchy to get just the process name, then take the first part before any underscore (if applicable), and convert it to lowercase.
  This is what determines that the results of `CAT_CAT` get published to `${params.outdir}/cat/`.
- **`mode:`** Controls how files are published (copy, symlink, etc.).
  This is configurable via the `params.publish_dir_mode` parameter.
- **`saveAs:`** Filters which files to publish.
  This example excludes `versions.yml` files by returning `null` for them, preventing them from being published.

This provides a consistent logic for organizing outputs.

The output looks even better when all the modules in a pipeline adopt this convention, so feel free to go delete the `publishDir` directives from the other modules in your pipeline.
This default will be applied even to modules that we didn't explicitly modify to follow nf-core guidelines.

That being said, you may decide you want to organize your inputs differently, and the good news is that it's easy to do so.

#### 1.5.3. Override the default

To override the default `publishDir` directive, you can simply add your own directives to the `conf/modules.config` file.

For example, you could override the default for a single process using the `withName:` selector, as in this example where we add a custom `publishDir` directive for the 'COWPY' process.

```groovy title="core-hello/conf/modules.config" linenums="13" hl_lines="8-10"
process {
    publishDir = [
        path: { "${params.outdir}/${task.process.tokenize(':')[-1].tokenize('_')[0].toLowerCase()}" },
        mode: params.publish_dir_mode,
        saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
    ]

    withName: 'COWPY' {
        ext.args = { "-c ${params.character}" }
        publishDir = [
            path: 'my_custom_results'
        ]
    }
}
```

We're not actually going to make that change, but feel free to play with this and see what logic you can implement.

The point is that this system allows gives you the best of both worlds: consistency by default and the flexibility to customize the configuration on demand.

To summarize, you get:

- **Single source of truth**: All publishing configuration lives in `modules.config`
- **Useful default**: Processes work out-of-the-box without per-module configuration
- **Easy customization**: Override publishing behavior in config, not in module code
- **Portable modules**: Modules don't hardcode output locations

This completes the set of nf-core module features you should absolutely learn to use, but there are others which you can read about in the [nf-core modules specifications](https://nf-co.re/docs/guidelines/components/modules).

### Takeaway

You now know how to adapt local modules to follow nf-core conventions:

- Design your modules to accept and propagate metadata tuples;
- Use `ext.args` to keep module interfaces minimal and portable;
- Use `ext.prefix` for configurable, standardized output file naming;
- Adopt the default centralized `publishDir` directive for a consistent results directory structure.

### What's next?

Learn how to use nf-core's built-in template-based tools to create modules the easy way.

---

## 2. Create a module with the nf-core tooling

Now that you've learned the nf-core module patterns by applying them manually, let's look at how you would create modules in practice.

### 2.1. Generate a module scaffold from a template

Similar to what exists for creating pipelines, the nf-core project provides tooling to generate properly structured modules based on a template, with all these patterns built in from the start.

#### 2.1.1. Run the module creation command

The `nf-core modules create` command generates a module template that already follows all the conventions you've learned.

Let's create a new version of the `COWPY` module with a minimal template by running this command:

```bash
nf-core modules create --empty-template COWPY
```

The `--empty-template` flag creates a clean starter template without extra code, making it easier to see the essential structure.

The command runs interactively, guiding you through the setup.
It automatically looks up tool information from package repositories like Bioconda and bio.tools to pre-populate metadata.

You'll be prompted for several configuration options:

- **Author information**: Your GitHub username for attribution
- **Resource label**: A predefined set of computational requirements.
  The nf-core project provides standard labels like `process_single` for lightweight tools and `process_high` for demanding ones.
  These labels help manage resource allocation across different execution environments.
- **Metadata requirement**: Whether the module needs sample-specific information via a `meta` map (usually yes for data processing modules).

The tool handles the complexity of finding package information and setting up the structure, allowing you to focus on implementing the tool's specific logic.

#### 2.1.2. Examine the module scaffold

The tool creates a complete module structure in `modules/local/` (or `modules/nf-core/` if you're in the nf-core/modules repository):

??? abstract "Directory contents"

    ```console
    modules/local/cowpy
    ├── environment.yml
    ├── main.nf
    ├── meta.yml
    └── tests
        └── main.nf.test
    ```

Each file serves a specific purpose:

- **`main.nf`**: Process definition with all the nf-core patterns built in
- **`meta.yml`**: Module documentation describing inputs, outputs, and the tool
- **`environment.yml`**: Conda environment specification for dependencies
- **`tests/main.nf.test`**: nf-test test cases to validate the module works

!!! tip "Learn more about testing"

    The generated test file uses nf-test, a testing framework for Nextflow pipelines and modules. To learn how to write and run these tests, see the [nf-test side quest](../side_quests/nf_test.md).

The generated `main.nf` includes all the patterns you just learned, plus some additional features:

```groovy title="modules/local/cowpy/main.nf" hl_lines="11 21 22"
process COWPY {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/YOUR-TOOL-HERE':
        'biocontainers/YOUR-TOOL-HERE' }"

    input:
    tuple val(meta), path(input)        // Pattern 1: Metadata tuples ✓

    output:
    tuple val(meta), path("*"), emit: output
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''                // Pattern 2: ext.args ✓
    def prefix = task.ext.prefix ?: "${meta.id}"  // Pattern 3: ext.prefix ✓

    """
    // Add your tool command here

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        COWPY: \$(cowpy --version)
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    echo $args
    touch ${prefix}.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        COWPY: \$(cowpy --version)
    END_VERSIONS
    """
}
```

Notice how all the patterns you applied manually above are already there!

The template also includes several additional nf-core conventions.
Some of these work out of the box, while others are placeholders we'll need to fill in, as described below.

**Features that work as-is:**

- **`tag "$meta.id"`**: Adds sample ID to process names in logs for easier tracking
- **`label 'process_single'`**: Resource label for configuring CPU/memory requirements
- **`when:` block**: Allows conditional execution via `task.ext.when` configuration

These features are already functional and make modules more maintainable.

**Placeholders we'll customize below:**

- **`input:` and `output:` blocks**: Generic declarations we'll update to match our tool
- **`script:` block**: Contains a comment where we'll add the `cowpy` command
- **`stub:` block**: Template we'll update to produce the correct outputs
- **Container and environment**: Placeholders we'll fill with package information

The next sections walk through completing these customizations.

### 2.2. Set up the container and conda environment

The nf-core guidelines require that we specify both a container and a Conda environment as part of the module.

#### 2.2.1. Container

For the container, you can use [Seqera Containers](https://seqera.io/containers/) to automatically build a container from any Conda package, including conda-forge packages.
In this case we are using the same prebuilt container as previously.

The default code offers to toggle between Docker and Singularity, but we're going to simplify that line and just specify the Docker container we got from Seqera Containers above.

=== "After"

```groovy title="modules/local/cowpy/main.nf" linenums="3" hl_lines="6"
process COWPY {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273"
```

=== "Before"

```groovy title="modules/local/cowpy/main.nf" linenums="3" hl_lines="6"
process COWPY {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/YOUR-TOOL-HERE':
        'biocontainers/YOUR-TOOL-HERE' }"
```

#### 2.2.2. Conda environment

For the Conda environment, the module code specifies `conda "${moduleDir}/environment.yml"` which means that it should be configured in the `environment.yml` file.

The module creation tool warned us that it couldn't find the `cowpy` package in Bioconda (the primary channel for bioinformatics tools).
However, `cowpy` is available in conda-forge, so you can complete the `environment.yml` like this:

=== "After"

    ```yaml title="modules/local/cowpy/environment.yml"  linenums="1" hl_lines="1 3 5"
    name: COWPY
    channels:
      - conda-forge
    dependencies:
      - cowpy=1.1.5
    ```

=== "Before"

    ```yaml title="modules/local/cowpy/environment.yml" linenums="1"
    ---
    # yaml-language-server: $schema=https://raw.githubusercontent.com/nf-core/modules/master/modules/environment-schema.json
    channels:
      - conda-forge
      - bioconda
    dependencies:
      # TODO nf-core: List required Conda package(s).
      #               Software MUST be pinned to channel (i.e. "bioconda"), version (i.e. "1.10").
      #               For Conda, the build (i.e. "h9402c20_2") must be EXCLUDED to support installation on different operating systems.
      - "YOUR-TOOL-HERE"
    ```

For submission to nf-core, we would have to follow the defaults more closely, but for our own use we can simplify the code in this way.

!!! tip "Bioconda vs conda-forge packages"

    - **Bioconda packages**: Automatically get BioContainers built, providing ready-to-use containers
    - **conda-forge packages**: Can use Seqera Containers to build containers on-demand from the Conda recipe

    Most bioinformatics tools are in Bioconda, but for conda-forge tools, Seqera Containers provides an easy solution for containerization.

### 2.3. Plug in the `COWPY` logic

Now let's update the code elements that are specific to what the `COWPY` process does: the inputs and outputs, and the script block.

#### 2.3.1. Inputs and outputs

The generated template includes generic input and output declarations that you'll need to customize for your specific tool.
Looking back at our manual `COWPY` module from section 1, we can use that as a guide.

Update the input and output blocks:

=== "After"

    ```groovy title="modules/local/cowpy/main.nf" linenums="8" hl_lines="2 5"
    input:
    tuple val(meta), path(input_file)

    output:
    tuple val(meta), path("${prefix}.txt"), emit: cowpy_output
    path "versions.yml"           , emit: versions
    ```

=== "Before"

    ```groovy title="modules/local/cowpy/main.nf" linenums="8" hl_lines="2 5"
    input:
    tuple val(meta), path(input)

    output:
    tuple val(meta), path("*"), emit: output
    path "versions.yml"           , emit: versions
    ```

This specifies:

- The input file parameter name (`input_file` instead of generic `input`)
- The output filename using the configurable prefix pattern (`${prefix}.txt` instead of wildcard `*`)
- A descriptive emit name (`cowpy_output` instead of generic `output`)

If you're using the Nextflow language server to validate syntax, the `${prefix}` part will be flagged as an error at this stage because we haven't added it to the script block yet.
Let's get to that now.

#### 2.3.2. The script block

The template provides a comment placeholder in the script block where you should add the actual tool command.

Based on the module we wrote manually earlier, we should make the following edits:

=== "After"

    ```groovy title="modules/local/cowpy/main.nf" linenums="15" hl_lines="3 6"
    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"

    """
    cat $input_file | cowpy $args > ${prefix}.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        COWPY: \$(cowpy --version)
    END_VERSIONS
    """
    ```

=== "Before"

    ```groovy title="modules/local/cowpy/main.nf" linenums="15" hl_lines="6"
    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    // Add your tool command here

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        COWPY: \$(cowpy --version)
    END_VERSIONS
    """
    ```

Key changes:

- Change `def prefix` to just `prefix` (without `def`) to make it accessible in the output block
- Replace the comment with the actual `cowpy` command that uses both `$args` and `${prefix}.txt`

Note that if we hadn't already done the work of adding the `ext.args` and `ext.prefix` configuration for the `COWPY` process to the `modules.config` file, we would need to do that now.

#### 2.3.3. Implementing the stub block

In the Nextflow context, a [stub](https://www.nextflow.io/docs/latest/process.html#stub) block allows you to define a lightweight, dummy script used for rapid prototyping and testing of a pipeline's logic without executing the actual command.

<!-- TODO (future) This is super glossed over but should really be explained or at least link out to an explanation about stubs (the reference doc isn't terribly helpful either). Right now this is likely to be mostly meaningless to anyone who doesn't already know about stubs. -->

Don't worry too much if this seems mysterious; we include this for completeness but you can also just delete the stub section if you don't want to deal with it, as it's completely optional.

=== "After"

    ```groovy title="modules/local/cowpy/main.nf" linenums="27" hl_lines="3 6"
    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch ${prefix}.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        COWPY: \$(cowpy --version)
    END_VERSIONS
    """
    ```

=== "Before"

    ```groovy title="modules/local/cowpy/main.nf" linenums="27" hl_lines="3 6"
    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    echo $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        COWPY: \$(cowpy --version)
    END_VERSIONS
    """
    ```

Key changes:

- Change `def prefix` to just `prefix` to match the script block
- Remove the `echo $args` line (which was just template placeholder code)
- The stub creates an empty `${prefix}.txt` file matching what the script block produces

This allows you to test workflow logic and file handling without waiting for the actual tool to run.

Once you've completed the environment setup (section 2.2), inputs/outputs (section 2.3.1), script block (section 2.3.2), and stub block (section 2.3.3), the module is ready to test!

### 2.4. Swap in the new `COWPY` module and run the pipeline

All we need to do to try out this new version of the `COWPY` module is to switch the import statement in the `hello.nf` workflow file to point to the new file.

=== "After"

    ```groovy title="workflows/hello.nf" linenums="1" hl_lines="10"
    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */
    include { paramsSummaryMap       } from 'plugin/nf-schema'
    include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
    include { sayHello               } from '../modules/local/sayHello.nf'
    include { convertToUpper         } from '../modules/local/convertToUpper.nf'
    include { COWPY                  } from '../modules/local/cowpy/main.nf'
    include { CAT_CAT                } from '../modules/nf-core/cat/cat/main'
    ```

=== "Before"

    ```groovy title="modules/local/cowpy/main.nf" linenums="1" hl_lines="10"
    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */
    include { paramsSummaryMap       } from 'plugin/nf-schema'
    include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
    include { sayHello               } from '../modules/local/sayHello.nf'
    include { convertToUpper         } from '../modules/local/convertToUpper.nf'
    include { COWPY                  } from '../modules/local/cowpy.nf'
    include { CAT_CAT                } from '../modules/nf-core/cat/cat/main'
    ```

Let's run the pipeline to test it.

```bash
nextflow run . --outdir core-hello-results -profile test,docker --validate_params false
```

??? success "Command output"

    ```console hl_lines="33"
      N E X T F L O W   ~  version 25.04.3

    Launching `./main.nf` [prickly_neumann] DSL2 - revision: b9e9b3b8de

    Input/output options
      input                     : /workspaces/training/hello-nf-core/core-hello/assets/greetings.csv
      outdir                    : core-hello-results

    Institutional config options
      config_profile_name       : Test profile
      config_profile_description: Minimal test dataset to check pipeline function

    Generic options
      validate_params           : false
      trace_report_suffix       : 2025-12-27_08-23-51

    Core Nextflow options
      runName                   : prickly_neumann
      containerEngine           : docker
      launchDir                 : /workspaces/training/hello-nf-core/core-hello
      workDir                   : /workspaces/training/hello-nf-core/core-hello/work
      projectDir                : /workspaces/training/hello-nf-core/core-hello
      userName                  : root
      profile                   : test,docker
      configFiles               : /workspaces/training/hello-nf-core/core-hello/nextflow.config

    !! Only displaying parameters that differ from the pipeline defaults !!
    ------------------------------------------------------
    executor >  local (8)
    [e9/008ede] CORE_HELLO:HELLO:sayHello (3)       [100%] 3 of 3 ✔
    [f0/d70cfe] CORE_HELLO:HELLO:convertToUpper (3) [100%] 3 of 3 ✔
    [be/0ecc58] CORE_HELLO:HELLO:CAT_CAT (test)     [100%] 1 of 1 ✔
    [11/8e082f] CORE_HELLO:HELLO:COWPY (test)       [100%] 1 of 1 ✔
    -[core/hello] Pipeline completed successfully-
    ```

This produces the same results as previously.

### Takeaway

You now know how to use the built-in nf-core tooling to create modules efficiently using templates rather than writing everything from scratch.

### What's next?

Learn what are the benefits of contributing modules to nf-core and what are the main steps and requirements involved.

---

## 3. Contributing modules back to nf-core

The [nf-core/modules](https://github.com/nf-core/modules) repository welcomes contributions of well-tested, standardized modules.

### 3.1. Why contribute?

Contributing your modules to nf-core:

- Makes your tools available to the entire nf-core community through the modules catalog at [nf-co.re/modules](https://nf-co.re/modules)
- Ensures ongoing community maintenance and improvements
- Provides quality assurance through code review and automated testing
- Gives your work visibility and recognition

### 3.2. Contributor's checklist

To contribute a module to nf-core, you will need to go through the following steps:

1. Check if it already exists at [nf-co.re/modules](https://nf-co.re/modules)
2. Fork the [nf-core/modules](https://github.com/nf-core/modules) repository
3. Use `nf-core modules create` to generate the template
4. Fill in the module logic and tests
5. Test with `nf-core modules test tool/subtool`
6. Lint with `nf-core modules lint tool/subtool`
7. Submit a pull request

For detailed instructions, see the [nf-core components tutorial](https://nf-co.re/docs/tutorials/nf-core_components/components).

### 3.3. Resources

- **Components tutorial**: [Complete guide to creating and contributing modules](https://nf-co.re/docs/tutorials/nf-core_components/components)
- **Module specifications**: [Technical requirements and guidelines](https://nf-co.re/docs/guidelines/components/modules)
- **Community support**: [nf-core Slack](https://nf-co.re/join) - Join the `#modules` channel

### Takeaway

You now know how to create nf-core modules! You learned the four key patterns that make modules portable and maintainable:

- **Metadata tuples** propagate metadata through the workflow
- **`ext.args`** simplifies module interfaces by handling optional arguments via configuration
- **`ext.prefix`** standardizes output file naming
- **Centralized publishing** via `publishDir` configured in `modules.config` rather than hardcoded in modules

By transforming `COWPY` step-by-step, you developed a deep understanding of these patterns, making you equipped to work with, debug, and create nf-core modules.
In practice, you'll use `nf-core modules create` to generate properly structured modules with these patterns built in from the start.

Finally, you learned how to contribute modules to the nf-core community, making tools available to researchers worldwide while benefiting from ongoing community maintenance.

### What's next?

When you're ready, continue to [Part 5: Input validation](./05_input_validation.md) to learn how to add schema-based input validation to your pipeline.
