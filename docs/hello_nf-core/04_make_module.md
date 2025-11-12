# Part 4: Make an nf-core module

In this fourth part of the Hello nf-core training course, we show you how to create an nf-core module by applying the key conventions that make modules portable and maintainable.

The nf-core project provides a command (`nf-core modules create`) that generates properly structured module templates automatically, similar to what we used for the workflow in Part 2.
However, for teaching purposes, we're going to start by doing it manually: transforming the local `cowpy` module in your `core-hello` pipeline into an nf-core-style module step-by-step.
After that, we'll show you how to use the template-based module creation to work more efficiently in the future.

!!! note

    This section assumes you have completed [Part 3: Use an nf-core module](./03_use_module.md) and have integrated the `CAT_CAT` module into your pipeline.

    If you didn't complete Part 3 or want to start fresh for this section, you can use the `core-hello-part3` solution as your starting point:

    ```bash
    cp -r hello-nf-core/solutions/core-hello-part3 core-hello
    cd core-hello
    ```

    This gives you a pipeline with the `CAT_CAT` module already integrated.

---

## 1. Transform `cowpy` into an nf-core module

In this section, we'll apply nf-core conventions to the local `cowpy` module in your `core-hello` pipeline, transforming it into a module that follows community standards.

We'll apply the following nf-core conventions incrementally:

1. **Update `cowpy` to use metadata tuples** to propagate sample metadata through the workflow.
2. **Centralize tool argument configuration with `ext.args`** to increase module versatility while keeping the interface minimal.
3. **Standardize output naming with `ext.prefix`** to promote consistency.
4. **Centralize the publishing configuration** to promote consistency.

<!-- TODO: any additional comment? -->

!!! tip "Working directory"

    Make sure you're in the `core-hello` directory (your pipeline root) for all the commands and file edits in this section.

    ```bash
    cd core-hello
    ```

### 1.1. Update `cowpy` to use metadata tuples

In the current version of the `core-hello` pipeline, we're extracting the file from `CAT_CAT`'s output tuple to pass to `cowpy`.

<!-- TODO: add a diagram to recap current state -->

It would be better to have `cowpy` accept metadata tuples directly, allowing metadata to flow on through the workflow.
To the end, we'll need to make the following changes:

1. Update the input and output definitions
2. Update the process call in the workflow
3. Update the emit block in the workflow

Once we've done all that, we'll run the pipeline to test that everything still works as before.

<!-- TODO: outline steps -->

#### 1.1.1. Update the input and output definitions

Let's get started!
Open the `cowpy.nf` module file (under `core-hello/modules/local/`) and modify it to accept metadata tuples as shown below.

=== "After"

    ```groovy title="core-hello/modules/local/cowpy.nf" linenums="1" hl_lines="12 16"
    #!/usr/bin/env nextflow

    // Generate ASCII art with cowpy (https://github.com/jeffbuttars/cowpy)
    process cowpy {

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

=== "Before"

    ```groovy title="core-hello/modules/local/cowpy.nf" linenums="1" hl_lines="12 16"
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

As you can see, we changed both the **main input** and the **output** to a tuple that follows the `tuple val(meta), path(input_file)` pattern introduced in Part 3 of this training.
For the output, we also took this opportunity to add `emit: cowpy_output` in order to give the output channel a descriptive name.

Now that we've changed what the process expects, we need to update what we provide to it in the process call.

#### 1.1.2. Update the process call in the workflow

The good news is that this change simplifies the process call.
Now that the output of `CAT_CAT` and the input of `cowpy` are the same 'shape', i.e. they both consist of a `tuple val(meta), path(input_file)` structure, we can simply connect them directly instead of having to extract the file explicitly from the output of the `CAT_CAT` process.

Open the `hello.nf` workflow file (under `core-hello/workflows/`) and update the call to `cowpy` as shown below.

=== "After"

    ```groovy title="core-hello/workflows/hello.nf" linenums="43" hl_lines="2"
        // generate ASCII art of the greetings with cowpy
        cowpy(CAT_CAT.out.file_out, params.character)
    ```

=== "Before"

    ```groovy title="core-hello/workflows/hello.nf" linenums="43" hl_lines="5"
        // extract the file from the tuple since cowpy doesn't use metadata yet
        ch_for_cowpy = CAT_CAT.out.file_out.map{ meta, file -> file }

        // generate ASCII art of the greetings with cowpy
        cowpy(ch_for_cowpy, params.character)
    ```

You can see that we no longer need to construct the `ch_for_cowpy` channel, so that line (and its comment line) can be deleted entirely.

#### 1.1.3. Update the emit block in the workflow

Since `cowpy` now emits a named output, `cowpy_output`, we can update the `hello.nf` workflow's `emit:` block to use that.

=== "After"

    ```groovy title="core-hello/workflows/hello.nf" linenums="60" hl_lines="2"
        emit:
        cowpy_hellos   = cowpy.out.cowpy_output
        versions       = ch_versions                 // channel: [ path(versions.yml) ]
    ```

=== "Before"

    ```groovy title="core-hello/workflows/hello.nf" linenums="60" hl_lines="2"
        emit:
        cowpy_hellos   = cowpy.out
        versions       = ch_versions                 // channel: [ path(versions.yml) ]
    ```

This is technically not required, but it's good practice to refer to named outputs whenever possible.

#### 1.1.4. Run the pipeline to test it

Let's run the workflow to test that everything is working correctly after these changes.

```bash
nextflow run . --outdir core-hello-results -profile test,docker --validate_params false
```

The pipeline should run successfully, with metadata now flowing from `CAT_CAT` through `cowpy`:

```console title="Output (excerpt)"
executor >  local (8)
[b2/4cf633] CORE_HELLO:HELLO:sayHello (2)       [100%] 3 of 3 ✔
[ed/ef4d69] CORE_HELLO:HELLO:convertToUpper (3) [100%] 3 of 3 ✔
[2d/32c93e] CORE_HELLO:HELLO:CAT_CAT (test)     [100%] 1 of 1 ✔
[da/6f3246] CORE_HELLO:HELLO:cowpy              [100%] 1 of 1 ✔
-[core/hello] Pipeline completed successfully-
```

That completes what we needed to do to make `cowpy` handle metadata tuples.
Now let's look at what else we can do to take advantage of nf-core module patterns.

### 1.2. Centralize tool argument configuration with `ext.args`

In its current state, the `cowpy` process expects to receive a value for the `character` parameter.
As a result, we have to provide a value every time we call the process, even if we'd be happy with the defaults set by the tool.
For `cowpy` this is admittedly not a big problem, but for tools with many optional parameters, it can get quite cumbersome.

The nf-core project recommends using a Nextflow feature called `ext.args`, which makes it possible to manage tool arguments more conveniently. <!-- TODO: add URL to doc -->
Specifically, nf-core modules use a special configuration variable called `task.ext.args`.

Instead of declaring process inputs for every tool option, you write the module to reference `task.ext.args` in its command line.
Then it's just a matter of adding the arguments and values you want to use in the `modules.config` file, which consolidates configuration details for all modules.
Nextflow will add those arguments with their values into the tool command line at runtime.

Let's apply it to the `cowpy` module.
We're going to need to make the following changes:

1. Update the `cowpy` module
2. Add the character parameter to `ext.args`
3. Update the `hello.nf` workflow

Once we've done all that, we'll run the pipeline to test that everything still works as before.

#### 1.2.1. Update the `cowpy` module

Let's get started!
Open the `cowpy.nf` module file (under `core-hello/modules/local/`) and modify it to use `ext.args` as shown below.

=== "After"

    ```groovy title="modules/local/cowpy.nf" linenums="1" hl_lines="16 18"
    #!/usr/bin/env nextflow

    // Generate ASCII art with cowpy (https://github.com/jeffbuttars/cowpy)
    process cowpy {

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

    ```groovy title="core-hello/modules/local/cowpy.nf" linenums="1" hl_lines="6 13"
    #!/usr/bin/env nextflow

    // Generate ASCII art with cowpy (https://github.com/jeffbuttars/cowpy)
    process cowpy {

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
   Going forward, we'll supply that argument via the `task.ext.args` configuration (see next step).

2. **In the `script:` block, we added the line `def args = task.ext.args ?: ''`.**
   That line uses the `?:` operator to determine the value of the `args` variable: the content of `task.ext.args` if it is not empty, or an empty string if it is.

3. **In the command line, we replaced `-c "$character"` with `$args`.**
   This is where Nextflow will inject any tool arguments set in `task.ext.args`.

The module interface is now simpler - it only accepts the essential metadata and file inputs. By removing the local `publishDir`, we follow the nf-core convention of centralizing all publishing configuration in `modules.config`.

!!! note

    The `?:` operator is often called the 'Elvis operator' because it looks like a sideways Elvis Presley face, with the `?` character symbolizing the wave in his hair.

#### 1.2.2. Add the `character` parameter to `ext.args`

Now that we've taken the `character` declaration out of the module, we've got to add it to `ext.args` in the `modules.config` configuration file.

Specifically, we're going to add this little chunk of code to the `process {}` block:

```groovy title="Configuration syntax"
withName: 'cowpy' {
    ext.args = { "-c ${params.character}" }
}
```

The `withName:` syntax assigns this configuration to the `cowpy` process only, and `ext.args = { "-c ${params.character}" }` simply composes a string that will include the value of the `character` parameter.
Note the use of curly braces, which tell Nextflow to evaluate the value of the parameter at runtime.

Makes sense? Let's add it in.

Open `conf/modules.config` and add the configuration code inside the `process {}` block as shown below.

=== "After"

    ```groovy title="core-hello/conf/modules.config" linenums="13" hl_lines="6-8"
    process {
        publishDir = [
            path: { "${params.outdir}/${task.process.tokenize(':')[-1].tokenize('_')[0].toLowerCase()}" },
        ]

        withName: 'cowpy' {
            ext.args = { "-c ${params.character}" }
        }
    }
    ```

=== "Before"

    ```groovy title="core-hello/conf/modules.config" linenums="13"
    process {
        publishDir = [
            path: { "${params.outdir}/${task.process.tokenize(':')[-1].tokenize('_')[0].toLowerCase()}" },
        ]
    }
    ```

Hopefully you can imagine having all the modules in a pipeline have their `ext.args` specified in this file, with the following benefits:

- The **module interface stays simple** - It only accepts the essential metadata and file inputs
- The **pipeline still exposes `params.character`** - End-users can still configure it as before
- The **module is now portable** - It can be reused in other pipelines without expecting a specific parameter name
- The configuration is **centralized** in `modules.config`, keeping workflow logic clean

By using the `modules.config` file as the place where all pipelines centralize per-module configuration, we make our modules more reusable across different pipelines.

#### 1.2.3. Update the `hello.nf` workflow

Since the cowpy module no longer requires the `character` parameter as an input, we need to update the workflow call.

Open `workflows/hello.nf` and update the cowpy call:

=== "After"

    ```groovy title="core-hello/workflows/hello.nf" linenums="39" hl_lines="2"
        // generate ASCII art of the greetings with cowpy
        cowpy(CAT_CAT.out.file_out)
    ```

=== "Before"

    ```groovy title="core-hello/workflows/hello.nf" linenums="39" hl_lines="2"
        // generate ASCII art of the greetings with cowpy
        cowpy(CAT_CAT.out.file_out, params.character)
    ```

The workflow code is now cleaner - we don't need to pass `params.character` directly to the process. The module interface is kept minimal, making it more portable, while the pipeline still provides the explicit option through configuration.

#### 1.2.4. Run the pipeline to test it

Test that the workflow still works with the ext.args configuration. Let's specify a different character to verify the configuration is working (using `kosh`, one of the more... enigmatic options):

```bash
nextflow run . --outdir core-hello-results -profile test,docker --validate_params false --character kosh
```

The pipeline should run successfully. In the output, look for the cowpy process execution line which will show something like:

```console title="Output (excerpt)"
[bd/0abaf8] CORE_HELLO:HELLO:cowpy              [100%] 1 of 1 ✔
```

Let's verify that the `ext.args` configuration worked by checking the output. Use the task hash (the `bd/0abaf8` part) to look at the output file:

```bash
cat work/bd/0abaf8*/cowpy-test.txt
```

```console title="Output"
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

You should see the ASCII art displayed with the kosh character, confirming that the `ext.args` configuration worked!

!!! note "Optional: Inspect the command file"

    If you want to see exactly how the configuration was applied, you can inspect the `.command.sh` file:

    ```bash
    cat work/bd/0abaf8*/.command.sh
    ```

    You'll see the cowpy command with the `-c kosh` argument:

    ```console
    #!/usr/bin/env bash
    ...
    cat test.txt | cowpy -c kosh > cowpy-test.txt
    ```

    This shows that the `.command.sh` file was generated correctly based on the `ext.args` configuration.

<!-- TODO: commentary -->

It may seem unnecessary for a simple tool like `cowpy`, but it can make a big difference for data analysis tools that have a lot of optional arguments.
This approach keeps the module interface focused on essential data (files, metadata, and any mandatory per-sample parameters), while options that control the behavior of the tool are handled separately through configuration.

To summarize the benefits:

- **Clean interface**: The module focuses on essential data inputs (metadata and files)
- **Flexibility**: Users can specify tool arguments via configuration, including sample-specific values
- **Consistency**: All nf-core modules follow this pattern
- **Portability**: Modules can be reused without hardcoded tool options
- **No workflow changes**: Adding or changing tool options doesn't require updating workflow code

!!! note "ext.args can do more"

    The `ext.args` system has powerful additional capabilities not covered here, including switching argument values dynamically based on metadata. See the [nf-core module specifications](https://nf-co.re/docs/guidelines/components/modules) for more details.

### 1.3. Standardize output naming with `ext.prefix`

There's one more nf-core pattern we can apply: using `ext.prefix` for configurable output file naming.

Currently, the `cowpy` module includes a `publishDir` directive, making publishing decisions at the module level. We can't control where outputs go at the workflow level - each module makes its own publishing decisions.

The `task.ext.prefix` pattern is another nf-core convention for standardizing output file naming across modules while keeping it configurable.

Benefits:

- **Standardized naming**: Output files are typically named using sample IDs from metadata
- **Configurable**: Users can override the default naming if needed
- **Consistent**: All nf-core modules follow this pattern
- **Predictable**: Easy to know what output files will be called

#### 1.3.1. Update the module

Let's update the cowpy module to use `ext.prefix` for output file naming.

Open `modules/local/cowpy.nf` and change as follows:

=== "After"

    ```groovy title="core-hello/modules/local/cowpy.nf" linenums="1" hl_lines="13 17 19"
    #!/usr/bin/env nextflow

    // Generate ASCII art with cowpy (https://github.com/jeffbuttars/cowpy)
    process cowpy {

        container 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'
        conda 'conda-forge::cowpy==1.1.5'

        input:
            tuple val(meta), path(input_file)

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

    ```groovy title="core-hello/modules/local/cowpy.nf" linenums="1" hl_lines="13 18"
    #!/usr/bin/env nextflow

    // Generate ASCII art with cowpy (https://github.com/jeffbuttars/cowpy)
    process cowpy {

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

Key changes:

1. **Added ext.prefix**: `prefix = task.ext.prefix ?: "${meta.id}"` provides a configurable prefix with a sensible default (the sample ID)
2. **Updated output**: Changed from hardcoded `cowpy-${input_file}` to `${prefix}.txt`
3. **Updated command**: Uses the configured prefix for the output filename

Note that the local `publishDir` has already been removed in the previous step, so we're continuing with the centralized configuration approach.

#### 1.3.2. Configure ext.prefix

To maintain the same output file naming as before (`cowpy-<id>.txt`), we can configure `ext.prefix` in modules.config.

Update `conf/modules.config`:

=== "After"

    ```groovy title="core-hello/conf/modules.config" linenums="21" hl_lines="3"
        withName: 'cowpy' {
            ext.args = { "-c ${params.character}" }
            ext.prefix = { "cowpy-${meta.id}" }
        }
    ```

=== "Before"

    ```groovy title="core-hello/conf/modules.config" linenums="21"
        withName: 'cowpy' {
            ext.args = { "-c ${params.character}" }
        }
    ```

Note that we use a closure (`{ "cowpy-${meta.id}" }`) which has access to `meta` because it's evaluated in the context of the process execution.

!!! note

    The `ext.prefix` closure has access to `meta` because the configuration is evaluated in the context of the process execution, where metadata is available.

#### 1.3.3. Test and verify

Test the workflow once more:

```bash
rm -rf core-hello-results
nextflow run . --outdir core-hello-results -profile test,docker --validate_params false
```

Check the outputs:

```bash
ls core-hello-results/cowpy/
```

You should see the cowpy output files with the same naming as before: `cowpy-test.txt` (based on the batch name). This demonstrates how `ext.prefix` allows you to maintain your preferred naming convention while keeping the module interface flexible.

If you wanted to change the naming (for example, to just `test.txt`), you would only need to modify the `ext.prefix` configuration - no changes to the module or workflow code would be required.

### 1.4. Centralize the publishing configuration

For output publishing, nf-core pipelines centralize control at the workflow level by configuring `publishDir` in `conf/modules.config` rather than in individual modules.

Currently, our `cowpy` module has `publishDir 'results', mode: 'copy'` which hardcodes the output location.
In nf-core pipelines, publishing is instead configured in `conf/modules.config`.

#### 1.4.1. Update the module

<!-- TODO: show removing the publishDir -->

        publishDir 'results', mode: 'copy'

#### 1.4.2. Run the pipeline to see what happens

<!-- TODO: show where the results end up -->
<!-- TODO: then show how that comes from the default configuration -->

The nf-core template includes a default `publishDir` configuration that applies to all processes:

```groovy
process {
    publishDir = [
        path: { "${params.outdir}/${task.process.tokenize(':')[-1].tokenize('_')[0].toLowerCase()}" },
        mode: params.publish_dir_mode,
        saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
    ]
}
```

This looks complicated, but it breaks down into three parts:

- **path**: Determines the output directory based on the process name. When processes run, their full name includes the workflow hierarchy (like `CORE_HELLO:HELLO:SAMTOOLS_SORT`). The `tokenize` operations strip away that hierarchy to get just the process name, then take the first part before any underscore, and convert it to lowercase. So `SAMTOOLS_SORT` would publish to `${params.outdir}/samtools/`.
- **mode**: Controls how files are published (copy, symlink, etc.), configurable via the `params.publish_dir_mode` parameter.
- **saveAs**: Filters which files to publish. This example excludes `versions.yml` files by returning `null` for them, preventing them from being published.

Individual processes can override this default using `withName:` blocks in the same config file.

#### 1.4.3. Override the default

<!-- TODO: show how to override -->

<!-- TODO: wrap-up commentary -->

Benefits of this approach:

- **Single source of truth**: All publishing configuration lives in `modules.config`
- **Useful default**: Processes work out-of-the-box without per-module configuration
- **Easy customization**: Override publishing behavior in config, not in module code
- **Portable modules**: Modules don't hardcode output locations

For more details, see the [nf-core modules specifications](https://nf-co.re/docs/guidelines/components/modules).

### Takeaway

You now know how to adapt local modules to follow nf-core conventions:

- Update modules to accept and propagate metadata tuples
- Use `ext.args` to keep module interfaces minimal and portable
- Use `ext.prefix` for configurable, standardized output file naming
- Configure process-specific parameters through `modules.config`

### What's next?

<!-- TODO: update this -->

---

## 2. Generate modules with nf-core tools

Now that you've learned the nf-core module patterns by applying them manually, let's look at how you'd create modules in practice.
The nf-core project provides the `nf-core modules create` command that generates properly structured module templates with all these patterns built in from the start.

### 2.1. Using nf-core modules create

The `nf-core modules create` command generates a module template that already follows all the conventions you've learned.

For example, to create the `cowpy` module with a minimal template:

```bash
nf-core modules create --empty-template cowpy
```

The `--empty-template` flag creates a clean starter template without extra code, making it easier to see the essential structure.

The command runs interactively, guiding you through the setup. It automatically looks up tool information from package repositories like Bioconda and bio.tools to pre-populate metadata.

You'll be prompted for several configuration options:

- **Author information**: Your GitHub username for attribution
- **Resource label**: A predefined set of computational requirements. nf-core provides standard labels like `process_single` for lightweight tools and `process_high` for demanding ones. These labels help manage resource allocation across different execution environments.
- **Metadata requirement**: Whether the module needs sample-specific information via a `meta` map (usually yes for data processing modules)

The tool handles the complexity of finding package information and setting up the structure, allowing you to focus on implementing the tool's specific logic.

### 2.2. What gets generated

The tool creates a complete module structure in `modules/local/` (or `modules/nf-core/` if you're in the nf-core/modules repository):

??? example "Directory contents"

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

    The generated test file uses nf-test, a testing framework for Nextflow pipelines and modules. To learn how to write and run these tests, see the [nf-test side quest](../../side_quests/nf_test/).

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
    def args = task.ext.args ?: ''              // Pattern 2: ext.args ✓
    def prefix = task.ext.prefix ?: "${meta.id}"  // Pattern 3: ext.prefix ✓

    """
    // Add your tool command here

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cowpy: \$(cowpy --version)
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
        cowpy: \$(cowpy --version)
    END_VERSIONS
    """
}
```

Notice how all three patterns you applied manually are already there!
The template also includes several additional nf-core conventions.
Some of these work out of the box, while others are placeholders we'll need to fill in:

**Features that work as-is:**

- **`tag "$meta.id"`**: Adds sample ID to process names in logs for easier tracking
- **`label 'process_single'`**: Resource label for configuring CPU/memory requirements
- **`when:` block**: Allows conditional execution via `task.ext.when` configuration

These features are already functional and make modules more maintainable.

**Placeholders we'll customize below:**

- **`input:` and `output:` blocks**: Generic declarations we'll update to match our tool
- **`script:` block**: Contains a comment where we'll add the cowpy command
- **`stub:` block**: Template we'll update to produce the correct outputs
- **Container and environment**: Placeholders we'll fill with package information

The next sections walk through completing these customizations.

### 2.3. Completing the environment and container setup

In the case of cowpy, the tool warned that it couldn't find the package in Bioconda (the primary channel for bioinformatics tools).
However, cowpy is available in conda-forge, so you would complete the `environment.yml` like this:

```yaml title="modules/local/cowpy/environment.yml"
name: cowpy
channels:
  - conda-forge
dependencies:
  - cowpy=1.1.5
```

For the container, you can use [Seqera Containers](https://seqera.io/containers/) to automatically build a container from any Conda package, including conda-forge packages:

```groovy
container "community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273"
```

!!! tip "Bioconda vs conda-forge packages"

    - **Bioconda packages**: Automatically get BioContainers built, providing ready-to-use containers
    - **conda-forge packages**: Can use Seqera Containers to build containers on-demand from the Conda recipe

    Most bioinformatics tools are in Bioconda, but for conda-forge tools, Seqera Containers provides an easy solution for containerization.

### 2.4. Defining inputs and outputs

The generated template includes generic input and output declarations that you'll need to customize for your specific tool.
Looking back at our manual cowpy module from section 1, we can use that as a guide.

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

### 2.5. Writing the script block

The template provides a comment placeholder where you add the actual tool command.
We can reference our manual module from section 1.3.2 for the command logic:

=== "After"

    ```groovy title="modules/local/cowpy/main.nf" linenums="15" hl_lines="3 6"
    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"

    """
    cat $input_file | cowpy $args > ${prefix}.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cowpy: \$(cowpy --version)
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
        cowpy: \$(cowpy --version)
    END_VERSIONS
    """
    ```

Key changes:

- Change `def prefix` to just `prefix` (without `def`) so it's accessible in the output block
- Replace the comment with the actual cowpy command that uses both `$args` and `${prefix}.txt`

### 2.6. Implementing the stub block

The stub block provides a fast mock implementation for testing pipeline logic without running the actual tool.
It must produce the same output files as the script block:

=== "After"

    ```groovy title="modules/local/cowpy/main.nf" linenums="27" hl_lines="3 6"
    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch ${prefix}.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cowpy: \$(cowpy --version)
    END_VERSIONS
    """
    ```

=== "Before"

    ```groovy title="modules/local/cowpy/main.nf" linenums="27" hl_lines="5-6"
    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    echo $args
    touch ${prefix}.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cowpy: \$(cowpy --version)
    END_VERSIONS
    """
    ```

Key changes:

- Change `def prefix` to just `prefix` to match the script block
- Remove the `echo $args` line (which was just template placeholder code)
- The stub creates an empty `${prefix}.txt` file matching what the script block produces

This allows you to test workflow logic and file handling without waiting for the actual tool to run.

Once you've completed the environment setup (section 2.1.2), inputs/outputs (section 2.1.3), script block (section 2.1.4), and stub block (section 2.1.5), the module is ready to test!

### Takeaway

<!-- TODO: update this -->

### What's next?

<!-- TODO: update this -->

---

## 3. Contributing modules back to nf-core

The [nf-core/modules](https://github.com/nf-core/modules) repository welcomes contributions of well-tested, standardized modules.

### 3.1. Why contribute?

Contributing your modules to nf-core:

- Makes your tools available to the entire nf-core community through the modules catalog at [nf-co.re/modules](https://nf-co.re/modules)
- Ensures ongoing community maintenance and improvements
- Provides quality assurance through code review and automated testing
- Gives your work visibility and recognition

### 3.2. Contributing workflow

To contribute a module to nf-core:

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

## Takeaway

You now know how to create nf-core modules! You learned the four key patterns that make modules portable and maintainable:

- **Metadata tuples** track sample information through the workflow
- **`ext.args`** simplifies module interfaces by handling optional arguments via configuration
- **`ext.prefix`** standardizes output file naming
- **Centralized publishing** via `publishDir` configured in `modules.config` rather than hardcoded in modules

By transforming `cowpy` step-by-step, you developed a deep understanding of these patterns, making you equipped to work with, debug, and create nf-core modules.
In practice, you'll use `nf-core modules create` to generate properly structured modules with these patterns built in from the start.

Finally, you learned how to contribute modules to the nf-core community, making tools available to researchers worldwide while benefiting from ongoing community maintenance.

## What's next?

Continue to [Part 5: Input validation](./05_input_validation.md) to learn how to add schema-based input validation to your pipeline, or explore other nf-core modules you might add to enhance your pipeline further.
