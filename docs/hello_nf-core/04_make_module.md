# Part 4: Make an nf-core module

In this fourth part of the Hello nf-core training course, we show you how to create an nf-core module by learning the key conventions that make modules portable and maintainable.

The nf-core project provides a command (`nf-core modules create`) that generates properly structured module templates automatically.
However, for teaching purposes, we're going to **learn by doing**: transforming the local `cowpy` module in your `core-hello` pipeline into an nf-core-style module step-by-step.
This hands-on approach will help you understand the patterns deeply, making you better equipped to work with nf-core modules in practice.

We'll apply three essential nf-core patterns incrementally:

1. **Metadata tuples**: Accept and propagate sample metadata through the workflow
2. **`ext.args`**: Simplify the module interface by handling optional arguments via configuration
3. **`ext.prefix`**: Standardize output file naming with configurable prefixes

Once you understand these patterns, we'll show you how to use the official nf-core tooling to create modules efficiently.

!!! note

    This section assumes you have completed [Part 3: Use an nf-core module](./03_use_module.md) and have integrated the `CAT_CAT` module into your pipeline.

    If you didn't complete Part 3 or want to start fresh for this section, you can use the `core-hello-part3` solution as your starting point:

    ```bash
    cp -r hello-nf-core/solutions/core-hello-part3 core-hello
    cd core-hello
    ```

    This gives you a pipeline with the `CAT_CAT` module already integrated.

---

## 1. Transform cowpy into an nf-core module

In this section, we'll apply nf-core conventions to the local `cowpy` module in your `core-hello` pipeline, transforming it into a module that follows community standards.

!!! tip "Working directory"

    Make sure you're in the `core-hello` directory (your pipeline root) for all the commands and file edits in this section.

    ```bash
    cd core-hello
    ```

### 1.1. Update cowpy to use metadata tuples

Currently, we're extracting the file from `CAT_CAT`'s output tuple to pass to `cowpy`. It would be better to have `cowpy` accept metadata tuples directly, allowing metadata to flow through the entire workflow.

Open `core-hello/modules/local/cowpy.nf` and modify it to accept metadata tuples:

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

Key changes:

1. **Input**: Changed from `path input_file` to `tuple val(meta), path(input_file)` to accept metadata
2. **Output**: Changed to emit a tuple with metadata: `tuple val(meta), path("cowpy-${input_file}"), emit: cowpy_output`
3. **Named emit**: Added `emit: cowpy_output` to give the output channel a descriptive name

Now update the workflow to pass the tuple directly instead of extracting the file. Open `core-hello/workflows/hello.nf`:

=== "After"

    ```groovy title="core-hello/workflows/hello.nf" linenums="43" hl_lines="2"
        // generate ASCII art of the greetings with cowpy
        cowpy(CAT_CAT.out.file_out, params.character)
    ```

=== "Before"

    ```groovy title="core-hello/workflows/hello.nf" linenums="43" hl_lines="2-4"
        // generate ASCII art of the greetings with cowpy
        // Extract the file from the tuple since cowpy doesn't use metadata yet
        ch_for_cowpy = CAT_CAT.out.file_out.map{ meta, file -> file }
        cowpy(ch_for_cowpy, params.character)
    ```

Also update the emit block to use the named emit:

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

Test the workflow to ensure metadata flows through correctly:

```bash
nextflow run . --outdir core-hello-results -profile test,docker --validate_params false
```

The pipeline should run successfully with metadata now flowing from `CAT_CAT` through `cowpy`:

```console title="Output (excerpt)"
executor >  local (8)
[b2/4cf633] CORE_HELLO:HELLO:sayHello (2)       [100%] 3 of 3 ✔
[ed/ef4d69] CORE_HELLO:HELLO:convertToUpper (3) [100%] 3 of 3 ✔
[2d/32c93e] CORE_HELLO:HELLO:CAT_CAT (test)     [100%] 1 of 1 ✔
[da/6f3246] CORE_HELLO:HELLO:cowpy              [100%] 1 of 1 ✔
-[core/hello] Pipeline completed successfully-
```

### 1.2. Simplify the interface with ext.args

Now let's address another nf-core pattern: simplifying module interfaces by using `ext.args` for optional command-line arguments.

Currently, our `cowpy` module requires the `character` parameter to be passed as a separate input. While this works, nf-core modules use a different approach for **tool configuration arguments**: instead of adding input parameters for every tool option, they use `ext.args` to pass these via configuration. This keeps the module interface focused on essential data (files, metadata, and any mandatory per-sample parameters), while tool configuration options are handled through `ext.args`.

#### Understanding ext.args

The `task.ext.args` pattern is an nf-core convention for passing command-line arguments to tools through configuration rather than as process inputs. Instead of adding input parameters for tool options, nf-core modules accept arguments through the `ext.args` configuration directive.

!!! note "ext.args can be dynamic"

    While `ext.args` is configured outside the module, it can access metadata to provide sample-specific values using closures. For example: `ext.args = { meta.single_end ? '--single' : '--paired' }`

Benefits of this approach:

- **Clean interface**: The module focuses on essential data inputs (metadata and files)
- **Flexibility**: Users can specify tool arguments via configuration, including sample-specific values
- **Consistency**: All nf-core modules follow this pattern
- **Portability**: Modules can be reused without hardcoded tool options
- **No workflow changes**: Adding or changing tool options doesn't require updating workflow code

#### Centralized publishing configuration

Before we update the module to use `ext.args`, let's address an important nf-core convention: **modules should not contain hardcoded `publishDir` directives**.

Currently, our `cowpy` module has `publishDir 'results', mode: 'copy'` which hardcodes the output location.
In nf-core pipelines, publishing is instead configured in `conf/modules.config`.

The nf-core template includes a **default publishDir configuration** that applies to all processes:

```groovy
process {
    publishDir = [
        path: { "${params.outdir}/${task.process.tokenize(':')[-1].tokenize('_')[0].toLowerCase()}" },
        mode: params.publish_dir_mode,
        saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
    ]
}
```

This default automatically publishes outputs to `${params.outdir}/<process_name>/` for every process.
Individual processes can customize their publishing using `withName:` blocks in the same config file.

Benefits of this approach:

- **Single source of truth**: All publishing configuration lives in `modules.config`
- **Useful default**: Processes work out-of-the-box without per-module configuration
- **Easy customization**: Override publishing behavior in config, not in module code
- **Portable modules**: Modules don't hardcode output locations

#### Update the module

Now let's update the cowpy module to use `ext.args` and remove the local `publishDir`.

Open `modules/local/cowpy.nf`:

=== "After"

    ```groovy title="modules/local/cowpy.nf" linenums="1" hl_lines="16 18"
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

Key changes:

1. **Removed character input**: The module no longer requires `character` as a separate input
2. **Removed local publishDir**: Deleted the `publishDir 'results', mode: 'copy'` directive to rely on centralized configuration
3. **Added ext.args**: The line `def args = task.ext.args ?: ''` uses the Elvis operator (`?:`) to provide an empty string as default if `task.ext.args` is not set
4. **Updated command**: Changed from hardcoded `-c "$character"` to using the configurable `$args`

The module interface is now simpler - it only accepts the essential metadata and file inputs. By removing the local `publishDir`, we follow the nf-core convention of centralizing all publishing configuration in `modules.config`.

#### Configure ext.args

Now we need to configure the `ext.args` to pass the character option. This allows us to keep the module interface simple while still providing the character option at the pipeline level.

Open `conf/modules.config` and add the cowpy configuration:

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

This configuration passes the `params.character` value to cowpy's `-c` flag. Note that we use a closure (`{ "-c ${params.character}" }`) to allow the parameter to be evaluated at runtime.

Key points:

- The **module interface stays simple** - it only accepts the essential metadata and file inputs
- The **pipeline still exposes `params.character`** - users can configure it as before
- The **module is now portable** - it can be reused in other pipelines without expecting a specific parameter name
- Configuration is **centralized** in `modules.config`, keeping workflow logic clean

!!! note

    The `modules.config` file is where nf-core pipelines centralize per-module configuration. This separation of concerns makes modules more reusable across different pipelines.

#### Update the workflow

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

#### Test

Test that the workflow still works with the ext.args configuration. Let's specify a different character to verify the configuration is working (using `kosh`, one of the more... enigmatic options):

```bash
nextflow run . --outdir core-hello-results -profile test,docker --validate_params false --character kosh
```

The pipeline should run successfully. In the output, look for the cowpy process execution line which will show something like:

```console title="Output (excerpt)"
[bd/0abaf8] CORE_HELLO:HELLO:cowpy              [100%] 1 of 1 ✔```

Now let's verify that the `ext.args` configuration actually passed the character argument to the cowpy command. Use the task hash (the `bd/0abaf8` part) to inspect the `.command.sh` file in the work directory:

```bash
cat work/bd/0abaf8*/.command.sh
```

You should see the cowpy command with the `-c cow` argument:

```console title="Output"
#!/usr/bin/env bash
...
cat test.txt | cowpy -c kosh > cowpy-test.txt
```

This confirms that `task.ext.args` successfully passed the character parameter through the configuration rather than requiring it as a process input.

We can also check the output:

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


### 1.3. Add configurable output naming with ext.prefix

There's one more nf-core pattern we can apply: using `ext.prefix` for configurable output file naming.

#### Understanding ext.prefix

The `task.ext.prefix` pattern is another nf-core convention for standardizing output file naming across modules while keeping it configurable.

Benefits:

- **Standardized naming**: Output files are typically named using sample IDs from metadata
- **Configurable**: Users can override the default naming if needed
- **Consistent**: All nf-core modules follow this pattern
- **Predictable**: Easy to know what output files will be called

#### Update the module

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

#### Configure ext.prefix

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

#### Test and verify

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

### Takeaway

You now know how to adapt local modules to follow nf-core conventions:

- Update modules to accept and propagate metadata tuples
- Use `ext.args` to keep module interfaces minimal and portable
- Use `ext.prefix` for configurable, standardized output file naming
- Configure process-specific parameters through `modules.config`

### What's next?

Clean up by optionally removing the now-unused local module.

---

### 1.4. Optional: Clean up unused local modules

Now that we're using the nf-core `cat/cat` module, the local `collectGreetings` module is no longer needed.

Remove or comment out the import line for `collectGreetings` in `workflows/hello.nf`:

```groovy title="core-hello/workflows/hello.nf" linenums="10"
include { sayHello               } from '../modules/local/sayHello.nf'
include { convertToUpper         } from '../modules/local/convertToUpper.nf'
// include { collectGreetings    } from '../modules/local/collectGreetings.nf'  // No longer needed
include { cowpy                  } from '../modules/local/cowpy.nf'
include { CAT_CAT                } from '../modules/nf-core/cat/cat/main'
```

You can optionally delete the `collectGreetings.nf` file:

```bash
rm modules/local/collectGreetings.nf
```

However, you might want to keep it as a reference for understanding the differences between local and nf-core modules.

---

## 2. Use nf-core tooling to create modules

Now that you understand the nf-core module patterns by applying them manually, let's look at how you'd create modules in practice.
The nf-core project provides the `nf-core modules create` command that generates properly structured module templates with all these patterns built in from the start.

### 2.1. Using nf-core modules create

The `nf-core modules create` command generates a module template that already follows all the conventions you've learned.

For example, to create the `cowpy` module with a minimal template:

```bash
nf-core modules create --empty-template cowpy
```

The `--empty-template` flag creates a clean starter template without extra boilerplate, making it easier to see the essential structure.

The command runs interactively, guiding you through the setup. It automatically looks up tool information from package repositories like Bioconda and bio.tools to pre-populate metadata.

You'll be prompted for:

- **Author information**: Your GitHub username for attribution
- **Resource label**: The computational requirements (e.g., `process_single` for lightweight tools, `process_high` for demanding ones)
- **Metadata requirement**: Whether the module needs sample-specific information via a `meta` map (usually yes for data processing modules)

The tool handles the complexity of finding package information and setting up the structure, allowing you to focus on implementing the tool's specific logic.

#### What gets generated

The tool creates a complete module structure in `modules/local/` (or `modules/nf-core/` if you're in the nf-core/modules repository):

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

```groovy title="modules/local/cowpy/main.nf" hl_lines="12 21 22"
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
The template also includes several additional nf-core conventions that we didn't cover in the hands-on section:

- **`tag "$meta.id"`**: Adds sample ID to process names in logs for easier tracking
- **`label 'process_single'`**: Resource label for configuring CPU/memory requirements
- **`when:` block**: Allows conditional execution via `task.ext.when` configuration
- **`stub:` block**: Provides a fast mock implementation for testing pipeline logic without running the actual tool
- **`versions.yml` output**: Captures software version information for reproducibility
- **Container placeholders**: Template structure for Singularity and Docker containers

These additional conventions make modules more maintainable and provide better visibility into pipeline execution.

#### Completing the environment and container setup

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

Once you've completed the environment setup and filled in the command logic, the module is ready to test!

### 2.2. Contributing modules back to nf-core

The [nf-core/modules](https://github.com/nf-core/modules) repository welcomes contributions of well-tested, standardized modules.

#### Why contribute?

Contributing your modules to nf-core:

- Makes your tools available to the entire nf-core community through the modules catalog at [nf-co.re/modules](https://nf-co.re/modules)
- Ensures ongoing community maintenance and improvements
- Provides quality assurance through code review and automated testing
- Gives your work visibility and recognition

#### Contributing workflow

To contribute a module to nf-core:

1. Check if it already exists at [nf-co.re/modules](https://nf-co.re/modules)
2. Fork the [nf-core/modules](https://github.com/nf-core/modules) repository
3. Use `nf-core modules create` to generate the template
4. Fill in the module logic and tests
5. Test with `nf-core modules test tool/subtool`
6. Lint with `nf-core modules lint tool/subtool`
7. Submit a pull request

For detailed instructions, see the [nf-core components tutorial](https://nf-co.re/docs/tutorials/nf-core_components/components).

#### Resources

- **Components tutorial**: [Complete guide to creating and contributing modules](https://nf-co.re/docs/tutorials/nf-core_components/components)
- **Module specifications**: [Technical requirements and guidelines](https://nf-co.re/docs/guidelines/components/modules)
- **Community support**: [nf-core Slack](https://nf-co.re/join) - Join the `#modules` channel

---

## Takeaway

You now know how to create nf-core modules! You learned the three key patterns that make modules portable and maintainable:

- **Metadata tuples** track sample information through the workflow
- **`ext.args`** simplifies module interfaces by handling optional arguments via configuration
- **`ext.prefix`** standardizes output file naming
- **Centralized configuration** in `modules.config` keeps modules reusable

By transforming `cowpy` step-by-step, you developed a deep understanding of these patterns—making you equipped to work with, debug, and create nf-core modules.
In practice, you'll use `nf-core modules create` to generate properly structured modules with these patterns built in from the start.

Finally, you learned how to contribute modules to the nf-core community, making tools available to researchers worldwide while benefiting from ongoing community maintenance.

## What's next?

Continue to [Part 5: Input validation](./05_input_validation.md) to learn how to add schema-based input validation to your pipeline, or explore other nf-core modules you might add to enhance your pipeline further.
