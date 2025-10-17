# Part 3: Add an existing nf-core module

In this third part of the Hello nf-core training course, we show you how to add an existing nf-core module to your pipeline.

One of the great advantages of nf-core pipelines is the ability to leverage pre-built, tested modules from the nf-core/modules repository. Rather than writing every process from scratch, you can install and use community-maintained modules that follow best practices.

In this section, we'll replace the custom `collectGreetings` module with the `cat/cat` module from nf-core/modules.

!!! note

    This section assumes you have completed [Part 2: Rewrite Hello for nf-core](./02_rewrite_hello.md) and have a working `core-hello` pipeline.

---

## 1. Find and explore the cat/cat module

The `collectGreetings` process in our pipeline uses the Unix `cat` command to concatenate multiple greeting files into one. This is a perfect use case for the nf-core `cat/cat` module, which is designed specifically for concatenating files.

### 1.1. Browse available modules on the nf-core website

The nf-core project maintains a centralized catalog of modules at [https://nf-co.re/modules](https://nf-co.re/modules).

Navigate to the modules page in your web browser and use the search bar to search for "cat_cat".

You should see `cat/cat` in the search results. Click on it to view the module documentation.

The module page shows:

- A description: "A module for concatenation of gzipped or uncompressed files"
- Installation command: `nf-core modules install cat/cat`
- Input and output channel structure
- Available parameters

### 1.2. List available modules from the command line

You can also search for modules directly from the command line using nf-core tools.

```bash
nf-core modules list remote
```

This will display a list of all available modules in the nf-core/modules repository. You can scroll through or pipe to `grep` to find specific modules:

```bash
nf-core modules list remote | grep cat
```

```console title="Output"
cat/cat
```

### 1.3. Get detailed information about the module

To see detailed information about a specific module, use the `info` command:

```bash
nf-core modules info cat/cat
```

This displays documentation about the module, including its inputs, outputs, and basic usage information.

### Takeaway

You now know how to find and explore available nf-core modules using both the website and command-line tools.

### What's next?

Learn how to install the module in your pipeline.

---

## 2. Install and import the module

Now that we've identified the `cat/cat` module as a suitable replacement for our custom `collectGreetings` process, let's install it in our pipeline.

### 2.1. Install the cat/cat module

!!! note

    Make sure you are in the `core-hello` directory (your pipeline root) in your terminal before running the module installation command.

Navigate to your pipeline directory and run the installation command:

```bash
cd core-hello
nf-core modules install cat/cat
```

The tool will prompt you to confirm the installation. Press Enter to accept the default options.

```console title="Output"
INFO     Installing 'cat/cat'
INFO     Include statement: include { CAT_CAT } from '../modules/nf-core/cat/cat/main'
```

The command automatically:

- Downloads the module files to `modules/nf-core/cat/cat/`
- Updates `modules.json` to track the installed module
- Provides you with the correct `include` statement to use in your workflow

### 2.2. Verify the module installation

Let's check that the module was installed correctly:

```bash
tree modules/nf-core/cat
```

```console title="Output"
modules/nf-core/cat
└── cat
    ├── environment.yml
    ├── main.nf
    ├── meta.yml
    └── tests
        ├── main.nf.test
        ├── main.nf.test.snap
        ├── nextflow.config
        └── tags.yml
```

You can also verify the installation by listing locally installed modules:

```bash
nf-core modules list local
```

```console title="Output"
INFO     Modules installed in '.':

┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━┳━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
┃ Module Name                ┃ Repository                  ┃
┡━━━━━━━━━━━━━━━━━━━━━━━━━━━━╇━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┩
│ cat/cat                    │ nf-core/modules             │
└────────────────────────────┴─────────────────────────────┘
```

### 2.3. Add the import statement to your workflow

Open [core-hello/workflows/hello.nf](core-hello/workflows/hello.nf) and add the `include` statement for the `CAT_CAT` module in the imports section.

The nf-core convention is to use uppercase for module names when importing them.

=== "After"

    ```groovy title="core-hello/workflows/hello.nf" linenums="1" hl_lines="12"
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
    include { CAT_CAT                } from '../modules/nf-core/cat/cat/main'
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
    include { sayHello               } from '../modules/local/sayHello.nf'
    include { convertToUpper         } from '../modules/local/convertToUpper.nf'
    include { collectGreetings       } from '../modules/local/collectGreetings.nf'
    include { cowpy                  } from '../modules/local/cowpy.nf'
    ```

Note how the path for the nf-core module differs from the local modules:

- **nf-core module**: `'../modules/nf-core/cat/cat/main'` (includes the tool name twice and references `main.nf`)
- **Local module**: `'../modules/local/collectGreetings.nf'` (single file reference)

### Takeaway

You know how to install nf-core modules using the command-line tools and add the appropriate import statements to your workflow.

### What's next?

Learn how to use the module in your workflow.

---

## 3. Wire up the module to the workflow

Now we need to replace the call to `collectGreetings` with a call to `CAT_CAT`, adapting the inputs and outputs to match the module's interface.

### 3.1. Examine the cat/cat module interface

Let's look at the `cat/cat` module's main.nf file to understand its interface:

```bash
head -30 modules/nf-core/cat/cat/main.nf
```

The key parts of the module are:

```groovy title="modules/nf-core/cat/cat/main.nf (excerpt)" linenums="1"
process CAT_CAT {
    tag "$meta.id"
    label 'process_single'

    input:
    tuple val(meta), path(files_in)

    output:
    tuple val(meta), path("${prefix}"), emit: file_out
    path "versions.yml"           , emit: versions
```

The module expects:

- **Input**: A tuple containing metadata (`meta`) and input file(s) (`files_in`)
- **Output**: A tuple containing metadata and the concatenated output file, plus a versions file

### 3.2. Compare with collectGreetings interface

Our custom `collectGreetings` module has a simpler interface:

```groovy title="modules/local/collectGreetings.nf (excerpt)" linenums="1"
process collectGreetings {

    publishDir 'results', mode: 'copy'

    input:
        path input_files
        val batch_name

    output:
        path "COLLECTED-${batch_name}-output.txt" , emit: outfile
        val count_greetings , emit: count
```

The main differences are:

- `CAT_CAT` requires a metadata map, while `collectGreetings` doesn't
- `CAT_CAT` outputs a tuple, while `collectGreetings` outputs a simple path
- `CAT_CAT` requires a filename prefix via the `meta.id` field

### 3.3. Understanding metadata maps

You've just seen that `CAT_CAT` expects inputs and outputs structured as tuples with metadata:

```groovy
input:
tuple val(meta), path(files_in)

output:
tuple val(meta), path("${prefix}"), emit: file_out
```

This pattern is standard across all nf-core modules. The metadata map (commonly called `meta`) is a Groovy-style map containing information about a sample or dataset, with `id` being the required field used for naming outputs and tracking samples.

Why use metadata maps?

- **Sample tracking**: Keep sample information with data throughout the workflow
- **Standardization**: All nf-core modules follow this pattern
- **Flexibility**: Easy to add custom metadata fields
- **Output naming**: Consistent file naming based on sample IDs

!!! note "Learn more about metadata"

    For a comprehensive introduction to working with metadata in Nextflow workflows, including how to read metadata from samplesheets and use it to customize processing, see the [Metadata in workflows](../side_quests/metadata) side quest.

Rather than extracting the file from the tuple and losing the metadata when passing it to `cowpy`, let's update the `cowpy` module to accept and use metadata maps. This follows nf-core best practices and makes our module more consistent with community standards.

### 3.4. Update the cowpy module to use metadata

Let's start by updating the cowpy module to accept and propagate metadata, following the same pattern we saw in `CAT_CAT`.

Open [core-hello/modules/local/cowpy.nf](core-hello/modules/local/cowpy.nf) and modify it to accept metadata tuples:

=== "After"

    ```groovy title="core-hello/modules/local/cowpy.nf" linenums="1" hl_lines="11 12 15 16"
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

Key changes:

1. **Input**: Changed the first input from `path input_file` to `tuple val(meta), path(input_file)` to accept metadata
2. **Output**: Changed to emit a tuple with metadata: `tuple val(meta), path("cowpy-${input_file}"), emit: cowpy_output`
3. **Named emit**: Added `emit: cowpy_output` to give the output channel a descriptive name

Now our `cowpy` module follows the same metadata pattern as `CAT_CAT`, allowing metadata to flow through our entire workflow.

### 3.5. Adapt the workflow to use CAT_CAT and updated cowpy

Now we need to modify our workflow code to:

1. Create a metadata map with an appropriate ID
2. Combine the metadata with the collected files into a tuple
3. Call `CAT_CAT` instead of `collectGreetings`
4. Pass the tuple output from `CAT_CAT` directly to `cowpy` (metadata flows through!)

Open [core-hello/workflows/hello.nf](core-hello/workflows/hello.nf) and modify the workflow logic in the `main` block:

=== "After"

    ```groovy title="core-hello/workflows/hello.nf" linenums="26" hl_lines="7-15"
        // emit a greeting
        sayHello(ch_samplesheet)

        // convert the greeting to uppercase
        convertToUpper(sayHello.out)

        // collect all the greetings into one file using nf-core cat/cat module
        // create metadata map with batch name as the ID
        def meta = [ id: params.batch ]
        ch_for_cat = convertToUpper.out.collect().map { files -> tuple(meta, files) }

        CAT_CAT(ch_for_cat)

        // generate ASCII art of the greetings with cowpy
        cowpy(CAT_CAT.out.file_out, params.character)
    ```

=== "Before"

    ```groovy title="core-hello/workflows/hello.nf" linenums="26" hl_lines="7-14"
        // emit a greeting
        sayHello(ch_samplesheet)

        // convert the greeting to uppercase
        convertToUpper(sayHello.out)

        // collect all the greetings into one file
        collectGreetings(convertToUpper.out.collect(), params.batch)

        // emit a message about the size of the batch
        collectGreetings.out.count.view { "There were $it greetings in this batch" }

        // generate ASCII art of the greetings with cowpy
        cowpy(collectGreetings.out.outfile, params.character)
    ```

Let's break down what we changed:

1. **Created metadata**: `def meta = [ id: params.batch ]` creates a Groovy-style map with an `id` field set to our batch name
2. **Created a tuple channel**: `ch_for_cat = convertToUpper.out.collect().map { files -> tuple(meta, files) }` combines the metadata and collected files into the tuple format expected by `CAT_CAT`
3. **Called CAT_CAT**: Replaced `collectGreetings(...)` with `CAT_CAT(ch_for_cat)`
4. **Passed tuple to cowpy**: Since we updated `cowpy` to accept metadata tuples, we can now pass `CAT_CAT.out.file_out` directly. The metadata flows through!
5. **Removed count view**: The `cat/cat` module doesn't emit a count, so we removed that line

Notice how the metadata now flows seamlessly from `CAT_CAT` to `cowpy` without needing to extract and repackage it.

!!! note

    We removed the `collectGreetings.out.count.view { ... }` line because the nf-core `cat/cat` module doesn't provide a count of files. If you want to keep this functionality, you would need to count the files before calling `CAT_CAT`.

### 3.6. Update the emit block

Update the `emit` block to use the named emit channel from our updated cowpy module:

=== "After"

    ```groovy title="core-hello/workflows/hello.nf" linenums="52" hl_lines="2"
        emit:
        cowpy_hellos   = cowpy.out.cowpy_output
        versions       = ch_versions                 // channel: [ path(versions.yml) ]
    ```

=== "Before"

    ```groovy title="core-hello/workflows/hello.nf" linenums="52"
        emit:
        cowpy_hellos   = cowpy.out
        versions       = ch_versions                 // channel: [ path(versions.yml) ]
    ```

Since we added a named emit (`emit: cowpy_output`) to the cowpy module, we now reference it explicitly as `cowpy.out.cowpy_output`.

### 3.7. Test the workflow with metadata

Let's test that the workflow works with metadata tuples:

```bash
nextflow run core-hello --outdir core-hello-results -profile test,docker --validate_params false
```

```console title="Output"
 N E X T F L O W   ~  version 24.10.4

Launching `core-hello/main.nf` [curious_davinci] DSL2 - revision: c31b966b36

Input/output options
  input                     : core-hello/assets/greetings.csv
  outdir                    : core-hello-results

Institutional config options
  config_profile_name       : Test profile
  config_profile_description: Minimal test dataset to check pipeline function

Generic options
  validate_params           : false

Core Nextflow options
  runName                   : curious_davinci
  containerEngine           : docker
  profile                   : test,docker

!! Only displaying parameters that differ from the pipeline defaults !!
------------------------------------------------------
executor >  local (7)
[a1/2f8d9c] CORE_HELLO:HELLO:sayHello (1)       | 3 of 3 ✔
[e2/9a8b3d] CORE_HELLO:HELLO:convertToUpper (2) | 3 of 3 ✔
[c4/7e1b2a] CORE_HELLO:HELLO:CAT_CAT             | 1 of 1 ✔
[f5/3d9c8b] CORE_HELLO:HELLO:cowpy              | 1 of 1 ✔
-[core/hello] Pipeline completed successfully-
```

Notice that `CAT_CAT` now appears in the process execution list instead of `collectGreetings`.

### 3.8. Verify the outputs

Check that the outputs look correct:

```bash
ls results/
```

You should see the cowpy output files. The pipeline now successfully passes metadata through from `CAT_CAT` to `cowpy`.

At this point, we've successfully integrated metadata tuples into our workflow. But there's more we can do to follow nf-core conventions even more closely. Let's continue refining the module.

---

## 4. Simplify the module interface with ext.args

Now that we have metadata working, let's address another nf-core pattern: simplifying module interfaces by using `ext.args` for optional command-line arguments.

Currently, our `cowpy` module still requires the `character` parameter to be passed as a separate input. While this works, nf-core modules follow a convention of keeping interfaces minimal - only essential inputs (metadata and files) should be declared. Optional tool arguments are instead passed via configuration.

### 4.1. Understanding ext.args

The `task.ext.args` pattern is an nf-core convention for passing optional command-line arguments to tools. Instead of adding multiple input parameters for every possible tool option, nf-core modules accept optional arguments through the `ext.args` configuration directive.

Benefits of this approach:

- **Minimal interface**: The module only requires essential inputs (metadata and files)
- **Flexibility**: Users can specify any tool arguments via configuration
- **Consistency**: All nf-core modules follow this pattern
- **Portability**: Modules can be reused in other pipelines without expecting specific parameter names
- **No workflow changes**: Adding new tool options doesn't require updating workflow code

### 4.2. Update cowpy to use ext.args

Let's update the cowpy module to use `ext.args` instead of the `character` input parameter.

Open [core-hello/modules/local/cowpy.nf](core-hello/modules/local/cowpy.nf):

=== "After"

    ```groovy title="core-hello/modules/local/cowpy.nf" linenums="1" hl_lines="11 14 16 17"
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

    ```groovy title="core-hello/modules/local/cowpy.nf" linenums="1"
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
2. **Added ext.args**: The line `def args = task.ext.args ?: ''` uses the Elvis operator (`?:`) to provide an empty string as default if `task.ext.args` is not set
3. **Updated command**: Changed from hardcoded `-c "$character"` to using the configurable `$args`

The module interface is now simpler - it only accepts the essential metadata and file inputs.

### 4.3. Configure ext.args in modules.config

Now we need to configure the `ext.args` to pass the character option. This allows us to keep the module interface simple while still providing the character option at the pipeline level.

Open [core-hello/conf/modules.config](core-hello/conf/modules.config) and add the cowpy configuration:

=== "After"

    ```groovy title="core-hello/conf/modules.config" linenums="1" hl_lines="10-12"
    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        Config file for defining DSL2 per module options and publishing paths
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    process {
        publishDir = [
            path: { "${params.outdir}/${task.process.tokenize(':')[-1].tokenize('_')[0].toLowerCase()}" },
        ]

        withName: 'CORE_HELLO:HELLO:cowpy' {
            ext.args = { "-c ${params.character}" }
        }
    }
    ```

=== "Before"

    ```groovy title="core-hello/conf/modules.config" linenums="1"
    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        Config file for defining DSL2 per module options and publishing paths
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

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

### 4.4. Update the workflow call

Since the cowpy module no longer requires the `character` parameter as an input, we need to update the workflow call.

Open [core-hello/workflows/hello.nf](core-hello/workflows/hello.nf) and update the cowpy call:

=== "After"

    ```groovy title="core-hello/workflows/hello.nf" linenums="39" hl_lines="2"
        // generate ASCII art of the greetings with cowpy
        cowpy(CAT_CAT.out.file_out)
    ```

=== "Before"

    ```groovy title="core-hello/workflows/hello.nf" linenums="39"
        // generate ASCII art of the greetings with cowpy
        cowpy(CAT_CAT.out.file_out, params.character)
    ```

The workflow code is now cleaner - we don't need to pass `params.character` directly to the process. The module interface is kept minimal, making it more portable, while the pipeline still provides the explicit option through configuration.

### 4.5. Test the workflow with ext.args

Test that the workflow still works with the ext.args configuration:

```bash
nextflow run core-hello --outdir core-hello-results -profile test,docker --validate_params false
```

The pipeline should run successfully, producing the same cowpy output as before.

### Takeaway

You now understand how to use `ext.args` to keep module interfaces minimal while maintaining flexibility through configuration.

---

## 5. Add configurable output naming with ext.prefix

There's one more nf-core pattern we can apply: using `ext.prefix` for configurable output file naming.

### 5.1. Understanding ext.prefix

The `task.ext.prefix` pattern is another nf-core convention for standardizing output file naming across modules while keeping it configurable.

Benefits:

- **Standardized naming**: Output files are typically named using sample IDs from metadata
- **Configurable**: Users can override the default naming if needed
- **Consistent**: All nf-core modules follow this pattern
- **Predictable**: Easy to know what output files will be called

### 5.2. Update cowpy to use ext.prefix

Let's update the cowpy module to use `ext.prefix` for output file naming.

Open [core-hello/modules/local/cowpy.nf](core-hello/modules/local/cowpy.nf):

=== "After"

    ```groovy title="core-hello/modules/local/cowpy.nf" linenums="1" hl_lines="15 17 18 20"
    #!/usr/bin/env nextflow

    // Generate ASCII art with cowpy (https://github.com/jeffbuttars/cowpy)
    process cowpy {

        publishDir 'results', mode: 'copy'

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

    ```groovy title="core-hello/modules/local/cowpy.nf" linenums="1"
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

Key changes:

1. **Added ext.prefix**: `prefix = task.ext.prefix ?: "${meta.id}"` provides a configurable prefix with a sensible default (the sample ID)
2. **Updated output**: Changed from hardcoded `cowpy-${input_file}` to `${prefix}.txt`
3. **Updated command**: Uses the configured prefix for the output filename

### 5.3. Configure ext.prefix to maintain original naming

To maintain the same output file naming as before (`cowpy-<id>.txt`), we can configure `ext.prefix` in modules.config.

Update [core-hello/conf/modules.config](core-hello/conf/modules.config):

=== "After"

    ```groovy title="core-hello/conf/modules.config" linenums="11" hl_lines="3"
        withName: 'CORE_HELLO:HELLO:cowpy' {
            ext.args = { "-c ${params.character}" }
            ext.prefix = { "cowpy-${meta.id}" }
        }
    ```

=== "Before"

    ```groovy title="core-hello/conf/modules.config" linenums="11"
        withName: 'CORE_HELLO:HELLO:cowpy' {
            ext.args = { "-c ${params.character}" }
        }
    ```

Note that we use a closure (`{ "cowpy-${meta.id}" }`) which has access to `meta` because it's evaluated in the context of the process execution.

!!! note

    The `ext.prefix` closure has access to `meta` because the configuration is evaluated in the context of the process execution, where metadata is available.

### 5.4. Test and verify

Test the workflow once more:

```bash
nextflow run core-hello --outdir core-hello-results -profile test,docker --validate_params false
```

Check the outputs:

```bash
ls results/
```

You should see the cowpy output files with the same naming as before: `cowpy-test.txt` (based on the batch name). This demonstrates how `ext.prefix` allows you to maintain your preferred naming convention while keeping the module interface flexible.

If you wanted to change the naming (for example, to just `test.txt`), you would only need to modify the `ext.prefix` configuration - no changes to the module or workflow code would be required.

### Takeaway

You now know how to:

- Adapt your workflow to use an nf-core module
- Understand and work with metadata maps
- Update local modules to follow nf-core patterns with metadata tuples
- Create metadata structures and pass them through your workflow
- Use the `ext.args` pattern to keep module interfaces simple and portable
- Use the `ext.prefix` pattern for configurable, standardized output file naming
- Configure process-specific parameters through `modules.config`
- Use named emit channels for clearer output references

### What's next?

Clean up by optionally removing the now-unused local module.

---

## 6. Optional: Clean up unused local modules

Now that we're using the nf-core `cat/cat` module, the local `collectGreetings` module is no longer needed.

### 6.1. Remove the collectGreetings import

Remove or comment out the import line for `collectGreetings`:

```groovy title="core-hello/workflows/hello.nf" linenums="10"
include { sayHello               } from '../modules/local/sayHello.nf'
include { convertToUpper         } from '../modules/local/convertToUpper.nf'
// include { collectGreetings    } from '../modules/local/collectGreetings.nf'  // No longer needed
include { cowpy                  } from '../modules/local/cowpy.nf'
include { CAT_CAT                } from '../modules/nf-core/cat/cat/main'
```

### 6.2. Optionally remove the module file

You can optionally delete the `collectGreetings.nf` file:

```bash
rm modules/local/collectGreetings.nf
```

However, you might want to keep it as a reference for understanding the differences between local and nf-core modules.

### Takeaway

You know how to replace custom local modules with nf-core modules and clean up unused code.

### What's next?

Continue to [Part 4: Input validation](./04_input_validation.md) to learn how to add schema-based input validation to your pipeline, or explore other nf-core modules you might add to enhance your pipeline further.
