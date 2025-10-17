# Part 3: Add an existing nf-core module

In this third part of the Hello nf-core training course, we show you how to add an existing nf-core module to your pipeline and adapt your local modules to follow nf-core conventions.

One of the great advantages of nf-core pipelines is the ability to leverage pre-built, tested modules from the nf-core/modules repository. Rather than writing every process from scratch, you can install and use community-maintained modules that follow best practices.

In this section, we'll replace the custom `collectGreetings` module with the `cat/cat` module from nf-core/modules, then progressively adapt our `cowpy` module to follow nf-core patterns.

!!! note

    This section assumes you have completed [Part 2: Rewrite Hello for nf-core](./02_rewrite_hello.md) and have a working `core-hello` pipeline.

---

## 1. Use an nf-core module

First, let's learn how to find, install, and use an existing nf-core module in our pipeline.

The `collectGreetings` process in our pipeline uses the Unix `cat` command to concatenate multiple greeting files into one. This is a perfect use case for the nf-core `cat/cat` module, which is designed specifically for concatenating files.

!!! note "Module naming convention"

    nf-core modules follow the naming convention `software/command`. The `cat/cat` module wraps the `cat` command from the `cat` software package. Other examples include `fastqc/fastqc` (FastQC software, fastqc command) or `samtools/view` (samtools software, view command).

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

### 1.4. Install and verify the cat/cat module

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

### 1.5. Add the import statement to your workflow

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

### 1.6. Examine the cat/cat module interface

Let's look at the `cat/cat` module's main.nf file to understand its interface:

```bash
head -30 modules/nf-core/cat/cat/main.nf
```

The key parts of the module are:

```groovy title="modules/nf-core/cat/cat/main.nf (excerpt)" linenums="1" hl_lines="6 9"
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

### 1.7. Compare with collectGreetings interface

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

### 1.8. Understanding metadata maps

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

For now, we'll pass the output from `CAT_CAT` to `cowpy` with the character parameter. In the next section, we'll adapt `cowpy` to follow nf-core conventions.

### 1.9. Wire up CAT_CAT in the workflow

Now we need to modify our workflow code to use `CAT_CAT` instead of `collectGreetings`. Since `CAT_CAT` requires metadata tuples, we need to create a metadata map and combine it with our files.

Open [core-hello/workflows/hello.nf](core-hello/workflows/hello.nf) and modify the workflow logic in the `main` block:

=== "After"

    ```groovy title="core-hello/workflows/hello.nf" linenums="26" hl_lines="7-17"
        // emit a greeting
        sayHello(ch_samplesheet)

        // convert the greeting to uppercase
        convertToUpper(sayHello.out)

        // collect all the greetings into one file using nf-core cat/cat module
        // create metadata map with batch name as the ID
        def cat_meta = [ id: params.batch ]
        ch_for_cat = convertToUpper.out.collect().map { files -> tuple(cat_meta, files) }

        CAT_CAT(ch_for_cat)

        // generate ASCII art of the greetings with cowpy
        // Extract the file from the tuple since cowpy doesn't use metadata yet
        ch_for_cowpy = CAT_CAT.out.file_out.map{ meta, file -> file }
        cowpy(ch_for_cowpy, params.character)
    ```

=== "Before"

    ```groovy title="core-hello/workflows/hello.nf" linenums="26" hl_lines="7-15"
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
4. **Extracted file for cowpy**: Since `cowpy` doesn't accept metadata tuples yet, we extract just the file from the tuple using `.map{ meta, file -> file }` and store it in `ch_for_cowpy`
5. **Called cowpy**: Pass the extracted file channel to `cowpy` along with the character parameter
6. **Removed count view**: The `cat/cat` module doesn't emit a count, so we removed that line

!!! note

    We removed the `collectGreetings.out.count.view { ... }` line because the nf-core `cat/cat` module doesn't provide a count of files. If you want to keep this functionality, you would need to count the files before calling `CAT_CAT`.

    Notice that we're extracting just the file from `CAT_CAT`'s output tuple to pass to `cowpy`. In the next section, we'll update `cowpy` to work with metadata tuples directly.

### 1.10. Test the workflow

Let's test that the workflow works with the `cat/cat` module:

```bash
nextflow run . --outdir core-hello-results -profile test,docker --validate_params false
```

```console title="Output"
 N E X T F L O W   ~  version 25.04.3

Launching `./main.nf` [extravagant_volhard] DSL2 - revision: 6aa79210e6

Input/output options
  input                     : /workspaces/training/hello-nf-core/nf-core-hello/assets/greetings.csv
  outdir                    : core-hello-results

Institutional config options
  config_profile_name       : Test profile
  config_profile_description: Minimal test dataset to check pipeline function

Generic options
  validate_params           : false
  trace_report_suffix       : 2025-10-17_19-51-31

Core Nextflow options
  runName                   : extravagant_volhard
  containerEngine           : docker
  launchDir                 : /workspaces/training/hello-nf-core/nf-core-hello
  workDir                   : /workspaces/training/hello-nf-core/nf-core-hello/work
  projectDir                : /workspaces/training/hello-nf-core/nf-core-hello
  userName                  : root
  profile                   : test,docker
  configFiles               : /workspaces/training/hello-nf-core/nf-core-hello/nextflow.config

!! Only displaying parameters that differ from the pipeline defaults !!
------------------------------------------------------
executor >  local (8)
[60/3ac109] NFCORE_HELLO:HELLO:sayHello (3)       [100%] 3 of 3 ✔
[58/073077] NFCORE_HELLO:HELLO:convertToUpper (3) [100%] 3 of 3 ✔
[00/4f3d32] NFCORE_HELLO:HELLO:CAT_CAT (test)     [100%] 1 of 1 ✔
[98/afab8b] NFCORE_HELLO:HELLO:cowpy              [100%] 1 of 1 ✔
-[nf-core/hello] Pipeline completed successfully-
```

Notice that `CAT_CAT` now appears in the process execution list instead of `collectGreetings`.

### Takeaway

You now know how to:

- Find and install nf-core modules
- Understand metadata maps and why nf-core modules use them
- Create metadata structures to pass to nf-core modules
- Wire up nf-core modules in your workflow

### What's next?

Adapt your local modules to follow nf-core conventions.

---

## 2. Adapt local modules to nf-core conventions

Now that we've successfully integrated the nf-core `CAT_CAT` module, let's adapt our local `cowpy` module to follow nf-core conventions. We'll do this incrementally, introducing one pattern at a time:

1. First, we'll update `cowpy` to accept and propagate metadata tuples
2. Then, we'll simplify its interface using `ext.args`
3. Finally, we'll add configurable output naming with `ext.prefix`

### 2.1. Update cowpy to use metadata tuples

Currently, we're extracting the file from `CAT_CAT`'s output tuple to pass to `cowpy`. It would be better to have `cowpy` accept metadata tuples directly, allowing metadata to flow through the entire workflow.

Open [core-hello/modules/local/cowpy.nf](core-hello/modules/local/cowpy.nf) and modify it to accept metadata tuples:

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

Now update the workflow to pass the tuple directly instead of extracting the file. Open [core-hello/workflows/hello.nf](core-hello/workflows/hello.nf):

=== "After"

    ```groovy title="core-hello/workflows/hello.nf" linenums="39" hl_lines="2"
        // generate ASCII art of the greetings with cowpy
        cowpy(CAT_CAT.out.file_out, params.character)
    ```

=== "Before"

    ```groovy title="core-hello/workflows/hello.nf" linenums="39" hl_lines="2-4"
        // generate ASCII art of the greetings with cowpy
        // Extract the file from the tuple since cowpy doesn't use metadata yet
        ch_for_cowpy = CAT_CAT.out.file_out.map{ meta, file -> file }
        cowpy(ch_for_cowpy, params.character)
    ```

Also update the emit block to use the named emit:

=== "After"

    ```groovy title="core-hello/workflows/hello.nf" linenums="58" hl_lines="2"
        emit:
        cowpy_hellos   = cowpy.out.cowpy_output
        versions       = ch_versions                 // channel: [ path(versions.yml) ]
    ```

=== "Before"

    ```groovy title="core-hello/workflows/hello.nf" linenums="58" hl_lines="2"
        emit:
        cowpy_hellos   = cowpy.out
        versions       = ch_versions                 // channel: [ path(versions.yml) ]
    ```

Test the workflow to ensure metadata flows through correctly:

```bash
nextflow run . --outdir core-hello-results -profile test,docker --validate_params false
```

The pipeline should run successfully with metadata now flowing from `CAT_CAT` through `cowpy`.

### 2.2. Simplify the interface with ext.args

Now let's address another nf-core pattern: simplifying module interfaces by using `ext.args` for optional command-line arguments.

Currently, our `cowpy` module requires the `character` parameter to be passed as a separate input. While this works, nf-core modules follow a convention of keeping interfaces minimal - only essential inputs (metadata and files) should be declared. Optional tool arguments are instead passed via configuration.

#### Understanding ext.args

The `task.ext.args` pattern is an nf-core convention for passing optional command-line arguments to tools. Instead of adding multiple input parameters for every possible tool option, nf-core modules accept optional arguments through the `ext.args` configuration directive.

Benefits of this approach:

- **Minimal interface**: The module only requires essential inputs (metadata and files)
- **Flexibility**: Users can specify any tool arguments via configuration
- **Consistency**: All nf-core modules follow this pattern
- **Portability**: Modules can be reused in other pipelines without expecting specific parameter names
- **No workflow changes**: Adding new tool options doesn't require updating workflow code

#### Update the module

Let's update the cowpy module to use `ext.args` instead of the `character` input parameter.

Open [core-hello/modules/local/cowpy.nf](core-hello/modules/local/cowpy.nf):

=== "After"

    ```groovy title="core-hello/modules/local/cowpy.nf" linenums="1" hl_lines="18 20"
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

    ```groovy title="core-hello/modules/local/cowpy.nf" linenums="1" hl_lines="13"
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

#### Configure ext.args

Now we need to configure the `ext.args` to pass the character option. This allows us to keep the module interface simple while still providing the character option at the pipeline level.

Open [core-hello/conf/modules.config](core-hello/conf/modules.config) and add the cowpy configuration:

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

Open [core-hello/workflows/hello.nf](core-hello/workflows/hello.nf) and update the cowpy call:

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

Test that the workflow still works with the ext.args configuration. Let's specify a different character to verify the configuration is working:

```bash
nextflow run . --outdir core-hello-results -profile test,docker --validate_params false --character cow
```

The pipeline should run successfully. In the output, look for the cowpy process execution line which will show something like:

```console title="Output (excerpt)"
[f3/abc123] process > CORE_HELLO:HELLO:cowpy [100%] 1 of 1 ✔
```

Now let's verify that the `ext.args` configuration actually passed the character argument to the cowpy command. Use the task hash (the `f3/abc123` part) to inspect the `.command.sh` file in the work directory:

```bash
cat work/f3/abc123*/command.sh
```

You should see the cowpy command with the `-c cow` argument:

```console title="Output"
#!/usr/bin/env bash
...
cat test.txt | cowpy -c cow > cowpy-test.txt
```

This confirms that `task.ext.args` successfully passed the character parameter through the configuration rather than requiring it as a process input.

### 2.3. Add configurable output naming with ext.prefix

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

Open [core-hello/modules/local/cowpy.nf](core-hello/modules/local/cowpy.nf):

=== "After"

    ```groovy title="core-hello/modules/local/cowpy.nf" linenums="1" hl_lines="15 19 21"
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

    ```groovy title="core-hello/modules/local/cowpy.nf" linenums="1" hl_lines="15 20"
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

#### Configure ext.prefix

To maintain the same output file naming as before (`cowpy-<id>.txt`), we can configure `ext.prefix` in modules.config.

Update [core-hello/conf/modules.config](core-hello/conf/modules.config):

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
nextflow run . --outdir core-hello-results -profile test,docker --validate_params false
```

Check the outputs:

```bash
ls results/
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

### 2.4. Optional: Clean up unused local modules

Now that we're using the nf-core `cat/cat` module, the local `collectGreetings` module is no longer needed.

Remove or comment out the import line for `collectGreetings` in [core-hello/workflows/hello.nf](core-hello/workflows/hello.nf):

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

## 3. Contributing modules back to nf-core

Now that you've learned how to create modules following nf-core conventions, you might develop modules that could benefit the wider community. The [nf-core/modules](https://github.com/nf-core/modules) repository welcomes contributions of well-tested, standardized modules.

### Why contribute?

Contributing your modules to nf-core:

- Makes your tools available to the entire nf-core community through the modules catalog at [nf-co.re/modules](https://nf-co.re/modules)
- Ensures ongoing community maintenance and improvements
- Provides quality assurance through code review and automated testing
- Gives your work visibility and recognition

### Getting started

Before contributing a new module:

1. Check if it already exists at [nf-co.re/modules](https://nf-co.re/modules) or in [open PRs](https://github.com/nf-core/modules/pulls)
2. Create an issue on [nf-core/modules](https://github.com/nf-core/modules/issues) to notify the community of your plans
3. Use `nf-core modules create <tool>/<subtool>` in the nf-core/modules repository to generate the module structure
4. Follow the patterns you've learned: metadata tuples, `ext.args`, `ext.prefix`, and comprehensive testing
5. Submit a pull request and request review from `@nf-core/modules-team`

### Resources

- **Comprehensive guide**: [nf-core components tutorial](https://nf-co.re/docs/tutorials/nf-core_components/components)
- **Module specifications**: [Module guidelines](https://nf-co.re/docs/guidelines/components/modules)
- **Community support**: Join the `#modules` channel on [nf-core Slack](https://nf-co.re/join)

Contributing to nf-core is a rewarding way to give back to the community while ensuring your tools follow best practices and reach researchers worldwide.

---

## Takeaway

You now know how to leverage pre-built nf-core modules in your pipeline and adapt your local modules to follow nf-core conventions. You learned how to use metadata tuples to track sample information through the workflow, simplify module interfaces with `ext.args` for configurable arguments, and use `ext.prefix` for standardized output file naming. These patterns keep modules portable and reusable while centralizing configuration in `modules.config`, making your pipeline more maintainable and consistent with nf-core best practices.

Finally, you learned how to contribute your modules back to the nf-core community, making them available to researchers worldwide and ensuring they benefit from ongoing community maintenance.

## What's next?

Continue to [Part 4: Input validation](./04_input_validation.md) to learn how to add schema-based input validation to your pipeline, or explore other nf-core modules you might add to enhance your pipeline further.
