# Part 3: Use an nf-core module

In this third part of the Hello nf-core training course, we show you how to find, install, and use an existing nf-core module in your pipeline.

One of the great advantages of nf-core pipelines is the ability to leverage pre-built, tested modules from the [nf-core/modules](https://github.com/nf-core/modules) repository. Rather than writing every process from scratch, you can install and use community-maintained modules that follow best practices. You can browse available modules at [nf-co.re/modules](https://nf-co.re/modules).

In this section, we'll replace the custom `collectGreetings` module with the `cat/cat` module from nf-core/modules.

!!! note

    This section assumes you have completed [Part 2: Rewrite Hello for nf-core](./02_rewrite_hello.md) and have a working `core-hello` pipeline.

    If you didn't complete Part 2 or want to start fresh for this section, you can use the `core-hello-part2` solution as your starting point:

    ```bash
    cp -r hello-nf-core/solutions/core-hello-part2 core-hello
    cd core-hello
    ```

    This gives you a fully functional nf-core pipeline ready for adding modules.

---

## 1. Use an nf-core module

First, let's learn how to find, install, and use an existing nf-core module in our pipeline.

The `collectGreetings` process in our pipeline uses the Unix `cat` command to concatenate multiple greeting files into one. This is a perfect use case for the nf-core `cat/cat` module, which is designed specifically for concatenating files.

!!! note "Module naming convention"

    nf-core modules follow the naming convention `software/command` when a tool provides multiple commands, like `samtools/view` (samtools package, view command) or `gatk/haplotypecaller` (GATK package, HaplotypeCaller command). For tools that provide only one main command, modules use a single level like `fastqc` or `multiqc`. The `cat/cat` naming reflects the organizational structure in the modules repository.

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
nf-core modules list remote | grep 'cat/cat'
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

```groovy title="modules/local/collectGreetings.nf (excerpt)" linenums="1" hl_lines="6-7 10"
process collectGreetings {

    publishDir 'results', mode: 'copy'

    input:
        path input_files
        val batch_name

    output:
        path "COLLECTED-${batch_name}-output.txt" , emit: outfile
```

The main differences are:

- `CAT_CAT` requires a metadata map (`tuple val(meta), path(files_in)`), while `collectGreetings` takes separate `path` and `val` inputs
- `CAT_CAT` outputs a tuple with metadata, while `collectGreetings` outputs a simple path
- `CAT_CAT` uses `meta.id` for the filename prefix, while `collectGreetings` uses the `batch_name` parameter

### 1.8. Understanding metadata maps

You've just seen that `CAT_CAT` expects inputs and outputs structured as tuples with metadata:

```groovy
input:
tuple val(meta), path(files_in)

output:
tuple val(meta), path("${prefix}"), emit: file_out
```

This pattern is standard across all nf-core modules. The metadata map (commonly called `meta`) is a Groovy-style map containing information about a sample or dataset, with `id` being the required field used for naming outputs and tracking samples.

For example, a typical metadata map might look like:

```groovy
[id: 'sample1', single_end: false, strandedness: 'forward']
```

In this tutorial, we use a simple metadata map with just the batch name:

```groovy
[id: 'test']
```

Why use metadata maps?

- **Sample tracking**: Keep sample information with data throughout the workflow
- **Standardization**: All nf-core modules follow this pattern
- **Flexibility**: Easy to add custom metadata fields
- **Output naming**: Consistent file naming based on sample IDs

!!! note "Learn more about metadata"

    For a comprehensive introduction to working with metadata in Nextflow workflows, including how to read metadata from samplesheets and use it to customize processing, see the [Metadata in workflows](../side_quests/metadata) side quest.

For now, we'll pass the output from `CAT_CAT` to `cowpy` with the character parameter. In the next section, we'll adapt `cowpy` to follow nf-core conventions.

### 1.9. Wire up CAT_CAT in the workflow

Now we need to modify our workflow code to use `CAT_CAT` instead of `collectGreetings`. Since `CAT_CAT` requires metadata tuples, we'll do this in several steps to make it clear how to work with metadata.

Open [core-hello/workflows/hello.nf](core-hello/workflows/hello.nf) and make the following changes to the workflow logic in the `main` block.

#### Step 1: Create a metadata map

First, we need to create a metadata map for `CAT_CAT`. Remember that nf-core modules require metadata with at least an `id` field.

Add these lines after the `convertToUpper` call, removing the `collectGreetings` call:

=== "After"

    ```groovy title="core-hello/workflows/hello.nf" linenums="26" hl_lines="7-8"
        // emit a greeting
        sayHello(ch_samplesheet)

        // convert the greeting to uppercase
        convertToUpper(sayHello.out)

        // create metadata map with batch name as the ID
        def cat_meta = [ id: params.batch ]

        // generate ASCII art of the greetings with cowpy
        cowpy(collectGreetings.out.outfile, params.character)
    ```

=== "Before"

    ```groovy title="core-hello/workflows/hello.nf" linenums="26" hl_lines="7-8"
        // emit a greeting
        sayHello(ch_samplesheet)

        // convert the greeting to uppercase
        convertToUpper(sayHello.out)

        // collect all the greetings into one file
        collectGreetings(convertToUpper.out.collect(), params.batch)

        // generate ASCII art of the greetings with cowpy
        cowpy(collectGreetings.out.outfile, params.character)
    ```

This creates a simple metadata map where the `id` is set to our batch name (which will be "test" when using the test profile).

#### Step 2: Create a channel with metadata tuples

Next, transform the channel of files into a channel of tuples containing metadata and files:

=== "After"

    ```groovy title="core-hello/workflows/hello.nf" linenums="26" hl_lines="8-10"
        // emit a greeting
        sayHello(ch_samplesheet)

        // convert the greeting to uppercase
        convertToUpper(sayHello.out)

        // create metadata map with batch name as the ID
        def cat_meta = [ id: params.batch ]
        // create a channel with metadata and files in tuple format
        ch_for_cat = convertToUpper.out.collect().map { files -> tuple(cat_meta, files) }

        // generate ASCII art of the greetings with cowpy
        cowpy(collectGreetings.out.outfile, params.character)
    ```

=== "Before"

    ```groovy title="core-hello/workflows/hello.nf" linenums="26" hl_lines="8"
        // emit a greeting
        sayHello(ch_samplesheet)

        // convert the greeting to uppercase
        convertToUpper(sayHello.out)

        // create metadata map with batch name as the ID
        def cat_meta = [ id: params.batch ]

        // generate ASCII art of the greetings with cowpy
        cowpy(collectGreetings.out.outfile, params.character)
    ```

This line does two things:

- `.collect()` gathers all files from the `convertToUpper` output into a single list
- `.map { files -> tuple(cat_meta, files) }` creates a tuple of `[metadata, files]` in the format `CAT_CAT` expects

#### Step 3: Call CAT_CAT

Now call `CAT_CAT` with the properly formatted channel:

=== "After"

    ```groovy title="core-hello/workflows/hello.nf" linenums="26" hl_lines="11-12"
        // emit a greeting
        sayHello(ch_samplesheet)

        // convert the greeting to uppercase
        convertToUpper(sayHello.out)

        // create metadata map with batch name as the ID
        def cat_meta = [ id: params.batch ]
        ch_for_cat = convertToUpper.out.collect().map { files -> tuple(cat_meta, files) }

        // concatenate files using the nf-core cat/cat module
        CAT_CAT(ch_for_cat)

        // generate ASCII art of the greetings with cowpy
        cowpy(collectGreetings.out.outfile, params.character)
    ```

=== "Before"

    ```groovy title="core-hello/workflows/hello.nf" linenums="26" hl_lines="10"
        // emit a greeting
        sayHello(ch_samplesheet)

        // convert the greeting to uppercase
        convertToUpper(sayHello.out)

        // create metadata map with batch name as the ID
        def cat_meta = [ id: params.batch ]
        ch_for_cat = convertToUpper.out.collect().map { files -> tuple(cat_meta, files) }

        // generate ASCII art of the greetings with cowpy
        cowpy(collectGreetings.out.outfile, params.character)
    ```

#### Step 4: Update cowpy to use CAT_CAT output

Finally, update the `cowpy` call to use the output from `CAT_CAT`. Since `cowpy` doesn't accept metadata tuples yet (we'll fix this in the next section), we need to extract just the file:

=== "After"

    ```groovy title="core-hello/workflows/hello.nf" linenums="26" hl_lines="13-17"
        // emit a greeting
        sayHello(ch_samplesheet)

        // convert the greeting to uppercase
        convertToUpper(sayHello.out)

        // create metadata map with batch name as the ID
        def cat_meta = [ id: params.batch ]
        ch_for_cat = convertToUpper.out.collect().map { files -> tuple(cat_meta, files) }

        CAT_CAT(ch_for_cat)

        // generate ASCII art of the greetings with cowpy
        // extract the file from the tuple since cowpy doesn't use metadata yet
        ch_for_cowpy = CAT_CAT.out.file_out.map{ meta, file -> file }
        cowpy(ch_for_cowpy, params.character)
    ```

=== "Before"

    ```groovy title="core-hello/workflows/hello.nf" linenums="26" hl_lines="13-14"
        // emit a greeting
        sayHello(ch_samplesheet)

        // convert the greeting to uppercase
        convertToUpper(sayHello.out)

        // create metadata map with batch name as the ID
        def cat_meta = [ id: params.batch ]
        ch_for_cat = convertToUpper.out.collect().map { files -> tuple(cat_meta, files) }

        CAT_CAT(ch_for_cat)

        // generate ASCII art of the greetings with cowpy
        cowpy(collectGreetings.out.outfile, params.character)
    ```

The `.map{ meta, file -> file }` operation extracts just the file from the `[metadata, file]` tuple that `CAT_CAT` outputs.

!!! note

    We're extracting just the file from `CAT_CAT`'s output tuple to pass to `cowpy`. In the next section, we'll update `cowpy` to work with metadata tuples directly, so this extraction step won't be necessary.

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

In [Part 4](./04_adapt_module.md), we'll adapt your local `cowpy` module to follow nf-core conventions.
