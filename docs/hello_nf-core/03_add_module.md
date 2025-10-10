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

Navigate to the modules page in your web browser and use the search bar to search for "cat".

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

From your `core-hello` directory, run the following command:

```bash
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

### 3.3. Adapt the workflow to use CAT_CAT

We need to modify our workflow code to:

1. Create a metadata map with an appropriate ID
2. Combine the metadata with the collected files into a tuple
3. Call `CAT_CAT` instead of `collectGreetings`
4. Adapt downstream processes to handle the tuple output

Open [core-hello/workflows/hello.nf](core-hello/workflows/hello.nf) and modify the workflow logic in the `main` block:

=== "After"

    ```groovy title="core-hello/workflows/hello.nf" linenums="26" hl_lines="11-14"
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
        cowpy(CAT_CAT.out.file_out.map{ meta, file -> file }, params.character)
    ```

=== "Before"

    ```groovy title="core-hello/workflows/hello.nf" linenums="26"
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

1. **Created metadata**: `def meta = [ id: params.batch ]` creates a map with an ID field set to our batch name
2. **Created a tuple channel**: `ch_for_cat = convertToUpper.out.collect().map { files -> tuple(meta, files) }` combines the metadata and collected files into the tuple format expected by `CAT_CAT`
3. **Called CAT_CAT**: Replaced `collectGreetings(...)` with `CAT_CAT(ch_for_cat)`
4. **Extracted file from tuple**: Modified the cowpy call to extract just the file from the output tuple using `.map{ meta, file -> file }`
5. **Removed count view**: The `cat/cat` module doesn't emit a count, so we removed that line

!!! note

    We removed the `collectGreetings.out.count.view { ... }` line because the nf-core `cat/cat` module doesn't provide a count of files. If you want to keep this functionality, you would need to count the files before calling `CAT_CAT`.

### 3.4. Update the emit block

Update the `emit` block to reflect the new output:

=== "After"

    ```groovy title="core-hello/workflows/hello.nf" linenums="52" hl_lines="2"
        emit:
        cowpy_hellos   = cowpy.out
        versions       = ch_versions                 // channel: [ path(versions.yml) ]
    ```

=== "Before"

    ```groovy title="core-hello/workflows/hello.nf" linenums="52"
        emit:
        cowpy_hellos   = cowpy.out
        versions       = ch_versions                 // channel: [ path(versions.yml) ]
    ```

In this case, the emit block doesn't need to change because we're still emitting the cowpy output.

### 3.5. Test the updated workflow

Let's test that our workflow still works with the nf-core module:

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

### 3.6. Verify the outputs

Check that the outputs look correct:

```bash
ls results/
```

You should still see the concatenated and cowpy output files, though the naming may be slightly different since `CAT_CAT` uses the metadata ID for the output filename.

### Takeaway

You know how to adapt your workflow to use an nf-core module, including creating the appropriate metadata structures and handling tuple-based inputs and outputs.

### What's next?

Clean up by optionally removing the now-unused local module.

---

## 4. Optional: Clean up unused local modules

Now that we're using the nf-core `cat/cat` module, the local `collectGreetings` module is no longer needed.

### 4.1. Remove the collectGreetings import

Remove or comment out the import line for `collectGreetings`:

```groovy title="core-hello/workflows/hello.nf" linenums="10"
include { sayHello               } from '../modules/local/sayHello.nf'
include { convertToUpper         } from '../modules/local/convertToUpper.nf'
// include { collectGreetings    } from '../modules/local/collectGreetings.nf'  // No longer needed
include { cowpy                  } from '../modules/local/cowpy.nf'
include { CAT_CAT                } from '../modules/nf-core/cat/cat/main'
```

### 4.2. Optionally remove the module file

You can optionally delete the `collectGreetings.nf` file:

```bash
rm modules/local/collectGreetings.nf
```

However, you might want to keep it as a reference for understanding the differences between local and nf-core modules.

### Takeaway

You know how to replace custom local modules with nf-core modules and clean up unused code.

### What's next?

Continue to [Part 4: Input validation](./04_input_validation.md) to learn how to add schema-based input validation to your pipeline, or explore other nf-core modules you might add to enhance your pipeline further.
