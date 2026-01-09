# Workflow output definitions

Nextflow provides a powerful new syntax for managing workflow outputs that centralizes output configuration and enables automatic generation of output manifests.

### Learning goals

In this side quest, you'll learn how to use the workflow output definition syntax, which provides a cleaner alternative to the traditional `publishDir` directive.

By the end of this side quest, you'll be able to:

- Understand the limitations of the traditional `publishDir` approach
- Use the `publish:` section in workflows to declare outputs
- Configure the `output {}` block to organize published files
- Use dynamic paths based on metadata to organize outputs
- Generate index files (manifests) documenting your outputs

These skills will help you build workflows with cleaner output management and better documentation of results.

### Prerequisites

Before taking on this side quest, you should:

- Have completed the [Hello Nextflow](../../hello_nextflow/) tutorial or equivalent beginner's course.
- Be comfortable using basic Nextflow concepts (processes, channels, operators)

---

## 0. Get started

#### Open the training codespace

If you haven't yet done so, make sure to open the training environment as described in the [Environment Setup](../envsetup/index.md).

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

#### Move into the project directory

Let's move into the directory where the files for this tutorial are located.

```bash
cd side-quests/workflow_outputs
```

You can set VSCode to focus on this directory:

```bash
code .
```

#### Review the materials

You'll find a simple workflow file called `main.nf`, a `modules` directory containing a module file with two processes, and a `greetings.csv` file containing sample data.

```console title="Directory contents"
.
├── greetings.csv
├── main.nf
├── modules
│   └── greetings.nf
└── nextflow.config
```

This directory contains a simple greeting pipeline similar to what you built in Hello Nextflow.
The CSV file contains greetings in different languages that we'll process and organize by language.

#### Review the assignment

Your challenge is to refactor this workflow to use the new workflow output definition syntax instead of `publishDir`, organizing outputs by language and generating index files that document what was produced.

#### Readiness checklist

Think you're ready to dive in?

- [ ] I understand the goal of this course and its prerequisites
- [ ] My codespace is up and running
- [ ] I've set my working directory appropriately
- [ ] I understand the assignment

If you can check all the boxes, you're good to go.

---

## 1. The traditional approach: publishDir

### 1.1. Review the current workflow

Let's start by examining how the current workflow uses `publishDir` to manage outputs.

Take a look at the modules file:

```groovy title="modules/greetings.nf" linenums="1"
#!/usr/bin/env nextflow

/*
 * Create a greeting file from input text
 */
process SAY_HELLO {

    publishDir 'results/greetings', mode: 'copy'

    input:
        tuple val(greeting), val(language)

    output:
        tuple val(greeting), val(language), path("${greeting}-output.txt")

    script:
    """
    echo '$greeting' > '${greeting}-output.txt'
    """
}

/*
 * Convert greeting to uppercase
 */
process CONVERT_TO_UPPER {

    publishDir 'results/uppercase', mode: 'copy'

    input:
        tuple val(greeting), val(language), path(input_file)

    output:
        tuple val(greeting), val(language), path("UPPER-${input_file}")

    script:
    """
    cat '$input_file' | tr '[a-z]' '[A-Z]' > 'UPPER-${input_file}'
    """
}
```

Each process has its own `publishDir` directive that specifies where outputs should be copied.
The processes pass metadata (greeting and language) through their channels, which we'll use later for organizing outputs.

Now look at the main workflow file:

```groovy title="main.nf" linenums="1"
#!/usr/bin/env nextflow

// Enable workflow output definition syntax (required for Nextflow < 25.10)
nextflow.preview.output = true

/*
 * Pipeline parameters
 */
params.input = 'greetings.csv'

// Include modules
include { SAY_HELLO } from './modules/greetings.nf'
include { CONVERT_TO_UPPER } from './modules/greetings.nf'

workflow {

    // Create a channel from the CSV file with metadata
    greeting_ch = channel.fromPath(params.input)
                        .splitCsv(header: true)
                        .map { row -> [row.greeting, row.language] }

    // Create greeting files
    SAY_HELLO(greeting_ch)

    // Convert to uppercase
    CONVERT_TO_UPPER(SAY_HELLO.out)
}
```

The workflow parses the CSV file and extracts both the greeting text and language into a tuple, which flows through both processes.

### 1.2. Run the workflow

Let's run the workflow and see how outputs are organized:

```bash
nextflow run main.nf
```

??? example "Output"

    ```console
    N E X T F L O W   ~  version 25.04.3

    Launching `main.nf` [friendly_ride] DSL2 - revision: abc123

    executor >  local (10)
    [12/abc123] SAY_HELLO (1)        | 5 of 5 ✔
    [34/def456] CONVERT_TO_UPPER (1) | 5 of 5 ✔
    ```

Check the results directory:

```bash
ls -la results/
```

You should see two subdirectories: `greetings/` and `uppercase/`, each containing the relevant output files.

### 1.3. Limitations of publishDir

While `publishDir` works well for simple cases, it has some limitations as workflows grow more complex:

1. **Scattered configuration**: Whether you define `publishDir` in processes or via configuration selectors like `withName` and `withLabel`, the publish settings for each process are specified separately. To understand the full output structure of a large workflow, you need to piece together information from many different places.

2. **No automatic record of outputs**: When a workflow completes, you're left with a directory of files but no machine-readable summary of what was produced. If you want a CSV or JSON file listing all outputs with their associated metadata (often called a "manifest"), you have to build that yourself.

3. **Repetitive patterns**: If multiple processes need similar publish configurations, you end up repeating yourself—either in process definitions or across many configuration selectors.

4. **Coupling between processes and output structure**: The output organization is tied to process-level configuration. If you want to reorganize outputs (say, grouping by sample instead of by process), you need to update the configuration for each affected process.

The workflow output definition syntax addresses these limitations by centralizing output configuration in one place.

### Takeaway

The traditional `publishDir` approach works but scatters output configuration across process definitions and doesn't provide automatic documentation of outputs.

### What's next?

In the next section, we'll introduce the workflow output definition syntax that centralizes output management.

---

## 2. Introducing workflow outputs

### 2.1. A different approach to publishing

With `publishDir`, each process is responsible for publishing its own outputs.
The process definition (or its configuration) specifies where files should go.
This means publishing logic is distributed across your workflow—every process that produces user-facing outputs needs its own publish configuration.

The workflow output definition syntax takes a different approach: **publishing is handled at the workflow level, not the process level**.

Instead of processes deciding where to put their outputs, you:

1. Write processes that simply emit their outputs to channels (no `publishDir`)
2. Declare in your workflow which channels contain outputs worth publishing
3. Configure in a separate `output {}` block how those outputs should be organized

This separation means processes focus purely on computation, while output organization is managed in one central place.

### 2.2. The syntax

The workflow output definition syntax uses two constructs:

1. A `publish:` section inside your workflow that declares which channels to publish
2. An `output {}` block that configures how those outputs are organized

!!! note "Feature flag"

    This training environment uses an older version of Nextflow that requires a feature flag for workflow outputs.
    The starter script already includes this line:

    ```groovy
    nextflow.preview.output = true
    ```

    In Nextflow 25.10 and later, this flag is no longer needed and can be removed.

Let's modify our workflow to use this new syntax.

### 2.3. Remove publishDir from processes

Edit `modules/greetings.nf` to remove the `publishDir` directives from both processes:

=== "After"

    ```groovy title="modules/greetings.nf" linenums="1" hl_lines="5-6"
    #!/usr/bin/env nextflow

    /*
     * Create a greeting file from input text
     * Note: No publishDir - outputs managed by workflow output block
     */
    process SAY_HELLO {

        input:
            tuple val(greeting), val(language)

        output:
            tuple val(greeting), val(language), path("${greeting}-output.txt")

        script:
        """
        echo '$greeting' > '${greeting}-output.txt'
        """
    }

    /*
     * Convert greeting to uppercase
     * Note: No publishDir - outputs managed by workflow output block
     */
    process CONVERT_TO_UPPER {

        input:
            tuple val(greeting), val(language), path(input_file)

        output:
            tuple val(greeting), val(language), path("UPPER-${input_file}")

        script:
        """
        cat '$input_file' | tr '[a-z]' '[A-Z]' > 'UPPER-${input_file}'
        """
    }
    ```

=== "Before"

    ```groovy title="modules/greetings.nf" linenums="1" hl_lines="8 27"
    #!/usr/bin/env nextflow

    /*
     * Create a greeting file from input text
     */
    process SAY_HELLO {

        publishDir 'results/greetings', mode: 'copy'

        input:
            tuple val(greeting), val(language)

        output:
            tuple val(greeting), val(language), path("${greeting}-output.txt")

        script:
        """
        echo '$greeting' > '${greeting}-output.txt'
        """
    }

    /*
     * Convert greeting to uppercase
     */
    process CONVERT_TO_UPPER {

        publishDir 'results/uppercase', mode: 'copy'

        input:
            tuple val(greeting), val(language), path(input_file)

        output:
            tuple val(greeting), val(language), path("UPPER-${input_file}")

        script:
        """
        cat '$input_file' | tr '[a-z]' '[A-Z]' > 'UPPER-${input_file}'
        """
    }
    ```

The only change is removing the `publishDir` directives.
The inputs, outputs, and script remain exactly the same.

### 2.4. Add the publish section to main.nf

Now update `main.nf` to add the `publish:` section inside the workflow:

=== "After"

    ```groovy title="main.nf" linenums="1" hl_lines="17 29-31"
    #!/usr/bin/env nextflow

    // Enable workflow output definition syntax (required for Nextflow < 25.10)
    nextflow.preview.output = true

    /*
     * Pipeline parameters
     */
    params.input = 'greetings.csv'

    // Include modules
    include { SAY_HELLO } from './modules/greetings.nf'
    include { CONVERT_TO_UPPER } from './modules/greetings.nf'

    workflow {

        main:
        // Create a channel from the CSV file with metadata
        greeting_ch = channel.fromPath(params.input)
                            .splitCsv(header: true)
                            .map { row -> [row.greeting, row.language] }

        // Create greeting files
        SAY_HELLO(greeting_ch)

        // Convert to uppercase
        CONVERT_TO_UPPER(SAY_HELLO.out)

        publish:
        greetings = SAY_HELLO.out
        uppercase = CONVERT_TO_UPPER.out
    }
    ```

=== "Before"

    ```groovy title="main.nf" linenums="1"
    #!/usr/bin/env nextflow

    // Enable workflow output definition syntax (required for Nextflow < 25.10)
    nextflow.preview.output = true

    /*
     * Pipeline parameters
     */
    params.input = 'greetings.csv'

    // Include modules
    include { SAY_HELLO } from './modules/greetings.nf'
    include { CONVERT_TO_UPPER } from './modules/greetings.nf'

    workflow {

        // Create a channel from the CSV file with metadata
        greeting_ch = channel.fromPath(params.input)
                            .splitCsv(header: true)
                            .map { row -> [row.greeting, row.language] }

        // Create greeting files
        SAY_HELLO(greeting_ch)

        // Convert to uppercase
        CONVERT_TO_UPPER(SAY_HELLO.out)
    }
    ```

When using a `publish:` section, the workflow content must be placed in a `main:` block.
The `publish:` section then declares named outputs—here `greetings` and `uppercase`—each assigned to an output channel from a process.

### 2.5. Add the output block

Now add the `output {}` block after the workflow to configure how outputs are organized:

=== "After"

    ```groovy title="main.nf" linenums="1" hl_lines="34-47"
    #!/usr/bin/env nextflow

    // Enable workflow output definition syntax (required for Nextflow < 25.10)
    nextflow.preview.output = true

    /*
     * Pipeline parameters
     */
    params.input = 'greetings.csv'

    // Include modules
    include { SAY_HELLO } from './modules/greetings.nf'
    include { CONVERT_TO_UPPER } from './modules/greetings.nf'

    workflow {

        main:
        // Create a channel from the CSV file with metadata
        greeting_ch = channel.fromPath(params.input)
                            .splitCsv(header: true)
                            .map { row -> [row.greeting, row.language] }

        // Create greeting files
        SAY_HELLO(greeting_ch)

        // Convert to uppercase
        CONVERT_TO_UPPER(SAY_HELLO.out)

        publish:
        greetings = SAY_HELLO.out
        uppercase = CONVERT_TO_UPPER.out
    }

    /*
     * Output block defines how published outputs are organized
     */
    output {
        greetings {
            mode 'copy'
            path 'greetings'
        }

        uppercase {
            mode 'copy'
            path 'uppercase'
        }
    }
    ```

=== "Before"

    ```groovy title="main.nf" linenums="1"
    #!/usr/bin/env nextflow

    // Enable workflow output definition syntax (required for Nextflow < 25.10)
    nextflow.preview.output = true

    /*
     * Pipeline parameters
     */
    params.input = 'greetings.csv'

    // Include modules
    include { SAY_HELLO } from './modules/greetings.nf'
    include { CONVERT_TO_UPPER } from './modules/greetings.nf'

    workflow {

        main:
        // Create a channel from the CSV file with metadata
        greeting_ch = channel.fromPath(params.input)
                            .splitCsv(header: true)
                            .map { row -> [row.greeting, row.language] }

        // Create greeting files
        SAY_HELLO(greeting_ch)

        // Convert to uppercase
        CONVERT_TO_UPPER(SAY_HELLO.out)

        publish:
        greetings = SAY_HELLO.out
        uppercase = CONVERT_TO_UPPER.out
    }
    ```

The `output {}` block:

- Configures each named output with its subdirectory path
- Sets the publish mode to `copy` for each output
- The base output directory defaults to `results` (override with `-output-dir`)

### 2.6. Run the updated workflow

Clean up previous results and run:

```bash
rm -rf results work .nextflow*
nextflow run main.nf
```

The outputs should be organized the same way as before, but now the configuration is centralized in one place.

### Takeaway

The workflow output definition syntax separates output configuration from process definitions:

- The `publish:` section declares which channels to publish
- The `output {}` block configures paths and options

### What's next?

In the next section, we'll use dynamic paths to organize outputs by metadata.

---

## 3. Dynamic publish paths

### 3.1. Organizing by metadata

One powerful feature of workflow outputs is the ability to use closures to dynamically determine output paths based on the data itself.

Since our outputs include language metadata, we can organize files by language.

Update the `output {}` block to use dynamic paths:

=== "After"

    ```groovy title="main.nf" linenums="34" hl_lines="7 12"
    /*
     * Output block defines how published outputs are organized
     */
    output {
        greetings {
            mode 'copy'
            path { greeting, language, file -> "greetings/${language}" }
        }

        uppercase {
            mode 'copy'
            path { greeting, language, file -> "uppercase/${language}" }
        }
    }
    ```

=== "Before"

    ```groovy title="main.nf" linenums="34" hl_lines="7 12"
    /*
     * Output block defines how published outputs are organized
     */
    output {
        greetings {
            mode 'copy'
            path 'greetings'
        }

        uppercase {
            mode 'copy'
            path 'uppercase'
        }
    }
    ```

The closure receives the elements of the output tuple (greeting, language, file) and returns the subdirectory path.

### 3.2. Run with dynamic paths

Clean up and run:

```bash
rm -rf results work .nextflow*
nextflow run main.nf
```

Now check the results:

```bash
find results -type f
```

??? example "Output"

    ```console
    results/greetings/English/Hello-output.txt
    results/greetings/French/Bonjour-output.txt
    results/greetings/Spanish/Holà-output.txt
    results/greetings/Italian/Ciao-output.txt
    results/greetings/German/Hallo-output.txt
    results/uppercase/English/UPPER-Hello-output.txt
    results/uppercase/French/UPPER-Bonjour-output.txt
    results/uppercase/Spanish/UPPER-Holà-output.txt
    results/uppercase/Italian/UPPER-Ciao-output.txt
    results/uppercase/German/UPPER-Hallo-output.txt
    ```

The outputs are now organized by language, making it easy to find results for specific languages.

### 3.3. Override the output directory

You can override the output directory from the command line:

```bash
nextflow run main.nf -output-dir my_results
```

This creates outputs in `my_results/` instead of `results/`.

### Takeaway

Dynamic paths let you organize outputs based on metadata:

- Use closures to compute paths from output tuple elements
- The `-output-dir` flag overrides the base directory

### What's next?

In the next section, we'll add index files to document our outputs.

---

## 4. Index files

### 4.1. Generating output manifests

Index files are CSV, JSON, or YAML manifests that document what outputs were produced.
They're useful for:

- Downstream pipelines that need to consume your outputs
- Documentation of what was generated
- Quality control and auditing

### 4.2. Add index file configuration

Update the `output {}` block to generate index files:

=== "After"

    ```groovy title="main.nf" linenums="34" hl_lines="8-10 16-18"
    /*
     * Output block defines how published outputs are organized
     */
    output {
        greetings {
            mode 'copy'
            path { greeting, language, file -> "greetings/${language}" }
            index {
                path 'greetings/index.csv'
            }
        }

        uppercase {
            mode 'copy'
            path { greeting, language, file -> "uppercase/${language}" }
            index {
                path 'uppercase/index.csv'
            }
        }
    }
    ```

=== "Before"

    ```groovy title="main.nf" linenums="34"
    /*
     * Output block defines how published outputs are organized
     */
    output {
        greetings {
            mode 'copy'
            path { greeting, language, file -> "greetings/${language}" }
        }

        uppercase {
            mode 'copy'
            path { greeting, language, file -> "uppercase/${language}" }
        }
    }
    ```

### 4.3. Run and check the index files

Clean up and run:

```bash
rm -rf results work .nextflow*
nextflow run main.nf
```

View the generated index file:

```bash
cat results/greetings/index.csv
```

??? example "Output"

    ```csv
    "greeting","language","file"
    "Hello","English","greetings/English/Hello-output.txt"
    "Bonjour","French","greetings/French/Bonjour-output.txt"
    "Holà","Spanish","greetings/Spanish/Holà-output.txt"
    "Ciao","Italian","greetings/Italian/Ciao-output.txt"
    "Hallo","German","greetings/German/Hallo-output.txt"
    ```

The index file contains all the metadata from the output tuples plus the path to each file.

### 4.4. Using JSON format

You can also generate JSON index files:

```groovy title="main.nf"
greetings {
    path { greeting, language, file -> "greetings/${language}" }
    index {
        path 'greetings/index.json'
    }
}
```

!!! tip "Use maps for cleaner index files"

    For better control over field names in index files, use maps instead of tuples in your outputs:

    ```groovy
    output:
        tuple path("${greeting}-output.txt"), val([greeting: greeting, language: language])
    ```

    This gives you named fields in the index file rather than positional names like "v0", "v1".

### Takeaway

Index files provide automatic documentation of workflow outputs:

- Generated in CSV, JSON, or YAML format
- Include all metadata from output tuples
- Useful for downstream pipelines and auditing

### What's next?

Let's summarize what we've learned.

---

## Summary

### Key patterns

**Workflow output definition syntax:**

```groovy
workflow {
    // ... process calls ...

    publish:
    output_name = PROCESS.out
}

output {
    output_name {
        mode 'copy'
        path { metadata -> "subdir/${metadata}" }
        index {
            path 'index.csv'
        }
    }
}
```

**Key benefits:**

- Centralized output configuration
- Dynamic paths based on metadata
- Automatic index file generation
- Override output directory with `-output-dir`

### When to use workflow outputs vs publishDir

| Use Case                                        | Approach             |
| ----------------------------------------------- | -------------------- |
| Simple pipelines with few outputs               | `publishDir` is fine |
| Complex output organization                     | Workflow outputs     |
| Need output manifests                           | Workflow outputs     |
| Multiple processes publishing to same structure | Workflow outputs     |
| Quick prototyping                               | `publishDir`         |

### Additional resources

- [Nextflow workflow outputs tutorial](https://www.nextflow.io/docs/latest/tutorials/workflow-outputs.html)
- [Output block reference](https://www.nextflow.io/docs/latest/reference/workflow-output-def.html)

---

## What's next?

Congratulations on completing this side quest!
You've learned how to use workflow output definitions to organize your pipeline outputs and generate documentation.

Return to the [Side Quests](./index.md) menu to continue your training journey.
