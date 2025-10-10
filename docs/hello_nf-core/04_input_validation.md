# Part 4: Input validation

In this fourth part of the Hello nf-core training course, we show you how to use the nf-schema plugin to validate inputs.

One of the key features of nf-core pipelines is robust input validation. The nf-schema plugin provides automated validation of pipeline parameters and sample sheets against JSON schemas, helping catch errors early and provide clear feedback to users.

In this section, we'll add proper input validation to our pipeline using the nf-schema plugin.

!!! note

    This section assumes you have completed [Part 2: Rewrite Hello for nf-core](./02_rewrite_hello.md) and have a working `core-hello` pipeline.

---

## 1. Understand schema-based validation

Before we implement validation, let's understand how nf-schema works and what schema files do.

### 1.1. What is nf-schema?

The [nf-schema plugin](https://nextflow-io.github.io/nf-schema/latest/) is a Nextflow plugin that provides functionality for:

- **Parameter validation**: Validating pipeline parameters against a JSON schema
- **Sample sheet validation**: Validating input sample sheets and converting them to channels
- **Help text generation**: Automatically generating help text from schema definitions
- **Parameter summary**: Displaying which parameters differ from defaults

nf-schema is the successor to the deprecated nf-validation plugin and uses standard JSON Schema Draft 2020-12 for validation.

### 1.2. The two schema files

An nf-core pipeline typically uses two schema files:

1. **`nextflow_schema.json`**: Defines and validates pipeline parameters (e.g., `--input`, `--outdir`, `--batch`)
2. **`assets/schema_input.json`**: Defines the structure and validates the contents of input sample sheets

Both files use [JSON Schema](https://json-schema.org/) format, a standard for describing and validating JSON data structures.

### 1.3. Examine an existing schema

Let's look at the `schema_input.json` file that was created with our pipeline template:

```bash
cat core-hello/assets/schema_input.json
```

```json title="core-hello/assets/schema_input.json" linenums="1"
{
    "$schema": "https://json-schema.org/draft/2020-12/schema",
    "$id": "https://raw.githubusercontent.com/core/hello/main/assets/schema_input.json",
    "title": "core/hello pipeline - params.input schema",
    "description": "Schema for the file provided with params.input",
    "type": "array",
    "items": {
        "type": "object",
        "properties": {
            "sample": {
                "type": "string",
                "pattern": "^\\S+$",
                "errorMessage": "Sample name must be provided and cannot contain spaces",
                "meta": ["id"]
            },
            "fastq_1": {
                "type": "string",
                "format": "file-path",
                "exists": true,
                "pattern": "^\\S+\\.f(ast)?q\\.gz$",
                "errorMessage": "FastQ file for reads 1 must be provided, cannot contain spaces and must have extension '.fq.gz' or '.fastq.gz'"
            },
            "fastq_2": {
                "type": "string",
                "format": "file-path",
                "exists": true,
                "pattern": "^\\S+\\.f(ast)?q\\.gz$",
                "errorMessage": "FastQ file for reads 2 cannot contain spaces and must have extension '.fq.gz' or '.fastq.gz'"
            }
        },
        "required": ["sample", "fastq_1"]
    }
}
```

This schema defines:

- The input is an array of objects (one per sample)
- Each object has fields: `sample`, `fastq_1`, and `fastq_2`
- Field types and validation rules (e.g., file paths must exist and match certain patterns)
- Which fields are required
- Custom error messages for validation failures

This schema is appropriate for FASTQ sequencing data but doesn't match our simple greetings CSV format.

### Takeaway

You now know what nf-schema is, what it does, and how schema files define validation rules for pipeline inputs.

### What's next?

Create a custom schema file for our greetings input format.

---

## 2. Create a schema for the greetings input

Our `greetings.csv` file is very simple - it just contains one greeting per line. Let's create a schema that validates this format.

### 2.1. Understand the greetings.csv format

Let's remind ourselves what our input looks like:

```bash
cat core-hello/assets/greetings.csv
```

```csv title="core-hello/assets/greetings.csv"
Hello
Bonjour
Holà
```

This is a simple CSV with:
- One column (no header)
- One greeting per line
- Text strings with no special format requirements

### 2.2. Design the schema structure

For our use case, we want to:

1. Accept CSV input with one column
2. Treat each row as a greeting string
3. Ensure greetings are not empty
4. Optionally, ensure no whitespace-only entries

We'll structure this as an array of objects, where each object has a `greeting` field.

### 2.3. Create the schema file

Replace the contents of [core-hello/assets/schema_input.json](core-hello/assets/schema_input.json) with the following:

```json title="core-hello/assets/schema_input.json" linenums="1"
{
    "$schema": "https://json-schema.org/draft/2020-12/schema",
    "$id": "https://raw.githubusercontent.com/core/hello/main/assets/schema_input.json",
    "title": "core/hello pipeline - params.input schema",
    "description": "Schema for the greetings file provided with params.input",
    "type": "array",
    "items": {
        "type": "object",
        "properties": {
            "greeting": {
                "type": "string",
                "pattern": "^\\S.*$",
                "errorMessage": "Greeting must be provided and cannot be empty or start with whitespace"
            }
        },
        "required": ["greeting"]
    }
}
```

Let's break down the key parts:

- **`type: "array"`**: The input is parsed as an array (list) of items
- **`items.type: "object"`**: Each item in the array is an object
- **`properties.greeting`**: Defines a field called `greeting`
  - **`type: "string"`**: Must be a text string
  - **`pattern: "^\\S.*$"`**: Must start with a non-whitespace character (but can contain spaces after that)
  - **`errorMessage`**: Custom error message shown if validation fails
- **`required: ["greeting"]`**: The `greeting` field is mandatory

### 2.4. Understand the samplesheet structure expectation

When nf-schema reads a CSV file, it expects either:

1. A **headerless CSV**: Each row becomes an object with fields named by column position (column1, column2, etc.)
2. A **CSV with headers**: The first row defines field names

For our simple case, we'll add a header to our greetings file to match the schema.

### 2.5. Update the greetings.csv file

Add a header line to the greetings file:

=== "After"

    ```csv title="core-hello/assets/greetings.csv" linenums="1"
    greeting
    Hello
    Bonjour
    Holà
    ```

=== "Before"

    ```csv title="core-hello/assets/greetings.csv" linenums="1"
    Hello
    Bonjour
    Holà
    ```

Now the CSV file has a header that matches the field name in our schema.

### Takeaway

You know how to create a JSON schema that defines validation rules for a simple sample sheet input.

### What's next?

Implement the validation in the pipeline code.

---

## 3. Implement input validation with samplesheetToList

Now we need to replace our simple CSV parsing with nf-schema's `samplesheetToList` function, which validates and converts the sample sheet.

### 3.1. Understand samplesheetToList

The `samplesheetToList` function:

1. Reads the input sample sheet (CSV, TSV, JSON, or YAML)
2. Validates it against the provided JSON schema
3. Returns a Groovy list where each entry corresponds to a row
4. Throws helpful error messages if validation fails

Basic usage:

```groovy
include { samplesheetToList } from 'plugin/nf-schema'

def input_list = samplesheetToList(params.input, "assets/schema_input.json")
```

### 3.2. Update the input handling code

Open [core-hello/subworkflows/local/utils_nfcore_hello_pipeline/main.nf](core-hello/subworkflows/local/utils_nfcore_hello_pipeline/main.nf) and locate the section where we create the input channel (around line 64).

We need to:

1. Import the `samplesheetToList` function
2. Use it to validate and parse the input
3. Extract just the greeting strings for our workflow

First, add the import at the top of the file:

=== "After"

    ```groovy title="core-hello/subworkflows/local/utils_nfcore_hello_pipeline/main.nf" linenums="1" hl_lines="13"
    //
    // Subworkflow with functionality that may be useful for any Nextflow pipeline
    //

    import org.yaml.snakeyaml.Yaml
    import groovy.json.JsonOutput

    include { UTILS_NFSCHEMA_PLUGIN  } from '../../nf-core/utils_nfschema_plugin'
    include { UTILS_NEXTFLOW_PIPELINE } from '../../nf-core/utils_nextflow_pipeline'
    include { UTILS_NFCORE_PIPELINE   } from '../../nf-core/utils_nfcore_pipeline'

    include { paramsSummaryMap          } from 'plugin/nf-schema'
    include { samplesheetToList         } from 'plugin/nf-schema'
    ```

=== "Before"

    ```groovy title="core-hello/subworkflows/local/utils_nfcore_hello_pipeline/main.nf" linenums="1"
    //
    // Subworkflow with functionality that may be useful for any Nextflow pipeline
    //

    import org.yaml.snakeyaml.Yaml
    import groovy.json.JsonOutput

    include { UTILS_NFSCHEMA_PLUGIN  } from '../../nf-core/utils_nfschema_plugin'
    include { UTILS_NEXTFLOW_PIPELINE } from '../../nf-core/utils_nextflow_pipeline'
    include { UTILS_NFCORE_PIPELINE   } from '../../nf-core/utils_nfcore_pipeline'

    include { paramsSummaryMap          } from 'plugin/nf-schema'
    ```

Now update the channel creation code:

=== "After"

    ```groovy title="core-hello/subworkflows/local/utils_nfcore_hello_pipeline/main.nf" linenums="64" hl_lines="4-8"
        //
        // Create channel from input file provided through params.input
        //
        ch_samplesheet = Channel.fromList(samplesheetToList(params.input, "${projectDir}/assets/schema_input.json"))
            .map { row ->
                // Extract just the greeting string from each row
                row.greeting
            }

        emit:
        samplesheet = ch_samplesheet
        versions    = ch_versions
    ```

=== "Before"

    ```groovy title="core-hello/subworkflows/local/utils_nfcore_hello_pipeline/main.nf" linenums="64"
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

Let's break down what changed:

1. **`samplesheetToList(params.input, "${projectDir}/assets/schema_input.json")`**: Validates the input file against our schema and returns a list
2. **`Channel.fromList(...)`**: Converts the list into a Nextflow channel
3. **`.map { row -> row.greeting }`**: Extracts just the greeting string from each validated row

### 3.3. Enable parameter validation

The pipeline template already includes code to validate parameters, but we had it disabled with `--validate_params false`. Now that we have proper schema validation set up, let's enable it.

The validation is controlled by the `params.validate_params` parameter in [core-hello/nextflow.config](core-hello/nextflow.config). Let's check that it's set to `true` (the default):

```groovy title="core-hello/nextflow.config" linenums="35"
    validate_params            = true
```

This should already be the default. The validation is performed by the `UTILS_NFSCHEMA_PLUGIN` subworkflow, which is called during pipeline initialization.

### Takeaway

You know how to use the `samplesheetToList` function to validate and parse input sample sheets using JSON schemas.

### What's next?

Test that the validation works correctly.

---

## 4. Test input validation

Let's verify that our validation works by testing both valid and invalid inputs.

### 4.1. Test with valid input

First, let's confirm the pipeline still runs successfully with valid input:

```bash
nextflow run core-hello --outdir core-hello-results -profile test,docker
```

Note that we no longer need `--validate_params false`!

```console title="Output"
 N E X T F L O W   ~  version 24.10.4

Launching `core-hello/main.nf` [serene_volta] DSL2 - revision: c31b966b36

Input/output options
  input                     : core-hello/assets/greetings.csv
  outdir                    : core-hello-results

Institutional config options
  config_profile_name       : Test profile
  config_profile_description: Minimal test dataset to check pipeline function

Core Nextflow options
  runName                   : serene_volta
  containerEngine           : docker
  profile                   : test,docker

!! Only displaying parameters that differ from the pipeline defaults !!
------------------------------------------------------
executor >  local (7)
[5a/1c3d8b] CORE_HELLO:HELLO:sayHello (1)       | 3 of 3 ✔
[2b/9f4a2c] CORE_HELLO:HELLO:convertToUpper (3) | 3 of 3 ✔
[8c/3e7d1a] CORE_HELLO:HELLO:collectGreetings   | 1 of 1 ✔
[4d/6b2f8e] CORE_HELLO:HELLO:cowpy              | 1 of 1 ✔
-[core/hello] Pipeline completed successfully-
```

Great! The pipeline runs successfully and validation passes silently.

### 4.2. Test with invalid input (empty greeting)

Now let's test that validation catches errors. Create a test file with an invalid entry:

```bash
cat > /tmp/invalid_greetings.csv << 'EOF'
greeting
Hello

Holà
EOF
```

This file has an empty second row (just whitespace), which should fail our validation rule.

Try running the pipeline with this invalid input:

```bash
nextflow run core-hello --input /tmp/invalid_greetings.csv --outdir test-results -profile docker
```

```console title="Output"
 N E X T F L O W   ~  version 24.10.4

Launching `core-hello/main.nf` [silly_cuvier] DSL2 - revision: c31b966b36

ERROR ~ Validation of '/tmp/invalid_greetings.csv' file failed!

 -- Check '/tmp/invalid_greetings.csv' --
   -> Entry 2: Greeting must be provided and cannot be empty or start with whitespace

 -- Check '.nextflow.log' file for details
```

Perfect! The validation caught the error and provided a clear, helpful error message pointing to:

- Which file failed validation
- Which entry (row 2) has the problem
- What the specific problem is

### 4.3. Test with missing required field

Let's create another invalid file, this time missing the header:

```bash
cat > /tmp/no_header.csv << 'EOF'
Hello
Bonjour
Holà
EOF
```

Try running with this file:

```bash
nextflow run core-hello --input /tmp/no_header.csv --outdir test-results -profile docker
```

The validation will fail because the first row is treated as a header, but then the actual data rows don't have a `greeting` field defined.

### 4.4. Understanding validation benefits

Input validation provides several important benefits:

1. **Early error detection**: Problems are caught before any processing begins
2. **Clear error messages**: Users get specific, actionable feedback about what's wrong
3. **Reduced debugging time**: No need to trace cryptic errors through the pipeline
4. **Documentation**: The schema serves as documentation of expected input format
5. **Type safety**: Ensures data types match expectations throughout the pipeline

### Takeaway

You know how to test input validation and understand the benefits it provides for pipeline usability and robustness.

### What's next?

Explore additional validation features and parameter schema validation.

---

## 5. Optional: Explore parameter validation

In addition to sample sheet validation, nf-schema also validates pipeline parameters against `nextflow_schema.json`.

### 5.1. Examine the parameter schema

Let's look at a section of the parameter schema:

```bash
grep -A 10 '"batch"' core-hello/nextflow_schema.json
```

You might notice that the `batch` parameter isn't defined yet in the schema. Let's add it.

### 5.2. Add the batch parameter to the schema

The parameter schema can be edited manually, but nf-core provides a tool to help:

```bash
nf-core pipelines schema build
```

This launches an interactive tool that helps you add and configure parameters. However, for our simple case, we can edit the JSON directly.

Open [core-hello/nextflow_schema.json](core-hello/nextflow_schema.json) and find the `"input_output_options"` section. Add the `batch` parameter:

```json title="core-hello/nextflow_schema.json (excerpt)" linenums="30" hl_lines="17-22"
"input_output_options": {
    "title": "Input/output options",
    "type": "object",
    "fa_icon": "fas fa-terminal",
    "description": "Define where the pipeline should find input data and save output data.",
    "required": ["input", "outdir"],
    "properties": {
        "input": {
            "type": "string",
            "format": "file-path",
            "exists": true,
            "schema": "assets/schema_input.json",
            "mimetype": "text/csv",
            "pattern": "^\\S+\\.csv$",
            "description": "Path to comma-separated file containing greetings.",
            "fa_icon": "fas fa-file-csv"
        },
        "batch": {
            "type": "string",
            "default": "batch-01",
            "description": "Name for this batch of greetings",
            "fa_icon": "fas fa-tag"
        },
        "outdir": {
            "type": "string",
            "format": "directory-path",
            "description": "The output directory where the results will be saved.",
            "fa_icon": "fas fa-folder-open"
        }
    }
}
```

### 5.3. Test parameter validation

Now try running with an invalid parameter type:

```bash
nextflow run core-hello --batch 12345 --outdir test-results -profile test,docker
```

The pipeline should run fine because `12345` will be converted to the string `"12345"`. Parameter validation is more useful for catching missing required parameters or invalid file paths.

Try running without the required `input` parameter:

```bash
nextflow run core-hello --outdir test-results -profile docker
```

```console title="Output"
ERROR ~ Validation of pipeline parameters failed!

 -- Check '.nextflow.log' file for details
 * --input: required property is missing
```

Excellent! The validation catches missing required parameters.

### Takeaway

You know how parameter validation works and how to add parameter definitions to the schema.

### What's next?

Continue to the next section of the training, or explore the [nf-schema documentation](https://nextflow-io.github.io/nf-schema/latest/) to learn about advanced features like metadata maps, multi-file samples, and custom validation rules.
