# Part 2: Single-sample processing

In Part 1, you tested the {TOOL_A} and {TOOL_B} commands manually in their respective containers.
Now we're going to wrap those same commands into a Nextflow workflow.

## Assignment

In this part of the course, we're going to develop a workflow that does the following:

1. {PROCESS_1_DESCRIPTION}
2. {PROCESS_2_DESCRIPTION}

<figure class="excalidraw">
--8<-- "docs/en/docs/nf4_science/{DOMAIN_DIR}/img/{DIAGRAM_1}.svg"
</figure>

This replicates the steps from Part 1, where you ran these commands manually in their containers.

As a starting point, we provide you with a workflow file, `{DOMAIN_DIR}.nf`, that outlines the main parts of the workflow, as well as two module files, {TOOL_A_MODULE}.nf and {TOOL_B_MODULE}.nf, that outline the structure of the modules.
These files are not functional; their purpose is just to serve as scaffolds for you to fill in with the interesting parts of the code.

## Lesson plan

In order to make the development process more educational, we've broken this down into {N} steps:

1. **Write a single-stage workflow that runs {TOOL_A_ACTION}.**
   This covers creating a module, importing it, and calling it in a workflow.
2. **Add a second process to run {TOOL_B_ACTION}.**
   This introduces chaining process outputs to inputs and handling accessory files.
3. **Adapt the workflow to run on a batch of samples.**
   This covers parallel execution and introduces tuples to keep associated files together.
4. **Make the workflow accept a samplesheet containing a batch of input files.**
   This demonstrates a common pattern for providing inputs in bulk.

Each step focuses on a specific aspect of workflow development.

---

## 1. Write a single-stage workflow that runs {TOOL_A_ACTION}

This first step focuses on the basics: loading {PRIMARY_INPUT_TYPE} and {TOOL_A_OUTPUT_DESCRIPTION}.

Recall the `{TOOL_A_COMMAND_NAME}` command from [Part 1](01_method.md):

```bash
{TOOL_A_COMMAND_SUMMARY}
```

The command takes {INPUT_DESCRIPTION} and produces {OUTPUT_DESCRIPTION}.
The container URI was `{TOOL_A_CONTAINER_URI}`.

We're going to take this information and wrap it in Nextflow in three stages:

1. Set up the input
2. Write the process and call it in the workflow
3. Configure the output handling

### 1.1. Set up the input

We need to declare an input parameter, create a test profile to provide a convenient default value, and create an input channel.

#### 1.1.1. Add an input parameter declaration

In the main workflow file `{DOMAIN_DIR}.nf`, under the `Pipeline parameters` section, declare a CLI parameter called `{PRIMARY_PARAM_NAME}`.

=== "After"

    ```groovy title="{DOMAIN_DIR}.nf" linenums="5" hl_lines="4-7"
    /*
     * Pipeline parameters
     */
    params {
        // Primary input
        {PRIMARY_PARAM_NAME}: Path
    }
    ```

=== "Before"

    ```groovy title="{DOMAIN_DIR}.nf" linenums="5"
    /*
     * Pipeline parameters
     */

    // Primary input
    ```

That sets up the CLI parameter, but we don't want to type out the file path every time we run the workflow during development.
There are multiple options for providing a default value; here we use a test profile.

#### 1.1.2. Create a test profile with a default value in `nextflow.config`

A test profile provides convenient default values for trying out a workflow without specifying inputs on the command line.
This is a common convention in the Nextflow ecosystem (see [Hello Config](../../hello_nextflow/06_hello_config.md) for more detail).

Add a `profiles` block to `nextflow.config` with a `test` profile that sets the `{PRIMARY_PARAM_NAME}` parameter to one of the test {PRIMARY_INPUT_TYPE}s.

=== "After"

    ```groovy title="nextflow.config" linenums="1" hl_lines="3-7"
    docker.enabled = true

    profiles {
        test {
            params.{PRIMARY_PARAM_NAME} = "${projectDir}/data/{TEST_INPUT_PATH}"
        }
    }
    ```

=== "Before"

    ```groovy title="nextflow.config" linenums="1"
    docker.enabled = true
    ```

Here, we're using `${projectDir}`, a built-in Nextflow variable that points to the directory where the workflow script is located.
This makes it easy to reference data files and other resources without hardcoding absolute paths.

#### 1.1.3. Set up the input channel

{INPUT_CHANNEL_INSTRUCTIONS}

### 1.2. Write the {TOOL_A_NAME} module

{MODULE_INSTRUCTIONS_COVERING: container, input, output, script}

### 1.3. Import and call the module in the workflow

{IMPORT_AND_CALL_INSTRUCTIONS}

### 1.4. Run the workflow

At this point, we have a one-step workflow that should be fully functional.

We can run it with `-profile test` to use the default value set up in the test profile and avoid having to write the path on the command line.

```bash
nextflow run {DOMAIN_DIR}.nf -profile test
```

??? success "Command output"

    ```console
    N E X T F L O W   ~  version 25.10.2

    ┃ Launching `{DOMAIN_DIR}.nf` [{RUN_NAME}] DSL2 - revision: {HASH}

    executor >  local (1)
    [{HASH}] {PROCESS_A_NAME} (1) | 1 of 1 ✔
    ```

{VERIFY_OUTPUT_INSTRUCTIONS}

### Takeaway

You know how to create a module containing a process, import it into a workflow, call it with an input channel, and publish the results.

### What's next?

Add a second process to chain additional analysis steps.

---

## 2. Add a second process to run {TOOL_B_ACTION}

{DESCRIBE_WHAT_THIS_STEP_ADDS}

Recall the `{TOOL_B_COMMAND_NAME}` command from [Part 1](01_method.md):

```bash
{TOOL_B_COMMAND_SUMMARY}
```

### 2.1. Write the {TOOL_B_NAME} module

{MODULE_INSTRUCTIONS}

### 2.2. Import and call the module in the workflow

{IMPORT_AND_CALL_INSTRUCTIONS_CHAINING_OUTPUTS}

### 2.3. Run the workflow

```bash
nextflow run {DOMAIN_DIR}.nf -profile test
```

??? success "Command output"

    ```console
    {EXPECTED_OUTPUT}
    ```

{VERIFY_OUTPUT_INSTRUCTIONS}

### Takeaway

You know how to chain process outputs to inputs and handle accessory files in the workflow.

### What's next?

Scale the workflow to process multiple samples in parallel.

---

## 3. Adapt the workflow to run on a batch of samples

So far, the workflow processes a single sample.
To handle multiple samples, we need to modify how inputs are provided and leverage Nextflow's dataflow paradigm to parallelize execution.

### 3.1. {BATCH_INPUT_INSTRUCTIONS}

{INSTRUCTIONS_FOR_SWITCHING_TO_SAMPLESHEET_OR_BATCH_INPUT}

### 3.2. Run the workflow on multiple samples

```bash
nextflow run {DOMAIN_DIR}.nf -profile test
```

??? success "Command output"

    ```console
    {EXPECTED_OUTPUT_MULTIPLE_SAMPLES}
    ```

### Takeaway

You know how to leverage Nextflow's dataflow paradigm to parallelize per-sample processing across multiple input samples.

### What's next?

In [Part 3](03_multi_sample.md), you will add multi-sample aggregation to combine results across all samples.
