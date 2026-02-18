# Part 1: Plugin Basics

In this section, you'll learn how plugins extend Nextflow, then try three different plugins to see them in action.

---

## 1. How plugins work

Plugins extend Nextflow through several types of extension:

| Extension Type    | What it does                                 | Example                      |
| ----------------- | -------------------------------------------- | ---------------------------- |
| Functions         | Add custom functions callable from workflows | `samplesheetToList()`        |
| Workflow monitors | Respond to events like task completion       | Custom logging, Slack alerts |
| Executors         | Add task execution backends                  | AWS Batch, Kubernetes        |
| Filesystems       | Add storage backends                         | S3, Azure Blob               |

Functions and workflow monitors (called "trace observers" in the Nextflow API) are the most common types for plugin authors.
Executors and filesystems are typically created by platform vendors.

The next three exercises show you function plugins and an observer plugin, so you can see both types in action.

---

## 2. Use a function plugin: nf-hello

The [nf-hello](https://github.com/nextflow-io/nf-hello) plugin provides a `randomString` function that generates random strings.
In this exercise, you'll start with a workflow that defines this function locally, then replace it with the plugin version.
This demonstrates the core plugin workflow: declare a plugin, import its functions, and use them like any other Nextflow function.

### 2.1. See the starting point

The `random_id_example.nf` file contains a workflow with a local `randomString` function:

```bash
cat random_id_example.nf
```

```groovy title="Output"
#!/usr/bin/env nextflow

/**
 * Generate a random alphanumeric string
 */
def randomString(int length) {
    def chars = ('a'..'z') + ('A'..'Z') + ('0'..'9')
    def random = new Random()
    return (1..length).collect { chars[random.nextInt(chars.size())] }.join()
}

workflow {
    // Generate random IDs for each sample
    Channel.of('sample_A', 'sample_B', 'sample_C')
        .map { sample -> "${sample}_${randomString(8)}" }
        .view()
}
```

The workflow creates a channel of three sample names and appends a random string to each one.
The `randomString` function is defined locally at the top of the file.

Run it:

```bash
nextflow run random_id_example.nf
```

```console title="Output"
sample_A_Pmlkc9S0
sample_B_ErYxlPsF
sample_C_AJvSgOVt
```

(Your random strings will differ.)

This works, but the `randomString` function is defined inside this script.
If you wanted to use it in another pipeline, you'd have to copy it.
A plugin solves this by packaging the function so any pipeline can import it.

### 2.2. Configure the plugin

Add the plugin to your `nextflow.config`:

```groovy title="nextflow.config"
// Configuration for plugin development exercises
plugins {
    id 'nf-hello@0.5.0'
}
```

Plugins are declared in `nextflow.config` using the `plugins {}` block.
Nextflow automatically downloads them from the registry at runtime.

### 2.3. Use the plugin function

Update `random_id_example.nf` to import `randomString` from the plugin instead of defining it locally:

```groovy title="random_id_example.nf"
#!/usr/bin/env nextflow

include { randomString } from 'plugin/nf-hello'

workflow {
    // Generate random IDs for each sample
    Channel.of('sample_A', 'sample_B', 'sample_C')
        .map { sample -> "${sample}_${randomString(8)}" }
        .view()
}
```

Import plugin functions with `include { function } from 'plugin/plugin-id'`.
This is the same `include` syntax used for Nextflow modules, with a `plugin/` prefix.
You can see the [source code for `randomString`](https://github.com/nextflow-io/nf-hello/blob/e67bddebfa589c7ae51f41bf780c92068dc09e93/plugins/nf-hello/src/main/nextflow/hello/HelloExtension.groovy#L110) in the nf-hello repository on GitHub.

### 2.4. Run it

```bash
nextflow run random_id_example.nf
```

```console title="Output"
Pipeline is starting! 🚀
sample_A_xcwzhtbm
sample_B_yqurtfsq
sample_C_lpxepimu
Pipeline complete! 👋
```

(Your random strings will be different.
The plugin version generates lowercase strings, while the local version used mixed case and digits.)

The first time you use a plugin, Nextflow downloads it automatically from the registry.
After that, any pipeline that declares `nf-hello@0.5.0` gets the exact same `randomString` function without needing to copy code between projects.

---

## 3. Use a function plugin: nf-schema

The [nf-schema](https://github.com/nextflow-io/nf-schema) plugin is one of the most widely-used Nextflow plugins.
It provides `samplesheetToList`, a function that parses CSV/TSV files using a JSON schema that defines the expected columns and types.

In `main.nf`, the pipeline reads `greetings.csv` using `splitCsv` and a manual `map`:

```groovy
greeting_ch = channel.fromPath(params.input)
                    .splitCsv(header: true)
                    .map { row -> row.greeting }
```

The nf-schema plugin can replace this with validated, schema-driven parsing.
A JSON schema file (`greetings_schema.json`) is provided in the exercise directory.

### 3.1. Look at the schema

```bash
cat greetings_schema.json
```

```json title="Output"
{
  "$schema": "https://json-schema.org/draft/2020-12/schema",
  "type": "array",
  "items": {
    "type": "object",
    "properties": {
      "greeting": {
        "type": "string",
        "description": "The greeting text"
      },
      "language": {
        "type": "string",
        "description": "The language of the greeting"
      }
    },
    "required": ["greeting"]
  }
}
```

The schema defines two columns (`greeting` and `language`) and marks `greeting` as required.
If someone passes a CSV missing the `greeting` column, nf-schema catches the error before the pipeline runs.

### 3.2. Add nf-schema to the config

Update `nextflow.config` to include both plugins:

=== "After"

    ```groovy title="nextflow.config" hl_lines="3"
    plugins {
        id 'nf-hello@0.5.0'
        id 'nf-schema@2.6.1'
    }
    ```

=== "Before"

    ```groovy title="nextflow.config"
    plugins {
        id 'nf-hello@0.5.0'
    }
    ```

### 3.3. See the starting point

The `schema_example.nf` file reads `greetings.csv` using the same `splitCsv` + `map` pattern from `main.nf`:

```bash
cat schema_example.nf
```

```groovy title="Output"
#!/usr/bin/env nextflow

params.input = 'greetings.csv'

workflow {
    Channel.fromPath(params.input)
        .splitCsv(header: true)
        .map { row -> row.greeting }
        .view { greeting -> "Greeting: $greeting" }
}
```

### 3.4. Update it to use samplesheetToList

Replace the `splitCsv` + `map` pattern with `samplesheetToList`:

=== "After"

    ```groovy title="schema_example.nf" hl_lines="3 7-8"
    #!/usr/bin/env nextflow

    include { samplesheetToList } from 'plugin/nf-schema'

    params.input = 'greetings.csv'

    workflow {
        Channel.fromList(samplesheetToList(params.input, 'greetings_schema.json'))
            .map { row -> row[0] }
            .view { greeting -> "Greeting: $greeting" }
    }
    ```

=== "Before"

    ```groovy title="schema_example.nf"
    #!/usr/bin/env nextflow

    params.input = 'greetings.csv'

    workflow {
        Channel.fromPath(params.input)
            .splitCsv(header: true)
            .map { row -> row.greeting }
            .view { greeting -> "Greeting: $greeting" }
    }
    ```

Instead of `splitCsv` and a manual `map` to extract fields, `samplesheetToList` parses the CSV according to the schema.
Each row becomes a list of values in column order, so `row[0]` is the greeting and `row[1]` is the language.

### 3.5. Run it

```bash
nextflow run schema_example.nf
```

```console title="Output"
Greeting: Hello
Greeting: Bonjour
Greeting: Holà
Greeting: Ciao
Greeting: Hallo
```

The result is the same, but the schema adds validation.
In real pipelines, `samplesheetToList` is used to parse complex sample sheets with many columns, where manual `splitCsv` + `map` would be error-prone.

---

## 4. Use an observer plugin: nf-co2footprint

Not all plugins provide functions to import.
The [nf-co2footprint](https://github.com/nextflow-io/nf-co2footprint) plugin uses a **trace observer** to monitor your pipeline's resource usage and estimate its carbon footprint.
You don't need to change any pipeline code; just add it to the config.

### 4.1. Add nf-co2footprint to the config

Update `nextflow.config`:

=== "After"

    ```groovy title="nextflow.config" hl_lines="4"
    plugins {
        id 'nf-hello@0.5.0'
        id 'nf-schema@2.6.1'
        id 'nf-co2footprint'
    }
    ```

=== "Before"

    ```groovy title="nextflow.config"
    plugins {
        id 'nf-hello@0.5.0'
        id 'nf-schema@2.6.1'
    }
    ```

### 4.2. Run the pipeline

```bash
nextflow run main.nf
```

You'll see additional log messages from the plugin during execution, including some warnings.
These are normal; the plugin is estimating resource usage with limited information from this small example.

At the end, look for a line like:

```console title="Output (partial)"
🌱 The workflow run used 163.4 uWh of electricity, resulting in the release of 78.43 ug of CO₂ equivalents into the atmosphere.
```

(Your numbers will differ.)

### 4.3. View the report

The plugin generates output files in your working directory:

```bash
ls co2footprint_*
```

```console title="Output"
co2footprint_report_<timestamp>.html
co2footprint_summary_<timestamp>.txt
co2footprint_trace_<timestamp>.txt
```

Look at the summary:

```bash
cat co2footprint_summary_*.txt
```

The summary converts the CO2 estimate into relatable comparisons, like the equivalent distance driven by car or the time a tree would need to sequester the same amount.

This plugin works entirely through the observer mechanism.
It hooks into workflow lifecycle events to collect resource metrics, then generates its report when the pipeline completes.
No `include` statement is needed because it doesn't provide functions; it runs automatically once declared in the config.

### 4.4. Clean up

Remove the nf-schema and nf-co2footprint plugins from `nextflow.config` before continuing (they add noise to the output in later exercises):

=== "After"

    ```groovy title="nextflow.config"
    plugins {
        id 'nf-hello@0.5.0'
    }
    ```

=== "Before"

    ```groovy title="nextflow.config"
    plugins {
        id 'nf-hello@0.5.0'
        id 'nf-schema@2.6.1'
        id 'nf-co2footprint'
    }
    ```

Also clean up the generated files:

```bash
rm -f co2footprint_*
```

---

## 5. Discovering plugins

The [Nextflow Plugin Registry](https://registry.nextflow.io/) is the central hub for finding available plugins.

![The nf-hello plugin page on registry.nextflow.io](img/plugin-registry-nf-hello.png)

Each plugin page shows its description, available versions, installation instructions, and links to documentation.

---

## Takeaway

You used three different plugins:

- **nf-hello**: A function plugin providing `randomString`, imported with `include`
- **nf-schema**: A function plugin providing `samplesheetToList` for schema-validated CSV parsing
- **nf-co2footprint**: An observer plugin that monitors resource usage automatically, with no `include` needed

Key patterns:

- Plugins are declared in `nextflow.config` with `plugins { id 'plugin-name@version' }`
- Function plugins require `include { function } from 'plugin/plugin-id'`
- Observer plugins work automatically once declared in the config
- The [Nextflow Plugin Registry](https://registry.nextflow.io/) lists available plugins

---

## What's next?

The following sections show you how to build your own plugin.
If you're not interested in plugin development, you can stop here or skip ahead to the [Summary](summary.md).

[Continue to Part 2 :material-arrow-right:](02_create_project.md){ .md-button .md-button--primary }
