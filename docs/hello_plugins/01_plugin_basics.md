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

The next exercises show you function plugins and an observer plugin, so you can see both types in action.

---

## 2. Use function plugins

Function plugins add callable functions that you import into your workflows.
You'll try two: nf-hello (a simple example) and nf-schema (a widely-used real-world plugin).
Both exercises modify the same `hello.nf` pipeline, so you can see how plugins enhance an existing workflow.

### 2.1. nf-hello: replace hand-written code

The [nf-hello](https://github.com/nextflow-io/nf-hello) plugin provides a `randomString` function that generates random strings.
The pipeline already defines its own inline version of this function, which you'll replace with the one from the plugin.

#### 2.1.1. See the starting point

Look at the pipeline:

```bash
cat hello.nf
```

```groovy title="Output"
#!/usr/bin/env nextflow

params.input = 'greetings.csv'

/**
 * Generate a random alphanumeric string
 */
def randomString(int length) {
    def chars = ('a'..'z') + ('A'..'Z') + ('0'..'9')
    def random = new Random()
    return (1..length).collect { chars[random.nextInt(chars.size())] }.join()
}

process SAY_HELLO {
    input:
        val greeting
    output:
        stdout
    script:
    """
    echo '$greeting'
    """
}

workflow {
    greeting_ch = channel.fromPath(params.input)
                        .splitCsv(header: true)
                        .map { row -> "${row.greeting}_${randomString(8)}" }
    SAY_HELLO(greeting_ch)
    SAY_HELLO.out.view { result -> "Output: ${result.trim()}" }
}
```

The pipeline defines its own `randomString` function inline, then uses it to append a random ID to each greeting.

Run it:

```bash
nextflow run hello.nf
```

```console title="Output"
Output: Hello_aBcDeFgH
Output: Bonjour_xYzWvUtS
Output: Holà_qRsPdMnK
Output: Ciao_jLhGfEcB
Output: Hallo_tNwOiAuR
```

Your output order and random strings will differ, and if you run the script again you'll get a different set of random greetings.

#### 2.1.2. Configure the plugin

Now let's replace the inline function with one from a plugin. Add this plugin to your `nextflow.config`:

```groovy title="nextflow.config"
// Configuration for plugin development exercises
plugins {
    id 'nf-hello@0.5.0'
}
```

Plugins are declared in `nextflow.config` using the `plugins {}` block.
Nextflow automatically downloads them from the [Nextflow Plugin Registry](https://registry.nextflow.io/), a central repository of community and official plugins.

#### 2.1.3. Use the plugin function

Replace the inline `randomString` function with the plugin version:

=== "After"

    ```groovy title="hello.nf" hl_lines="3"
    #!/usr/bin/env nextflow

    include { randomString } from 'plugin/nf-hello'

    params.input = 'greetings.csv'

    process SAY_HELLO {
        input:
            val greeting
        output:
            stdout
        script:
        """
        echo '$greeting'
        """
    }

    workflow {
        greeting_ch = channel.fromPath(params.input)
                            .splitCsv(header: true)
                            .map { row -> "${row.greeting}_${randomString(8)}" }
        SAY_HELLO(greeting_ch)
        SAY_HELLO.out.view { result -> "Output: ${result.trim()}" }
    }
    ```

=== "Before"

    ```groovy title="hello.nf" hl_lines="5-12"
    #!/usr/bin/env nextflow

    params.input = 'greetings.csv'

    /**
     * Generate a random alphanumeric string
     */
    def randomString(int length) {
        def chars = ('a'..'z') + ('A'..'Z') + ('0'..'9')
        def random = new Random()
        return (1..length).collect { chars[random.nextInt(chars.size())] }.join()
    }

    process SAY_HELLO {
        input:
            val greeting
        output:
            stdout
        script:
        """
        echo '$greeting'
        """
    }

    workflow {
        greeting_ch = channel.fromPath(params.input)
                            .splitCsv(header: true)
                            .map { row -> "${row.greeting}_${randomString(8)}" }
        SAY_HELLO(greeting_ch)
        SAY_HELLO.out.view { result -> "Output: ${result.trim()}" }
    }
    ```

The `include` statement replaces 7 lines of hand-written code with a single import from a tested, versioned plugin.
The syntax `#!groovy include { function } from 'plugin/plugin-id'` is the same `include` used for Nextflow modules, with a `plugin/` prefix.
You can see the [source code for `randomString`](https://github.com/nextflow-io/nf-hello/blob/e67bddebfa589c7ae51f41bf780c92068dc09e93/plugins/nf-hello/src/main/nextflow/hello/HelloExtension.groovy#L110) in the nf-hello repository on GitHub.

#### 2.1.4. Run it

```bash
nextflow run hello.nf
```

```console title="Output"
Pipeline is starting! 🚀
Output: Hello_yqvtclcc
Output: Bonjour_vwwpyzcs
Output: Holà_wrghmgab
Output: Ciao_noniajuy
Output: Hallo_tvrtuxtp
Pipeline complete! 👋
```

(Your random strings will differ.)

The output still has random suffixes, but now `randomString` comes from the nf-hello plugin instead of inline code.
The "Pipeline is starting!" and "Pipeline complete!" messages come from nf-hello's built-in observer, which demonstrates that a single plugin can provide both functions and observers.

Nextflow downloads plugins automatically the first time they're used, so any pipeline that declares `nf-hello@0.5.0` gets the exact same tested `randomString` function without copying code between projects.

You've now seen the three steps for using a function plugin: declare it in `nextflow.config`, import the function with `include`, and call it in your workflow.
The next exercise applies these same steps to a real-world plugin.

### 2.2. nf-schema: validated CSV parsing

The [nf-schema](https://github.com/nextflow-io/nf-schema) plugin is one of the most widely-used Nextflow plugins.
It provides `samplesheetToList`, a function that parses CSV/TSV files using a JSON schema that defines the expected columns and types.

The pipeline currently reads `greetings.csv` using `splitCsv` and a manual `map`, but nf-schema can replace this with validated, schema-driven parsing.
A JSON schema file (`greetings_schema.json`) is already provided in the exercise directory.

#### 2.2.1. Look at the schema

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

#### 2.2.2. Add nf-schema to the config

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

#### 2.2.3. Update hello.nf to use samplesheetToList

Replace the `splitCsv` input with `samplesheetToList`:

=== "After"

    ```groovy title="hello.nf" hl_lines="4 20 21"
    #!/usr/bin/env nextflow

    include { randomString } from 'plugin/nf-hello'
    include { samplesheetToList } from 'plugin/nf-schema'

    params.input = 'greetings.csv'

    process SAY_HELLO {
        input:
            val greeting
        output:
            stdout
        script:
        """
        echo '$greeting'
        """
    }

    workflow {
        greeting_ch = Channel.fromList(samplesheetToList(params.input, 'greetings_schema.json'))
                            .map { row -> "${row[0]}_${randomString(8)}" }
        SAY_HELLO(greeting_ch)
        SAY_HELLO.out.view { result -> "Output: ${result.trim()}" }
    }
    ```

=== "Before"

    ```groovy title="hello.nf" hl_lines="19 20 21"
    #!/usr/bin/env nextflow

    include { randomString } from 'plugin/nf-hello'

    params.input = 'greetings.csv'

    process SAY_HELLO {
        input:
            val greeting
        output:
            stdout
        script:
        """
        echo '$greeting'
        """
    }

    workflow {
        greeting_ch = channel.fromPath(params.input)
                            .splitCsv(header: true)
                            .map { row -> "${row.greeting}_${randomString(8)}" }
        SAY_HELLO(greeting_ch)
        SAY_HELLO.out.view { result -> "Output: ${result.trim()}" }
    }
    ```

Instead of `splitCsv` and a manual `map` to extract fields, `samplesheetToList` parses the CSV according to the schema.
Each row becomes a list of values in column order, so `row[0]` is the greeting and `row[1]` is the language.

#### 2.2.4. Run it

```bash
nextflow run hello.nf
```

```console title="Output"
Pipeline is starting! 🚀
Output: Hello_diozjdwm
Output: Bonjour_speathmm
Output: Holà_dllxnzap
Output: Ciao_wzueddzc
Output: Hallo_hsxwrjbh
Pipeline complete! 👋
```

(Your random strings will differ.)

The output is the same, but now the schema validates the CSV structure before the pipeline runs.
In real pipelines with complex sample sheets and many columns, this kind of validation prevents errors that manual `splitCsv` + `map` would miss.

Both nf-hello and nf-schema are function plugins: they provide functions that you import with `include` and call in your workflow code.
The next exercise shows a different type of plugin that works without any `include` statements at all.

---

## 3. Use an observer plugin: nf-co2footprint

Not all plugins provide functions to import.
The [nf-co2footprint](https://github.com/nextflow-io/nf-co2footprint) plugin uses a **trace observer** to monitor your pipeline's resource usage and estimate its carbon footprint.
You don't need to change any pipeline code; just add it to the config.

### 3.1. Add nf-co2footprint to the config

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

### 3.2. Run the pipeline

```bash
nextflow run hello.nf
```

You'll see additional log messages and warnings from the plugin during execution, which are normal; the plugin is estimating resource usage with limited information from this small example.

At the end, look for a line like:

```console title="Output (partial)"
🌱 The workflow run used 163.4 uWh of electricity, resulting in the release of 78.43 ug of CO₂ equivalents into the atmosphere.
```

(Your numbers will differ.)

### 3.3. View the report

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

This plugin works entirely through the observer mechanism, hooking into workflow lifecycle events to collect resource metrics and generate its report when the pipeline completes.
No `include` statement is needed because it doesn't provide functions; it runs automatically once declared in the config.

You've now tried function plugins (imported with `include`) and an observer plugin (activated through config alone).
These are the two most common extension types, but as the table in section 1 shows, plugins can also add executors and filesystems.

---

## 4. Discovering plugins

The [Nextflow Plugin Registry](https://registry.nextflow.io/) is the central hub for finding available plugins.

![The nf-hello plugin page on registry.nextflow.io](img/plugin-registry-nf-hello.png)

Each plugin page shows its description, available versions, installation instructions, and links to documentation.

---

## 5. Prepare for plugin development

The following sections (Parts 2-6) use a separate pipeline file, `greet.nf`, which relies on nf-schema but not nf-hello or nf-co2footprint.

Update `nextflow.config` to keep only nf-schema:

```groovy title="nextflow.config"
// Configuration for plugin development exercises
plugins {
    id 'nf-schema@2.6.1'
}
```

Remove the co2footprint output files:

```bash
rm -f co2footprint_*
```

The `hello.nf` file retains your Part 1 work for reference; going forward, you'll work with `greet.nf`.

---

## Takeaway

You used three different plugins:

- **nf-hello**: A function plugin providing `randomString`, imported with `include`
- **nf-schema**: A function plugin providing `samplesheetToList` for schema-validated CSV parsing
- **nf-co2footprint**: An observer plugin that monitors resource usage automatically, with no `include` needed

Key patterns:

- Plugins are declared in `nextflow.config` with `#!groovy plugins { id 'plugin-name@version' }`
- Function plugins require `#!groovy include { function } from 'plugin/plugin-id'`
- Observer plugins work automatically once declared in the config
- The [Nextflow Plugin Registry](https://registry.nextflow.io/) lists available plugins

---

## What's next?

The following sections show you how to build your own plugin.
If you're not interested in plugin development, you can stop here or skip ahead to the [Summary](summary.md).

[Continue to Part 2 :material-arrow-right:](02_create_project.md){ .md-button .md-button--primary }
