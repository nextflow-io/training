# Part 1: Plugin Basics

In this section, you'll use an existing plugin in a Nextflow workflow, then learn how plugins extend Nextflow.

---

## 1. Use a plugin in a workflow

The [nf-hello](https://github.com/nextflow-io/nf-hello) plugin provides a `randomString` function that generates random strings.
In this exercise, you'll start with a workflow that defines this function locally, then replace it with the plugin version.
This demonstrates the core plugin workflow: declare a plugin, import its functions, and use them like any other Nextflow function.

### 1.1. See the starting point

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

### 1.2. Configure the plugin

Add the plugin to your `nextflow.config`:

```groovy title="nextflow.config"
// Configuration for plugin development exercises
plugins {
    id 'nf-hello@0.5.0'
}
```

Plugins are declared in `nextflow.config` using the `plugins {}` block.
Nextflow automatically downloads them from the registry at runtime.

### 1.3. Use the plugin function

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

### 1.4. Run it

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

## 2. How plugins work

Plugins extend Nextflow through several types of extension:

| Extension Type    | What it does                                 | Example                      |
| ----------------- | -------------------------------------------- | ---------------------------- |
| Functions         | Add custom functions callable from workflows | `samplesheetToList()`        |
| Workflow monitors | Respond to events like task completion       | Custom logging, Slack alerts |
| Executors         | Add task execution backends                  | AWS Batch, Kubernetes        |
| Filesystems       | Add storage backends                         | S3, Azure Blob               |

Functions and workflow monitors (called "trace observers" in the Nextflow API) are the most common types for plugin authors.
Executors and filesystems are typically created by platform vendors.

---

## 3. Discovering plugins

The [Nextflow Plugin Registry](https://registry.nextflow.io/) is the central hub for finding available plugins.

![The nf-hello plugin page on registry.nextflow.io](img/plugin-registry-nf-hello.png)

Each plugin page shows its description, available versions, installation instructions, and links to documentation.

---

## Takeaway

You learned that:

- Plugins are declared in `nextflow.config` with `plugins { id 'plugin-name@version' }`
- Import plugin functions with `include { function } from 'plugin/plugin-id'`
- The [Nextflow Plugin Registry](https://registry.nextflow.io/) lists available plugins
- Plugins extend Nextflow through functions, workflow monitors, executors, and filesystems

---

## What's next?

The following sections show you how to build your own plugin.
If you're not interested in plugin development, you can stop here or skip ahead to the [Summary](summary.md).

[Continue to Part 2 :material-arrow-right:](02_create_project.md){ .md-button .md-button--primary }
