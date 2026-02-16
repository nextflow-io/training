# Part 1: Plugin Basics

In this section, you'll use an existing plugin in a Nextflow workflow, then learn how plugins extend Nextflow.

---

## 1. Use a plugin in a workflow

The [nf-hello](https://github.com/nextflow-io/nf-hello) plugin provides a `randomString` function that generates random strings.
You'll replace a local function with the plugin version to see how plugins work in practice.

### 1.1. See the starting point

The `random_id_example.nf` file contains a workflow with a local `randomString` function:

```bash
cat random_id_example.nf
```

Run it to see the output:

```bash
nextflow run random_id_example.nf
```

This works, but the function is trapped in this file.
The plugin version can be shared across any pipeline.

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

The first run downloads the plugin automatically.
Any pipeline using `nf-hello@0.5.0` gets the exact same `randomString` function.

The development burden is on the plugin developer, not the pipeline developer.
Nextflow handles installing and updating plugins on your behalf.

---

## 2. How plugins work

Plugins extend Nextflow through well-defined extension points:

| Extension Type  | Purpose                                  | Example                 |
| --------------- | ---------------------------------------- | ----------------------- |
| Functions       | Custom functions callable from workflows | `samplesheetToList()`   |
| Executors       | Custom task execution backends           | AWS Batch, Kubernetes   |
| Filesystems     | Custom storage backends                  | S3, Azure Blob          |
| Trace Observers | Monitor workflow execution               | Custom logging, metrics |

Functions and trace observers are the most common types for plugin authors.
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
- Plugins extend Nextflow through functions, observers, executors, and filesystems

---

## What's next?

The following sections show you how to build your own plugin.
If you're not interested in plugin development, you can stop here or skip ahead to the [Summary](summary.md).

[Continue to Part 2 :material-arrow-right:](02_create_project.md){ .md-button .md-button--primary }
