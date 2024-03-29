# Anatomy of the nf-core template

The nf-core pipeline template comes packed with a lot of files and folders.

Fortunately, while the template can feel overwhelming, a complete understanding isn't required to get started.

![Anatomy of the nf-core template](img/4_template.png)

!!! note "Viewing dot files"

    Dot files are not shown using the `ls` command. You can use `ll` or `ls -la` to view all files. Alternatively, `tree -a` command can be used to show the branching structure of the pipeline repository.

## Template tour

The most important folder and files in the nf-core pipeline template are outlined below.

### Workflows, subworkflows, and modules

The nf-core pipeline template has a `main.nf` file that calls a `<workflow>.nf` file from the `workflows` folder. The `<workflow>.nf` file is where the pipeline is built and is used to bring everything else together. Instead of having one large monolithic script, it is broken up into a combination of subworkflows and modules.

A subworkflow is a groups of modules that are used in combination with each other and have a common purpose. Subworkflows improve pipeline readability and help with the reuse of modules within a pipeline. Subworkflow are typically stored as either `nf-core` or `local` subworkflows within the `subworkflows` folder. An nf-core subworkflow is developed by the community and are shared in the nf-core subworkflows GitHub repository while local subworkflows are pipeline specific and are not shared in the nf-core subworkflows repository.

A module is a wrapper for a process, the basic processing primitive to execute a user script. Like subworkflows, nf-core modules are developed by the community and are shared in the nf-core modules GitHub repository while local modules are pipeline specific are not shared in the nf-core subworkflows repository.

<figure class="excalidraw">
--8<-- "docs/nf_templates/img/4_nested.excalidraw.svg"
</figure>

To enable reporting and reproducibility, subworkflows and modules from the nf-core repository are tracked using hashes in the `modules.json` file in the main pipeline repository.

### Configuration files

The nf-core pipeline template utilizes Nextflows flexible customization and has a series of abstracted configuration files throughout the template. Many of which can be flexibly applied using profiles.  

In the template, the `nextflow.config` file is a central configuration file that is used to set pipeline defaults.

There are a series of additional config files that stored in the `conf/` folder which are added to the `nextflow.config` file using `include` statements.

```console "nextflow.config"
// Load base.config by default for all pipelines
includeConfig 'conf/base.config'
```

These configuration files are either included by default (e.g., `base.config` and `modules.config`) or are included as a part of the profiles scope and can be flexibly applied at runtime (e.g., `test.config` and `test_full.config`).

### `nf-core.yml`

The `nf-core.yml` file

### `nextflow_schema.json`

The `nextflow_schema.json` file 

## Changes to the template structure

Occasionally, the structure of the nf-core pipeline template is updated during a new release of the nf-core tooling. Most of the time these changes are minor. However, sometimes, larger structural changes are adopted to align with changes in the wider ecosystem.

For example, as of nf-core tools 2.13, the groovy code that lived in the `lib` folder has been moved to `subworkflows/`. Moving this code has made it easier to find, modify, and test the code. Importantly, it's modular and is paving the way for a more flexible template in the future.

!!! note

    The `TEMPLATE` branch is essential for adopting these changes.
