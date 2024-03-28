# Creating a pipeline using the nf-core template

The nf-core pipeline template is a standardized framework designed to streamline the development of Nextflow-based bioinformatics pipelines.

The pipeline template can be generated using nf-core tools and the `nf-core create` command:

```bash
nf-core create
```

Although you can provide options on the command line, it’s easiest to use the interactive prompts.

```bash
? Workflow name <pipeline name>
? Description <pipeline description>
? Author <your name>
? Do you want to customize which parts of the template are used? (y/N) n
```

If you do not customize the template, all features of the template will be included by default.

## Customizing the template

There is also flexibility for which template areas you include in your template.

The following template areas can be skipped:

-   **GitHub hosting:** Files required for GitHub hosting of the pipeline. E.g., `.github/` and `.gitignore`.
-   **GitHub CI:** Files required for GitHub continuous integration tests. E.g., `.github/workflows/`.
-   **GitHub badges:** GitHub badges in the `README.md` file.
-   **iGenomes config:** Pipeline options related to iGenomes. E.g., `conf/igenomes.config`.
-   **nf-core/configs:** Repository options that integrate nf-core config profiles.

If you choose to modify the template, the nf-core tooling will provide a series of interactive prompts help guide your choices.

```bash
? Workflow name mypipeline
? Description My pipeline
? Author Chris
? Do you want to customize which parts of the template are used? (y/N) y
? Pipeline prefix myorg
Skip template areas?
   ○ GitHub hosting
   ○ GitHub CI
   ● GitHub badges
   ○ iGenomes config
   ○ nf-core/configs
```

!!! note

    GitHub badges were skipped in the example above.

!!! question "Exercise"

    Create a pipeline named `mypipeline` using the `nf-core create` command. Do **not** modify the template.
