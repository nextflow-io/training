# Adopting nf-core pipeline template

Workflow management systems, such as Nextflow, provide a powerful way to develop pipelines. However, they do not provide standards for how pipelines should be developed. This gap has led to the establishment of pipeline registries, such as [nf-core](https://nf-co.re/), with implementation guidelines that provide standards for pipeline development.

nf-core developed a base template for pipelines that comes pre-packed with features that support the adoption of good coding practices. It receives regular updates that can be semi-automatically synchronised. The template establishes a pipeline file structure, with code, documentation, and continuous integration (CI) tests.

Pipelines are tested using automated linting tests for code formatting and syntax, and flags outdated code. To help guarantee deployability, pipelines can be configured with a minimal test dataset that is used to test the pipeline for every pull request, and options to test a full-size dataset for each stable release. The template plugs into the nf-core tooling, offering a full suite of commands that support pipeline development.

As the Nextflow ecosystem evolves, so does the nf-core template. New nf-core tool releases frequently come with changes that support the adoption of new features (e.g., new plugins and operators). Leveraging the nf-core template enables the maintenance and grow of pipelines within the wider ecosystem.

This workshop will explore the nf-core template and how it supports the adoption of best practises.

## Gitpod

This workshop can be completed using the [![Nextflow Training GitPod](https://img.shields.io/badge/Gitpod-%20Open%20in%20Gitpod-908a85?logo=gitpod)](https://gitpod.io/#https://github.com/nextflow-io/training), which provides a virtual machine with everything you need.

!!! question "Exercise"

    Use the following command to switch to the empty `nf-templates` folder:

    ```bash
    cd /workspace/gitpod/nf-templates
    ```
