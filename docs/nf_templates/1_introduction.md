# Adopting the nf-core pipeline template

Workflow management systems, such as Nextflow, provide a powerful way to develop pipelines. However, they do not provide standards for how pipelines should be developed. This gap has led to the establishment of pipeline registries, such as [nf-core](https://nf-co.re/), with implementation guidelines that provide standards for pipeline development.

nf-core developed a base template for pipelines that comes pre-packed with features that support the adoption of good coding practices. It receives regular updates that can be semi-automatically synchronised. The template establishes a pipeline file structure, with code, documentation, and continuous integration (CI) tests. The template plugs into the nf-core tooling, offering a full suite of commands that support pipeline development.

As the Nextflow ecosystem evolves, so does the nf-core template. New nf-core tool releases frequently come with changes that support the adoption of new features (e.g., new plugins and operators). Leveraging the nf-core template enables the maintenance and grow of pipelines within the wider ecosystem.

This workshop will explore the nf-core template and how it supports the adoption of best practises.

## Audience & prerequisites

Please note that this is not a beginners workshop and some familiarity with Nextflow, the command line, and common file formats is assumed.

**Prerequisites**

-   A GitHub account
-   Experience with command line
-   Experience writing pipelines with Nextflow

## Gitpod

This workshop can be completed using the [Nextflow Training Gitpod Environment](https://gitpod.io/#https://github.com/nextflow-io/training).

Gitpod offers a virtual machine with everything already set up for you, accessible from your web browser or built into your code editor (e.g., VSCode). To start, click on the button below.

[![Open in GitPod](https://img.shields.io/badge/Gitpod-%20Open%20in%20Gitpod-908a85?logo=gitpod)](https://gitpod.io/#https://github.com/nextflow-io/training)

!!! question "Exercise"

    Open the [Gitpod training environment](https://gitpod.io/#https://github.com/nextflow-io/training) and use the following command to switch to the empty `nf-templates` folder:

    ```bash
    cd /workspace/gitpod/docs/nf_templates
    ```
