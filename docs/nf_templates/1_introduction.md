# Adopting the nf-core pipeline template

Workflow management systems, such as Nextflow, provide a powerful way to develop pipelines. However, they do not provide standards for how pipelines should be developed. This gap has led to the establishment of pipeline registries, such as [nf-core](https://nf-co.re/), with implementation guidelines that provide standards for pipeline development.

nf-core developed a base template for pipelines that comes pre-packed with features that support the adoption of good coding practices. The template establishes a pipeline file structure, with code, documentation, and continuous integration (CI) tests. The template plugs into the nf-core tooling, offering a full suite of commands that support pipeline usage and development.

As the Nextflow ecosystem evolves, so does the base template. New nf-core tool releases frequently come with changes that support the adoption of new features (e.g., new plugins and operators). Updates can be semi-automatically synchronised using the nf-core tooling. Leveraging the nf-core template enables the maintenance and growth of pipelines within the wider ecosystem. 

This workshop will explore the nf-core template and how it supports the adoption of best practises.

## Learning objectives

By the end of this workshop you will be able to:

- Discuss the benefits of using the nf-core pipeline template
- Identify the features of nf-core pipeline template that enable best practises
- Utilize nf-core tooling to accelerate your pipeline development

## Audience & prerequisites

Please note that this is not a beginners workshop and some familiarity with Nextflow, the command line, and common file formats is assumed.

**Prerequisites**

-   A GitHub account
-   Experience with command line
-   Experience writing pipelines with Nextflow
