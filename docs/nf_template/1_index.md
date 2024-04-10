# Leveraging the nf-core template

Nextflow provide a powerful way to develop pipelines. However, it does not provide standards for how pipelines should be developed. This gap has led to the establishment of pipeline registries, such as [nf-core](https://nf-co.re/), with tools and implementation guidelines that provide support and standards for pipeline development.

To enable standardisation and the adoption of good coding practices nf-core developed a base template for pipelines. The template establishes a pipeline file structure, with code, documentation, and continuous integration (CI) tests. The template plugs into the nf-core tooling, offering a full suite of commands that support pipeline usage and development.

Leveraging the nf-core template enables accelerated maintenance and development of pipelines.

## Learning objectives

In this workshop you will create a pipeline using the nf-core template and utilize the nf-core tooling to build a pipeline using best practices.

By the end of this workshop you will be able to:

-   Discuss the benefits of using the nf-core pipeline template
-   Identify the features of nf-core pipeline template that enable best practises
-   Utilize nf-core tooling to accelerate your pipeline development

## Audience & prerequisites

Please note that this is **not** a beginners workshop and familiarity with Nextflow, the command line, and common file formats is assumed.

**Prerequisites**

-   A GitHub account
-   Experience with command line
-   Experience writing pipelines with Nextflow
