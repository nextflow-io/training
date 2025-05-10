# Part 1: Run nf-core/demo

In this first part of the Hello nf-core training course, we show you how to find and try out an nf-core pipeline, understand how the code is organized, and recognize how it differs from plain Nextflow code as shown in [Hello Nextflow](../hello_nextflow/index.md).

We are going to use a pipeline called nf-core/demo that is maintained by the nf-core project as part of its inventory of pipelines for demonstrating code structure and tool operations.

---

## 0. Warmup

Before we go looking for the pipeline, let's create a project directory where we're going to do the work.

Make sure you are in the `hello-nf-core/` directory as instructed in the [Orientation](./00_orientation.md), then create the directory as follows:

```bash
mkdir demo_run/
cd demo_run
```

<!-- TODO: Add commentary about what will be stored/created where -->

---

## 1. Find and retrieve the nf-core/demo pipeline

Let's start by locating the nf-core/demo pipeline on the project website at [nf-co.re](https://nf-co.re), which centralizes all information such as: general documentation and help articles, documentation for each of the pipelines, blog posts, event announcements and so forth.

### 1.1. Find the pipeline on the website

In your web browser, go to https://nf-co.re/pipelines/ and type `demo` in the search bar.

![search results](./img/search-results.png)

Click on the pipeline name, `demo`, to access the pipeline details page.

Each released pipeline has a dedicated page that includes the following documentation sections:

- **Introduction:** An introduction and overview of the pipeline
- **Usage:** Descriptions of how to execute the pipeline
- **Parameters:** Grouped pipeline parameters with descriptions
- **Output:** Descriptions and examples of the expected output files
- **Results:** Example output files generated from the full test dataset
- **Releases & Statistics:** Pipeline version history and statistics

Whenever you are considering adopting a new pipeline, you should read the pipeline documentation carefully first to understand what it does and how it should be configured before attempting to run it.

Have a look now and see if you can find out:

- which tools the pipeline will run (Check the tab: `Introduction`)
- which inputs and parameters the pipeline accepts or requires (Check the tab: `Parameters`)
- what are the outputs produced by the pipeline (Check the tab: `Output`)

  The `Introduction` tab provides an overview of the pipeline, including a visual representation (called a subway map) and a list of tools that are run as part of the pipeline.

  ![pipeline subway map](./img/nf-core-demo-subway-cropped.png)

  1. Read QC (FASTQC)
  2. Adapter and quality trimming (SEQTK_TRIM)
  3. Present QC for raw reads (MULTIQC)

  The documentation also provides an example input file (see below) and an example command line.

  ```bash
  nextflow run nf-core/demo \
   -profile <docker/singularity/.../institute> \
   --input samplesheet.csv \
   --outdir <OUTDIR>
  ```

You'll notice that the example command does NOT specify a workflow file, just the reference to the pipeline repository, `nf-core/demo`.

When invoked this way, Nextflow will assume that the code is organized in a certain way.
Let's retrieve the code so we can examine this structure.

### 1.2. Retrieve the pipeline code

Once we've determined the pipeline appears to be suitable for our purposes, we're going to want to try it out.
Fortunately Nextflow makes it easy to retrieve pipeline from correctly-formatted repositories without having to download anything manually.

Return to your terminal and run the following:

```bash
nextflow pull nf-core/demo
```

Nextflow will `pull` the pipeline code, meaning it will download the full repository to your local drive.

```console title="Output"
Checking nf-core/demo ...
 downloaded from https://github.com/nf-core/demo.git - revision: 04060b4644 [master]
```

By default, the files will be saved to `$HOME/.nextflow/assets`.
This can be customized using the `NXF_ASSETS` environment variable; see the configuration documentation at https://www.nextflow.io/docs/latest/config.html.

!!! note
To be clear, you can do this with any Nextflow pipeline that is appropriately set up in GitHub, not just nf-core pipelines.
However nf-core is the largest open-source collection of Nextflow pipelines.

### Takeaway

You now know how to find a pipeline via the nf-core website and retrieve a local copy of the source code.

### What's next?

Explore how the code is organized and why.

---

## 2. Examine the pipeline code structure

Now that we've retrieved the pipeline's source code, let's have a look at how it's organized.

```bash
tree -L 1 $HOME/.nextflow/assets/nf-core/demo
```

```console title="Output"
/root/.nextflow/assets/nf-core/demo
├── assets
├── CHANGELOG.md
├── CITATIONS.md
├── CODE_OF_CONDUCT.md
├── conf
├── docs
├── LICENSE
├── main.nf
├── modules
├── modules.json
├── nextflow_schema.json
├── nextflow.config
├── README.md
├── subworkflows
├── tower.yml
└── workflows
```

### 2.1. Examine the pipeline code structure

Now that we've retrieved the pipeline's source code, let's have a look at how it's organized.

You can view the files

### Takeaway

You know how to [...].

### What's next?

Learn how to [...].

---

## 3. Run the nf-core/demo pipeline

TODO: instructions

### Takeaway

You know how to [...].

### What's next?

[...].

---

!!! note
nf-core tools are pre-installed for you in our training environment.
If you are using a different environment, you need to install the nf-core tools package as described here: https://nf-co.re/docs/nf-core-tools/installation.
