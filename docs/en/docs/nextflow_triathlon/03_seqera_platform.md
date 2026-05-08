# Part 3: Run on Seqera Platform

In this third and final part of the Nextflow Triathlon, you will use Seqera Platform to launch, monitor, and manage Nextflow pipelines through a web interface.

We will start by running nf-core/rnaseq, a production-scale RNA sequencing pipeline, to see what Platform can do with a real workload.
We will then explore the `tw` CLI to drive Platform from the terminal.
Finally, we will add the nf-core/demo pipeline directly from its GitHub repository — demonstrating how to bring any correctly-configured pipeline onto the Platform.

!!! tip

    You will need access to a Seqera Platform workspace for this part.
    Your instructor will provide login details and workspace information.
    The compute environment has already been configured for you.

---

## 1. Launch nf-core/rnaseq

nf-core/rnaseq is a community pipeline for bulk RNA sequencing analysis.
It runs quality control, adapter trimming, alignment, and quantification — the core steps of most RNA-seq experiments.

In this section, you will add the pipeline to your workspace, prepare input data, launch a run, and monitor its execution.

### 1.1. Add the pipeline to your workspace

Seqera Platform provides a curated collection of nf-core and other community pipelines that can be added to your workspace in a few clicks.

1. In the left sidebar, click **Pipelines**, then click **Add pipeline**.
2. Select **Seqera Pipelines** to browse the community collection.
3. Search for `rnaseq` and select **nf-core/rnaseq**.
4. Choose your workspace from the dropdown and click **Add pipeline**.

The pipeline is now listed in your workspace's **Pipelines** panel and is ready to launch.

### 1.2. Explore the input data

The test dataset for this tutorial consists of paired-end RNA-seq reads from human cell lines, hosted in cloud storage.

Use **Data Explorer** to browse the input files before launching:

1. In the left sidebar, click **Data Explorer**.
2. Navigate to the bucket your instructor pointed you to.
3. Browse the sample FASTQ files and the samplesheet CSV.

The samplesheet is the key input to nf-core/rnaseq.
It maps sample identifiers to their FASTQ files and specifies strandedness:

```csv title="samplesheet.csv"
sample,fastq_1,fastq_2,strandedness
GM12878_rep1,s3://...GM12878_rep1_R1.fastq.gz,s3://...GM12878_rep1_R2.fastq.gz,auto
GM12878_rep2,s3://...GM12878_rep2_R1.fastq.gz,s3://...GM12878_rep2_R2.fastq.gz,auto
K562_rep1,s3://...K562_rep1_R1.fastq.gz,s3://...K562_rep1_R2.fastq.gz,auto
K562_rep2,s3://...K562_rep2_R1.fastq.gz,s3://...K562_rep2_R2.fastq.gz,auto
```

### 1.3. Launch the pipeline

1. In the **Pipelines** panel, click **Launch** next to nf-core/rnaseq.

   The launch form has three sections: **General config**, **Run parameters**, and **Advanced settings**.

2. In **General config**:

   - **Revision**: select the latest stable release (e.g. `3.18.0`)
   - **Config profiles**: add `test` if using the test dataset, or leave blank for your own data
   - **Compute environment**: select the pre-configured environment your instructor specified

3. In **Run parameters**, set:

   - **input**: path to your samplesheet (e.g. `s3://your-bucket/samplesheet.csv`) or select it from Data Explorer
   - **outdir**: path to the output directory (e.g. `s3://your-bucket/rnaseq-results`)
   - **genome**: `GRCh38` (or as specified by your instructor)

4. In **Advanced settings**:

   - **Nextflow version**: set this to match the version required by the pipeline (check the pipeline's `tower.yml` or ask your instructor)
   - **Environment variables**: add `NXF_SYNTAX_PARSER=v1` — nf-core pipelines use the v1 Nextflow syntax and will fail without this on Nextflow 25.10 and later

5. Click **Launch** to start the run.

### 1.4. Monitor execution

After launching, you will be taken to the **Runs** panel for your pipeline.

The run view shows:

- **Status**: current state of the run (submitted, running, succeeded, failed)
- **Command line**: the exact `nextflow run` command that Platform constructed and submitted
- **Parameters**: all parameter values used for this run
- **Tasks**: a table of every process call, with status, duration, and resource usage

Click on any task row to inspect its execution details, including:

- The `.command.sh` script that was run
- stdout and stderr logs
- CPU, memory, and I/O metrics

The **Reports** tab will show a MultiQC report once the run completes, aggregating quality control metrics across all samples.

!!! note

    This is a production-scale run — it may take longer than the exercises in Parts 1 and 2.
    Your instructor will let you know whether to wait for completion or move on while it runs in the background.

### Takeaway

You know how to add a pipeline to a Seqera Platform workspace, configure and launch a run, and monitor execution at scale.

### What's next?

Learn how to do all of this from the command line using the `tw` CLI.

---

## 2. Launch pipelines from the command line

In the run view, click the **Command line** tab.
You will see the exact `nextflow run` command that Platform constructed and submitted on your behalf — the same kind of command you have been running manually throughout this course.

Platform does not replace Nextflow; it orchestrates it.
Everything you can do through the web interface, you can also do from the terminal using the `tw` CLI, the command-line tool for interacting with the Platform API.
This is useful for automating launches from scripts or CI/CD pipelines.

### 2.1. Get an access token

The `tw` CLI authenticates with Platform using a personal access token.

1. In the Platform web interface, click your avatar in the top-right corner and select **Your tokens**.
2. Click **Add token**, give it a name (e.g. `training`), and click **Add**.
3. Copy the token value — it will only be shown once.

### 2.2. Configure the CLI

The `tw` CLI is pre-installed in the training environment.
Set your access token as an environment variable:

```bash
export TOWER_ACCESS_TOKEN=<your-token>
```

Verify the connection:

```bash
tw info
```

??? success "Command output"

    ```console
    System info:
      Tower API endpoint  : https://api.cloud.seqera.io
      Tower API version   : 1.25.0
      Tower version       : 24.2.0
      Authenticated as    : Your Name (your@email.com)
    ```

### 2.3. Explore your workspace from the CLI

List the workspaces you have access to:

```bash
tw workspaces list
```

??? success "Command output"

    ```console
    Available workspaces:
    ID              | Name            | Full Name                | Visibility
    --------------- | --------------- | ------------------------ | ----------
    <workspace-id>  | my-workspace    | my-org/my-workspace      | PRIVATE
    ```

View the runs in your workspace, including the nf-core/rnaseq run you just launched:

```bash
tw runs list --workspace <org>/<workspace>
```

??? success "Command output"

    ```console
    Pipeline runs for my-org/my-workspace workspace:
    ID        | Status   | Name          | Pipeline               | Run name
    --------- | -------- | ------------- | ---------------------- | --------
    <run-id>  | RUNNING  | ...           | nf-core/rnaseq         | happy_curie
    ```

The same run you are monitoring in the web interface is visible here.

### 2.4. Launch nf-core/rnaseq from the CLI

Launch nf-core/rnaseq using the same parameters as section 1, this time from the terminal:

```bash
tw launch nf-core/rnaseq \
    --workspace <org>/<workspace> \
    --compute-env <compute-env-name> \
    --revision 3.18.0 \
    --env NXF_SYNTAX_PARSER=v1 \
    --params '{"input": "s3://your-bucket/samplesheet.csv", "outdir": "s3://your-bucket/rnaseq-results-cli", "genome": "GRCh38"}'
```

Replace `<org>/<workspace>`, `<compute-env-name>`, and the parameter values with those your instructor provided.

??? success "Command output"

    ```console
    Launching pipeline nf-core/rnaseq [3.18.0]
    Run name: focused_einstein
    https://cloud.seqera.io/orgs/my-org/workspaces/my-workspace/watch/<run-id>
    ```

The CLI prints a direct link to the run.
Open it in your browser and confirm the run appears in the **Runs** panel with a **submitted** or **running** status.

You do not need to wait for this run to complete — once you can see it running, you have confirmed that the CLI and the web interface are two views onto the same Platform.
Move on to section 3 while it runs in the background.

### Takeaway

You know how to authenticate the `tw` CLI, inspect your workspace, and launch pipelines from the terminal with the same options available in the web interface.

### What's next?

Add nf-core/demo from its GitHub repository to your workspace via the web interface, with launch defaults configured so the pipeline is ready for your whole team to use.

---

## 3. Add a pipeline from GitHub and run it

Not every pipeline is available in the Seqera Pipelines catalog.
Any Nextflow pipeline hosted in a GitHub repository can be added to the Platform directly — as long as the repository has a `main.nf` entry point and a `nextflow.config` at its root.

nf-core/demo is a good example to practice with: you already ran it from the command line in Part 2, so you know what it does and what outputs to expect.

### 3.1. What makes a repository compatible

For the Platform to launch a pipeline from GitHub, the repository needs:

- **`main.nf`** (or a named entry point) at the repository root
- **`nextflow.config`** at the repository root
- A valid GitHub URL that the Platform can access (public repo, or private with credentials configured)

nf-core pipelines always meet these requirements.
The same is true for well-structured community pipelines, and for your own pipelines if you follow the same conventions.

### 3.2. Add nf-core/demo to your workspace

1. In the left sidebar, click **Pipelines**, then click **Add pipeline**.
2. Select **GitHub repository** instead of Seqera Pipelines.
3. Enter the repository URL: `https://github.com/nf-core/demo`
4. Select a **Revision** (branch, tag, or commit) — use `master` or the latest release tag.
5. Choose your workspace and compute environment.
6. Click **Add pipeline**.

The pipeline will appear in your workspace's **Pipelines** panel.

### 3.3. Configure launch defaults

Before launching, set default values that will pre-fill the launch form for future runs.
This is useful when you want to standardize how a pipeline is run within your team.

1. Click the pipeline name to open its settings, then click **Edit**.
2. Under **Pre-run script** or **Run parameters**, set:
   - **Config profiles**: `docker,test`
   - **outdir**: `demo-results`
3. Save the changes.

### 3.4. Launch and monitor the run

1. Click **Launch** next to nf-core/demo.
2. Confirm or adjust the pre-filled parameters and click **Launch**.
3. Watch the run progress in the **Runs** panel.

Because you are using the `test` profile, this run completes quickly.
Once it finishes, click into the run to explore the task table and execution reports — the same views you saw for the rnaseq run, at a smaller scale.

### Takeaway

You know how to add any GitHub-hosted Nextflow pipeline to Seqera Platform, configure its launch defaults, and run it.

---

## Summary

In this part you learned to:

- Add a pipeline from the curated Seqera Pipelines catalog, configure it, and launch a production-scale run
- Prepare input data using Data Explorer
- Monitor task-level execution and review pipeline reports
- Authenticate the `tw` CLI and drive Platform from the terminal
- Add a pipeline directly from a GitHub repository and configure its launch defaults for team use
