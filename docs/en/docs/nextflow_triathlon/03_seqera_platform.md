# Part 3: Run on Seqera Cloud

In this third and final part of the Nextflow Triathlon, you will use Seqera to launch, monitor, and manage Nextflow pipelines through a web interface.

We will start by setting up access to Seqera and exploring what's available.
We will then launch the same nf-core/rnaseq pipeline we ran in Part 2, first from the web interface, then again from the terminal.
Finally, we will add the nf-core/demo pipeline directly from its GitHub repository, demonstrating how to bring any correctly-configured pipeline into a Seqera workspace.

---

## 1. Get started with Seqera

Seqera provides a comprehensive platform for launching, monitoring, and managing Nextflow pipelines.
This section walks you through signing up and getting oriented before running your first pipeline.

### 1.1. Sign up for a free account

Go to [cloud.seqera.io](https://cloud.seqera.io) and create a free account.
You can sign up using your email address, GitHub, or Google credentials.

A free account gives you:

- **Personal workspace**: your own space to add pipelines, configure compute environments, and manage runs
- **Access to the Community Showcase**: a curated collection of nf-core and community pipelines with pre-configured settings and example run data

See the [Seqera Platform documentation](https://docs.seqera.io) for a full overview of account tiers and available features.

### 1.2. Explore the Community Showcase

Before launching your own pipelines, take a few minutes to explore the Community Showcase.
It gives you a realistic preview of what the Platform looks like with real pipelines and data.

1. Log in at [cloud.seqera.io](https://cloud.seqera.io).
2. In the left sidebar, click **Showcase**.
3. Browse the available pipelines — you will recognize several nf-core pipelines from Part 2.
4. Click a pipeline to view its configuration and launch settings.
5. Click **Runs** to explore example run histories, including task-level details and reports from previous executions.

This is a read-only view, but it shows you how the interface works before you run anything yourself.

### 1.3. Access a workspace with compute

Launching pipelines requires a workspace with a configured compute environment.
Seqera Platform supports two ways to provide compute:

- **Connect your own infrastructure**: AWS, Azure, Google Cloud, and HPC schedulers (SLURM, LSF, PBS, and others).
  See the [compute environments documentation](https://docs.seqera.io) for setup guides.
- **Seqera Compute**: a managed service that provides pre-provisioned compute environments on AWS, for a fee, with no cloud account setup required.
  You can activate it directly from your workspace settings.

**Group training:**
If you are attending a group training session, you may have been added to an organization and workspace that already has compute configured.
Your instructor will give you the organization name, workspace name, and any other details you need.

**Working independently:**
If you are working through this training by yourself, you will need to set up a compute environment in your personal workspace using one of the options above.
Free credits to try out Seqera Compute are [available on request](https://seqera.io/platform/compute/).

!!! note

    The rest of this part assumes you have access to a workspace with a configured compute environment.
    If you are in a group training session, your instructor will confirm which workspace and compute environment to use.

### Takeaway

You have a Seqera Platform account, you've explored the Community Showcase, and you're able to access a workspace with compute.

### What's next?

Launch a production-scale RNA-seq pipeline from the Seqera Cloud web interface.

---

## 2. Launch nf-core/rnaseq from the web interface

The nf-core/rnaseq pipeline is a community-curated pipeline for bulk RNA sequence data analysis.
It runs quality control, adapter trimming, alignment, and quantification — the core steps of most RNA-seq experiments.

In this section, you will add the pipeline to your workspace, launch a run, and monitor its execution.

### 2.1. Add the pipeline to your workspace

Seqera Platform provides a curated collection of nf-core and other community pipelines that can be added to your workspace in a few clicks.

1. In the left sidebar, click **Pipelines**, then click **Add pipeline**.
2. Select **Seqera Pipelines** to browse the community collection.
3. Search for `rnaseq` and select **nf-core/rnaseq**.
4. Choose your workspace from the dropdown and click **Add pipeline**.

The pipeline is now listed in your workspace's **Pipelines** panel and is ready to launch.

### 2.2. Launch the pipeline

In the **Pipelines** panel, click **Launch** next to nf-core/rnaseq.
The pipeline is already configured with the `test` profile, so the input data, output directory, and genome reference are pre-filled.
Select your compute environment, then click **Launch**.

!!! info

    The launch form has three sections. For reference, here is what each contains:

    **General config**: pipeline revision, config profiles (e.g. `test`), and compute environment.

    **Run parameters**: `input` (path to a samplesheet), `outdir` (output directory), `genome`, and other pipeline-specific parameters.

    **Advanced settings**: Nextflow version, and environment variables such as `NXF_SYNTAX_PARSER=v1`, which nf-core pipelines require on Nextflow 25.10 and later.

### 2.3. Monitor execution

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

### Takeaway

You know how to add a pipeline to a Seqera Cloud workspace, configure and launch a run, and monitor execution at scale.

### What's next?

Learn how to do all of this from the command line using the `tw` CLI.

---

## 3. Launch pipelines from the command line

In the run view, click the **Command line** tab.
You will see the exact `nextflow run` command that Platform constructed and submitted on your behalf — the same kind of command you have been running manually throughout this course.

Platform does not replace Nextflow; it orchestrates it.
Everything you can do through the web interface, you can also do from a terminal using the `tw` CLI, the command-line tool for interacting with the Platform API.
This is useful for automating launches from scripts or CI/CD pipelines.

We're going to do this now from the Github codespace we used for the first two parts of this training.

### 3.1. Get an access token

The `tw` CLI authenticates with Seqera using a personal access token.

1. In the Seqera web interface, click your avatar in the top-right corner and select **Your tokens**.
2. Click **Add token**, give it a name (e.g. `training`), and click **Add**.
3. Copy the token value — it will only be shown once.
   If you lose it, you can always invalidate it and generate another one.

### 3.2. Install the tw CLI

Run the following commands in your Codespace terminal to download and install the `tw` binary:

```bash
curl -fsSL https://github.com/seqeralabs/tower-cli/releases/latest/download/tw-linux-x86_64 -o tw
chmod +x tw
sudo mv tw /usr/local/bin/
```

Verify the installation:

```bash
tw --version
```

### 3.3. Configure the CLI

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
        Details
    -------------------------+-----------------------------
     Tower API endpoint      | https://api.cloud.seqera.io
     Tower API version       | 1.150.0
     Tower version           | 26.1.0-cycle54
     CLI version             | 0.30.0 (fde9dec)
     CLI minimum API version | 1.148.0
     Authenticated user      | <your-name>>

    System health status
    ---------------------------------------+----
     Remote API server connection check    | OK
     Tower API version check               | OK
     Authentication API credential's token | OK
    ```

[TODO: add a brief recap sentence]

### 3.4. Explore your workspace from the CLI

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

### 3.5. Launch nf-core/rnaseq from the CLI

The pipeline you added to your workspace in section 2.1 is available by name in the CLI.
Launch it with the `test` profile:

```bash
tw launch nf-core-rnaseq \
    --workspace <org>/<workspace> \
    --compute-env <compute-env-name> \
    -p test
```

Replace `<org>/<workspace>` and `<compute-env-name>` with the values your instructor provided.

??? success "Command output"

    ```console
    Launching pipeline nf-core-rnaseq
    Run name: focused_einstein
    https://cloud.seqera.io/orgs/my-org/workspaces/my-workspace/watch/<run-id>
    ```

Open the link in your browser and confirm the run appears in the **Runs** panel.

Once you can see it running, you have confirmed that the CLI and the web interface are two views onto the same workspace.

!!! note

    You can also pass a full GitHub URL directly to `tw launch` without adding the pipeline first.
    However, adding it to the workspace with `tw pipelines add` is generally better: it saves the pipeline configuration for future runs, makes it available by name, and makes it visible to all workspace members in the Launchpad.
    Section 4 walks through that workflow for nf-core/demo.

### Takeaway

You know how to authenticate the `tw` CLI, inspect your workspace, and launch a saved pipeline from the terminal.

### What's next?

Add a new pipeline to your workspace from the command line and launch it.

---

## 4. Add a new pipeline and run it

Any Nextflow pipeline on GitHub can be added to your workspace with `tw pipelines add`, as long as it has a `main.nf` entry point and a `nextflow.config` at its root.
nf-core/demo is a good example to practice with: you already ran it in Part 2, so you know what it does and what to expect.

### 4.1. Add nf-core/demo to your workspace

```bash
tw pipelines add \
    --workspace <org>/<workspace> \
    --name nf-core-demo \
    https://github.com/nf-core/demo
```

??? success "Command output"

    ```console
    New pipeline 'nf-core-demo' added at <org>/<workspace> workspace
    ```

### 4.2. Verify it appears in the Launchpad

```bash
tw pipelines list --workspace <org>/<workspace>
```

??? success "Command output"

    ```console
    Pipelines at <org>/<workspace>:
    ID   | Name            | Repository
    ---- | --------------- | -----------------------------------------
    ...  | nf-core-demo    | https://github.com/nf-core/demo
    ...  | nf-core-rnaseq  | ...
    ```

Open your workspace in the browser and click **Launchpad** to confirm nf-core/demo now appears alongside nf-core/rnaseq.

!!! tip

    You can also add pipelines via the web interface: in the left sidebar, click **Pipelines**, then **Add pipeline**, and select **GitHub repository**.
    The web form lets you set a compute environment and configure default launch parameters at the same time.

### 4.3. Launch nf-core/demo

```bash
tw launch nf-core-demo \
    --workspace <org>/<workspace> \
    --compute-env <compute-env-name> \
    -p test
```

??? success "Command output"

    ```console
    Launching pipeline nf-core-demo
    Run name: compassionate_darwin
    https://cloud.seqera.io/orgs/my-org/workspaces/my-workspace/watch/<run-id>
    ```

Because you are using the `test` profile, this run completes quickly.
Once it finishes, click into the run to explore the task table and execution reports — the same views you saw for the rnaseq run, at a smaller scale.

### Takeaway

You know how to add any GitHub-hosted Nextflow pipeline to your workspace using the CLI and launch it.

---

## Summary

In this part you learned to:

- Sign up for a Seqera account and explore the Community Showcase
- Add a pipeline from the curated catalog, launch a production-scale run, and monitor execution
- Authenticate the `tw` CLI and launch a saved pipeline from the terminal
- Add a new pipeline from GitHub using the CLI and verify it appears in the Launchpad
