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

See the [Seqera documentation](https://docs.seqera.io) for a full overview of account tiers and available features.

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

Seqera supports two ways to provide compute:

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

You have a Seqera account, you've explored the Community Showcase, and you're able to access a workspace with compute.

### What's next?

Launch a production-scale RNA-seq pipeline from the Seqera Cloud web interface.

---

## 2. Launch nf-core/rnaseq from the web interface

As we've already covered, the nf-core/rnaseq pipeline is a community-curated pipeline for bulk RNA sequence data analysis.

In this section, you will add the pipeline to your workspace, launch a run, and monitor its execution.

### 2.1. Add the pipeline to your workspace

Conveniently, nf-core/rnaseq is part of a curated collection of pipelines that can be added to your workspace in a few clicks through the Seqera Pipelines service.

_We'll show you how to add your own pipelines later in this lesson._

1. Navigate to [**Seqera Pipelines**](https://seqera.io/pipelines) to browse the community collection.
2. Search for `rnaseq` and select **nf-core/rnaseq**.
3. Click **Launch Pipeline** or scroll to the bottom of the page to the **Launch Pipeline** section.
4. Make sure you are logged in and select the appropriate values from the **Organizations**, **Workspace** and **Compute Environment** dropdown menus.
   **Tip for groups:** If you are using a shared workspace, add a unique identifier (such as your username) to the pipeline name.
5. Click **Add pipeline to your Seqera account**

A box will appear showing the message: **Pipeline added: View Pipeline**.
Clicking the link will take you to the pipeline entry in your launchpad.

The pipeline is now listed in your workspace's **Pipelines** panel and is ready to launch.

### 2.2. Launch the pipeline

Click the pipeline's **Launch** button, either in the **Launchpad** panel or on the pipeline details page.
This opens the configuration interface.

The pipeline is already configured with the `test` profile, so the input data, output directory, and genome reference are pre-filled.
You can ignore the rest of the parameters and advanced settings for now.

Click the blue **Launch** button to actually start the run.

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

This will take a while to run, so we'll continue on for now and we'll circle back later to look at outputs and so on.

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

### 3.1. Install the tw CLI

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

??? success "Command output"

    ```console
    tw version 0.30.0 (fde9dec)
    ```

The `tw` CLI is installed and ready to configure.

### 3.2. Get an access token

The `tw` CLI authenticates with Seqera using a personal access token.

1. In the Seqera web interface, click your avatar in the top-right corner and select **Your tokens**.
2. Click **Add token**, give it a name (e.g. `training`), and click **Add**.
3. Copy the token value — it will only be shown once.
   If you don't save it somewhere right away, you will need to generate another one.

### 3.3. Configure the CLI

For convenience, we're going to set up a configuration file containing the
access token you just generated and the workspace identifier.

We made a stub for you called `.seqera_config` under the `triathlon` directory in your training environment (in Codespaces).
Open the file in the editor and set the two variables:

- **`TOWER_ACCESS_TOKEN`**: the token you generated in section 3.2
- **`TOWER_WORKSPACE_ID`**: your workspace in `org-name/workspace-name` format

Once the values are filled in, load the config:

```bash
source ../.seqera_config
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
     Authenticated user      | <your-name>

    System health status
    ---------------------------------------+----
     Remote API server connection check    | OK
     Tower API version check               | OK
     Authentication API credential's token | OK
    ```

The `tw` CLI is now authenticated and connected to your Seqera account.
Run `source ../.seqera_config` at the start of each Codespace session to reload the config.

!!! tip

    If your workspace has no primary compute environment set, you can add `export TOWER_COMPUTE_ENV=<compute-env-name>` to your config file.
    Any config value can be overridden on the command line by passing the flag explicitly (e.g. `--compute-env other-env`).
    See the [tw CLI reference](https://docs.seqera.io/platform/latest/cli/reference) for the full list of options and environment variables.

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
tw runs list
```

??? success "Command output"

    ```console
    Pipeline runs for my-org/my-workspace workspace:
    ID        | Status   | Name          | Pipeline               | Run name
    --------- | -------- | ------------- | ---------------------- | --------
    <run-id>  | RUNNING  | ...           | nf-core/rnaseq         | happy_curie
    ```

The same run you are monitoring in the web interface is visible here.

!!! note

    Because `TOWER_WORKSPACE_ID` is set in `.seqera_config`, you can omit `--workspace` from all `tw` commands.
    Without the config, you would pass it explicitly:

    ```bash
    tw runs list --workspace <org>/<workspace>
    ```

Everything visible in the web interface is accessible from the CLI.

### 3.5. Launch nf-core/rnaseq from the CLI

The pipeline you added to your workspace in section 2.1 is available by name in the CLI.
Launch it with the `test` profile:

```bash
tw launch nf-core-rnaseq -p test
```

??? success "Command output"

    ```console
    Launching pipeline nf-core-rnaseq
    Run name: focused_einstein
    https://cloud.seqera.io/orgs/my-org/workspaces/my-workspace/watch/<run-id>
    ```

Open the link in your browser and confirm the run appears in the **Runs** panel.

Once you can see it running, you have confirmed that the CLI and the web interface are two views onto the same workspace.

!!! note

    You can also pass a full GitHub URL directly to `tw launch` without adding the pipeline to a workspace first.
    However, adding the pipeline explicitly before launching it is generally better: it saves the pipeline configuration for future runs, makes it available by name, and makes it visible to all workspace members in the Launchpad.

    It is possible to add a pipeline to a workspace directly from the command line using `tw`.
    Section 4 shows how to do this with the nf-core/demo pipeline.

### Takeaway

You know how to authenticate the `tw` CLI, inspect your workspace, and launch a saved pipeline from the terminal.

### What's next?

Add a new pipeline to your workspace from the command line and launch it.

---

## 4. Add a new pipeline and run it

Any Nextflow pipeline on GitHub can be added to your workspace with `tw pipelines add`, as long as it has a `main.nf` entry point and a `nextflow.config` at its root.
nf-core/demo is a good example to practice with: you already ran it in Part 2, so you know what it does and what to expect.

### 4.1. Add nf-core/demo to your workspace

Run the following command to register the pipeline in your workspace:

```bash
tw pipelines add \
    --name nf-core-demo \
    https://github.com/nf-core/demo
```

??? success "Command output"

    ```console
    New pipeline 'nf-core-demo' added at my-org/my-workspace workspace
    ```

The pipeline is now registered and will appear in the Launchpad.

### 4.2. Verify it appears in the Launchpad

List the pipelines in your workspace to confirm it was added:

```bash
tw pipelines list
```

??? success "Command output"

    ```console
    Pipelines at my-org/my-workspace:
    ID   | Name            | Repository
    ---- | --------------- | -----------------------------------------
    ...  | nf-core-demo    | https://github.com/nf-core/demo
    ...  | nf-core-rnaseq  | ...
    ```

Open your workspace in the browser and click **Launchpad** to confirm nf-core/demo now appears alongside nf-core/rnaseq.

!!! tip

    You can also add pipelines via the web interface: in the left sidebar, click **Launchpad**, then **Add pipeline**, and fill out the form accordingly.

### 4.3. Launch nf-core/demo

Launch the pipeline with the `test` profile:

```bash
tw launch nf-core-demo -p test
```

??? success "Command output"

    ```console
    Launching pipeline nf-core-demo
    Run name: compassionate_darwin
    https://cloud.seqera.io/orgs/my-org/workspaces/my-workspace/watch/<run-id>
    ```

Because you are using the `test` profile, this run should complete quickly.
Once it finishes, click into the run to explore the task table and execution reports, which are equivalent to the views you saw for the rnaseq run, at a smaller scale.

### Takeaway

You know how to add any GitHub-hosted Nextflow pipeline to your workspace using the CLI and launch it.

---

## Summary

In this part you learned to:

- Sign up for a Seqera account and explore the Community Showcase
- Add a pipeline from the curated catalog, launch a production-scale run, and monitor execution
- Authenticate the `tw` CLI and launch a saved pipeline from the terminal
- Add a new pipeline from GitHub using the CLI and verify it appears in the Launchpad
