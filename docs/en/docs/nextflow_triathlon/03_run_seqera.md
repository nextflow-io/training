# Part 3: Run on Seqera Cloud

In this third and final part of the Nextflow Triathlon, you will use Seqera to launch, monitor, and manage Nextflow pipelines through a web interface.

We will start by setting up access to Seqera and exploring what's available in the platform out of the box.
We will then launch the same nf-core/rnaseq pipeline we ran in Part 2, first from the web interface, then again from the terminal.
Finally, we will add the nf-core/demo pipeline directly from its GitHub repository and launch it two ways: from the web interface by filling in parameters manually, and from the CLI using a parameter file.

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

The pipeline is now listed in your workspace's **Launchpad** panel and is ready to launch.

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

You know how to add a pipeline to a Seqera workspace, configure and launch a run, and monitor execution at scale.

### What's next?

Drive the same workspace from the command line using the `seqera` CLI and `nextflow launch`.

---

## 3. Launch pipelines from the command line

In the run view, click the **Command line** tab.
You will see the exact `nextflow run` command that Platform constructed and submitted on your behalf — the same kind of command you have been running manually throughout this course.

Platform does not replace Nextflow; it orchestrates it.
Everything you can do through the web interface, you can also do from a terminal.
You will use two commands:

- `seqera` for authentication and workspace context. It logs you in via your browser and remembers your default org and workspace.
- `nextflow launch` for actually submitting pipelines to Platform. It is the same `nextflow` binary you used in Part 2, with `launch` swapped in for `run` to push the workload to your workspace instead of running it locally.

This is the same mental model as Part 2: you stay in the `nextflow` CLI you already know.

!!! note "Why not `tw`?"

    Earlier versions of this training used the legacy `tw` CLI (Tower CLI).
    It is being replaced by the `seqera` CLI for authentication and by `nextflow launch` for pipeline submission.
    The new flow removes manual token handling and a parallel launch command.
    The same flags you learned for `nextflow run` work for `nextflow launch`.

We will run these commands in the GitHub Codespace we used for the first two parts of this training.

### 3.1. Install the `seqera` CLI

Run the install script in your Codespace terminal:

```bash
curl -fsSL https://ai.seqera.io/install | bash
```

Verify the installation:

```bash
seqera version
```

The `seqera` CLI is installed and ready to authenticate.

### 3.2. Log in

Run:

```bash
seqera login
```

The command opens your browser and walks you through Seqera's OAuth flow.
After you approve, the CLI receives a refresh token and stores it in your operating system's keychain (or a file fallback on headless systems).
You do not have to copy or paste anything.

This persists across Codespace sessions; you only need to log in once per environment.

!!! tip "Headless or scripted environments"

    If you need to authenticate non-interactively (CI, automation, scripts that can't open a browser), set `SEQERA_ACCESS_TOKEN` instead.
    Generate a personal access token from your avatar menu → **Your tokens** in the web interface and export it:

    ```bash
    export SEQERA_ACCESS_TOKEN=<your-token>
    ```

### 3.3. Select your org and workspace

`seqera login` authenticates you, but you also need to tell the CLI which organization and workspace to operate on.

List the organizations you have access to:

```bash
seqera org list
```

Pick the one you're using for the training and set it as the default:

```bash
seqera org select <org-name>
```

Set the workspace within that org:

```bash
seqera workspace select <workspace-name>
```

Verify everything is wired up:

```bash
seqera info
```

??? success "Command output"

    ```console
     Platform API     | https://api.cloud.seqera.io
     Authenticated as | <your-name>
     Organization     | <org-name>
     Workspace        | <workspace-name>
     Default compute  | <compute-env-name>
     Status           | OK
    ```

Once your org and workspace are set, every `seqera` and `nextflow launch` command operates against that workspace.
You can always override per command with `--workspace <org>/<workspace>` if you need to.

!!! tip "Machine-readable output"

    Pass `-o json` (or `--output json`) to get a structured response.
    Useful if you're driving the CLI from a script or a coding agent:

    ```bash
    seqera info -o json
    ```

### 3.4. Explore your workspace from the CLI

List the runs in your workspace, including the nf-core/rnaseq run you launched from the web interface in section 2:

```bash
seqera runs list
```

??? success "Command output"

    ```console
    Pipeline runs in <org>/<workspace>:
    ID        | Status   | Pipeline               | Run name
    --------- | -------- | ---------------------- | --------------
    <run-id>  | RUNNING  | nf-core/rnaseq         | happy_curie
    ```

You can filter the listing:

```bash
seqera runs list --status RUNNING --limit 5
```

Dig into a specific run by ID:

```bash
seqera runs show <run-id>
```

Everything visible in the web interface's **Runs** panel is also reachable from the CLI.

### 3.5. Launch nf-core/rnaseq from the CLI

You don't need to do anything special to launch a pipeline that's already in your Launchpad.
You can refer to it by its short name, exactly the way `nextflow run` accepts an nf-core pipeline name:

```bash
nextflow launch nf-core/rnaseq -profile test
```

??? success "Command output"

    ```console
    Launching `nf-core/rnaseq` on Seqera Platform
    Run name: focused_einstein
    https://cloud.seqera.io/orgs/<org>/workspaces/<workspace>/watch/<run-id>
    ```

Open the link in your browser and confirm the run appears in the **Runs** panel.

Once you can see it running, you have confirmed that the CLI and the web interface are two views onto the same workspace.

!!! tip "Same flags as `nextflow run`"

    `nextflow launch` uses the same flag grammar as `nextflow run`: `-profile`, `-revision`, `-params-file`, `-resume`, `-w`, plus `--<param>=<value>` for arbitrary pipeline parameters.
    The only difference is that `launch` submits to Seqera Platform instead of running locally.
    Useful Platform-specific flags: `-workspace <org>/<ws>`, `-compute-env <name>`, `-name <mnemonic>`.

### Takeaway

You know how to log in with `seqera`, pin your default org and workspace, list runs, and launch a saved pipeline from the terminal with `nextflow launch`.

### What's next?

Add a new pipeline to your workspace from the command line and launch it with inline parameters.

---

## 4. Add a new pipeline and run it

Any Nextflow pipeline on GitHub can be launched from your workspace, as long as it has a `main.nf` entry point and a `nextflow.config` at its root.
nf-core/demo is a good example to practice with: you already ran it in Part 2, so you know what it does and what to expect.

You don't actually need to register a pipeline to launch it — `nextflow launch` accepts a pipeline name or URL directly.
Registering is still useful when you want the pipeline to appear in the **Launchpad** so other workspace members can find and run it.
We'll do both: register it (section 4.1), then launch it from the web (4.3) and from the CLI (4.4).

### 4.1. Add nf-core/demo to your workspace

Register the pipeline in your workspace:

```bash
seqera pipelines add \
    --name nf-core-demo \
    https://github.com/nf-core/demo
```

??? success "Command output"

    ```console
    Pipeline 'nf-core-demo' added to <org>/<workspace>
    ```

The pipeline is now registered and will appear in the Launchpad.

### 4.2. Verify it appears in the Launchpad

List the pipelines in your workspace to confirm it was added:

```bash
seqera pipelines list
```

??? success "Command output"

    ```console
    Pipelines in <org>/<workspace>:
    Name            | Repository
    --------------- | -----------------------------------------
    nf-core-demo    | https://github.com/nf-core/demo
    nf-core-rnaseq  | ...
    ```

Open your workspace in the browser and click **Launchpad** to confirm nf-core/demo now appears alongside nf-core/rnaseq.

!!! tip

    You can also add pipelines via the web interface: in the left sidebar, click **Launchpad**, then **Add pipeline**, and fill out the form accordingly.

Click the **Launch** button on the nf-core/demo entry to open its launch form.
You will see that the `input` and `outdir` parameters are highlighted in red — they are required fields with no default values, because `seqera pipelines add` registered only the pipeline source without pre-configuring any parameters.
The next two sections walk through how to provide those values: first through the web form, then from the command line.

### 4.3. Launch nf-core/demo from the web interface

With the launch form open, fill in the two required parameters.

For `input`, enter the test samplesheet URL from the nf-core/demo test profile.
You can find it in `conf/test.config` inside the pipeline repository, which you examined in Part 2, section 2.1:

```
https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/samplesheet/samplesheet_test_illumina_amplicon.csv
```

For `outdir`, enter a cloud storage path where the pipeline can write its results.
Use the bucket configured for your workspace, with a subdirectory to keep runs organized:

```
s3://my-bucket/demo-results
```

Once both fields are filled, click the blue **Launch** button.

The run appears in the **Runs** panel and should complete in a few minutes on the test dataset.
Click into the run to explore the task table and any execution reports.

### 4.4. Launch nf-core/demo from the CLI

`nextflow launch` accepts pipeline parameters the same way `nextflow run` does: as `--<param>=<value>` flags directly on the command line.
You don't need to create a parameter file just to set one or two values.

```bash
nextflow launch nf-core/demo \
    -profile test \
    --outdir s3://my-bucket/demo-results-cli
```

??? success "Command output"

    ```console
    Launching `nf-core/demo` on Seqera Platform
    Run name: focused_feynman
    https://cloud.seqera.io/orgs/<org>/workspaces/<workspace>/watch/<run-id>
    ```

Open the link to confirm the run appears in the **Runs** panel.

The `test` profile provides the `input` samplesheet, and the `--outdir` flag supplies the output path.

!!! tip "When to use a parameter file"

    Inline `--<param>` flags work well for one or two values.
    For pipelines with many parameters, or when you want a reproducible record of exactly what was launched, put them in a YAML or JSON file and pass it with `-params-file`:

    ```yaml title="params.yaml"
    outdir: "s3://my-bucket/demo-results-cli"
    ```

    ```bash
    nextflow launch nf-core/demo -profile test -params-file params.yaml
    ```

    You can also bake defaults into the registered pipeline so they're applied to every launch:

    ```bash
    seqera pipelines add \
      --name nf-core-demo-defaults \
      --description "Demo pipeline with defaults for testing" \
      --params-file params.yaml \
      --profile test \
      https://github.com/nf-core/demo
    ```

### Takeaway

You know how to add any GitHub-hosted Nextflow pipeline to your workspace and launch it, both from the web interface by filling in parameters manually, and from the CLI with `nextflow launch` using inline parameter flags or a parameter file.

---

## Summary

In this part you learned to:

- Sign up for a Seqera account and explore the Community Showcase
- Add a pipeline from the curated catalog, launch a production-scale run, and monitor execution
- Authenticate the `seqera` CLI with `seqera login` and pin a default org and workspace
- List and inspect runs from the terminal with `seqera runs list` / `seqera runs show`
- Launch a saved pipeline on Platform with `nextflow launch`, using the same flag grammar as `nextflow run`
- Add a new pipeline from GitHub with `seqera pipelines add` and verify it appears in the Launchpad
- Launch a pipeline from the web interface by filling in required parameters manually
- Launch a pipeline from the CLI with inline `--<param>` flags or a parameter file
