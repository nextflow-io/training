# Part 1: Launch pipelines from the web interface

In this part of the Run with Seqera training course, you will set up access to Seqera Platform and launch a production-scale pipeline from the web interface.

Make sure your working directory is set to `seqera-run/` as instructed on the [Getting started](./00_orientation.md) page.

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
3. Browse the available pipelines — you will recognize several nf-core pipelines from the Run nf-core course.
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

    The rest of this course assumes you have access to a workspace with a configured compute environment.
    If you are in a group training session, your instructor will confirm which workspace and compute environment to use.

### Takeaway

You have a Seqera account, you've explored the Community Showcase, and you're able to access a workspace with compute.

### What's next?

Launch a production-scale RNA-seq pipeline from the Seqera Cloud web interface.

---

## 2. Launch nf-core/rnaseq from the web interface

As covered in Run nf-core, the nf-core/rnaseq pipeline is a community-curated pipeline for bulk RNA sequence data analysis.

In this section, you will add the pipeline to your workspace, launch a run, and monitor its execution.

### 2.1. Add the pipeline to your workspace

Conveniently, nf-core/rnaseq is part of a curated collection of pipelines that can be added to your workspace in a few clicks through the Seqera Pipelines service.

_We'll show you how to add your own pipelines later in this course._

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

This will take a while to run, so we'll continue on for now and circle back later to look at outputs and so on.

### Takeaway

You know how to add a pipeline to a Seqera workspace, configure and launch a run, and monitor execution at scale.

### What's next?

Head on to [Part 2](./02_launch_from_cli.md), where you'll learn how to do all of this from the command line using the `tw` CLI.

---

## Summary

In this part you learned to:

- Sign up for a Seqera account and explore the Community Showcase
- Add a pipeline from the curated catalog, launch a production-scale run, and monitor execution
