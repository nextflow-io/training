# Exploring the Platform

Use CoScientist to inspect and add Platform assets so you can see MCP doing real work against your workspace.

---

## 1. List what is already in the workspace.

Send the following prompt to get a snapshot of the current workspace contents:

```text
List the pipelines on the Launchpad and the datasets in this workspace.
```

??? example "What CoScientist typically does"

    It returns the current Launchpad pipelines and datasets it can see over MCP.
    The exact wording and the number of items listed will differ depending on the workspace state.

## 2. Add rnaseq-nf to the Launchpad.

Ask CoScientist to register a pipeline on your behalf:

```text
Add the pipeline https://github.com/nextflow-io/rnaseq-nf to the Launchpad in this workspace, using the available compute environment.
```

CoScientist will confirm the action or ask which compute environment to use if more than one is available.

!!! note "Checkpoint"

    In the Seqera Platform web app, open **Launchpad**: an entry for `rnaseq-nf` is now present.

<!-- TODO: screenshot: Launchpad showing the new rnaseq-nf pipeline -->

## 3. Inspect the compute environment.

Ask CoScientist to read the configuration for the compute environment attached to the pipeline:

```text
Show me the details of the compute environment this pipeline will run on.
```

CoScientist reads the compute-environment configuration over MCP and returns the key fields.

!!! note "Checkpoint"

    CoScientist reports the compute environment name and platform (for example AWS Batch) matching the training compute environment.

## 4. Reach all workspace assets.

The same MCP access that covers pipelines also covers datasets, reference genomes, and data links.
Ask CoScientist to list the datasets to confirm the breadth of its read access:

```text
List the datasets available in this workspace.
```

!!! note "Checkpoint"

    CoScientist returns the datasets registered in the workspace, or reports that none exist yet.

### Takeaway

CoScientist can read and modify Platform assets (Launchpad pipelines, compute environments, and datasets) directly through MCP tool calls.

### What's next?

In the next lesson, [start developing against rnaseq-nf](03_development.md).
