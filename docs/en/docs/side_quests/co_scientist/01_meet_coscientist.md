# Meet CoScientist

CoScientist is Seqera's AI assistant for bioinformatics.
Open the chat interface, connect it to your training workspace, and have a first conversation to see what it can do.
Then use it to inspect and add Platform assets, so you see how it reaches your Platform and acts on it.

---

## 1. Open the chat and connect to your workspace.

Navigate to [https://ai.seqera.io/chat](https://ai.seqera.io/chat) and sign in with your Seqera credentials.

<!-- TODO: verify the CoScientist chat URL is stable and publicly correct before publishing -->
<!-- TODO: screenshot: the CoScientist chat landing view after login -->

Use the workspace selector to choose the provided training workspace.

<!-- TODO: verify exact label of the workspace selector control in the UI -->
<!-- TODO: screenshot: workspace selector showing the training workspace -->

## 2. Confirm the connection.

Send the following prompt to confirm CoScientist can see your workspace:

```text
Which Seqera workspace am I connected to, and what compute environments are available here?
```

??? example "What CoScientist typically does"

    It names the connected workspace and lists the available compute environment(s).
    The exact wording will differ from run to run.

!!! note "Checkpoint"

    CoScientist correctly names the **training workspace** and lists at least one compute environment.
    If it cannot, the workspace connection is not set up; revisit step 1.

## 3. Ask CoScientist what it can do.

Get a baseline picture of the agent's capabilities in this workspace:

```text
What can you help me with in this workspace? What can you see and what actions can you take?
```

??? example "What CoScientist typically does"

    It describes that it can help develop pipelines, launch and monitor runs, browse data, and act on GitHub.
    The exact wording will differ from run to run.

## 4. How CoScientist reaches your assets.

CoScientist does not generate text and hand it back for you to act on.
It takes real actions on your behalf, using your authenticated Platform access.
It can list compute environments, launch a run, read run logs, and (once you connect GitHub in a later lesson) open a pull request.

!!! info "MCP and the Seqera Platform"

    CoScientist reaches the Platform through the **Seqera MCP server** (Model Context Protocol).
    MCP is the interface that exposes Platform operations as tools the agent can call: list compute environments, launch a run, read run logs, query datasets.
    The agent calls these tools with your own Platform token, so it acts with your permissions rather than a shared account.

## 5. Survey the workspace.

Ask CoScientist for a snapshot of what is already in the workspace:

```text
List the pipelines on the Launchpad and the datasets in this workspace.
```

??? example "What CoScientist typically does"

    It returns the current Launchpad pipelines and datasets it can see in the workspace.
    The exact wording and the number of items listed will differ depending on the workspace state.

The same access also covers reference genomes and data links.

!!! note "Checkpoint"

    CoScientist lists the workspace's pipelines and datasets, or reports that none exist yet.

## 6. Add rnaseq-nf to the Launchpad.

Ask CoScientist to register a pipeline on your behalf:

```text
Add the pipeline https://github.com/nextflow-io/rnaseq-nf to the Launchpad in this workspace, using the available compute environment.
```

CoScientist will confirm the action or ask which compute environment to use if more than one is available.

!!! note "Checkpoint"

    In the Seqera Platform web app, open **Launchpad**: an entry for `rnaseq-nf` is now present.

<!-- TODO: screenshot: Launchpad showing the new rnaseq-nf pipeline -->

## 7. Inspect the compute environment.

Ask CoScientist to read the configuration for the compute environment attached to the pipeline:

```text
Show me the details of the compute environment this pipeline will run on.
```

CoScientist reads the compute-environment configuration and returns the key fields.

!!! note "Checkpoint"

    CoScientist reports the compute environment name and platform (for example AWS Batch) matching the training compute environment.

### Takeaway

You connected CoScientist to your training workspace, had a first conversation, saw how it reaches the Platform through MCP, and used it to register and inspect a pipeline.
It acts on real Platform assets on your behalf, with your own permissions.

### What's next?

In the next lesson, [learn how to work effectively with the agent](02_working_with_the_agent.md) before you start changing code.
