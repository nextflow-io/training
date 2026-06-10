# Development

Connect your GitHub account, fork rnaseq-nf into your own namespace, introduce a change that will cause a runtime failure, and launch the pipeline.
The failure you create here is what you will diagnose in the next lesson.

---

## 1. Connect GitHub access.

CoScientist can act on GitHub repositories on your behalf once you have authorized the GitHub App integration.
Follow the authorization flow in the CoScientist interface to grant it access to your GitHub account.

<!-- TODO: screenshot: GitHub connection / authorization flow -->
<!-- TODO: verify the exact GitHub connection UI flow and labels -->

!!! note "Checkpoint"

    CoScientist reports that it can access your GitHub account.

## 2. Fork rnaseq-nf.

Ask CoScientist to create a personal fork of the upstream pipeline:

```text
Fork nextflow-io/rnaseq-nf into my GitHub account.
```

CoScientist calls the GitHub API over MCP and creates the fork in your account.

!!! note "Checkpoint"

    The repository `your-user/rnaseq-nf` exists on GitHub.

## 3. Make a change that we will run.

Ask CoScientist to introduce a small configuration change that will cause the pipeline to fail at runtime:

```text
In my fork, lower the memory available to the QUANT process to something very small (for example 1 MB) so we can see how the pipeline behaves under tight resources, and commit it.
```

!!! warning

    The change is deliberately a runtime failure, not a syntax error.
    A syntax error is caught before the pipeline starts: Nextflow parses the configuration and halts immediately, so no run record is created on the Platform.
    A resource constraint that is too low passes parsing but causes the process to fail mid-execution, producing a run record you can open, inspect, and debug.
    That failed run record is what the next lesson works with.

<!-- TODO: verify the exact break reproduces a FAILED run on the training compute environment, and record the precise resulting config (e.g. a withName: QUANT { memory = ... } directive) here -->

!!! note "Checkpoint"

    The change is committed to your fork and visible in the commit history.

## 4. Launch the pipeline.

Ask CoScientist to submit a run for your fork using the training compute environment:

```text
Launch my fork of rnaseq-nf on the Launchpad using the test profile and the available compute environment.
```

The run is expected to fail.
The next lesson uses that failure as the starting point for debugging with CoScientist.

!!! note "Checkpoint"

    A run for your fork appears in the Platform **Runs** list with status `submitted` or `running`.

---

> **Break point for live sessions:** this is a natural place to pause.

### Takeaway

You forked rnaseq-nf, committed a change, and launched the pipeline entirely through CoScientist — the run is now executing on the Platform.

### What's next?

In the next lesson, [debug the failed run](04_debugging.md).
