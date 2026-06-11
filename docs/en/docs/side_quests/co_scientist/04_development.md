# Development

Connect your GitHub account, fork rnaseq-nf into your own namespace, introduce a change that will cause a runtime failure, and launch the pipeline.
The failure you create here is what you will diagnose in the next lesson.

---

## 1. Connect GitHub access.

CoScientist can work with GitHub repositories on your behalf once it has access to your GitHub account.
In the web interface, connect your GitHub account so the agent can fork the pipeline, commit changes, and open pull requests.

<!-- TODO: screenshot: GitHub connection / authorization flow -->
<!-- TODO: verify the exact GitHub connection flow for the web interface. The Cloud docs do not document a GitHub connection mechanism; confirm whether it is a GitHub App, a token, or the agent running git/gh in the workspace. -->

!!! note "Checkpoint"

    CoScientist reports that it can access your GitHub account.

## 2. Fork rnaseq-nf.

Ask CoScientist to create a personal fork of the upstream pipeline:

```text
Fork nextflow-io/rnaseq-nf into my GitHub account.
```

CoScientist calls the GitHub API and creates the fork in your account.

!!! note "Checkpoint"

    The repository `your-user/rnaseq-nf` exists on GitHub.

## 3. Make a change that we will run.

Ask CoScientist to introduce a small configuration change that will cause the pipeline to fail at runtime:

```text
In my fork, add a config setting that gives the QUANT process only 10 MB of memory, so it fails at runtime under tight resources, and commit the change.
```

!!! warning

    The change is deliberately a runtime failure, not a syntax error.
    A syntax error is caught before the pipeline starts: Nextflow parses the configuration and halts immediately, so no run record is created on the Platform.
    A resource constraint that is too low passes parsing but causes the process to fail mid-execution, producing a run record you can open, inspect, and debug.
    That failed run record is what the next lesson works with.

<!-- TODO: verify the 10 MB QUANT limit reliably fails on the training compute environment (expected: salmon quant killed for exceeding memory, exit code 137). CoScientist should add `process { withName: QUANT { memory = '10.MB' } }` to nextflow.config; confirm the exact committed config and that the run reaches a failed status. -->

!!! note "Checkpoint"

    The change is committed to your fork and visible in the commit history.

Apply the habit from [Working with the agent](03_working_with_the_agent.md): open the commit and read the actual change rather than trusting the agent's summary.
Confirm it lowered the `QUANT` memory and left the rest of the configuration untouched.

## 4. Launch the pipeline.

Ask CoScientist to submit a run for your fork using the training compute environment:

```text
Launch my fork of rnaseq-nf on the Launchpad using the available compute environment.
```

The run is expected to fail.
The next lesson uses that failure as the starting point for debugging with CoScientist.

!!! note "Checkpoint"

    A run for your fork appears in the Platform **Runs** list with status `submitted` or `running`.

<!-- TODO: verify the exact casing the Platform uses for run statuses (submitted / running / succeeded / failed) and make all pages consistent with it -->

---

> **Break point for live sessions:** this is a natural place to pause.

### Takeaway

You forked rnaseq-nf, committed a change, and launched the pipeline entirely through CoScientist.
The run is now executing on the Platform.

### What's next?

In the next lesson, [debug the failed run](05_debugging.md).
