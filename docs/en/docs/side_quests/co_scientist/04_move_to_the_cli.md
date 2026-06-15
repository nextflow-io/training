# Move to the CLI

You developed the change and opened the pull request in the web chat.
You can keep working there, but running and iterating on tests is faster in your own terminal, so now you move to the `seqera ai` command-line agent.

---

## 1. Why move to the CLI.

Developing and testing a pipeline is better suited to your own terminal than the browser:

- Compute: the CLI uses your machine, which has more resources than the in-browser environment.
- Docker: it is available locally, so containerized processes and test runs work.
- Local access: the agent works against your own files, repos, editor, and personal credentials and config.
- Iteration: you edit, run, and see results in one terminal, without round-tripping through the web chat.

## 2. Install the CLI and get your fork.

Install the CLI (it needs Node.js 18 or later) and sign in:

```bash
npm install -g seqera
seqera login
```

`seqera login` opens a browser to authenticate against your Seqera account.
Start an interactive session, then bring the branch with your `fastp` change down to work on:

```bash
seqera ai
```

```text
Clone my fork of rnaseq-nf locally, on the branch with the fastp pull request, so I can run and test it here.
```

!!! note "Checkpoint"

    A local clone of your fork is on your machine, on the branch with the `fastp` change, and `seqera ai` is running in it.

## 3. Ask for an nf-test.

`rnaseq-nf` has no `nf-test` coverage.
If you are new to `nf-test`, see [Testing with nf-test](../nf_test/index.md) for background on snapshot assertions.

The simplest place to start is a pipeline-level test that runs the whole pipeline on its test data and checks the outputs.
Send CoScientist the following prompt:

```text
rnaseq-nf has no nf-test coverage. Add a pipeline-level nf-test that runs the whole pipeline on its test data and checks the outputs.
```

??? example "What CoScientist typically does"

    It scaffolds a `tests/` directory and a pipeline-level `.nf.test` that runs the pipeline end to end and asserts on the published outputs.
    The exact wording will differ from run to run.

For a reference of what the generated test looks like, see [`solutions/pipeline.nf.test`](solutions/pipeline.nf.test).

## 4. Snapshot-stable versus unstable output.

Not all output is safe to snapshot.
Snapshotting unstable output causes tests to fail on every run for reasons unrelated to correctness.
Use the table below to decide what to assert on.

| Output                                                | Snapshot? | Why                                                     |
| ----------------------------------------------------- | --------- | ------------------------------------------------------- |
| `quant_<id>/quant.sf` (Salmon per-transcript counts)  | Yes       | Deterministic counts for fixed inputs and version       |
| Output file and directory existence                   | Yes       | Structural, stable                                      |
| `quant_<id>/cmd_info.json`, `aux_info/meta_info.json` | No        | Embed the command line, Salmon version, and metadata    |
| `quant_<id>/logs/salmon_quant.log`                    | No        | Contains timestamps and runtimes                        |
| `multiqc_report.html`                                 | No        | Embeds a timestamp, tool versions, and an absolute path |
| FastQC `*_fastqc.zip` / `*_fastqc.html`               | No        | Zips embed timestamps and the FastQC version            |
| Any path containing the work-directory hash           | No        | Changes every run                                       |

!!! tip

    Salmon can introduce tiny nondeterminism across threads.
    Run it single-threaded (`--threads 1`) when you want a byte-stable `quant.sf` to snapshot.

Send CoScientist the following prompt to steer the assertion:

```text
Assert on the columns in quant.sf and that the expected output files exist. Do not snapshot the MultiQC HTML, the Salmon logs, cmd_info.json, or anything containing timestamps, versions, or work directory paths.
```

## 5. Run the test to green.

To run the test repeatedly until it passes, use **goal mode**:

```text
/goal run the pipeline nf-test and fix it until it passes
```

Goal mode keeps working toward an objective across several attempts and stops once the goal is met.
The agent runs `nf-test`, reads any failure, narrows the assertion or updates the snapshot, and re-runs, repeating until the test passes.
By default the CLI asks for your approval before it runs a command, so you see each `nf-test` invocation before it executes.

!!! note "Checkpoint"

    `nf-test test` reports a passing test for the pipeline.

## 6. Add a test workflow to CI.

Ask the agent to wire the test suite into continuous integration:

```text
Add a GitHub Actions workflow that installs nf-test and runs the test suite on every push and pull request.
```

A CI workflow runs the tests automatically on every change, so a regression is caught before it merges.

!!! note "Checkpoint"

    A workflow file exists under `.github/workflows/` that runs `nf-test`.

<!-- TODO: verify the exact generated GHA workflow (nf-test install action, Nextflow setup) against a real run -->

## 7. Add the tests and CI to your pull request.

The pull request you opened earlier contains the `fastp` step.
Ask the agent to commit the tests and the workflow to the same branch so the pull request picks them up:

```text
Commit the nf-test and the CI workflow to the branch of my open pull request and push.
```

!!! note "Checkpoint"

    The open pull request now also contains the tests and the CI workflow.

Open the pull request and confirm the diff contains the `fastp` step, the tests, and the workflow, and that it targets the branch you intend.

<!-- TODO: screenshot: the open PR on GitHub -->

### Takeaway

You moved to the CLI, added test coverage the pipeline lacked, chose what was safe to snapshot, and wired the tests into CI.
All of it joined the pull request you opened from the web chat.

### What's next?

[Package this into a reusable skill](05_build_a_skill.md).
