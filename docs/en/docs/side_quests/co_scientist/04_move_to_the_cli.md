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

## 2. Open your fork in the devcontainer.

You forked `rnaseq-nf` in the previous lesson, and CoScientist opened a pull request with the `fastp` step.
Open your fork in the devcontainer, which comes with the tools you need already installed.
Then switch to the pull request's branch so you can run and test the change:

```bash
git checkout <branch>   # CoScientist names it like seqera-co-scientist/add-fastp-trimming-...
```

<!-- TODO: confirm the exact devcontainer entry (Codespaces vs local reopen-in-container) once the rnaseq-nf devcontainer PR is merged -->

!!! note "Checkpoint"

    You are in the devcontainer for your fork, on the branch with the `fastp` change.

## 3. Install and start the CLI.

The devcontainer already includes the `seqera` CLI.
To install it yourself, or to work outside the devcontainer, run the install script:

```bash
curl -fsSL https://ai.seqera.io/install | bash
```

Then authenticate with a Seqera Platform access token.
Create one at [cloud.seqera.io/tokens](https://cloud.seqera.io/tokens) and export it before starting the CLI:

```bash
export TOWER_ACCESS_TOKEN=<your token>
```

The Seqera AI credits are charged on a **per-organization** basis, so you must select the right organization before using the

```bash
seqera org switch
```

The interactive prompt will tell you which org to select

```console
Available organizations:

  1. community (pro)
  2. Your-Org

Enter selection (1-2), or 'q' to cancel:
```

!!! note "Checkpoint"

    `seqera ai` is running in your local clone, configured to use your Seqera Platform workspace AI credits.

## 4. Ask for an nf-test and run it to green.

`rnaseq-nf` has no `nf-test` coverage.
If you are new to `nf-test`, see [Testing with nf-test](../nf_test/index.md) for background on snapshot assertions.

The simplest place to start is a pipeline-level test that runs the whole pipeline on its test data and checks the outputs.
Before writing the prompt, decide what is safe to assert on: not all output is stable, and snapshotting unstable output causes tests to fail on every run for reasons unrelated to correctness.

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

With those rules in mind, ask CoScientist to create the test, telling it what to snapshot and what to leave out:

```text
Add a pipeline-level nf-test that runs the whole pipeline on its test data. Assert on the quant.sf columns and that the expected output files exist; do not snapshot the MultiQC HTML, Salmon logs, cmd_info.json, or anything containing timestamps, versions, or work directory paths.
```

??? example "What CoScientist typically does"

    It scaffolds a `tests/` directory and a pipeline-level `.nf.test` that runs the pipeline end to end and asserts on the stable outputs.
    The exact wording will differ from run to run.

For a reference of what the generated test looks like, see [`solutions/pipeline.nf.test`](solutions/pipeline.nf.test).

Then use **goal mode** to run the test repeatedly until it passes, fixing the test along the way:

```text
/goal run the nf-tests until they pass and fix any issues in the tests
```

Goal mode keeps working toward an objective across several attempts and stops once the goal is met.
The agent runs `nf-test`, reads any failure, narrows the assertion or updates the snapshot, and re-runs, repeating until the test passes.
By default the CLI asks for your approval before it runs a command, so you see each `nf-test` invocation before it executes.

!!! note "Checkpoint"

    `nf-test test` reports a passing test for the pipeline.

## 5. Add the tests and CI to your pull request.

A CI workflow runs the tests automatically on every change, so a regression is caught before it merges.
The pull request you opened earlier contains the `fastp` step.
Ask the agent to create the workflow and commit it and the tests to the same branch, so the pull request picks them up, in one prompt:

```text
Add a GitHub Actions workflow that installs nf-test and runs the test suite on every push and pull request, then commit the workflow and the nf-test to the branch of my open pull request and push.
```

<!-- TODO: verify the exact generated GHA workflow (nf-test install action, Nextflow setup) against a real run -->

!!! note "Checkpoint"

    A workflow file exists under `.github/workflows/` that runs `nf-test`, and the open pull request now contains the `fastp` step, the tests, and the workflow.

Open the pull request and confirm the diff contains the `fastp` step, the tests, and the workflow, and that it targets the branch you intend.

### Takeaway

You moved to the CLI, added test coverage the pipeline lacked, chose what was safe to snapshot, and wired the tests into CI.
All of it joined the pull request you opened from the web chat.

### What's next?

[Package this into a reusable skill](05_build_a_skill.md).
