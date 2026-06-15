# Test, add CI, and open a PR

`rnaseq-nf` has no tests.
In this lesson you add `nf-test` coverage with the agent, decide what is safe to snapshot, wire the tests into CI, and open a pull request.
You are already in an interactive `seqera ai` session in your forked pipeline directory from the previous lesson.

---

## 1. Ask for an nf-test.

`rnaseq-nf` currently has no `nf-test` coverage.
If you are new to `nf-test`, see [Testing with nf-test](../nf_test/index.md) for background on snapshot assertions and process-level tests.

Send CoScientist the following prompt:

```text
rnaseq-nf has no nf-test coverage. Add an nf-test for the QUANT process.
```

??? example "What CoScientist typically does"

    It scaffolds a `tests/` directory and a `.nf.test` file for the `QUANT` process, with a first assertion on the process output.
    The exact wording will differ from run to run.

For a reference of what the generated test looks like, see [`solutions/quant.nf.test`](solutions/quant.nf.test).
Because `QUANT` needs a Salmon index, the test builds one in a `setup` block before running the process.

## 2. Snapshot-stable versus unstable output.

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

## 3. Run the test to green.

To run the test repeatedly until it passes, use **goal mode**:

```text
/goal run the nf-test for the QUANT process and fix it until it passes
```

Goal mode keeps working toward an objective across several attempts and stops once the goal is met.
The agent runs `nf-test`, reads any failure, narrows the assertion or updates the snapshot, and re-runs, repeating until the test passes.
By default the CLI asks for your approval before it runs a command, so you see each `nf-test` invocation before it executes.

!!! note "Checkpoint"

    `nf-test test` for the `QUANT` test reports a passing test.

## 4. Add a test workflow to CI.

Ask the agent to wire the test suite into continuous integration:

```text
Add a GitHub Actions workflow that installs nf-test and runs the test suite on every push and pull request.
```

A CI workflow runs the tests automatically on every change, so a regression is caught before it merges.

!!! note "Checkpoint"

    A workflow file exists under `.github/workflows/` that runs `nf-test`.

<!-- TODO: verify the exact generated GHA workflow (nf-test install action, Nextflow setup) against a real run -->

## 5. Open a pull request.

Ask the agent to propose the work upstream:

```text
Open a pull request from my fork back to nextflow-io/rnaseq-nf describing the fastp step, the tests, and the CI workflow.
```

!!! note

    In a real contribution the pull request targets the upstream repository.
    In a training setting, if upstream PRs are restricted, it is fine to open the PR against your own fork's default branch.
    The trainer can adjust the target as appropriate for the environment.

!!! note "Checkpoint"

    A pull request is open and visible on GitHub.

Open the pull request and confirm the diff contains the fastp step, the tests, and the workflow, and that it targets the branch you intend.

<!-- TODO: screenshot: the open PR on GitHub -->

### Takeaway

You added test coverage to a pipeline that had none, made a deliberate choice about what is safe to snapshot, and put the tests in CI.
You then proposed it all as a pull request.

### What's next?

[Package this into a reusable skill](05_build_a_skill.md).
