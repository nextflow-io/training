# Adding a test

`rnaseq-nf` has no `nf-test` coverage.
Ask CoScientist to add a test, make a deliberate choice about what is safe to snapshot, then switch to the `seqera ai` CLI to run the test to green.

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

## 2. Snapshot-stable versus unstable output.

Not all output is safe to snapshot.
Snapshotting unstable output causes tests to fail on every run for reasons unrelated to correctness.
Use the table below to decide what to assert on.

| Output                                      | Snapshot? | Why                                     |
| ------------------------------------------- | --------- | --------------------------------------- |
| `Salmon` per-transcript count columns       | Yes       | Deterministic given fixed inputs        |
| `FASTQC` pass/fail status                   | Yes       | Stable for fixed data                   |
| File existence and line counts              | Yes       | Structural, stable                      |
| `MultiQC` HTML report                       | No        | Embeds timestamps and software versions |
| `Salmon` `cmd_info.json` and log files      | No        | Contain timestamps and absolute paths   |
| Any path containing the work-directory hash | No        | Changes every run                       |
| Version strings                             | No        | Change with tool and container updates  |

Send CoScientist the following prompt to steer the assertion:

```text
Assert on the Salmon quantification columns and that the expected output files exist. Do not snapshot the MultiQC HTML or anything containing timestamps, versions, or work directory paths.
```

## 3. Switch to the CLI and run the test to green.

The same agent is available from the command line via `seqera ai`, which is useful for running and iterating on tests locally.
Install and authenticate:

```bash
npm install -g seqera
seqera login
```

<!-- TODO: verify the exact seqera ai CLI install, login, and invocation commands against the current product -->

Send the agent the following prompt from the CLI to drive the test to green:

```text
Run the nf-test for QUANT and fix it until it passes.
```

The agent runs `nf-test`, reads any failures, adjusts the assertions or snapshot, and re-runs.
It repeats this loop until the test reports a pass.
Failures caused by unstable output (timestamps, paths) are resolved by narrowing the assertion rather than updating the snapshot.

!!! note "Checkpoint"

    `nf-test test` for the `QUANT` test reports a passing test.

### Takeaway

You used CoScientist to add `nf-test` coverage to a pipeline that had none, made a deliberate choice about what is safe to snapshot, and drove the test to green from the CLI.

### What's next?

[Wrap up the side quest](next_steps.md).
