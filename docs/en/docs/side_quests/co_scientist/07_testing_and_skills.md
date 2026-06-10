# Testing, and building a reusable skill

First, write an `nf-test` for `rnaseq-nf` and learn what is safe to snapshot.
Then, package that testing workflow into a reusable skill you can invoke as a slash command any time.

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

## 4. Why turn this into a skill.

You just walked CoScientist through the same reasoning a careful tester applies every time: assert on stable output, exclude timestamps, version strings, and work-directory paths.
A reusable **skill** captures that reasoning once and exposes it as a slash command, so you or a teammate get the same disciplined `nf-test` every time without re-explaining the rules.
CoScientist exposes reusable skills as slash commands.

## 5. Author the skill.

A skill is a small file with a name, a description, and instructions the agent follows when the skill is invoked.
The file below defines the `write-nf-test` skill; save it as `write-nf-test.md` in your skills directory.

```markdown
---
name: write-nf-test
description: Generate an nf-test for a Nextflow process that asserts on stable output and excludes unstable content.
---

<!-- NOTE: illustrative skill file; verify CoScientist's actual skill format/location before use -->

# Write nf-test

When this skill is invoked with a process name, generate an `nf-test` for that process following the rules below.

## Steps

1. Scaffold an `nf-test` for the named process in the `tests/` directory.
   Use the `nextflow_process` block with `name`, `script`, and `process` fields.
   Provide a representative `input` in the `when` block.

2. Assert on deterministic output only:

   - `Salmon` per-transcript count columns (stable given fixed inputs)
   - `FASTQC` pass/fail status (stable for fixed data)
   - File existence and line counts (structural, stable)

3. Do NOT snapshot unstable content:

   - `MultiQC` HTML reports (embed timestamps and software versions)
   - `Salmon` `cmd_info.json` and log files (contain timestamps and absolute paths)
   - Any path or value containing the work-directory hash (changes every run)
   - Version strings (change with tool and container updates)

4. Run `nf-test test` on the generated file.
   If the test fails, read the failure message and narrow the assertion rather than updating the snapshot when the cause is unstable content.
   Repeat until the test reports a pass.
```

<!-- TODO: verify the exact CoScientist skill file format, frontmatter fields, and on-disk/in-workspace location against the product; this example mirrors the Claude Code skill shape -->

## 6. Install and invoke the skill.

Place the skill file in the location CoScientist reads for registered skills, then restart or reload the agent session so it picks up the new command.

<!-- TODO: verify how a skill is installed/registered and the exact slash-command invocation syntax -->

Invoke the skill as a slash command:

```text
/write-nf-test for the FASTQC process
```

!!! note "Checkpoint"

    Invoking the skill produces an nf-test that asserts on stable output and excludes unstable content, without you re-explaining the rules.

### Takeaway

You added `nf-test` coverage to a pipeline that had none, made a deliberate choice about what is safe to snapshot, and captured that testing discipline as a reusable skill exposed as a slash command.

### What's next?

[Wrap up the side quest](next_steps.md).
