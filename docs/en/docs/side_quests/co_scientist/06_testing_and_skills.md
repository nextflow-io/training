# Testing, and building a reusable skill

This lesson moves you from the web chat to the `seqera ai` command-line agent.
You start in the web chat to write an `nf-test` for `rnaseq-nf` and learn what is safe to snapshot.
Then you switch to the CLI to run that test until it passes, and finally package the testing workflow into a reusable skill you can invoke as a slash command any time.

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

## 3. Move to the Seqera AI CLI.

This is the milestone where you leave the web chat and move to the `seqera ai` command-line agent.
Sections 1 and 2 happened in the browser; everything from here runs in your own terminal.

You are about to run and iterate on a test, which is exactly what the CLI is built for:

- Compute: the CLI uses your machine, which has more resources than the in-browser environment.
- Docker: it is available locally, so containerized processes and `nf-test` runs work.
- Local access: the agent works against your own files, repos, editor, and personal credentials and config.
- Iteration: you edit, run, and see results in the same terminal, without round-tripping through the web chat.

Install the CLI (it needs Node.js 18 or later) and sign in:

```bash
npm install -g seqera
seqera login
```

`seqera login` opens a browser to authenticate against your Seqera account.
Then start an interactive session from your pipeline directory:

```bash
seqera ai
```

To run the test repeatedly until it passes, use **goal mode**.
Goal mode keeps working toward an objective across several attempts and stops once the goal is met:

```text
/goal run the nf-test for the QUANT process and fix it until it passes
```

The agent runs `nf-test`, reads any failure, narrows the assertion or updates the snapshot, and re-runs, repeating until the test passes.
By default the CLI asks for your approval before it runs a command, so you see each `nf-test` invocation before it executes.

!!! note "Checkpoint"

    `nf-test test` for the `QUANT` test reports a passing test.

## 4. Enhance CoScientist with your own skills.

CoScientist ships with built-in skills, and you can add your own to extend what it does.
A reusable **skill** captures a workflow once and exposes it as a slash command, so you or a teammate get the same result every time without re-explaining the steps.
Packaging the testing rules you just applied makes a good first skill: it produces a disciplined `nf-test` on demand.

## 5. Author the skill.

A CoScientist skill is a directory containing a `SKILL.md` file.
The file has YAML frontmatter with a `name` and a `description`, followed by a markdown body of instructions the agent reads when the skill runs.
CoScientist discovers skills from your project and user skill directories and surfaces each one as a slash command.

Create the skill in your project at `.agents/skills/write-nf-test/SKILL.md` with this content:

```markdown
---
name: write-nf-test
description: Generate an nf-test for a Nextflow process that asserts on stable output and excludes unstable content.
---

# Write nf-test

When this skill is invoked with a process name, generate an `nf-test` for that process following the rules below.

## Steps

1. Scaffold an `nf-test` for the named process under `tests/`.
   Use a `nextflow_process` block with the process `name`, `script`, and `process` fields, and provide a representative `input` in the `when` block.

2. Assert on deterministic output only:

   - The per-transcript count columns in `quant.sf`
   - The existence of the expected output files and directories

3. Do NOT snapshot unstable content:

   - MultiQC HTML reports (embedded timestamps and versions)
   - Salmon `cmd_info.json`, `meta_info.json`, and log files (timestamps, command line, version)
   - Any path or value containing the work-directory hash
   - Version strings

4. Run `nf-test test` on the generated file.
   When it fails, narrow the assertion rather than snapshotting unstable content, and re-run until it passes.
```

Keep skills small: CoScientist caps each discovered skill's context at 5 KB.
The `.agents/skills/` location follows the cross-agent Agent Skills convention, so the same skill works in other compatible agents.

## 6. Install and invoke the skill.

CoScientist reads the `.agents/skills/` directory in your project automatically.
Restart `seqera ai` so it discovers the new skill, then type `/` and confirm `write-nf-test` appears in the command palette.

Invoke it as a slash command, naming a different process:

```text
/write-nf-test for the FASTQC process
```

To make the skill available across all your projects rather than one, place the same directory under `~/.agents/skills/` instead.

!!! note "Checkpoint"

    Invoking the skill produces an nf-test that asserts on stable output and excludes unstable content, without you re-explaining the rules.

### Takeaway

You added `nf-test` coverage to a pipeline that had none, made a deliberate choice about what is safe to snapshot, and captured that testing discipline as a reusable skill exposed as a slash command.

### What's next?

[Wrap up the side quest](next_steps.md).
