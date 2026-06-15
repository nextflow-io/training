# Build a reusable skill

You just walked the agent through writing a disciplined `nf-test`, steering it toward stable output and away from unstable content.
Capture that workflow as a reusable skill so you get the same result every time without re-explaining the rules.

---

## 1. Enhance CoScientist with your own skills.

CoScientist ships with built-in skills, and you can add your own to extend what it does.
A reusable **skill** captures a workflow once and exposes it as a slash command, so you or a teammate get the same result every time without re-explaining the steps.
Packaging the testing rules you just applied makes a good first skill: it produces a disciplined `nf-test` on demand.

## 2. Author the skill.

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

## 3. Install and invoke the skill.

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

You captured the testing discipline as a reusable skill exposed as a slash command.
You and your teammates now get a consistent `nf-test` on demand, without re-explaining the rules each time.

### What's next?

[Wrap up the side quest](next_steps.md).
