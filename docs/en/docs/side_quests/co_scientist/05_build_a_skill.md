# Build a reusable skill

In the last lesson you hand-wrote a pipeline-level `nf-test` and learned which output is safe to snapshot.
Capture that discipline as a reusable skill, so CoScientist applies the same rules in every session without you re-explaining them.

---

## 1. Enhance CoScientist with your own skills.

CoScientist ships with built-in skills, and you can add your own.
A skill is a small instruction set that CoScientist loads as session context, so it follows your guidance automatically whenever you ask for related work.
Packaging the testing rules you just applied makes a good first skill.

## 2. Author the skill.

A skill lives in its own directory with a `SKILL.md` file: YAML frontmatter with a `name` and a `description` (both required, or the skill is skipped), followed by a markdown body of instructions.
CoScientist discovers skills from your project and user skill directories and sends them to the session as context.
A project directory like `.agents/skills/` takes priority, so a repository can override a global skill.

Create the directory and the empty skill file:

```bash
mkdir -p .agents/skills/write-nf-test
touch .agents/skills/write-nf-test/SKILL.md
```

Then open `.agents/skills/write-nf-test/SKILL.md` and add this content:

```markdown
---
name: write-nf-test
description: Generate an nf-test that asserts on stable output and excludes unstable content.
---

# Write nf-test

When asked to write an nf-test for a process or a pipeline, follow the rules below.

## Steps

1. Scaffold the test under `tests/`.
   For a process use a `nextflow_process` block; for a pipeline use a `nextflow_pipeline` block.
   Provide representative inputs or test data so the test runs end to end.

2. Assert on deterministic output only:

   - Numeric or tabular result files, such as a Salmon `quant.sf` table
   - The existence of the expected output files and directories

3. Do NOT snapshot unstable content:

   - Reports that embed timestamps or versions, such as MultiQC HTML
   - Log files, and files that record the command line or tool version (for example Salmon `cmd_info.json` and `meta_info.json`)
   - Any path or value containing the work-directory hash
   - Version strings

4. Run `nf-test test` on the generated file.
   When it fails, narrow the assertion rather than snapshotting unstable content, and re-run until it passes.
```

Keep skills small: CoScientist caps each discovered skill's context at 5 KB.
The `.agents/skills/` location follows the cross-agent Agent Skills convention, so the same skill works in other compatible agents.

## 3. Use the skill.

Restart `seqera ai` so it discovers the new skill and loads it into the session context.
Your own skills do not appear as slash commands; the `/` palette is reserved for built-in, backend-provided skills.
Instead, CoScientist now follows your skill whenever you ask for a test:

```text
Add an nf-test for the FASTQC process.
```

It applies the stable-versus-unstable rules from your skill automatically, without you restating them.
To make the skill available across all your projects rather than one, place the same directory under `~/.agents/skills/` instead.

!!! note "Checkpoint"

    With the skill in place, asking for a test produces one that asserts on stable output and excludes unstable content, without you re-explaining the rules.

### Takeaway

You captured the testing discipline as a reusable skill that CoScientist loads as session context.
You and your teammates now get a consistent `nf-test` on demand, without re-explaining the rules each time.

### What's next?

[Wrap up the side quest](next_steps.md).
