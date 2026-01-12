---
name: Run Tutorial Walkthrough
description: Walk through a training tutorial as a user would, progressively building examples exactly as instructed, running all commands, and verifying results against solutions. Use when testing tutorials end-to-end or validating that instructions are correct and complete.
---

# Run Tutorial Walkthrough

Walk through a training tutorial lesson as a learner would, progressively building code, running commands, and verifying the learning journey works end-to-end.

**This skill focuses on the unique value of simulating a learner's experience.** It automatically invokes other skills during the walkthrough:

- **Highlight validation** → Invokes `/validate` skill (includes highlight checking)
- **Pedagogical quality & lesson structure** → Invokes `/validate` skill
- **Testing final scripts** → Invokes `/test-example` skill for solution files

---

## Initial Setup Questions

**IMPORTANT**: Before starting the walkthrough, ask the user:

### 1. Which tutorial/lesson to run?

- A specific lesson file (e.g., `docs/hello_nextflow/01_hello_world.md`)
- A module to run all lessons (e.g., `hello_nextflow`)
- A side quest (e.g., `debugging`)

### 2. Which software environment to use?

- **Docker container (Recommended)** - Uses the training Docker image `ghcr.io/nextflow-io/training:latest`. This matches what learners use in Codespaces/Gitpod.
- **Local installation** - Uses locally installed Nextflow (requires matching tool versions)
- **Existing environment** - User has already set up an environment (just use it as-is)

**Why Docker is recommended**: The training materials are designed for the Codespaces/Gitpod environment which uses the `ghcr.io/nextflow-io/training:latest` image. This image has:

- Java 21 (required for plugin development)
- All dependencies pre-configured
- Nextflow version controlled via `NXF_VER` environment variable

---

## Docker Container Setup

### Determine NXF_VER from devcontainer.json

Before running Docker commands, read the Nextflow version from `.devcontainer/devcontainer.json`:

```bash
NXF_VER=$(grep -o '"NXF_VER":\s*"[^"]*"' .devcontainer/devcontainer.json | cut -d'"' -f4)
echo "Using NXF_VER=${NXF_VER}"
```

This ensures you use the same version learners get in Codespaces.

### Start a persistent container

First, clean up any existing container, then start a new one:

```bash
# Clean up any existing container
docker stop nf-training 2>/dev/null; docker rm nf-training 2>/dev/null

# Start fresh container with UTF-8 locale support
NXF_VER=$(grep -o '"NXF_VER":\s*"[^"]*"' .devcontainer/devcontainer.json | cut -d'"' -f4)
docker run -d --name nf-training \
  -e NXF_VER=${NXF_VER} \
  -e LANG=C.UTF-8 \
  -e LC_ALL=C.UTF-8 \
  -v ${PWD}:/workspaces/training \
  -w /workspaces/training \
  ghcr.io/nextflow-io/training:latest \
  sleep infinity
```

**Important**: The `LANG=C.UTF-8` and `LC_ALL=C.UTF-8` environment variables are critical for handling non-ASCII characters (like "Holà", "Grüß Gott") in file names and content. Without these, Nextflow may fail with "Malformed input or input contains unmappable characters" errors.

### Run commands in the container

Always include UTF-8 locale settings when executing commands:

```bash
docker exec -e LANG=C.UTF-8 -e LC_ALL=C.UTF-8 -w /workspaces/training/[working-dir] nf-training [command]
```

Example:

```bash
docker exec -e LANG=C.UTF-8 -e LC_ALL=C.UTF-8 -w /workspaces/training/side-quests/plugin_development nf-training nextflow run main.nf
```

### Clean up when done

```bash
docker stop nf-training && docker rm nf-training
```

---

## Working Directory Setup

| Tutorial Type        | Documentation                | Working Directory       | Solutions                         |
| -------------------- | ---------------------------- | ----------------------- | --------------------------------- |
| hello_nextflow       | `docs/hello_nextflow/`       | `hello-nextflow/`       | `hello-nextflow/solutions/`       |
| hello_nf-core        | `docs/hello_nf-core/`        | `hello-nf-core/`        | `hello-nf-core/solutions/`        |
| nf4_science/genomics | `docs/nf4_science/genomics/` | `nf4-science/genomics/` | `nf4-science/genomics/solutions/` |
| side_quests/\*       | `docs/side_quests/*.md`      | `side-quests/<name>/`   | `side-quests/solutions/<name>/`   |

---

## Walkthrough Process

### Phase 1: Preparation

1. **Read the lesson file** to understand the structure and sections
2. **Run `/validate` on the lesson file** - This checks:
   - Lesson structure and heading numbering
   - Code block `hl_lines` correctness
   - Pedagogical quality and clarity
   - Admonition syntax
3. **Identify starting files** - What files should already exist vs. what the user creates
4. **Prepare a clean working state** - Reset or backup existing files to avoid conflicts
5. **Start Docker container** (if using Docker)
6. **Verify prerequisites** - Check that required data files and configs exist

### Phase 2: Progressive Execution (Core Unique Value)

For each numbered section in the lesson (0, 1, 2, 3...):

#### 2.1 Read and Understand

- Read the section instructions
- Identify what code changes are expected (Before/After blocks, inline snippets)
- Note what commands should be run

#### 2.2 Verify Current State Matches "Before"

**This is the key unique check**: Before applying any code change:

1. Read the current file content
2. Compare with the "Before" block in the documentation
3. If they don't match, **STOP and report the discrepancy**
   - This catches issues where earlier steps didn't work correctly
   - This catches documentation that skipped steps

#### 2.3 Apply Code Changes Progressively

- **Don't copy the final solution** - Apply only the changes shown in this section
- Use the Edit tool to make incremental changes, mimicking how a learner would type
- If creating a new file, use the exact content shown at that point (not the final version)

#### 2.4 Run Commands Exactly As Shown

- Execute bash commands in the order they appear
- Use the exact syntax shown (don't "improve" or modify commands)
- Capture output for comparison

#### 2.5 Verify Outputs

- Compare console output with documented "Output" blocks
- **Acceptable differences** (don't flag):
  - Work directory hashes (random)
  - Run names (e.g., `[goofy_torvalds]`)
  - Timestamps
- **Flag as issues**:
  - Different text content
  - Missing or extra output lines
  - Failed commands
  - Missing expected files

### Phase 3: Final Verification

After completing all sections:

#### 3.1 Compare Final Code with Solutions

```bash
diff -u [built-file] [solution-file]
```

- Read the solution file from the solutions directory
- Report any differences between what was built and the solution
- Minor differences (whitespace, comments) may be acceptable

#### 3.2 Run `/test-example` on Solution Files

Invoke the `/test-example` skill on each solution file to verify:

- Fresh run works
- Resume functionality (processes should cache)
- Parameter handling (if applicable)
- Output matches documentation

#### 3.3 Cleanup (if walkthrough succeeded with no issues)

If the walkthrough completed successfully without issues, clean up the working directory:

```bash
# Reset modified files to their starting state
git checkout [working-directory]/

# Remove generated files and directories
rm -rf [working-directory]/work
rm -rf [working-directory]/.nextflow
rm -f [working-directory]/.nextflow.log*
rm -rf [working-directory]/[any-created-directories]  # e.g., nf-greeting for plugin tutorial
```

**Important**: Only clean up if no problems were encountered. If issues were found, leave the files in place for Phase 4.

### Phase 4: Propose Fixes and Create PR (if issues found)

If the walkthrough identified any fixable issues, offer to create a PR with the fixes.

#### 4.1 Identify Fixable Issues

Categorize issues into:

**Auto-fixable** (can fix programmatically):

- Extra/missing whitespace in code blocks
- Incorrect `hl_lines` values
- Heading numbering errors (use `--fix` flag)
- Minor formatting inconsistencies

**Requires manual review** (present to user for decision):

- Content accuracy issues
- Missing steps in documentation
- Incorrect command outputs
- Structural changes to lesson flow

#### 4.2 Present Proposed Changes to User

Before making any changes, clearly present:

1. **Summary of proposed fixes** - List each fix with:
   - File and line number
   - What the current content is
   - What the proposed fix would be
   - Why this fix is needed

2. **Ask for user approval** using AskUserQuestion:
   - "Do you want me to apply these fixes and create a PR?"
   - Options: "Yes, create PR", "Let me review/modify first", "No, skip PR"

3. **If user wants to review/modify**:
   - Apply fixes one at a time
   - After each fix, show the diff and ask if it's correct
   - Allow user to request modifications before proceeding

#### 4.3 Create the PR

Only after user approval:

1. **Create a new branch**:

   ```bash
   git checkout -b fix/[tutorial-name]-walkthrough-fixes
   ```

2. **Stage and commit changes**:

   ```bash
   git add [modified-files]
   git commit -m "Fix issues in [tutorial-name] tutorial

   - [list of fixes applied]

   Found during tutorial walkthrough testing.

   Co-Authored-By: Claude Opus 4.5 <noreply@anthropic.com>"
   ```

3. **Push and create PR**:

   ```bash
   git push -u origin fix/[tutorial-name]-walkthrough-fixes
   gh pr create --title "Fix issues in [tutorial-name] tutorial" --body "..."
   ```

4. **Return to original branch** and clean up working directory

#### 4.4 PR Body Format

```markdown
## Summary

Fixes issues found during tutorial walkthrough testing.

## Changes

- [File]: [description of fix]
- [File]: [description of fix]

## Issues Fixed

| Issue         | Severity      | Fix Applied        |
| ------------- | ------------- | ------------------ |
| [description] | Minor/Warning | [what was changed] |

## Testing

- [x] Tutorial walkthrough completed successfully after fixes
- [x] All commands execute as documented
- [x] Output matches documentation

🤖 Generated with [Claude Code](https://claude.ai/code)
```

---

## Output Format

```
# Tutorial Walkthrough: [lesson-name]

## Environment
- Mode: Docker (ghcr.io/nextflow-io/training:latest)
- Nextflow version: [from devcontainer.json]
- Working directory: /workspaces/training/[path]

## Section 0: [title]

### State Verification
- Current file matches "Before": ✓
[or]
- ⚠ State mismatch: Current file has X, "Before" shows Y

### Code Changes Applied
- [describe changes made]

### Commands Run
- `echo 'Hello World!'` → ✓ Output matches
- `nextflow run hello-world.nf` → ✓ Completed successfully

## Section 1: [title]
...

## Solution Comparison
- ✓ Final code matches solution
[or]
- ⚠ Differences:
  - Line 15: built has 'foo', solution has 'bar'

## Final Script Tests
- Fresh run: ✓
- Resume test: ✓ (all cached)
- Parameter test: ✓

## Validation Results (from /validate)

Report the actual results from the /validate skill:

- Heading numbering: ✓ 0 errors (or list specific issues with file:line)
- Code block highlights: ✓ All N blocks correct (or list incorrect hl_lines with what they highlight vs. what they should)
- TODO/FIXME comments: ✓ None found (or "Found N: list locations")
- Admonition syntax: ✓ All N properly formatted (or list malformed ones)
- Nextflow script conventions: ✓ All scripts have shebang and follow DSL2
- Orphaned files: ✓ All files referenced in mkdocs.yml

## Solution Tests (from /test-example)

Report the actual results from the /test-example skill:

- [script-name].nf:
  - Fresh run: ✓ Completed (N processes)
  - Resume test: ✓ All cached
  - Output verification: ✓ Files match documentation
- [other-script].nf: ✓ All tests passed

## Summary
- Sections completed: N/N
- State mismatches found: N
- Command failures: N
- Output discrepancies: N
- Solution match: ✓ / ⚠ (minor whitespace) / ✗ (significant differences)
- Validation issues: N
- Solution test failures: N

## Issues Found
[List any problems that would block or confuse a learner, organized by severity]

### Critical (blocks tutorial completion)
- [issue description]

### Warning (confusing but workable)
- [issue description]

### Minor (cosmetic or documentation-only)
- [issue description]

## Proposed Fixes
[If fixable issues were found, present them here before asking user]

| # | File | Line | Current | Proposed | Reason |
|---|------|------|---------|----------|--------|
| 1 | docs/side_quests/example.md | 123 | `hl_lines="1 11"` | `hl_lines="1"` | Line 11 is just a closing brace |
| 2 | ... | ... | ... | ... | ... |

**Ready to create PR?** [Ask user with AskUserQuestion before proceeding]
```

---

## When to Use Individual Skills Instead

This skill invokes `/validate` and `/test-example` automatically. Use individual skills when:

| Situation                                   | Use Instead           |
| ------------------------------------------- | --------------------- |
| Quick check of just highlights or structure | `/validate` alone     |
| Testing a script outside tutorial context   | `/test-example` alone |
| Finding TODO items across the codebase      | `/find-todos`         |

---

## Special Cases

### Exercises

When encountering `??? exercise` blocks:

1. Attempt the exercise based on instructions given
2. Check against `??? solution` block
3. Report if instructions are unclear or incomplete

### Branching/Optional Steps

Some tutorials have optional paths:

- Note which path was taken
- Ask user if they want to test alternative paths
- Document any untested branches

---

## Notes

- **Don't skip ahead** - The whole point is testing the progressive journey
- **State verification is critical** - This catches most tutorial bugs
- Reset files between full walkthrough runs to ensure clean state
- Save work directory paths for troubleshooting failed commands
