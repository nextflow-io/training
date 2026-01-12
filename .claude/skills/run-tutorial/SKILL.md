---
name: Run Tutorial Walkthrough
description: Walk through a training tutorial as a user would, progressively building examples exactly as instructed, running all commands, and verifying results against solutions. Use when testing tutorials end-to-end or validating that instructions are correct and complete.
---

# Run Tutorial Walkthrough

Walk through a training tutorial lesson as a learner would, progressively building code, running commands, and verifying the learning journey works end-to-end.

**This skill focuses on the unique value of simulating a learner's experience.** It automatically invokes other skills during the walkthrough:

- **Lesson validation** → Invokes `/validate` skill (includes highlight checking, structure, formatting)
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

### Docker-outside-of-Docker (DooD) for Containerized Processes

Many tutorials use Nextflow processes with `container` directives (e.g., FASTP, BWA, SAMTOOLS). In Codespaces/Gitpod, Docker-in-Docker is pre-configured. For **local testing**, you need Docker-outside-of-Docker (DooD) to run containerized processes.

**When DooD is needed**: Any tutorial where processes specify containers, including:

- `hello_nextflow` (later lessons with containers)
- `nf4_science/genomics` and other domain modules
- Side quests like `essential_scripting_patterns`, `metadata`, etc.

#### DooD Container Setup

The key differences from basic setup:

1. **Mount the Docker socket** to allow Nextflow to spawn sibling containers
2. **Use matching host paths** so work directories resolve correctly between containers

```bash
# Clean up any existing container
docker stop nf-training 2>/dev/null; docker rm nf-training 2>/dev/null

# Get NXF_VER and host path
NXF_VER=$(grep -o '"NXF_VER":\s*"[^"]*"' .devcontainer/devcontainer.json | cut -d'"' -f4)
HOST_PATH="${PWD}"

# Start container with DooD support
docker run -d --name nf-training \
  -e NXF_VER=${NXF_VER} \
  -e LANG=C.UTF-8 \
  -e LC_ALL=C.UTF-8 \
  -v /var/run/docker.sock:/var/run/docker.sock \
  -v "${HOST_PATH}:${HOST_PATH}" \
  -w "${HOST_PATH}" \
  ghcr.io/nextflow-io/training:latest \
  sleep infinity
```

**Critical**: The `-v "${HOST_PATH}:${HOST_PATH}"` mount ensures paths match between the training container and any containers Nextflow spawns. Without this, you'll see `.command.sh: No such file or directory` errors.

#### Running Commands with DooD

```bash
docker exec -e LANG=C.UTF-8 -e LC_ALL=C.UTF-8 -e USER=testuser \
  -w ${HOST_PATH}/[working-dir] \
  nf-training \
  nextflow run [script.nf] [options]
```

#### Apple Silicon (ARM) Macs: Platform Emulation

Most bioinformatics containers are built for x86_64/amd64. On ARM Macs, you need platform emulation:

1. Create a Nextflow config with platform option:

```bash
docker exec nf-training bash -c 'cat > /tmp/platform.config << EOF
docker.runOptions = "--platform linux/amd64"
EOF'
```

2. Include the config when running:

```bash
docker exec -e LANG=C.UTF-8 -e LC_ALL=C.UTF-8 -e USER=testuser \
  -w ${HOST_PATH}/side-quests/essential_scripting_patterns \
  nf-training \
  nextflow run main.nf -c /tmp/platform.config
```

**Note**: Platform emulation uses more memory. If you encounter OOM errors (exit code 137), increase Docker Desktop's memory allocation in Preferences → Resources.

#### Troubleshooting DooD

| Error                                    | Cause               | Solution                                             |
| ---------------------------------------- | ------------------- | ---------------------------------------------------- |
| `Cannot connect to Docker daemon`        | Socket not mounted  | Add `-v /var/run/docker.sock:/var/run/docker.sock`   |
| `.command.sh: No such file or directory` | Path mismatch       | Use matching paths: `-v "${HOST_PATH}:${HOST_PATH}"` |
| `exec format error`                      | ARM/x86 mismatch    | Add `--platform linux/amd64` to docker.runOptions    |
| Exit code 137 (OOM)                      | Insufficient memory | Increase Docker Desktop memory allocation            |

#### Symlink Trick for Codespaces Paths

Training materials often contain hardcoded paths like `/workspaces/training/...` which work in Codespaces/Gitpod but fail when testing locally with DooD (where paths are like `/Users/username/projects/training/...`).

**Solution**: Create a symlink inside the container so Codespaces paths resolve to your host paths:

```bash
# After starting the container, create the symlink
docker exec nf-training bash -c "rm -rf /workspaces/training && mkdir -p /workspaces && ln -sf ${HOST_PATH} /workspaces/training"

# Verify it works
docker exec nf-training ls /workspaces/training/
```

**Why this works**: Files like `sample_bams.txt` contain paths like `/workspaces/training/nf4-science/genomics/data/bam/reads_mother.bam`. The symlink makes these paths resolve correctly without modifying the training materials.

**Important**: The symlink must be recreated each time you restart the container.

#### Container Stability Issues

The `nf-training` container may stop unexpectedly during long tutorial walkthroughs, especially when running many containerized processes via DooD.

**Symptoms**:

- `Error: No such container: nf-training`
- Commands that worked before suddenly fail

**Recovery procedure**:

```bash
# 1. Check if container is still running
docker ps | grep nf-training

# 2. If not running, restart it
docker stop nf-training 2>/dev/null; docker rm nf-training 2>/dev/null

HOST_PATH="/path/to/training"  # Use your actual path
NXF_VER=$(grep -o '"NXF_VER":\s*"[^"]*"' .devcontainer/devcontainer.json | cut -d'"' -f4)

docker run -d --name nf-training \
  -e NXF_VER=${NXF_VER} \
  -e LANG=C.UTF-8 \
  -e LC_ALL=C.UTF-8 \
  -v /var/run/docker.sock:/var/run/docker.sock \
  -v "${HOST_PATH}:${HOST_PATH}" \
  -w "${HOST_PATH}" \
  ghcr.io/nextflow-io/training:latest \
  sleep infinity

# 3. Recreate the symlink (critical!)
docker exec nf-training bash -c "rm -rf /workspaces/training && mkdir -p /workspaces && ln -sf ${HOST_PATH} /workspaces/training"

# 4. Recreate platform config if needed (for ARM Macs)
docker exec nf-training bash -c 'cat > /tmp/platform.config << EOF
docker.runOptions = "--platform linux/amd64"
EOF'
```

**Tip**: When running long tutorials, periodically check that the container is still running before executing commands.

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

**CRITICAL**: Complete ALL steps (2.1-2.6) for EACH section before moving to the next section. Do NOT batch multiple sections together. The entire point of this skill is to catch issues that only appear when following the tutorial step-by-step.

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

#### 2.6 Confirm Section Complete Before Proceeding

**MANDATORY**: Before moving to the next section:

1. Verify the workflow/script ran successfully (exit code 0)
2. Confirm output matches documentation (within acceptable differences)
3. Only then proceed to section N+1

**Why this matters**: Skipping incremental testing defeats the purpose of this skill. A tutorial might work when you jump to the final solution but fail at intermediate steps - exactly the bugs learners encounter.

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
   - **Actual section heading** (read from the document, don't guess!)
   - What the current content is
   - What the proposed fix would be
   - Why this fix is needed

**CRITICAL**: Always verify section numbers by reading the document. Do NOT guess or infer section numbers. Search for the nearest `### N.N.` heading above the line you're referencing and use that exact text.

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

**IMPORTANT**: Before writing the PR body, read the document to verify the actual section headings for each fix. Do not guess section numbers.

```markdown
## Summary

Fixes issues found during tutorial walkthrough testing.

## Changes

| Location                           | Before    | After     | Reason      |
| ---------------------------------- | --------- | --------- | ----------- |
| Section X.Y (exact title from doc) | old value | new value | explanation |

## Testing

- [x] Tutorial walkthrough completed successfully after fixes
- [x] All commands execute as documented
- [x] Output matches documentation

🤖 Generated with [Claude Code](https://claude.ai/code)
```

---

## Output Format

**IMPORTANT**: When referencing sections in your report, always read the actual section heading from the document. Do not guess or infer section numbers - search for `### N.N.` headings near the relevant line numbers.

```
# Tutorial Walkthrough: [lesson-name]

## Environment
- Mode: Docker (ghcr.io/nextflow-io/training:latest)
- Nextflow version: [from devcontainer.json]
- Working directory: /workspaces/training/[path]

## Section 0: [exact title from document]

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

- **ONE SECTION AT A TIME** - Complete steps 2.1-2.6 for each section before moving to the next. Never batch multiple sections together.
- **Test after EVERY change** - Run the workflow after each section's code changes. This is non-negotiable.
- **Don't skip ahead** - The whole point is testing the progressive journey as a learner experiences it
- **State verification is critical** - This catches most tutorial bugs
- Reset files between full walkthrough runs to ensure clean state
- Save work directory paths for troubleshooting failed commands
- If a section doesn't include a "run this command" instruction but changes code, run the workflow anyway to verify the changes work
