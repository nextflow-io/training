---
name: Run Tutorial Walkthrough
description: Walk through a training tutorial as a user would, progressively building examples exactly as instructed, running all commands, and verifying results against solutions. Use when testing tutorials end-to-end or validating that instructions are correct and complete.
---

# Run Tutorial Walkthrough

Walk through a training tutorial lesson as a learner would, progressively building code, running commands, and verifying the learning journey works end-to-end.

**This skill focuses on the unique value of simulating a learner's experience.** It invokes other skills during the walkthrough:

- **Lesson validation** → Use `Skill` tool with `skill="validate"`
- **Testing final scripts** → Use `Skill` tool with `skill="test-example"`
- **Docker setup** → Use `Skill` tool with `skill="docker-setup"`

---

## Initial Setup Questions

**IMPORTANT**: Before starting the walkthrough, use `AskUserQuestion` to gather required information.

### Question 1: Which tutorial/lesson?

```
AskUserQuestion:
  question: "Which tutorial or lesson should I walk through?"
  header: "Tutorial"
  options:
    - label: "Specific lesson file"
      description: "e.g., docs/hello_nextflow/01_hello_world.md"
    - label: "Entire module"
      description: "e.g., hello_nextflow, nf4_science/genomics"
    - label: "Side quest"
      description: "e.g., debugging, metadata, plugin_development"
```

### Question 2: Which environment?

```
AskUserQuestion:
  question: "Which software environment should I use?"
  header: "Environment"
  options:
    - label: "Docker container (Recommended)"
      description: "Uses training Docker image - matches Codespaces/Gitpod"
    - label: "Existing environment"
      description: "Use current environment as-is (you've already set it up)"
```

**Why Docker is recommended**: The training materials are designed for Codespaces/Gitpod which uses `ghcr.io/nextflow-io/training:latest`. This image has Java 21, all dependencies, and controlled Nextflow version.

**If Docker selected**: Invoke `/docker-setup` skill using `Skill` tool to configure the container.

**If Existing environment selected**: Skip Docker setup. Verify Nextflow is installed: `nextflow -version`

---

## Working Directory Mapping

| Tutorial Type        | Documentation                | Working Directory       | Solutions                         |
| -------------------- | ---------------------------- | ----------------------- | --------------------------------- |
| hello_nextflow       | `docs/hello_nextflow/`       | `hello-nextflow/`       | `hello-nextflow/solutions/`       |
| hello_nf-core        | `docs/hello_nf-core/`        | `hello-nf-core/`        | `hello-nf-core/solutions/`        |
| nf4_science/genomics | `docs/nf4_science/genomics/` | `nf4-science/genomics/` | `nf4-science/genomics/solutions/` |
| side_quests/\*       | `docs/side_quests/*.md`      | `side-quests/<name>/`   | `side-quests/solutions/<name>/`   |

---

## Walkthrough Process

Use `TodoWrite` to track progress through these phases and sections.

### Phase 1: Preparation

1. **Read the lesson file** to understand the structure and sections
2. **Run validation** - Use `Skill` tool with `skill="validate"` on the lesson file
3. **Identify starting files** - What files should already exist vs. what the user creates
4. **Prepare a clean working state** - Reset or backup existing files to avoid conflicts
5. **Set up environment** - If Docker, invoke `/docker-setup` skill
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

1. Read the current file content using the `Read` tool
2. Compare with the "Before" block in the documentation
3. If they don't match, **STOP and report the discrepancy**
   - This catches issues where earlier steps didn't work correctly
   - This catches documentation that skipped steps

#### 2.3 Apply Code Changes Progressively

- **Don't copy the final solution** - Apply only the changes shown in this section
- Use the `Edit` tool to make incremental changes, mimicking how a learner would type
- If creating a new file, use the exact content shown at that point (not the final version)

#### 2.4 Run Commands Exactly As Shown

- Execute bash commands in the order they appear using the `Bash` tool
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
3. Update `TodoWrite` to mark section complete
4. Only then proceed to section N+1

**Why this matters**: Skipping incremental testing defeats the purpose of this skill. A tutorial might work when you jump to the final solution but fail at intermediate steps - exactly the bugs learners encounter.

### Phase 3: Final Verification

After completing all sections:

#### 3.1 Compare Final Code with Solutions

```bash
diff -u <built-file> <solution-file>
```

- Read the solution file from the solutions directory
- Report any differences between what was built and the solution
- Minor differences (whitespace, comments) may be acceptable

#### 3.2 Run Test on Solution Files

Use `Skill` tool with `skill="test-example"` on each solution file to verify:

- Fresh run works
- Resume functionality (processes should cache)
- Parameter handling (if applicable)
- Output matches documentation

#### 3.3 Cleanup (if walkthrough succeeded with no issues)

If the walkthrough completed successfully without issues, clean up:

```bash
# Determine working directory (replace with actual path)
WORKING_DIR="hello-nextflow"  # or side-quests/plugin_development, etc.

# Reset modified files to their starting state
git checkout "${WORKING_DIR}/"

# Remove generated files and directories
rm -rf "${WORKING_DIR}/work"
rm -rf "${WORKING_DIR}/.nextflow"
rm -f "${WORKING_DIR}/.nextflow.log"*
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

2. **Ask for user approval** using `AskUserQuestion`:

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
   git checkout -b fix/<tutorial-name>-walkthrough-fixes
   ```

2. **Stage and commit changes**:

   ```bash
   git add <modified-files>
   git commit -m "Fix issues in <tutorial-name> tutorial

   - <list of fixes applied>

   Found during tutorial walkthrough testing.

   Co-Authored-By: Claude <noreply@anthropic.com>"
   ```

3. **Push and create PR**:

   ```bash
   git push -u origin fix/<tutorial-name>-walkthrough-fixes
   gh pr create --title "Fix issues in <tutorial-name> tutorial" --body "..."
   ```

4. **Return to original branch** and clean up working directory

#### 4.4 PR Body Format

**IMPORTANT**: Before writing the PR body, read the document to verify the actual section headings for each fix.

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

Generated with [Claude Code](https://claude.ai/code)
```

---

## Output Format

**IMPORTANT**: When referencing sections in your report, always read the actual section heading from the document. Do not guess or infer section numbers.

```
# Tutorial Walkthrough: <lesson-name>

## Environment
- Mode: Docker / Local
- Nextflow version: <from devcontainer.json or local>
- Working directory: <path>

## Section 0: <exact title from document>

### State Verification
- Current file matches "Before": ✓
[or]
- State mismatch: Current file has X, "Before" shows Y

### Code Changes Applied
- <describe changes made>

### Commands Run
- `echo 'Hello World!'` → ✓ Output matches
- `nextflow run hello-world.nf` → ✓ Completed successfully

## Section 1: <title>
...

## Solution Comparison
- ✓ Final code matches solution
[or]
- Differences:
  - Line 15: built has 'foo', solution has 'bar'

## Final Script Tests
- Fresh run: ✓
- Resume test: ✓ (all cached)
- Parameter test: ✓

## Validation Results (from /validate)
- Heading numbering: ✓ 0 errors
- Code block highlights: ✓ All N blocks correct
- TODO/FIXME comments: ✓ None found
- Admonition syntax: ✓ All properly formatted

## Solution Tests (from /test-example)
- <script-name>.nf: ✓ All tests passed

## Summary
- Sections completed: N/N
- State mismatches found: N
- Command failures: N
- Output discrepancies: N
- Solution match: ✓ / (minor whitespace) / ✗ (significant)

## Issues Found

### Critical (blocks tutorial completion)
- <issue description>

### Warning (confusing but workable)
- <issue description>

### Minor (cosmetic or documentation-only)
- <issue description>

## Proposed Fixes
[If fixable issues were found]

| # | File | Line | Current | Proposed | Reason |
|---|------|------|---------|----------|--------|
| 1 | docs/side_quests/example.md | 123 | `hl_lines="1 11"` | `hl_lines="1"` | Line 11 is just a closing brace |

**Ready to create PR?** [Use AskUserQuestion before proceeding]
```

---

## When to Use Individual Skills Instead

This skill invokes `/validate` and `/test-example` automatically. Use individual skills when:

| Situation                                   | Use Instead           |
| ------------------------------------------- | --------------------- |
| Quick check of just highlights or structure | `/validate` alone     |
| Testing a script outside tutorial context   | `/test-example` alone |
| Finding TODO items across the codebase      | `/find-todos`         |
| Just need Docker container setup            | `/docker-setup` alone |

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

- **Use TodoWrite** - Track progress through phases and sections with the `TodoWrite` tool
- **ONE SECTION AT A TIME** - Complete steps 2.1-2.6 for each section before moving to the next
- **Test after EVERY change** - Run the workflow after each section's code changes
- **Don't skip ahead** - The whole point is testing the progressive journey as a learner experiences it
- **State verification is critical** - This catches most tutorial bugs
- Reset files between full walkthrough runs to ensure clean state
- Save work directory paths for troubleshooting failed commands
- If a section doesn't include a "run this command" instruction but changes code, run the workflow anyway to verify the changes work
