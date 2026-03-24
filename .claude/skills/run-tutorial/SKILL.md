---
name: Run Tutorial Walkthrough
description: Walk through a training tutorial as a user would, progressively building examples exactly as instructed, running all commands, and verifying results against solutions. Use when testing tutorials end-to-end or validating that instructions are correct and complete.
---

# Run Tutorial Walkthrough

Walk through a training tutorial lesson as a learner would, progressively building code, running commands, and verifying the learning journey works end-to-end.

---

## Operating Modes

This skill supports two distinct modes. **Ask the user which mode** if not clear from context:

### Mode A: Learner Simulation (Progressive Build)

Walk through as a learner would, editing starter files progressively. Use when:

- Testing that tutorial instructions are complete and correct
- Verifying the learning path works step-by-step

### Mode B: Documentation Verification

Run solution files and verify documentation matches actual outputs. Use when:

- Checking that code snippets, terminal outputs, and examples in docs are accurate
- Validating solution files work correctly

**Key difference:** In Mode A, you progressively edit starter files (then reset). In Mode B, you run solutions directly and check docs for accuracy.

---

## Starter Files

Starter files (e.g., `hello-nextflow/hello-config.nf`) represent the **end state of the previous chapter**. They are the starting point for learners beginning a new lesson.

### When to Edit Starter Files

| Scenario                                           | Action                  |
| -------------------------------------------------- | ----------------------- |
| Bug or syntax error inherited from previous chapter | ✅ Fix and commit        |
| Changes made while running through exercises       | ❌ Reset before committing |

### After Testing

If you modified starter files while walking through exercises, reset them:

```bash
git checkout hello-nextflow/*.nf
```

Only keep changes to:

- **Solution files** (`solutions/**/*`) - these should be correct and complete
- **Documentation** (`docs/**/*.md`) - fix code snippets, outputs, examples
- **Starter files with genuine bugs** - errors that exist before the lesson begins

---

## Skill Dependencies (MANDATORY)

This walkthrough MUST invoke other skills. **Do not skip these.**

| Checkpoint | Skill | Condition |
|------------|-------|-----------|
| Before Phase 2 | `/docker-setup` | If Docker environment selected |
| Before Phase 2 | `/validate` | Always - check lesson structure first |
| Phase 3 | `/test-example` | For each solution file |

---

## Initial Setup

Use `AskUserQuestion` to determine:

1. **Which tutorial/lesson** to walk through (specific file, module, or side quest)
2. **Environment**: Docker (recommended) or existing local environment

### Environment Setup

**If Docker selected:**

```
>>> STOP. INVOKE /docker-setup NOW.

Do not proceed until the container is confirmed running.
```

**If existing environment:** Verify Nextflow is installed with `nextflow -version`

### Strict Syntax Parser (v2)

**Always use the strict syntax parser** when running Nextflow commands:

```bash
NXF_SYNTAX_PARSER=v2 nextflow run ...
```

In Docker:

```bash
docker exec -e NXF_SYNTAX_PARSER=v2 nf-training nextflow run ...
```

This validates tutorials against the strict syntax that will become default in future Nextflow versions.

---

## Working Directory Mapping

See [../shared/repo-conventions.md](../shared/repo-conventions.md) for full directory mapping table.

Key pattern: documentation uses underscores (`hello_nextflow`), working directories use hyphens (`hello-nextflow`).

---

## Common Documentation Issues

When verifying documentation (especially in Mode B), watch for these frequent problems:

### Code Snippets

- **Non-existent properties/methods**: e.g., `.process` instead of `.name`
- **Parameter name mismatches**: e.g., `params.greeting` vs `params.input`
- **Missing path prefixes**: e.g., `'file.csv'` vs `'data/file.csv'`

### Configuration Examples

- **Profile parameter names** must match the actual params block
- **`nextflow config` output examples** must match what the config actually produces
- **Path values** in test profiles must be valid relative to the working directory

### Terminal Output Examples

- Version numbers (acceptable to differ from your run)
- Work directory hashes (acceptable to differ from your run)
- Process execution order (acceptable to differ)
- **Error messages** (must match if showing expected errors)

### Work Directory Hashes (within docs - CRITICAL)

While your actual hashes will differ from documented ones, check that docs are internally consistent:

- **Invalid hex**: Hashes like `[j6/abc123]` contain invalid characters (only 0-9, a-f allowed)
- **Copy-paste duplicates**: Multiple different `nextflow run` commands showing identical hashes (unless cached)
- **Prose mismatch**: Text referencing a hash that doesn't appear in nearby output
- **Cache confusion**: Non-cached processes (no "cached: N") should have unique hashes per run

### Before/After Code Blocks

- Verify "Before" matches the actual prior state
- Verify "After" produces working code
- Check that `hl_lines` highlight the correct lines (count from line 1 of snippet, not `linenums`)

---

## Walkthrough Process

Use `TodoWrite` to track progress through phases and sections.

### Phase 1: Preparation

1. **Read the lesson file** to understand the structure and sections
2. **Identify starting files** - What files should already exist vs. what the user creates
3. **Prepare a clean working state** - Reset or backup existing files to avoid conflicts
4. **Verify prerequisites** - Check that required data files and configs exist

#### CHECKPOINT: Validate Before Proceeding

```
>>> STOP. INVOKE /validate on the lesson file NOW.

Record any issues found. Do not proceed to Phase 2 until validation completes.
```

### Phase 2: Progressive Execution

Complete steps 2.1-2.6 for each section before moving to the next. Do not batch sections together.

For each numbered section in the lesson:

#### 2.1 Read and Understand

- Read the section instructions
- Identify what code changes are expected (Before/After blocks, inline snippets)
- Note what commands should be run

#### 2.2 Verify Current State Matches "Before"

Before applying any code change:

1. Read the current file content
2. Compare with the "Before" block in the documentation
3. If they don't match, stop and report the discrepancy

#### 2.3 Apply Code Changes Progressively

- **Don't copy the final solution** - Apply only the changes shown in this section
- Use the `Edit` tool to make incremental changes, mimicking how a learner would type
- If creating a new file, use the exact content shown at that point (not the final version)

#### 2.4 Run Commands Exactly As Shown

- Execute bash commands in the order they appear using the `Bash` tool
- Use the exact syntax shown (don't "improve" or modify commands)
- Capture output for comparison
- **If a section changes code but doesn't include a "run this" instruction, run the workflow anyway** to verify changes work

#### 2.5 Verify Outputs

- Compare console output with documented "Output" blocks
- See [references/acceptable-differences.md](references/acceptable-differences.md) for what to flag vs. ignore

#### 2.6 Confirm Section Complete

Before moving to the next section:

1. Verify the workflow/script ran successfully (exit code 0)
2. Confirm output matches documentation (within acceptable differences)
3. Update `TodoWrite` to mark section complete

### Phase 2B: Solution Testing (Mode B only)

For documentation verification mode, skip progressive editing of starter files. Instead:

#### 2B.1 Run Each Solution Directly

```bash
nextflow run solutions/<lesson>/script.nf -c solutions/<lesson>/nextflow.config
```

#### 2B.2 Test with Profiles (if applicable)

```bash
nextflow run solutions/<lesson>/script.nf -c solutions/<lesson>/nextflow.config -profile test
nextflow run solutions/<lesson>/script.nf -c solutions/<lesson>/nextflow.config -profile my_laptop,test
```

#### 2B.3 Verify Config Resolution

```bash
nextflow config -profile <profiles>
```

Compare output against documented `nextflow config` examples in the lesson.

#### 2B.4 Compare All Outputs Against Documentation

For each documented output block:

1. Find the corresponding section in the docs
2. Run the command/workflow
3. Compare actual vs documented output
4. Flag discrepancies (use acceptable-differences.md to filter noise)

#### 2B.5 Verify Hash Consistency (MANDATORY)

Run these checks on ALL documentation files in scope:

```bash
# 1. Check for invalid hex characters in hashes
grep -oE '\[[^]]+\]' docs/path/to/*.md | grep -E '\[[^0-9a-f/\]]' | head -20

# 2. Find duplicate hashes (may indicate copy-paste errors)
grep -oE '\[[a-z0-9]{2}/[a-z0-9]{6}\]' docs/path/to/*.md | sort | uniq -c | sort -rn | head -20
```

For any duplicates found, verify they are legitimate:
- Same output shown twice (acceptable)
- `-resume` showing cached results from earlier in same lesson (acceptable)
- Prose reference to nearby terminal output (acceptable)
- Different `nextflow run` commands without `-resume` (NOT acceptable - flag as issue)

See [references/acceptable-differences.md](references/acceptable-differences.md) for detailed hash validation rules.

### Phase 3: Final Verification

After completing all sections:

#### 3.1 Compare Final Code with Solutions

Read the solution file and compare with what was built. Minor differences (whitespace, comments) may be acceptable. Report any significant differences.

#### 3.2 Test Solution Files

```
>>> STOP. INVOKE /test-example on each solution file.

For each solution script, verify:
- Fresh run works
- Resume caches correctly
- Parameter handling (if applicable)
- Output matches documentation

Record results before proceeding.
```

#### 3.3 Cleanup

**Only if walkthrough succeeded with no issues:**

- Reset modified files: `git checkout "${WORKING_DIR}/"`
- Remove work directories and logs

**If issues found:** Leave files in place for Phase 4. Save work directory paths for troubleshooting.

### Phase 4: Propose Fixes and Create PR (if issues found)

If fixable issues were identified, follow [references/pr-workflow.md](references/pr-workflow.md) to create a PR.

Key points:
- Categorize as auto-fixable vs. requires review
- Present fixes to user with **actual section headings** (read from document - do not guess!)
- Get user approval via `AskUserQuestion` before making changes

---

## Output Format

Structure your report with these sections:

```
# Tutorial Walkthrough: <lesson-name>

## Environment
Mode, Nextflow version, working directory

## Section N: <exact title from document>
- State verification: matches "Before" ✓ or mismatch details
- Code changes applied
- Commands run with results

## Validation Results (from /validate)
## Solution Tests (from /test-example)
## Solution Comparison

## Summary
Sections completed, mismatches, failures, discrepancies

## Issues Found
Categorize as Critical / Warning / Minor

## Proposed Fixes (if any)
Table with file, line, current, proposed, reason
```

### Mode B Output Format

```
# Documentation Verification: <module-name>

## Environment
Container name, Nextflow version, working directory

## Solutions Tested
| Solution | Base Run | With Profiles | Issues |
|----------|----------|---------------|--------|
| 1-hello-world | ✓ | N/A | None |
| 6-hello-config | ✓ | ✓ test | Fixed: .process→.name |

## Hash Consistency Check
- Invalid hex characters found: [list or "None"]
- Duplicate hashes across different runs: [list or "None - all legitimate"]
- Prose/output mismatches: [list or "None"]

## Documentation Issues Found

### Critical (breaks functionality)
- File: `06_hello_config.md`, Lines 893-931
  - Issue: `.process` property doesn't exist, should be `.name`
  - Also affects: `solutions/6-hello-config/hello-config.nf`

### Minor (cosmetic/clarity)
- ...

## Fixes Applied
| File | Change | Reason |
|------|--------|--------|
| docs/.../06_hello_config.md | `.process` → `.name` | Property doesn't exist |
| solutions/.../nextflow.config | `params.greeting` → `params.input` | Wrong parameter name |
| docs/.../03_config.md | `[f0/35723c]` → `[2b/9a7d1e]` | Duplicate hash for different run |
```

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
