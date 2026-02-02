---
name: Validate Training Materials
description: Run comprehensive validation and review checks including heading numbering, TODO/FIXME comments, Nextflow script conventions, orphaned files, admonition syntax, lesson structure, formatting, content accuracy, and teaching effectiveness. Use when validating, reviewing, or checking training materials quality, lesson quality, or before committing changes.
---

# Validate Training Materials

Run comprehensive validation checks on training materials. Execute from repository root.

See [../shared/repo-conventions.md](../shared/repo-conventions.md) for directory mapping and file conventions.

## Skill Dependencies (MANDATORY)

This skill MUST invoke other skills during validation. **Do not skip these.**

| Task                         | Skill                | When   |
| ---------------------------- | -------------------- | ------ |
| Check code block highlights  | `/check-highlights`  | Always |
| Check inline code formatting | `/check-inline-code` | Always |

## Scope

Ask user to specify what to validate if not clear:

- Specific side quest (e.g., "debugging")
- Specific module (e.g., "hello_nextflow")
- Entire repository (only if explicitly requested)

## Determining Scope

Based on the user's request, determine the scope:

1. **Specific side quest** (e.g., "debugging"):
   - Documentation file: `docs/side_quests/{name}.md`
   - Example scripts: `side-quests/{name}/**/*.nf`
   - Solution files: `side-quests/solutions/{name}/**/*.nf`

2. **Specific module** (e.g., "hello_nextflow"):
   - Documentation files: `docs/{module}/**/*.md`
   - Example scripts: `{module}/**/*.nf` or `hello-nextflow/**/*.nf`
   - Solution files: `{module}/solutions/**/*.nf`

3. **Entire repository** (only if explicitly requested):
   - All files: `docs/**/*.md`, `**/*.nf`

## Tasks to Perform

Perform the following checks **only on files within the determined scope**:

1. **Check Heading Numbering**
   - Run: `uv run .github/check_headings.py [scoped-path]/**/*.md`
   - Report any heading numbering issues found
   - If errors exist, ask if user wants to auto-fix with `--fix` flag

2. **Find TODO/FIXME Comments**
   - Search markdown files in scope
   - Search Nextflow scripts in scope
   - Categorize by priority (high, medium, low)
   - Report files with most TODOs

3. **Check Nextflow Script Conventions**
   - Find .nf files in scope
   - Verify they start with `#!/usr/bin/env nextflow`
   - Check for DSL2 syntax
   - Report any that don't follow conventions

4. **Find Orphaned Files**
   - Check if main documentation file is referenced in `docs/en/mkdocs.yml`
   - Look for solution files without corresponding exercise documentation
   - Report any orphaned files within scope

5. **Verify Admonition Syntax**
   - Search for common admonition formatting errors in scoped files
   - Check for proper indentation (4 spaces)
   - Report any malformed admonitions

6. **Check Code Block Highlights**

   ```
   >>> STOP. INVOKE /check-highlights on the scoped files NOW.

   Record results before continuing.
   ```

7. **Check Inline Code Formatting**

   ```
   >>> STOP. INVOKE /check-inline-code on the scoped files NOW.

   Record results before continuing.
   ```

8. **Check Writing Style**
   - Search for LLM-style patterns that should be avoided:
     - `Let's` or `let's` at start of sentences
     - `Remember when` or `Remember that` callbacks
     - `Don't worry` reassurances
     - `worth mentioning` or `important to note` padding
     - Exclamation marks used for emphasis (not in code/output)
   - Check for em-dash elaborations (space-hyphen-space followed by lowercase) that could be periods
   - Flag any issues found for manual review

9. **Deep Lesson Review** (when reviewing a specific lesson file)

   If the user asks to review a specific lesson, use the checklist in [references/deep-review-checklist.md](references/deep-review-checklist.md).

## Output Format

Structure your report with these sections:

```
# Validation Report

## Heading Numbering
## TODO/FIXME Comments (count by priority, list top files)
## Nextflow Scripts (conventions check)
## Orphaned Files
## Admonitions
## Code Block Highlights (from /check-highlights)
## Inline Code Formatting (from /check-inline-code)
## Writing Style
## Summary
```

For deep lesson reviews, add sections for: Structure, Formatting, Content Accuracy, Teaching Effectiveness, Cross-References, Examples & Code, Writing Style, Overall Assessment, Positive Aspects.

**Important**: Always provide actionable next steps for any issues found. If no issues found, give clear confirmation.
