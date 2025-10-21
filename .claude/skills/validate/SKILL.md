---
name: Validate Training Materials
description: Run comprehensive validation and review checks including heading numbering, TODO/FIXME comments, Nextflow script conventions, orphaned files, admonition syntax, lesson structure, formatting, content accuracy, and teaching effectiveness. Use when validating, reviewing, or checking training materials quality, lesson quality, or before committing changes.
---

# Validate Training Materials

Run comprehensive validation checks on training materials to ensure quality and consistency. This includes both automated checks and deep lesson reviews.

## User Input

**IMPORTANT**: Ask the user to specify what to validate if not clear from the request:

- A specific side quest (e.g., "debugging", "metadata")
- A specific module (e.g., "hello_nextflow", "nf4_science/genomics")
- The entire repository (only if explicitly requested)

If the user provides a path or module name, **only validate that specific content**.

## Working Directory

**IMPORTANT**: All commands in this skill must be executed from the repository root directory.

- The repository root is the directory containing `mkdocs.yml`, `docs/`, and `.github/`
- Verify you are in the correct directory before running any commands (check for these files/folders)
- All file paths in this skill are relative to the repository root
- Do not change directories during skill execution
- Use paths relative to repository root only

## Determining Scope

Based on the user's request, determine the scope:

1. **Specific side quest** (e.g., "debugging"):

   - Documentation file: `docs/side_quests/{name}.md`
   - Example scripts: `side-quests/{name}/**/*.nf`
   - Solution files: `side-quests/solutions/{name}/**/*.nf`

2. **Specific module** (e.g., "hello_nextflow"):

   - Documentation files: `docs/{module}/**/*.md`
   - Example scripts: `{module}/**/*.nf` or `nf-training/**/*.nf`
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

   - Check if main documentation file is referenced in mkdocs.yml
   - Look for solution files without corresponding exercise documentation
   - Report any orphaned files within scope

5. **Verify Admonition Syntax**

   - Search for common admonition formatting errors in scoped files
   - Check for proper indentation (4 spaces)
   - Report any malformed admonitions

6. **Deep Lesson Review** (when reviewing a specific lesson file)

   If the user asks to review a specific lesson, perform this comprehensive checklist:

   a. **Structure Check**:

   - [ ] Proper heading numbering with trailing periods
   - [ ] Heading levels match numbering depth (## for 1., ### for 1.1.)
   - [ ] Each major section has "### Takeaway"
   - [ ] Each major section has "### What's next?"
   - [ ] Logical flow from simple to complex
   - [ ] Clear learning objectives stated or implied

   b. **Formatting Check**:

   - [ ] Code blocks have proper titles, line numbers, and highlighting
   - [ ] Console output properly formatted with `console title="Output"`
   - [ ] File paths use proper markdown formatting
   - [ ] Before/After comparisons use tabbed blocks
   - [ ] Admonitions properly formatted and indented
   - [ ] Each sentence on new line (for clean diffs)

   c. **Content Check**:

   - [ ] Technical accuracy of Nextflow code
   - [ ] Commands are correct and runnable
   - [ ] Parameter syntax correct (-- for pipeline, - for Nextflow)
   - [ ] Examples progress logically
   - [ ] Common pitfalls addressed
   - [ ] Edge cases explained

   d. **Teaching Effectiveness**:

   - [ ] Clear explanations for beginners
   - [ ] Concepts introduced before use
   - [ ] Examples are relevant and motivating
   - [ ] Exercises appropriate for skill level
   - [ ] Solutions available for exercises
   - [ ] Good use of tips and warnings

   e. **Cross-References**:

   - [ ] Links to related lessons work
   - [ ] References to files/scripts are correct
   - [ ] External links are valid
   - [ ] Prerequisites clearly stated

   f. **Examples & Code**:

   - [ ] All Nextflow examples are syntactically correct
   - [ ] Variable names are clear and consistent
   - [ ] Comments explain non-obvious code
   - [ ] Examples can be run as shown
   - [ ] Output examples match what code produces

## Output Format

Provide a structured report. For automated validation checks:

```
# Validation Report

## Heading Numbering
✓ All headings correctly numbered
[or list of issues with file:line]

## TODO/FIXME Comments
Found 15 total:
- High priority: 3
- Documentation: 8
- Code: 4

Top files:
- docs/side_quests/debugging.md (5 items)
- ...

## Nextflow Scripts
✓ All 23 scripts follow conventions
[or list of non-compliant scripts]

## Orphaned Files
Found 2 orphaned markdown files:
- docs/old/deprecated.md
- ...

## Admonitions
✓ All admonitions properly formatted
[or list of issues]

## Summary
[Overall assessment and recommended actions]
```

For deep lesson reviews, add:

```
## Lesson Review: [lesson-name]

### Structure
[Assessment of structure with specific issues/successes]

### Formatting
[Assessment of formatting with specific issues/successes]

### Content Accuracy
[Assessment of technical content]

### Teaching Effectiveness
[Assessment of pedagogical quality]

### Cross-References
[Assessment of links and references]

### Examples & Code
[Assessment of code quality]

### Overall Assessment
[Summary with severity-organized issues and recommendations]

### Positive Aspects
[Things worth preserving]
```

## Notes

- Use Grep and Glob tools for efficient searching
- Run heading check script directly via Bash
- Provide actionable next steps for any issues found
- If no issues found, give clear confirmation
