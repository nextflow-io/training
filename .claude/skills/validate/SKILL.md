---
name: Validate Training Materials
description: Run comprehensive validation and review checks including the site build and links validator, heading numbering, code block meta, TODO/FIXME comments, Nextflow script conventions, orphaned pages, MDX syntax, lesson structure, formatting, content accuracy, and teaching effectiveness. Use when validating, reviewing, or checking training materials quality, lesson quality, or before committing changes.
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

Based on the user's request, determine the scope. `<pages>` is `docs/src/content/docs`.

1. **Specific side quest** (e.g., "debugging"):
   - Documentation files: `<pages>/side_quests/{name}/*.mdx`
   - Example scripts: `side-quests/{name}/**/*.nf`
   - Solution files: `side-quests/solutions/{name}/**/*.nf` or `side-quests/{name}/solutions/**/*.nf`

2. **Specific module** (e.g., "hello_nextflow"):
   - Documentation files: `<pages>/{module}/**/*.mdx`
   - Example scripts: the module's working directory, e.g. `hello-nextflow/**/*.nf`
   - Solution files: `{working-dir}/solutions/**/*.nf`

3. **Entire repository** (only if explicitly requested):
   - All English pages: `<pages>/**/*.mdx`, excluding the translation folders `<pages>/{ca,de,es,fr,hi,it,ko,pl,pt,tr}/`
   - All scripts: `**/*.nf`

Never check or report on translated pages. They are generated.

## Tasks to Perform

Perform the following checks **only on files within the determined scope**, except the build, which always covers the whole site:

1. **Build the Site (with links validator)**
   - Run: `cd docs && npm ci` (first time only), then `npm run build`. It takes about 80 s and 6 GB RAM
   - The build fails on invalid MDX, invalid frontmatter (schema in `docs/src/content.config.ts`) and broken internal links or anchors. Report each error with its file and line
   - Search the build output for `[WARN] [starlight-codeblocks]` and `[WARN] [astro-expressive-code]`. CI fails on these (unknown directive, bad range)

2. **Check Heading Numbering**
   - Run: `uv run .github/check_headings.py [scoped-files]`, e.g. `uv run .github/check_headings.py docs/src/content/docs/hello_nextflow/*.mdx`
   - For a recursive scope: `find [scoped-dir] -name '*.mdx' -print0 | xargs -0 uv run .github/check_headings.py`
   - Report any heading numbering issues found
   - If errors exist, ask if user wants to auto-fix with `--fix` flag

3. **Check Code Block Meta**
   - Run: `node docs/scripts/check-fences.mjs [scoped-dir]`, e.g. `node docs/scripts/check-fences.mjs docs/src/content/docs/hello_nextflow`
   - Report errors (they fail CI) and summarise warnings

4. **Find TODO/FIXME Comments**
   - Search pages in scope (`TODO` usually sits in an MDX comment, `{/* TODO ... */}`)
   - Search Nextflow scripts in scope
   - Categorize by priority (high, medium, low)
   - Report files with most TODOs

5. **Check Nextflow Script Conventions**
   - Find .nf files in scope
   - Verify they start with `#!/usr/bin/env nextflow` (module files under `modules/` must not)
   - Check for DSL2 syntax
   - Report any that don't follow conventions

6. **Find Orphaned Files**
   - Check that every page in scope is listed in `docs/src/sidebar.mjs` (as a slug without `.mdx`, e.g. `"hello_nextflow/01_hello_world"`; a module index is listed by its directory slug)
   - Look for solution files without corresponding exercise documentation
   - Report any orphaned files within scope

7. **Check Page Syntax**
   - Every page has `title:` and `description:` frontmatter and no `# ` H1 in the body
   - Leftover MkDocs syntax that renders literally: `!!! `, `??? `, `=== "`, `hl_lines=`, `linenums=`, `--8<--`, `` `#!groovy ``. Replace with the MDX forms in CLAUDE.md
   - `<details>` blocks have a blank line after `<summary>` and before `</details>`; `<CodeWalkthrough>` has a blank line after the opening tag and before the closing tag
   - `:::` asides are closed; `:::code-tabs` contains only code blocks
   - Nextflow code blocks use ` ```nextflow `, not ` ```groovy ` (real Groovy stays `groovy`)
   - Command output sits in the `bash` block after `# [!output]`, not in a separate block

8. **Check Code Block Highlights**

   ```
   >>> STOP. INVOKE /check-highlights on the scoped files NOW.

   Record results before continuing.
   ```

9. **Check Inline Code Formatting**

   ```
   >>> STOP. INVOKE /check-inline-code on the scoped files NOW.

   Record results before continuing.
   ```

10. **Check Writing Style**
    - Search for LLM-style patterns that should be avoided:
      - `Remember when` or `Remember that` callbacks
      - `Don't worry` reassurances
      - Exclamation marks used for emphasis (not in code/output)
    - Check for em-dash elaborations (space-hyphen-space followed by lowercase) that could be periods
    - Flag any issues found for manual review
    - **Do NOT flag**: `Let's` / `let's`, `worth mentioning`, or `important to note` — these are acceptable in tutorial prose

11. **Deep Lesson Review** (when reviewing a specific lesson file)

    If the user asks to review a specific lesson, use the checklist in [references/deep-review-checklist.md](references/deep-review-checklist.md).

## Output Format

Structure your report with these sections:

```
# Validation Report

## Build and Links
## Heading Numbering
## Code Block Meta (from check-fences)
## TODO/FIXME Comments (count by priority, list top files)
## Nextflow Scripts (conventions check)
## Orphaned Files
## Page Syntax
## Code Block Highlights (from /check-highlights)
## Inline Code Formatting (from /check-inline-code)
## Writing Style
## Summary
```

For deep lesson reviews, add sections for: Structure, Formatting, Content Accuracy, Teaching Effectiveness, Cross-References, Examples & Code, Writing Style, Overall Assessment, Positive Aspects.

**Important**: Always provide actionable next steps for any issues found. If no issues found, give clear confirmation.
