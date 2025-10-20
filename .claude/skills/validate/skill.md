# Validate Training Materials

Run comprehensive validation checks on training materials to ensure quality and consistency.

## Tasks to Perform

1. **Check Heading Numbering**
   - Run: `uv run .github/check_headings.py docs/**/*.md`
   - Report any heading numbering issues found
   - If errors exist, ask if user wants to auto-fix with `--fix` flag

2. **Find TODO/FIXME Comments**
   - Search markdown files: `docs/**/*.md`
   - Search Nextflow scripts: `**/*.nf`
   - Search config files: `mkdocs.yml`, `.github/**/*.yml`
   - Categorize by priority (high, medium, low)
   - Report files with most TODOs

3. **Check Nextflow Script Conventions**
   - Find all .nf files
   - Verify they start with `#!/usr/bin/env nextflow`
   - Check for DSL2 syntax
   - Report any that don't follow conventions

4. **Find Orphaned Files**
   - Search for markdown files in docs/ not referenced in mkdocs.yml
   - Look for solution files without corresponding exercise documentation
   - Report any orphaned files

5. **Verify Admonition Syntax**
   - Search for common admonition formatting errors
   - Check for proper indentation (4 spaces)
   - Report any malformed admonitions

## Output Format

Provide a structured report:

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

## Notes

- Use Grep and Glob tools for efficient searching
- Run heading check script directly via Bash
- Provide actionable next steps for any issues found
- If no issues found, give clear confirmation
