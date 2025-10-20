---
description: Run all validation checks on training materials
---

Run validation checks on the training materials. Execute the following:

1. **Check heading numbering**:

   ```bash
   uv run .github/check_headings.py docs/**/*.md
   ```

   If errors found, ask user if they want to auto-fix:

   ```bash
   uv run .github/check_headings.py --fix docs/**/*.md
   ```

2. **Check for TODO/FIXME comments**:
   Search for TODO and FIXME comments in recently modified files and report them.

3. **Verify Nextflow scripts**:

   - Find all .nf files in the repository
   - Check that they start with `#!/usr/bin/env nextflow`
   - Check for DSL2 syntax
   - Suggest which scripts might need testing

4. **Check for orphaned files**:

   - Look for markdown files in docs/ that aren't referenced in mkdocs.yml
   - Look for solution files that don't have corresponding exercise documentation

5. **Verify admonition syntax**:

   - Check for common admonition formatting errors
   - Verify proper indentation

6. Provide a summary report with:
   - Number of issues found
   - Files that need attention
   - Suggested next steps
