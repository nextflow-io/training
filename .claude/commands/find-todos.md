---
description: Find all TODO and FIXME comments in the training materials
---

Search for TODO and FIXME comments in the training materials to identify work that needs to be done.

Search the following locations:
1. All markdown files in docs/
2. All Nextflow scripts (.nf files)
3. Configuration files (mkdocs.yml, .yml files in .github/)
4. Python scripts in .github/

For each TODO/FIXME found:
- Show the file path and line number
- Show the complete comment with context
- Categorize by priority if indicated (e.g., TODO(high), FIXME)

Organize results by:
1. **High Priority**: Issues marked as critical or blocking
2. **Documentation**: TODOs in markdown files
3. **Code**: TODOs in Nextflow or Python scripts
4. **Configuration**: TODOs in config files

Provide a summary with:
- Total count of TODOs/FIXMEs
- Files with the most items
- Suggestions for which to tackle first
- Any TODOs mentioned in CONTRIBUTING.md as known issues

Note: The repo maintainers recommend using the "Todo Tree" VSCode extension for ongoing TODO tracking.
