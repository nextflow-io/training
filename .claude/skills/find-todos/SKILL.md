---
name: Find TODO and FIXME Comments
description: Search markdown files, Nextflow scripts, and config files for TODO/FIXME comments, categorize by priority, and provide actionable recommendations. Use when you need to identify pending work or track technical debt.
---

# Find TODO and FIXME Comments

Search the training materials codebase for TODO and FIXME comments to identify pending work.

Execute from repository root. See [../shared/repo-conventions.md](../shared/repo-conventions.md) for directory structure.

## Tasks to Perform

1. **Search Markdown Files**
   - Use Grep to find TODO and FIXME in `docs/**/*.md`
   - Capture file path, line number, and full comment
   - Note context (which module/lesson)

2. **Search Nextflow Scripts**
   - Find TODO/FIXME in `**/*.nf` files
   - Distinguish between training examples and solution code
   - Note if in core examples or side-quests

3. **Search Configuration Files**
   - Check `docs/en/mkdocs.yml`
   - Check `.github/**/*.yml` and `.github/**/*.py`
   - Check `CONTRIBUTING.md`

4. **Categorize Results**
   - **High Priority**: Marked as FIXME, TODO(urgent), or blocking
   - **Documentation**: TODOs in markdown files
   - **Code**: TODOs in .nf or .py files
   - **Configuration**: TODOs in config files

5. **Check Known Issues**
   - Read CONTRIBUTING.md for any documented TODOs
   - Note which are tracked vs untracked

## Output Format

```
# TODO/FIXME Report

## Summary
Total items, High priority count, by category (Documentation/Code/Configuration)

## High Priority Items
[List FIXME and urgent items with file:line and context]

## Documentation TODOs
[Group by file with line numbers]

## Code TODOs
[Group by file with line numbers]

## Recommendations
1. Immediate attention: [blocking items]
2. Next sprint: [significant items]
3. Low priority: [nice-to-haves]

## Files with Most TODOs
[Top 3-5 files]
```

## Notes

- Use Grep with pattern `TODO|FIXME` (case insensitive)
- Distinguish between legitimate TODOs and example comments in training materials
- Provide actionable priorities, not just a list
