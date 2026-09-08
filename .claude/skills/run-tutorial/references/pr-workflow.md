# PR Workflow for Tutorial Fixes

When the walkthrough identifies fixable issues, follow this workflow.

## 1. Independent Verification (MANDATORY)

Before a proposed fix reaches step 2, it must be confirmed by a subagent that did not diagnose it - see [Delegation Pattern](../SKILL.md#delegation-pattern) in the main skill. Spawn a fresh subagent, hand it only the lesson section text and the proposed diff, and ask it a single question: does this diff resolve the discrepancy against what the section actually teaches?

- **Pass**: carry the fix into step 2.
- **Fail or uncertain**: send the fix back to diagnosis with the independent reviewer's objection. Do not weaken the objection or re-verify it yourself - get a second independent pass on the revised fix.

This catches a fix that suppresses a symptom (e.g., deleting a failing assertion, or loosening a code snippet until it stops looking wrong) without actually making the tutorial teach the right thing.

## 2. Categorize Issues

**Auto-fixable** (apply programmatically):
- Extra/missing whitespace in code blocks
- Incorrect `hl_lines` values
- Heading numbering errors (use `uv run .github/check_headings.py --fix`)
- Minor formatting inconsistencies

**Requires manual review** (present options to user):
- Content accuracy issues
- Missing steps in documentation
- Incorrect command outputs
- Structural changes to lesson flow

## 3. Present Fixes to User

Before making any changes, list each fix with:

| Field | Description |
|-------|-------------|
| File | Path to the file |
| Line | Line number |
| Section | **Actual heading from document** (do not guess!) |
| Current | What it says now |
| Proposed | What it should say |
| Reason | Why this fix is needed |

Use `AskUserQuestion` with options:
- "Yes, create PR"
- "Let me review/modify first"
- "No, skip PR"

## 4. Create Branch

```bash
git checkout -b fix/<tutorial-name>-walkthrough-fixes
```

Use descriptive branch names based on the tutorial being fixed.

## 5. Apply and Commit

Stage only the files you're fixing:

```bash
git add docs/path/to/file.md
git commit -m "Fix issues in <tutorial-name> tutorial

- <bullet point for each fix>

Found during tutorial walkthrough testing.

Co-Authored-By: Claude <noreply@anthropic.com>"
```

## 6. Push and Create PR

```bash
git push -u origin fix/<tutorial-name>-walkthrough-fixes
gh pr create --title "Fix issues in <tutorial-name> tutorial" --body "$(cat <<'EOF'
## Summary

Fixes issues found during tutorial walkthrough testing.

## Changes

| Location | Before | After | Reason |
|----------|--------|-------|--------|
| Section X.Y (exact title) | old | new | why |

## Testing

- [x] Tutorial walkthrough completed successfully after fixes
- [x] All commands execute as documented
- [x] Output matches documentation

Generated with [Claude Code](https://claude.ai/code)
EOF
)"
```

## 7. Cleanup

After PR is created:

```bash
git checkout <original-branch>
```

Leave working directory files in place until PR is merged, in case revisions are needed.
