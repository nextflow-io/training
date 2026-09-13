# PR Workflow for Tutorial Fixes

When the walkthrough identifies fixable issues, follow this workflow.

## 1. Independent Verification (MANDATORY for fixes that change lesson meaning)

This step applies to any fix that changes what the lesson teaches or claims - a new or reworded explanation, a behavior claim, pedagogical framing, or where an explanation is placed. A fix that's mechanically checkable instead (a `grep` for a banned pattern returns zero matches, an `hl_lines` count is arithmetically correct, a hash-consistency check passes) can be confirmed directly by rerunning that exact check - skip ahead to step 2 once you've done so.

For a fix that does qualify, it must be confirmed by a subagent that did not diagnose it - see [Delegation Pattern](../SKILL.md#delegation-pattern) in the main skill. Spawn a fresh subagent and hand it:

- the lesson section text and the proposed diff
- the original discrepancy that triggered the fix
- the observed command output or file contents that exposed it
- instructions to rerun or reproduce the smallest relevant check, not just read the diff
- instructions to confirm the fix's full scope: re-derive the complete set of instances of the underlying problem (e.g. grep the file for the pattern being fixed) so a fix that only addresses the specifically-cited instance isn't mistaken for a complete one

Ask it to report a verdict, the evidence it checked, and any remaining uncertainty - a bare pass/fail isn't enough to catch a fix that only looks right. Note that the fixer's own sanity check (confirming its change builds/runs) is a separate, necessary step but does not substitute for this independent pass - it shares whatever blind spot produced the original diagnosis.

- **Pass**: carry the fix into step 2.
- **Fail or uncertain**: send the fix back to diagnosis with the independent reviewer's objection. Do not weaken the objection or re-verify it yourself - get a second independent pass on the revised fix. If a fix fails twice, stop looping: surface the discrepancy and both objections to the user instead of attempting a third revision unsupervised.

This catches a fix that suppresses a symptom (e.g., deleting a failing assertion, or loosening a code snippet until it stops looking wrong) without actually making the tutorial teach the right thing, and a fix that's correct as far as it goes but leaves the same problem uncorrected elsewhere in the file.

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
