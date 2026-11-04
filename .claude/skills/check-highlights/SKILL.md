---
name: Check Code Block Highlights
description: Verify that line highlights in code blocks ({3,5-7} ranges on the fence line, [!code highlight] directives, meta on <Code>) mark the intended lines. Use when reviewing documentation with code examples, especially Before/After code walkthroughs, or when highlights seem incorrect.
---

# Check Code Block Highlights

Verify that line highlights in code blocks mark the intended lines.

See [../shared/repo-conventions.md](../shared/repo-conventions.md) for file conventions, and the `starlight-codeblocks` skill for the full directive syntax.

## Where highlights come from

- **Fence ranges**: `{3}`, `{2,5}`, `{4-6,14}` on the fence line, e.g. ` ```nextflow {4-6} title="hello.nf" `
- **`<Code>` includes**: the same range inside the `meta` prop, e.g. `<Code code={x} lang="csv" meta="{3-6} showLineNumbers" />`. Count against the included file
- **Directives**: `// [!code highlight]` (or `# [!code highlight]`) at the end of a line highlights that line. In a `bash` block with `# [!output]`, this is the only way to highlight output lines

## Critical Understanding

**Fence ranges are snippet-relative, NOT related to `startLineNumber`**:

- `startLineNumber=21` sets the _displayed_ starting line number
- `{3}` highlights the _3rd line of the snippet_ (which displays as line 23)
- These are completely independent attributes

Example:

````mdx
```nextflow {3} showLineNumbers startLineNumber=21
process FOO {      <- displayed as line 21 (snippet line 1)
    input:         <- displayed as line 22 (snippet line 2)
    val x          <- displayed as line 23 (snippet line 3) - HIGHLIGHTED
}
```
````

Two more counting rules:

- A line that holds only a directive comment (e.g. `// [!code highlight]` on its own line) is removed from the rendered block and does not count
- In a `bash` block with command output, ranges count the command lines only. Everything from `# [!output]` on is output, highlighted with `# [!code highlight]` at the end of the line

## How to Check

For each code block with a range or a highlight directive:

1. **Identify the code block** - Find fences with `{...}` after the language, `<Code ... meta="{...}">`, and lines ending in `[!code highlight]`. In `.mdx` pages under `docs/src/content/docs/` (skip the `<lang>/` translation folders)

2. **Count snippet lines** - Number each line of the code block starting from 1:

   ```
   1: #!/usr/bin/env nextflow
   2: (blank line)
   3: process FOO {
   4:     publishDir 'results'  <- if this should be highlighted, the range should include 4
   ...
   ```

3. **Verify intent** - For `<CodeWalkthrough>` Before/After pairs:

   - The `step="After"` block should highlight the NEW or CHANGED lines
   - The `step="Before"` block, if it has highlights, should mark the lines that WILL change (same content positions)
   - Both blocks should highlight corresponding lines showing the difference

4. **Check for common errors**:
   - Off-by-one errors (highlighting blank lines instead of code)
   - Confusing `startLineNumber` with the range (thinking `{21}` highlights displayed line 21)
   - Highlighting section headers (`input:`, `output:`) instead of the actual values
   - A fence range in a `bash` block that points into the output instead of a command line
   - Ranges past the end of the block. `node docs/scripts/check-fences.mjs docs/src/content/docs` reports these as errors

## Fixing Issues

When you find incorrect highlights:

1. Count the actual line positions in the snippet (1-indexed)
2. Identify which lines contain the meaningful changes
3. Update the range to reference those snippet line numbers (comma-separated, ranges with `-`: `{4-6,14}`)
4. Verify Before/After blocks highlight corresponding positions

## Output Format

For each file checked, report:

```
## [filename]

### Block at line [N]: [description]
- Current: {X,Y,Z}
- Lines X, Y, Z contain: [what's on those lines]
- Issue: [description of problem, if any]
- Suggested fix: {A,B,C} (highlighting [what those lines contain])

### Summary
- Blocks checked: N
- Issues found: M
- [List of fixes needed]
```

## Example Analysis

Given this code block:

````mdx
```nextflow {8,11,14} title="example.nf" showLineNumbers
#!/usr/bin/env nextflow       <- snippet line 1
                              <- snippet line 2 (blank)
process SAY_HELLO {           <- snippet line 3
                              <- snippet line 4 (blank)
    publishDir 'results'      <- snippet line 5
                              <- snippet line 6 (blank)
    input:                    <- snippet line 7
        val greeting          <- snippet line 8 - HIGHLIGHTED
                              <- snippet line 9 (blank)
    output:                   <- snippet line 10
        path "*.txt"          <- snippet line 11 - HIGHLIGHTED
                              <- snippet line 12 (blank)
    script:                   <- snippet line 13
        """echo hi"""         <- snippet line 14 - HIGHLIGHTED
}
```
````

Analysis: `{8,11,14}` correctly highlights the input value (line 8), output value (line 11), and script content (line 14).

## Notes

- Blank lines count! They're common sources of off-by-one errors
- In Before/After code walkthroughs, both blocks should highlight the same conceptual elements
- When in doubt, count every line including blanks and comments
