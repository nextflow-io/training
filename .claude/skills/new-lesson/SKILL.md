---
name: new-lesson
description: Create a new numbered lesson page within an existing training module. Use when adding a new topic like '05_error_handling.mdx' to an existing course.
---

Create a new lesson page within an existing training module.

See [../shared/repo-conventions.md](../shared/repo-conventions.md) for directory structure and file conventions, and the `starlight-codeblocks` skill for code block features.

**Note**: A "module" is a complete training course (like "Hello Nextflow" or "Side Quests"), while a "lesson" is a single numbered page within that module (like "01_hello_world.mdx").

Follow these steps:

1. Ask the user:

   - Which module to add the lesson to (e.g., hello_nextflow, side_quests, nf4_science/genomics)
   - Lesson number (e.g., 01, 02, 03)
   - Lesson title (e.g., "Hello World", "Working with Channels")

2. Create the page in `docs/src/content/docs/<module>/`, following the naming convention: `[number]_[title_with_underscores].mdx`

   - Example: `docs/src/content/docs/hello_nextflow/03_hello_workflow.mdx`

3. Use this lesson template structure:

````mdx
---
title: 'Part [N]: [Title]'
description: [One sentence summarising the lesson, used for search, share cards and llms.txt]
---

[Introduction paragraph explaining what this lesson covers]

---

## 1. [First Major Section]

[Content explaining the concept]

### 1.1. [Subsection]

[Detailed explanation with examples]

```bash
command example
# [!output]
expected output
```

### Takeaway

[Summary of what was learned in this section]

### What's next?

[Preview of next section]

---

## 2. [Second Major Section]

[Continue pattern...]

### 2.1. [Subsection]

[Content with code examples]

```nextflow {3,5} title="example.nf" showLineNumbers
#!/usr/bin/env nextflow

process EXAMPLE {
    // highlighted lines
}
```

### Takeaway

[Summary]

### What's next?

[Next steps or move to next lesson]
````

4. Include appropriate elements:

   - Code blocks with proper formatting (`title`, `showLineNumbers`, `startLineNumber=N`, highlights)
   - **For highlights (`{3,5-7}` on the fence)**: identify which lines you want highlighted, then count their position from line 1 of the snippet. Blank lines count. The numbers are completely independent of `startLineNumber`.
   - Before/After comparisons using a code walkthrough, Before first. Keep the blank lines after the opening tag and before the closing tag:

     ````mdx
     <CodeWalkthrough>

     ```nextflow title="example.nf" showLineNumbers step="Before"
     // broken code
     ```

     ```nextflow {5} title="example.nf" showLineNumbers step="After"
     // corrected code, with line 5 highlighted
     ```

     </CodeWalkthrough>
     ````

   - Command output inside the `bash` block, after `# [!output]`. Mark key output lines with `# [!code highlight]` at the end of the line; start error output with `# [!output stderr]`
   - Asides: `:::note`, `:::tip`, `:::caution`, `:::danger` (optional `[Title]`), closed by `:::`
   - Exercises and collapsible content: `<details class="exercise">` (see `/add-exercise`)
   - Links to other pages as relative `.mdx` links, e.g. `[Part 1](01_hello_world.mdx#11-examine-the-code)`
   - Clear explanations for beginners
   - No import lines: components are auto-imported

5. Add the page to `docs/src/sidebar.mjs`, in the module's `items` list in lesson order, as a slug without the extension (e.g. `"hello_nextflow/03_hello_workflow"`). The sidebar label comes from the page `title`.

6. Remind the user to:
   - Create corresponding Nextflow example scripts if needed
   - Add solution files for any exercises in `[module]/solutions/`
   - Test any Nextflow examples before committing
   - Run heading validation: `uv run .github/check_headings.py --fix docs/src/content/docs/<module>/*.mdx`
   - Preview locally with `/preview` to verify formatting
   - Translations are handled by the translation workflow, which maintainers run before a release
