---
description: Create a new lesson page within an existing training module
---

Create a new lesson page within an existing training module.

**Note**: A "module" is a complete training course (like "Hello Nextflow" or "Side Quests"), while a "lesson" is a single numbered page within that module (like "01_hello_world.md").

Follow these steps:

1. Ask the user:

   - Which module to add the lesson to (e.g., hello_nextflow, side_quests, nf4_science/genomics)
   - Lesson number (e.g., 01, 02, 03)
   - Lesson title (e.g., "Hello World", "Working with Channels")

2. Create the markdown file following the naming convention: `[number]_[title_with_underscores].md`

   - Example: `03_hello_workflow.md`

3. Use this lesson template structure:

````markdown
# Part [N]: [Title]

[Introduction paragraph explaining what this lesson covers]

---

## 1. [First Major Section]

[Content explaining the concept]

### 1.1. [Subsection]

[Detailed explanation with examples]

```bash
command example
```
````

```console title="Output"
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

```groovy title="example.nf" linenums="1" hl_lines="3 5"
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
   - Code blocks with proper formatting (linenums, titles, highlighting)
   - **For `hl_lines`**: Before writing this attribute, identify which lines you want highlighted, then count their position from line 1 of the snippet. Blank lines count. The `hl_lines` values are completely independent of `linenums`.
   - Before/After comparisons using tabbed blocks where relevant:
     ```markdown
     === "After"
         ```groovy title="example.nf" hl_lines="5" linenums="1"
         // corrected code
         ```
     === "Before"
         ```groovy title="example.nf" hl_lines="5" linenums="1"
         // broken code
         ```
     ```
   - Admonitions: `!!! note`, `!!! tip`, `!!! warning`, `??? exercise`
   - Console output examples
   - Clear explanations for beginners

5. Remind the user to:
   - Add the lesson to `mkdocs.yml` nav section in the correct module
   - Create corresponding Nextflow example scripts if needed
   - Add solution files for any exercises in `[module]/solutions/`
   - Test any Nextflow examples before committing
   - Run heading validation: `uv run .github/check_headings.py --fix docs/**/*.md`
   - Preview locally with `mkdocs serve` or Docker to verify formatting
````
