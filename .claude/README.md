# Claude Configuration for Nextflow Training Materials

This directory contains Claude AI assistant configuration to help developers create and maintain Nextflow training materials.

## Files Overview

### `.clinerules`

Core rules and conventions for this repository. Claude will automatically follow these guidelines when helping with:

- Nextflow script development
- Markdown documentation formatting
- Repository structure and organization
- Content style and best practices

### `commands/` Directory

Custom slash commands for common development tasks.

## Available Commands

### Content Creation

- **`/new-module`** - Scaffold a complete new training module with all standard files (index, orientation, lessons, survey, next_steps)
- **`/new-lesson`** - Create a new lesson page within an existing module following template structure
- **`/add-exercise`** - Add an exercise with solution to an existing lesson

### Quality Assurance

- **`/validate`** - Run all validation checks (heading numbering, TODO comments, Nextflow syntax, orphaned files)
- **`/review-lesson`** - Comprehensive review of a lesson for structure, formatting, content, and teaching effectiveness
- **`/test-example`** - Test a Nextflow script and verify it matches documentation

### Development Workflow

- **`/preview`** - Start local MkDocs preview server (Docker or Python)
- **`/update-nav`** - Update mkdocs.yml navigation after adding or reorganizing content
- **`/find-todos`** - Search for TODO and FIXME comments across the codebase

## Quick Start

### Creating New Content

```
User: /new-lesson
```

Claude will guide you through creating a properly formatted lesson page.

```
User: /new-module
```

Claude will scaffold a complete module with all necessary files.

### Quality Checks

Before committing, run validation:

```
User: /validate
```

To review a specific lesson:

```
User: /review-lesson
```

### Preview Changes

```
User: /preview
```

Claude will provide commands to start the local preview server.

## Repository Conventions

### Markdown Files

- Each sentence on new line (cleaner git diffs)
- Numbered headings with trailing periods (1., 1.1., 1.2.)
- Use admonitions for notes, tips, warnings, exercises
- Takeaway and What's next? sections at end of major sections

### Nextflow Scripts

- DSL2 syntax only
- UPPERCASE process names
- Shebang: `#!/usr/bin/env nextflow`
- Use params for configurable values
- Educational and simple examples

### Code Blocks

```groovy title="filename.nf" linenums="1" hl_lines="3"
#!/usr/bin/env nextflow
// Line 3 will be highlighted
params.example = 'value'
```

### Admonitions

```markdown
!!! note
Informational content

!!! tip
Helpful suggestions

!!! warning
Important warnings

??? exercise "Title"
Exercise content (collapsible)

??? solution
Solution content (collapsible)
```

## Development Workflow

1. **Create/Edit Content**

   - Use `/new-lesson` or `/new-module` for structure
   - Follow conventions in `.clinerules`
   - Add exercises with `/add-exercise`

2. **Preview Locally**

   - Use `/preview` to start server
   - View at http://127.0.0.1:8000/
   - Check formatting and navigation

3. **Validate**

   - Run `/validate` to check all files
   - Fix heading numbering if needed
   - Ensure no broken links

4. **Test Examples**

   - Use `/test-example` to verify Nextflow scripts work
   - Confirm outputs match documentation

5. **Update Navigation**

   - Use `/update-nav` to add new content to mkdocs.yml
   - Verify in preview

6. **Review**

   - Use `/review-lesson` for thorough quality check
   - Address any issues found

7. **Commit**
   - Write clear commit message
   - Push to fork
   - Open pull request

## Additional Resources

- **CONTRIBUTING.md** - Full contribution guidelines
- **Training Site** - https://training.nextflow.io
- **MkDocs Material** - https://squidfunk.github.io/mkdocs-material/
- **Nextflow Docs** - https://www.nextflow.io/docs/latest/

## Tips

- Use the Todo Tree VSCode extension to track TODO comments
- Install Prettier VSCode extension for auto-formatting
- Install Excalidraw VSCode extension for editing diagrams
- Test Nextflow examples before documenting them
- Each sentence on a new line makes git diffs cleaner
- Preview frequently to catch formatting issues early

## Getting Help

If you're unsure about conventions or best practices, ask Claude:

- "How should I format this code example?"
- "What admonition type should I use here?"
- "Is this Nextflow syntax correct for training?"
- "Review this section for consistency"

Claude has been configured with all the repository-specific knowledge needed to help you create high-quality training materials.
