# Claude Configuration for Nextflow Training Materials

This directory contains Claude AI assistant configuration to help developers create and maintain Nextflow training materials.

## Files Overview

### `CLAUDE.md` (in repository root)

Core rules and conventions for this repository. Claude automatically follows these guidelines for:

- Nextflow script development (DSL2, process naming, params)
- Page formatting (MDX frontmatter, heading numbering, asides, code blocks)
- Repository structure and organization
- Content style and best practices
- Common commands and gotchas

### `skills/` Directory

Skills for content creation, preview, validation and testing. Run one with `/<name>`, or ask Claude and it picks the skill when relevant.

`skills/shared/repo-conventions.md` holds the directory mapping and file conventions that every skill uses.

## Skill Types: Design Rationale

### Authoring skills (Interactive, Template-Based)

Used for **content creation** tasks that:

- Follow a template or standard structure
- Need user questions/answers for customization
- Are step-by-step guided workflows

### Checking skills (Autonomous, Multi-Step)

Used for **validation and testing** tasks that:

- Require searching across multiple files
- Need to run commands and analyze outputs
- Make decisions based on what they find
- Can work independently without constant user input

## Authoring Skills

### `/new-module`

Creates a complete training module with standard structure under `docs/src/content/docs/`:

- `index.mdx` (overview with the index page frontmatter)
- `00_orientation.mdx` (prerequisites and intro)
- Numbered lesson files
- `survey.mdx` and `next_steps.mdx`
- A sidebar group in `docs/src/sidebar.mjs`

### `/new-lesson`

Creates a new lesson page within an existing module with:

- `title` and `description` frontmatter
- Proper heading numbering (1., 1.1., etc.)
- Takeaway and What's next? sections
- Placeholder code blocks with proper formatting
- A sidebar entry in `docs/src/sidebar.mjs`

### `/add-exercise`

Adds an exercise with solution to an existing lesson:

- Uses `<details class="exercise">` and `<details class="solution">` blocks
- Creates corresponding solution files if needed
- Proper formatting and structure

### `/preview` and `/stop-preview`

Start and stop the Astro dev server (`npm run dev` in `docs/`) at http://localhost:4321/. Pages reload on save.

## Checking Skills

### `validate`

Runs comprehensive validation and review checks:

- Heading numbering validation (runs check_headings.py)
- Code block meta audit (runs check-fences.mjs)
- Site build with the links validator
- TODO/FIXME comment search and categorization
- Nextflow script convention checking
- Orphaned page detection (pages missing from the sidebar)
- Deep lesson review (structure, formatting, content accuracy, teaching effectiveness)

Outputs structured report with actionable recommendations.

### `check-highlights` and `check-inline-code`

Check line highlight ranges (`{3-5}`, `# [!code highlight]`) and inline code (`` `x{:nextflow}` ``). `validate` runs both.

### `test-example`

Tests a Nextflow script and verifies documentation accuracy:

- Runs the script and captures output
- Tests resume functionality
- Tests with different parameters
- Compares actual behavior with documented behavior
- Reports discrepancies and suggests fixes

### `run-tutorial` and `docker-setup`

Walk through a lesson as a learner would, in the training Docker image, and compare every command's output with the page.

### `find-todos`

Searches for TODO/FIXME comments in pages, Nextflow scripts and config files, and groups them by priority.

### `starlight-codeblocks`

Reference for every code block feature on the site. Vendored unchanged from the `starlight-codeblocks` npm package (`docs/node_modules/starlight-codeblocks/skills/`); copy it again after upgrading the package.

## Quick Examples

**Create content:**

```
/new-lesson
/add-exercise
/new-module
```

**Preview locally:**

```
/preview
```

**Validate and review:**

- "Can you validate the training materials?"
- "Please review this lesson for quality"
- "Test the example script"
- "Find all the TODOs in the codebase"

## Key Repository Conventions

### Pages

- `.mdx` files with `title:` and `description:` frontmatter, no H1
- Each sentence on new line (cleaner git diffs)
- Numbered headings with trailing periods: `## 1. Section`, `### 1.1. Subsection`
- Takeaway and What's next? sections at end of major sections
- Asides `:::note`, `:::tip`, `:::caution`; collapsibles `<details class="exercise">`, `<details class="solution">`

### Nextflow Scripts

- DSL2 syntax only
- UPPERCASE process names
- Always include shebang: `#!/usr/bin/env nextflow` in scripts with an entry workflow
- Do not include shebang in module files, i.e. `**/modules/*.nf`
- Use params for configurable values: `params.input = 'default'`

### Code Blocks

Include titles, line numbers and highlighting:

````mdx
```nextflow {3} title="example.nf" showLineNumbers
#!/usr/bin/env nextflow

params.greeting = 'Hello'  // highlighted
```
````

## Development Workflow

1. **Create/Edit Content**

   - Use `/new-lesson` or `/new-module` for structure
   - Follow conventions in `CLAUDE.md`
   - Add exercises with `/add-exercise`

2. **Review Quality**

   - Ask Claude to validate and review: "Can you review this lesson?"
   - The validate skill runs automated checks and deep lesson reviews

3. **Test Examples**

   - Ask Claude to test scripts: "Can you test this Nextflow example?"
   - The test-example skill verifies scripts and compares with documentation

4. **Update Navigation**

   - Add new pages to `docs/src/sidebar.mjs`
   - Or just ask Claude to update it

5. **Preview Locally**

   - Run `cd docs && npm run dev`, or use `/preview` (see CONTRIBUTING.md)
   - View at http://localhost:4321/

6. **Commit**
   - Run the heading check (CI runs it on every PR)
   - Write clear commit message
   - Push to fork and open PR

## Additional Tools

Consider also using:

- **Vale**: Linting for docs (can be added as separate PR/hook)
- **Todo Tree**: VSCode extension to track TODO comments
- **Excalidraw**: VSCode extension for editing diagrams

## Additional Resources

- **CONTRIBUTING.md** - Full contribution guidelines
- **TRANSLATING.md** - Translation guide (all translations are AI-generated)
- **Training Site** - https://training.nextflow.io
- **Nextflow Docs** - https://www.nextflow.io/docs/latest/
- **Starlight Docs** - https://starlight.astro.build/
- **starlight-codeblocks Docs** - https://ewels.github.io/starlight-codeblocks/
