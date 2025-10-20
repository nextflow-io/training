# Claude Configuration for Nextflow Training Materials

This directory contains Claude AI assistant configuration to help developers create and maintain Nextflow training materials.

## Files Overview

### `.clinerules`
Core rules and conventions for this repository. Claude automatically follows these guidelines for:
- Nextflow script development (DSL2, process naming, params)
- Markdown documentation formatting (heading numbering, admonitions, code blocks)
- Repository structure and organization
- Content style and best practices

### `commands/` Directory
Slash commands for interactive content creation tasks that use templates and require user input.

### `skills/` Directory
Skills for autonomous multi-step tasks like validation, testing, and analysis.

## Commands vs Skills: Design Rationale

### Commands (Interactive, Template-Based)
Commands are used for **content creation** tasks that:
- Follow a template or standard structure
- Need user questions/answers for customization
- Benefit from showing the full prompt/checklist to the user
- Are step-by-step guided workflows

### Skills (Autonomous, Multi-Step)
Skills are used for **validation and testing** tasks that:
- Require searching across multiple files
- Need to run commands and analyze outputs
- Make decisions based on what they find
- Can work independently without constant user input

## Available Commands

### `/new-module`
**Why a command**: Template-based scaffolding requiring user input about module name, location, and content.

Creates a complete training module with standard structure:
- index.md (overview)
- 00_orientation.md (prerequisites and intro)
- Numbered lesson files
- survey.md and next_steps.md
- solutions/ directory

### `/new-lesson`
**Why a command**: Guided workflow that asks questions and fills in a template structure.

Creates a new lesson page within an existing module with:
- Proper heading numbering (1., 1.1., etc.)
- Takeaway and What's next? sections
- Placeholder code blocks with proper formatting
- Standard structure

### `/add-exercise`
**Why a command**: Interactive insertion requiring decisions about placement, content, and difficulty level.

Adds an exercise with solution to an existing lesson:
- Uses `??? exercise` and `??? solution` admonitions
- Creates corresponding solution files if needed
- Proper formatting and structure

### `/review-lesson`
**Why a command**: The comprehensive checklist prompt is valuable for users to see what's being checked.

Comprehensive lesson review checking:
- Structure (heading numbering, Takeaway sections)
- Formatting (code blocks, admonitions, line-by-line sentences)
- Content accuracy (Nextflow syntax, command correctness)
- Teaching effectiveness (clear explanations, logical flow)

## Available Skills

### `validate`
**Why a skill**: Autonomous multi-step validation requiring searches, tool execution, and analysis.

Runs comprehensive validation checks:
- Heading numbering validation (runs check_headings.py)
- TODO/FIXME comment search and categorization
- Nextflow script convention checking
- Orphaned file detection
- Admonition syntax verification

Outputs structured report with actionable recommendations.

### `test-example`
**Why a skill**: Complex autonomous task requiring script execution, output verification, and comparison with docs.

Tests a Nextflow script and verifies documentation accuracy:
- Runs the script and captures output
- Tests resume functionality
- Tests with different parameters
- Compares actual behavior with documented behavior
- Reports discrepancies and suggests fixes

### `find-todos`
**Why a skill**: Search and analysis task that works autonomously across the codebase.

Searches for TODO/FIXME comments:
- Markdown files, Nextflow scripts, config files
- Categorizes by priority (high, medium, low)
- Groups by file and provides context
- Recommends prioritization

## Quick Examples

**Create content:**
```
/new-lesson
/add-exercise
/new-module
```

**Review quality:**
```
/review-lesson
```

**Run validation:**
```
@validate
@test-example
@find-todos
```

## Key Repository Conventions

### Markdown Files
- Each sentence on new line (cleaner git diffs)
- Numbered headings with trailing periods: `## 1. Section`, `### 1.1. Subsection`
- Takeaway and What's next? sections at end of major sections
- Use admonitions: `!!! note`, `!!! tip`, `!!! warning`, `??? exercise`, `??? solution`

### Nextflow Scripts
- DSL2 syntax only
- UPPERCASE process names
- Always include shebang: `#!/usr/bin/env nextflow`
- Use params for configurable values: `params.input = 'default'`

### Code Blocks
Include line numbers, titles, and highlighting:
```groovy title="example.nf" linenums="1" hl_lines="3"
#!/usr/bin/env nextflow

params.greeting = 'Hello'  // highlighted
```

## Hooks

### User Prompt Submit Hook
Automatically runs heading validation when you submit a prompt:
- Executes: `uv run .github/check_headings.py docs/**/*.md`
- Provides immediate feedback on heading numbering issues
- Non-blocking (allows operation to continue even if issues found)

## Development Workflow

1. **Create/Edit Content**
   - Use `/new-lesson` or `/new-module` for structure
   - Follow conventions in `.clinerules`
   - Add exercises with `/add-exercise`

2. **Review Quality**
   - Use `/review-lesson` for comprehensive checks
   - Use `@validate` skill for automated validation

3. **Test Examples**
   - Use `@test-example` skill to verify Nextflow scripts
   - Confirm outputs match documentation

4. **Update Navigation**
   - Manually edit mkdocs.yml to add new content
   - Or just ask Claude to update it

5. **Preview Locally**
   - Run: `mkdocs serve` or use Docker (see CONTRIBUTING.md)
   - View at http://127.0.0.1:8000/

6. **Commit**
   - Heading validation runs automatically
   - Write clear commit message
   - Push to fork and open PR

## Additional Tools

Consider also using:
- **Vale**: Linting for docs (can be added as separate PR/hook)
- **Prettier**: Markdown formatting (install VSCode extension)
- **Todo Tree**: VSCode extension to track TODO comments
- **Excalidraw**: VSCode extension for editing diagrams

## Additional Resources

- **CONTRIBUTING.md** - Full contribution guidelines
- **Training Site** - https://training.nextflow.io
- **Nextflow Docs** - https://www.nextflow.io/docs/latest/
