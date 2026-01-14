# Nextflow Training Materials Repository

This repository contains training materials for Nextflow, built with Material for MkDocs and published at https://training.nextflow.io.

## Quick Context

- **Purpose**: Educational materials teaching Nextflow workflow development
- **Build system**: MkDocs with Material theme and custom plugins
- **Languages**: Multiple (en, pt, es, fr, it, ko) - be aware when editing
- **Target audience**: Scientists and developers learning Nextflow

## Common Commands

### Preview locally

```bash
docker run --rm -it -p 8000:8000 -v ${PWD}:/docs ghcr.io/nextflow-io/training-mkdocs:latest
# View at http://0.0.0.0:8000/
```

### Validate heading numbering

```bash
uv run .github/check_headings.py docs/**/*.md
uv run .github/check_headings.py --fix docs/**/*.md  # auto-fix
```

### Format markdown

```bash
prettier --write docs/**/*.md
```

## Repository Structure

- `docs/` - All training content (markdown)
  - `hello_nextflow/` - Basic Nextflow introduction
  - `hello_nf-core/` - nf-core framework training
  - `nf4_science/` - Domain-specific training (genomics, RNAseq)
  - `side_quests/` - Advanced topics
- `nf-training/`, `nf4-science/`, etc. - Example Nextflow scripts for lessons
- `mkdocs.yml` - Site navigation and configuration
- `.github/check_headings.py` - Validates heading numbering

## Important Conventions

### Markdown Files

- **Heading numbering**: Must use `## 1.`, `### 1.1.` with trailing periods
- **One sentence per line**: For cleaner git diffs
- **Takeaway sections**: Each major section ends with "### Takeaway" and "### What's next?"
- **Admonitions**: Use `!!! note`, `!!! tip`, `!!! warning`, `??? exercise`, `??? solution`
- **Before/After comparisons**: Use tabbed blocks with highlighting to show code fixes (see debugging side quest):
  ```markdown
  === "After"
  `groovy title="fixed.nf" hl_lines="14" linenums="1"
      // corrected code with line 14 highlighted
      `
  === "Before"
  `groovy title="fixed.nf" hl_lines="14" linenums="1"
      // broken code with line 14 highlighted
      `
  ```
- **Code block line highlighting (`hl_lines`)**: This attribute specifies which lines to highlight, counted from the START of the code block (1-indexed). It is completely independent of `linenums`.

  **To set `hl_lines` correctly**: Before writing the `hl_lines` attribute, identify the lines you want highlighted, then count their position from line 1 of the snippet. Blank lines count. For example, given this snippet:

  ```
  #!/usr/bin/env nextflow    <- line 1
                             <- line 2 (blank)
  process FOO {              <- line 3
                             <- line 4 (blank)
      publishDir 'results'   <- line 5
                             <- line 6 (blank)
      input:                 <- line 7
          val x              <- line 8
  ```

  To highlight the `publishDir` line and the `val x` line, count their positions (5 and 8), then write: `hl_lines="5 8"`. The `linenums` attribute (which controls displayed line numbers) does not affect this counting.

  **Common errors**:

  - Skipping blank lines when counting (blank lines are lines too)
  - Assuming `hl_lines` numbers match `linenums` display numbers (they don't - always count from 1)
  - Highlighting structural keywords (`input:`, `output:`) instead of the meaningful content

### Nextflow Scripts

- **DSL2 only**: All examples use DSL2 syntax
- **UPPERCASE processes**: Process names like `PROCESS_NAME`
- **Shebang required**: `#!/usr/bin/env nextflow`
- **Educational focus**: Keep examples simple and well-commented
- **Dot notation for channels**: Use `.map{}`, `.view{}` instead of pipe operators

### Writing Style

Training materials should be clear, direct, and professional. Avoid patterns common in LLM-generated text:

**Tone:**

- Get to the point. Don't pad sentences with unnecessary elaboration
- Avoid casual interjections like "Let's", "Now let's", "Remember when we...?"
- Don't use exclamations for emphasis ("Much faster!", "So powerful!")
- Avoid phrases like "worth mentioning", "it's important to note", "don't worry"
- Use professional, neutral language rather than enthusiasm or reassurance

**Sentence structure:**

- Replace em-dash elaborations with periods or semicolons
  - Instead of: "This downloads the plugin - which happens automatically"
  - Write: "This downloads the plugin. It happens automatically."
- Use colons for list item explanations, not hyphens
  - Instead of: "**Option A** - does something useful"
  - Write: "**Option A**: does something useful"
- Keep explanations concise. One clear explanation is better than two different phrasings

**Content organization:**

- Don't explain the same concept two different ways. Pick one approach
- Avoid redundant sections (e.g., listing plugins when linking to a registry)
- Focus on practical, actionable information
- When teaching a concept, show one clear path rather than multiple alternatives

**Examples of what to avoid:**

```markdown
<!-- Avoid these patterns -->

Let's see how this works!
Remember when we created the config file earlier?
This is really powerful - it means you can do X, Y, and Z!
Don't worry - the tool handles this automatically.
It's worth mentioning that...

<!-- Better alternatives -->

Here's how this works.
The config file from section 2.1 defines...
This enables X, Y, and Z.
The tool handles this automatically.
Note that... (or just state the fact directly)
```

### Module Structure

Standard pattern for training modules:

```
module_name/
├── index.md (overview)
├── 00_orientation.md (introduction)
├── 01_topic.md (numbered lessons)
├── survey.md (feedback)
├── next_steps.md
└── solutions/ (working code)
```

## Gotchas

1. **Social cards are slow**: If preview is slow, disable with `CARDS=false`
2. **Multilingual builds**: Initial build takes minutes due to multiple languages
3. **Heading validation**: Auto-runs on commit via pre-commit hook
4. **Excalidraw diagrams**: Must use `.excalidraw.svg` extension for editability
5. **Module vs Lesson**: Module = entire course (e.g. "Hello Nextflow"), Lesson = single page (e.g. "01_hello_world.md")

## Testing Nextflow Examples

Always test examples before documenting:

```bash
cd [example-directory]
nextflow run example.nf
nextflow run example.nf -resume  # verify caching works
```

## Adding New Content

1. Create lesson/module using Claude commands: `/new-lesson`, `/new-module`
2. Add to `mkdocs.yml` nav section
3. Test locally with `/preview`
4. Validate by asking Claude to check quality (uses skills automatically)
5. Run heading check (happens automatically via hook)

## DO NOT

- Don't edit generated files in `.cache/` or `site/`
- Don't commit `work/` directories from Nextflow runs
- Don't skip heading validation - it will fail CI
- Don't forget to test Nextflow examples before documenting them
- Don't use hardcoded version numbers - they go stale

## Resources

- **Full contribution guide**: See CONTRIBUTING.md
- **Training site**: https://training.nextflow.io
- **Nextflow docs**: https://nextflow.io/docs/latest/
- **MkDocs Material**: https://squidfunk.github.io/mkdocs-material/
