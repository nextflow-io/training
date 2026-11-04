# Nextflow Training Materials Repository

This repository contains training materials for Nextflow, built with Astro Starlight and published at https://training.nextflow.io.

## Quick Context

- **Purpose**: Educational materials teaching Nextflow workflow development
- **Build system**: Astro 7 + Starlight, in `docs/`. Markdown is processed by Sätteri (remark/rehype plugins do not run). Code blocks use Expressive Code with the starlight-codeblocks plugin
- **Languages**: English plus ca, de, es, fr, hi, it, ko, pl, pt, tr. Translations are AI-generated, see TRANSLATING.md
- **Target audience**: Scientists and developers learning Nextflow
- **History**: The site moved from MkDocs Material in 2026. MkDocs syntax (`!!!`, `???`, `===` tabs, `hl_lines`, `--8<--`) no longer works

## Common Commands

### Preview locally

Node 22 or later:

```bash
cd docs
npm ci
npm run dev        # http://localhost:4321/, reloads on save
```

In an agent session `astro dev` starts in the background. Manage it with `npx astro dev status`, `npx astro dev logs` and `npx astro dev stop` (run in `docs/`).

Production-like build (about 80 s and 6 GB RAM; runs the links validator):

```bash
cd docs
BASE=/latest/ npm run build
BASE=/latest/ npm run preview   # http://localhost:4321/latest/
```

### Validate heading numbering

```bash
uv run .github/check_headings.py docs/src/content/docs/hello_nextflow/*.mdx
uv run .github/check_headings.py --fix docs/src/content/docs/hello_nextflow/*.mdx
find docs/src/content/docs -name '*.mdx' -print0 | xargs -0 uv run .github/check_headings.py  # all pages; translations are skipped
```

### Other checks

```bash
node docs/scripts/check-fences.mjs docs/src/content/docs   # code block meta; errors fail CI
nextflow lint .   # also runs in CI and posts results as a PR comment
prek run --all-files   # Prettier, editorconfig and whitespace hooks from prek.toml; also runs in CI
```

Prettier ignores `*.mdx`: it parses MDX 1 and can turn a code block into inline code. The build is the check: a page that is not valid MDX fails `npm run build`.

### Preview a release

`preview_release.py` serves the current branch at `https://training.nextflow.io/` as if it were a published version (useful for recording videos before a release ships). Run it without `sudo`; it asks for a password only to edit `/etc/hosts`. Details in CONTRIBUTING.md:

```bash
uv run ./preview_release.py --version 4.0
```

## Repository Structure

- `docs/` - The Astro app (`astro.config.mjs`, `ec.config.mjs`, `package.json`)
  - `src/content/docs/` - English pages (`.mdx`), the source of truth
  - `src/content/docs/<lang>/` - Translations, same relative paths. Generated, never edit by hand
  - `src/content/i18n/<lang>.json` - UI strings
  - `src/sidebar.mjs` - Site navigation. Every new page goes here
  - `src/redirects.mjs` - Redirects for moved pages
  - `src/content.config.ts` - Frontmatter schema, validated at build time
  - `src/components/` - Custom components (`CourseInfo`, `Banner`, ...); `src/announcement.json` sets the banner
- `_scripts/translate/` - Translation package; `languages/<lang>/` holds `llm-prompt.md` and `glossary.yml`
- `hello-nextflow/`, `nextflow-run/`, `nfcore-use/`, `nfcore-build/`, `config-exec/`, `seqera-scale/`, `nf4-science/`, `side-quests/` - Runnable Nextflow scripts paired with the matching lesson tree. Some lessons import files from these dirs (`?raw` imports), so changes to a script usually need a matching page edit (and vice versa).
- `.github/check_headings.py` - Validates heading numbering
- `.github/workflows/docs-build.yml` - PR build: heading check, fence check, `astro build` with links validator. PRs get a preview URL comment

## Important Conventions

### Pages

- Every page is `.mdx`. Frontmatter holds `title:` (no H1 in the body) and `description:` (one sentence, used for search, share cards and llms.txt)
- `tableOfContents: false` hides the table of contents
- Components (`Tabs`, `TabItem`, `Card`, `CardGrid`, `LinkCard`, `LinkButton`, `Code`, `CodeWalkthrough`, `Quiz`, `CourseInfo`) are auto-imported. Don't add import lines, except `?raw` code includes
- **Heading numbering**: Must use `## 1.`, `### 1.1.` with trailing periods
- **One sentence per line**: For cleaner git diffs
- **Takeaway sections**: Each major section ends with "### Takeaway" and "### What's next?"
- **Asides**: `:::note`, `:::tip`, `:::caution`, `:::danger`, optional title in brackets (`:::tip[Faster runs]`), closed by `:::`. Nest by giving the outer aside more colons (`::::note` around `:::tip`)
- **Collapsibles** (exercises, solutions, file contents): `<details class="TYPE">` (add `open` to start expanded), `<summary>Title</summary>`, a blank line, the body, a blank line, `</details>`. Types: `exercise`, `solution`, `abstract`, `example`, `info`, `note`, `tip`, `warning`, `success`, `failure`, `question`, `full-code`, `optional`, `courses`, `catalog`

### Code blocks

- Use `nextflow` for Nextflow scripts and config (add `nextflow="config"` to config snippets whose title doesn't end in `.config`); keep `groovy` for real Groovy (plugin classes, `build.gradle`). Fence meta: `{1,3-5}` highlights, `title="main.nf"`, `showLineNumbers`, `startLineNumber=21`
- **Line highlighting (`{...}`)**: Counts lines from the START of the code block (1-indexed), independent of `startLineNumber`.

  **To set it correctly**: identify the lines you want highlighted, then count their position from line 1 of the snippet. Blank lines count. For example, given this snippet:

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

  To highlight the `publishDir` line and the `val x` line, write `{5,8}`. `startLineNumber` only changes the displayed numbers.

  **Common errors**:

  - Skipping blank lines when counting (blank lines are lines too)
  - Assuming highlight numbers match the displayed line numbers (they don't - always count from 1)
  - Highlighting structural keywords (`input:`, `output:`) instead of the meaningful content

- **Before/After comparisons**: `<CodeWalkthrough>` with two blocks, `step="Before"` then `step="After"`, each keeping its title, line numbers and highlights. Blank line after the opening tag and before the closing tag. If a version needs prose or several blocks, use `<Tabs>` instead

  ````mdx
  <CodeWalkthrough>

  ```nextflow title="hello.nf" showLineNumbers startLineNumber=17 step="Before"
  workflow {
      sayHello()
  }
  ```

  ```nextflow {3-4} title="hello.nf" showLineNumbers startLineNumber=17 step="After"
  workflow {
      sayHello()
      publish:
      first_output = sayHello.out
  }
  ```

  </CodeWalkthrough>
  ````

- **Command output**: in the command's `bash` block, after a `# [!output]` line (readers click **Run code**; copy takes only the command). Fence ranges count command lines only; highlight output with `# [!code highlight]` at the line end. `# [!output stderr]` starts error output. Escape a literal `[!x]` in output as `[\!x]`. Over 60 output lines: add `runnable.outputDelay=0`

  ````mdx
  ```bash
  nextflow run hello.nf
  # [!output]
   N E X T F L O W   ~  version 26.04.4

  [65/7be2fa] sayHello | 1 of 1 ✔ # [!code highlight]
  ```
  ````

- **Annotations**: `// [!annotate] Note text` at the end of a code line (`# [!annotate]` in shell)
- **Code-only tabs**: `:::code-tabs`, then two or more blocks each with `label="..."`, then `:::`. Only code blocks inside. Tabs with prose: `<Tabs>` / `<TabItem label="...">`
- **Inline code highlighting**: suffix inside the backticks, `` `channel.of('a'){:nextflow}` ``. Only where syntax colour helps (braces, closures, interpolation)
- **Including a file from the repo**: after the frontmatter, `import cfg from '../../../../../nfcore-use/custom.config?raw';` (relative to the page), then `<Code code={cfg} lang="nextflow" title="custom.config" meta="{3-6} showLineNumbers" />`

The vendored `starlight-codeblocks` skill documents every code block feature.

### Links, images and components

- **Links**: relative links to the `.mdx` file, `[Part 1](../hello_nextflow/01_hello_world.mdx#11-examine-the-code)`. The links validator fails the build on a missing page or anchor. Set a heading's ID with `{#id}` at the end of the heading
- **Images**: Markdown images relative to the page. Diagrams are Excalidraw files saved as `img/name.excalidraw.svg` and embedded as `![Alt text](./img/name.excalidraw.svg)`; dark mode inverts them. Mermaid blocks do not render
- **Quizzes**: `<Quiz id="quiz-1">`, a blank line, the question, a `- [ ]` / `- [x]` list, a blank line, `</Quiz>`. IDs are unique per page
- **Other components**: `<LinkButton href="02_hello_channels.mdx" variant="primary" icon="right-arrow">Continue</LinkButton>`; `<CardGrid>` with `<Card title="...">`; captions as `<figure>` + `<figcaption>`; keys as `<kbd>Ctrl</kbd>+<kbd>C</kbd>`
- **MDX escaping in prose**: comments are `{/* ... */}` (HTML comments fail); write `\{` `\}` for braces and `&lt;` for a bare `<`; use `[https://x.org](https://x.org)` not `<https://x.org>`. Code spans and code blocks need no escaping

### Nextflow Scripts

- **DSL2 only**: All examples use DSL2 syntax
- **UPPERCASE processes**: Process names like `PROCESS_NAME`
- **Shebang required**: `#!/usr/bin/env nextflow` in scripts with an entry workflow
- **No Shebang** in module files, i.e. `**/modules/*.nf`
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

Standard pattern for training modules, under `docs/src/content/docs/`:

```
module_name/
├── index.mdx (overview)
├── 00_orientation.mdx (introduction)
├── 01_topic.mdx (numbered lessons)
├── survey.mdx (feedback)
├── next_steps.mdx
└── img/ (images and .excalidraw.svg diagrams)
```

Solutions live with the scripts, e.g. `hello-nextflow/solutions/`.

### Index page template

Course and module `index.mdx` pages render their info cards from frontmatter with `<CourseInfo frontmatter={frontmatter} />`, placed where the cards go:

```yaml
---
title: Course Title
description: One-sentence summary.
tableOfContents: false
page_type: index_page
index_type: course # or "module"; renders as a badge
additional_information:
  technical_requirements: true
  learning_objectives: [...]
  audience_prerequisites: [...]
  videos_playlist: <url> # mutually exclusive with `videos`
---
```

The schema is in `docs/src/content.config.ts`; the build fails on invalid fields. Full spec in CONTRIBUTING.md.

## Gotchas

1. **MDX is strict**: a stray `{`, `<` or HTML comment in prose fails the build with a file and line number
2. **Pages must be in the sidebar**: a page missing from `docs/src/sidebar.mjs` builds but is unreachable from navigation
3. **Build warnings fail CI**: starlight-codeblocks warnings (unknown directive, bad range) fail `docs-build.yml`. Read them after every build
4. **Excalidraw diagrams**: Must use `.excalidraw.svg` extension for editability
5. **Module vs Lesson**: Module = entire course (e.g. "Hello Nextflow"), Lesson = single page (e.g. "01_hello_world.mdx")

## Testing Nextflow Examples

Always test examples before documenting:

```bash
cd [example-directory]
nextflow run example.nf
nextflow run example.nf -resume  # verify caching works
```

## Adding New Content

1. Create lesson/module using Claude skills: `/new-lesson`, `/new-module`
2. Add the page to `docs/src/sidebar.mjs` (a slug without extension, e.g. `"hello_nextflow/07_hello_x"`)
3. Preview locally with `/preview`
4. Validate with `/validate`
5. Run the heading check (CI runs it too)

## Project-specific Claude skills

This repo ships custom skills tailored to training authoring. Prefer them over generic edits when the task fits:

- **Authoring**: `/new-module`, `/new-lesson`, `/add-exercise`
- **Preview & test**: `/preview`, `/stop-preview`, `/test-example`, `/run-tutorial`, `/docker-setup`
- **Review**: `/validate`, `/check-highlights`, `/check-inline-code`, `/find-todos`
- **Reference**: `starlight-codeblocks` (vendored from the npm package): every code block feature, its syntax and limits. Components are auto-imported here, so ignore its import instructions

## DO NOT

- Don't edit generated files in `docs/dist/`, `docs/.astro/` or `docs/node_modules/`
- Don't edit translated pages under `docs/src/content/docs/<lang>/`; fix the English page or the translation prompt
- Don't commit `work/` directories from Nextflow runs
- Don't skip heading validation - it will fail CI
- Don't forget to test Nextflow examples before documenting them
- Don't use hardcoded version numbers - they go stale

## Resources

- **Full contribution guide**: See CONTRIBUTING.md
- **Translation guide**: See TRANSLATING.md (all translations are AI-generated)
- **Training site**: https://training.nextflow.io
- **Nextflow docs**: https://nextflow.io/docs/latest/
- **Starlight**: https://starlight.astro.build/
- **starlight-codeblocks**: https://ewels.github.io/starlight-codeblocks/
