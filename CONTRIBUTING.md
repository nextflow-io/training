# Contributing

Table of contents:

- [Contributing](#contributing)
  - [Contribution model](#contribution-model)
  - [Installation](#installation)
  - [Previewing a release](#previewing-a-release)
  - [Translations](#translations)
  - [Announcement banner](#announcement-banner)
  - [Figures \& diagrams](#figures--diagrams)
  - [Content style and formatting](#content-style-and-formatting)
    - [Pages](#pages)
    - [Formatting / linting](#formatting--linting)
    - [Nextflow linting](#nextflow-linting)
    - [Headings CI tests](#headings-ci-tests)
    - [Code block checks](#code-block-checks)
    - [Asides and collapsibles](#asides-and-collapsibles)
    - [Code blocks](#code-blocks)
    - [Index page template](#index-page-template)
  - [TODO / FIXME](#todo--fixme)

## Contribution model

The typical workflow for contributing to training content is as follows:

1. Make a _fork_ of the GitHub repository to your own account
2. Develop locally (see below) and make your changes
3. Commit and push to your forked repository
4. Open a pull-request against the main repo, which can be reviewed and merged

The training website is built with [Astro Starlight](https://starlight.astro.build/).
Code blocks use [Expressive Code](https://expressive-code.com/) with the [starlight-codeblocks](https://ewels.github.io/starlight-codeblocks/) plugin.

Every pull request is built in CI, and a comment on the PR links to a preview of the built site.
Merges to `master` publish the `/development/` version of the site; releases publish a numbered version and `/latest/`.

To make changes, you should run the website locally so that you can preview changes.

## Installation

The site lives in `docs/` and needs [Node.js](https://nodejs.org/) 22 or later.

```bash
cd docs
npm ci
npm run dev
```

Open <http://localhost:4321/> in your web browser.
Pages reload when you save changes in your editor.
To work faster on English pages, run `TRANSLATIONS=false npm run dev`: it leaves out the 10 translations.

To check a change the way CI does, run a full build.
It builds every language and runs the links validator, so it takes about 80 seconds and 6 GB of RAM:

```bash
cd docs
BASE=/latest/ npm run build
BASE=/latest/ npm run preview
```

The preview serves the built site at <http://localhost:4321/latest/>.

## Previewing a release

The `preview_release.py` script serves the training site locally at `https://training.nextflow.io/` with the current branch appearing as a specified version release.
This is useful for recording videos or previewing how a release will look before it's published.

### How it works

1. Fetches the MkDocs-era released versions and `versions.json` from the `gh-pages` branch
2. Builds the current branch twice, as `/<version>/` and `/latest/`, as a release does (Astro writes the base path into every URL)
3. Updates `versions.json` to show your version as "latest"
4. Generates trusted TLS certificates using mkcert
5. Adds a temporary entry to `/etc/hosts` to redirect the domain locally
6. Serves the site using Caddy on port 443
7. Cleans up the hosts entry when you stop the server (Ctrl+C)

### First-time setup

Install the mkcert root CA (one-time, as your regular user):

```bash
mkcert -install
```

This adds a local certificate authority to your system keychain so browsers trust the generated certificates.
Restart your browser after running this command.

The script also needs Homebrew (it installs mkcert and Caddy), Node.js and an `upstream` remote pointing at `https://github.com/nextflow-io/training.git`.

### Usage

Run the script with [uv](https://docs.astral.sh/uv/), which handles its dependencies, and without `sudo`: it asks for your password only to edit `/etc/hosts`.

```bash
uv run ./preview_release.py --version 4.0

# Check current status
uv run ./preview_release.py status

# Clean up work directory on exit (default keeps it for faster restarts)
uv run ./preview_release.py --version 4.0 --clean
```

### Caching

The script caches the `gh-pages` content, both builds and the TLS certificates in `.preview-release/`.
Each build takes about a minute and is redone only when files in `docs/` change.
Use `--clean` to delete the work directory on exit, or remove `.preview-release/` to force a fresh build.

## Translations

Training materials are available in multiple languages.

> [!WARNING]
> All translations are **generated and maintained by AI**.
> Do not submit manual translations - they will be overwritten by automated updates.

### For Readers

If you find a translation error:

1. Check if it's a glossary issue (wrong term used consistently)
2. Open an issue describing: language, file, current text, expected text
3. A maintainer will update the translation prompt and regenerate

### For Contributors

To improve translations, update the LLM prompts - not the translated files directly:

- **Language-specific rules**: Edit `_scripts/translate/languages/<lang>/llm-prompt.md`
- **General formatting rules**: Edit `_scripts/general-llm-prompt.md`

Then use the GitHub Actions workflow to regenerate:

1. Go to **Actions** → **Translate** → **Run workflow**
2. Select the language and mode
3. The workflow will regenerate the translation with your prompt changes

Translation runs are triggered manually by maintainers, typically ahead of a new release.
Resulting translation PRs should be merged promptly to avoid them stacking up or becoming outdated.

### Why AI-Only?

Manual translations don't scale and aren't maintainable:

- They get overwritten when English content changes
- There's no way to track intentional vs accidental differences
- The same errors reappear in new content
- Quality varies between contributors

By using AI with carefully tuned prompts, we ensure:

- Consistent terminology across all pages
- Automatic updates when English changes
- Reproducible, auditable translations
- Easy fixes via prompt improvements

For detailed information, see [TRANSLATING.md](TRANSLATING.md).

## Announcement banner

The announcement banner at the top of every page is set in `docs/src/announcement.json`:

```json
{
  "active": true,
  "text": "Nextflow Training Week (2025-Q1):",
  "text_bold": "March 10-14, 2025",
  "url": "https://community.seqera.io/t/nextflow-training-week-2025-q1/1775/6",
  "url_text": "More info",
  "hide_after": "2025-03-15"
}
```

The banner shows while `active` is `true` and the date is before `hide_after`.
If you need more customisation, edit `docs/src/components/Banner.astro`.

## Figures & diagrams

Graphics should be drawn using [Excalidraw](https://excalidraw.com/).
Please use the [VSCode extension](https://marketplace.visualstudio.com/items?itemName=pomdtr.excalidraw-editor) and give files a `.excalidraw.svg` filename suffix.
Files will continue to be editable by others using this method.

Save diagrams in the `img/` directory next to the page and embed them as Markdown images:

```mdx
![The sayHello process](./img/hello_world.excalidraw.svg)
```

The site inverts `.excalidraw.svg` images in dark mode.
Mermaid code blocks do not render, so draw diagrams in Excalidraw.

## Content style and formatting

### Pages

Every page is an MDX file (`.mdx`) under `docs/src/content/docs/`.
English pages sit at the top level; translations sit in `docs/src/content/docs/<lang>/` and are generated.

Each page starts with frontmatter.
The `title` is the page heading, so the body has no `# H1`:

```mdx
---
title: "Part 1: Hello World"
description: One sentence that summarises the page, used for search, share cards and llms.txt.
---

Introduction paragraph.

## 1. First section
```

New pages must be added to the site navigation in `docs/src/sidebar.mjs`, as a slug without the extension (e.g. `"hello_nextflow/01_hello_world"`).

Starlight components (`Tabs`, `TabItem`, `Card`, `CardGrid`, `LinkCard`, `LinkButton`, `Code`), `CodeWalkthrough`, `Quiz` and `CourseInfo` are auto-imported, so pages need no import lines for them.

Link to other pages with relative links to the `.mdx` file, for example `[Part 1](../hello_nextflow/01_hello_world.mdx#11-examine-the-code)`.
The build checks every internal link and anchor, and fails on broken ones.

MDX is stricter than Markdown.
In prose, write `{/* comment */}` instead of HTML comments, escape braces as `\{` and `\}`, write a bare `<` as `&lt;`, and write links as `[https://example.com](https://example.com)` instead of `<https://example.com>`.
Code spans and code blocks need no escaping.

### Formatting / linting

Prettier ignores `.mdx` files: it parses them as MDX 1, and can turn a code block into inline code without failing the build.
The build is the check: a page that is not valid MDX fails `npm run build` with the file and line.

Prettier formats the other files, including the `.astro` components.
The checks run as git hooks with [prek](https://prek.j178.dev/), configured in `prek.toml`: run `prek install` once and they run on every commit (CI runs them too).
Editors use `docs/prettier.config.mjs` for files in `docs/`, once `npm ci` has run there.

For paragraphs containing multiple sentences, please start each sentence on a new line.
Doing that produces cleaner diffs for future code reviews.
The sentences will still be rendered as a continuous paragraph.
For an example, have a look at the code for this paragraph.

### Nextflow linting

This repository uses `nextflow lint` to check Nextflow scripts for errors and warnings.
The linter runs automatically in CI and posts results as a PR comment.

To run locally:

```bash
nextflow lint .
```

### Headings CI tests

This repository includes a Python tool to validate heading numbering consistency across training materials.

The `check_headings.py` script ensures:

- Sequential numbering at each level (1., 1.1., 1.2., etc.)
- Trailing periods after heading numbers
- Heading levels match numbering depth (## for 1., ### for 1.1.)

The easiest way to run it is [with `uv`](https://docs.astral.sh/uv/), which handles dependencies for you automatically:

```bash
# Check files for issues
uv run .github/check_headings.py docs/src/content/docs/hello_nextflow/*.mdx
```

```bash
# Auto-fix detected issues
uv run .github/check_headings.py --fix docs/src/content/docs/hello_nextflow/*.mdx
```

```bash
# Check every page (translations are skipped)
find docs/src/content/docs -name '*.mdx' -print0 | xargs -0 uv run .github/check_headings.py
```

The script runs in CI on every pull request and will cause a CI failure if any incorrect headings are found.

### Code block checks

`astro build` does not check code block attributes, so an off-by-one highlight range builds without error.
`check-fences.mjs` catches the mechanical problems:

```bash
node docs/scripts/check-fences.mjs docs/src/content/docs
```

Errors fail CI; warnings do not.
CI also fails when the build logs a code block warning, such as an unknown directive.

### Asides and collapsibles

Asides make a piece of content stand out:

```mdx
:::note
Text.
:::

:::tip[Custom title]
Text.
:::
```

The types are `note`, `tip`, `caution` and `danger`.
To nest asides, give the outer one more colons (`::::note` around a `:::tip`).

Collapsible blocks use HTML `<details>` with a class that sets the colour and icon.
Leave a blank line after `<summary>` and before `</details>`, so the body renders as Markdown:

```mdx
<details class="exercise">
<summary>Add a second greeting</summary>

Instructions.

<details class="solution">
<summary>Solution</summary>

Answer.

</details>

</details>
```

Add `open` (`<details class="tip" open>`) to show the content without a click.
The available classes are `exercise`, `solution`, `abstract`, `example`, `info`, `note`, `tip`, `warning`, `success`, `failure`, `question`, `full-code`, `optional`, `courses` and `catalog`.

### Code blocks

Code blocks use Expressive Code with the starlight-codeblocks plugin.
The full reference is the [starlight-codeblocks docs](https://ewels.github.io/starlight-codeblocks/); `.claude/skills/starlight-codeblocks/` has the same reference for AI agents.
The patterns this site uses:

- **Language**: `nextflow` for Nextflow scripts and config, `groovy` only for real Groovy (plugin classes, `build.gradle`). Add `nextflow="config"` to a config snippet whose title does not end in `.config`.
- **Attributes**: `title="main.nf"`, `showLineNumbers`, `startLineNumber=21` and line highlights such as `{1,3-5}`. Highlight numbers count from the first line of the block, blank lines included, whatever `startLineNumber` says.
- **Before/After**: a `<CodeWalkthrough>` with the Before block then the After block, marked `step="Before"` and `step="After"`. Leave a blank line after the opening tag and before the closing tag.
- **Command output**: put the output in the same `bash` block as the command, after a `# [!output]` line. Readers click **Run code** to see it, and the copy button copies only the command. Start error output with `# [!output stderr]`, and highlight an output line with `# [!code highlight]` at its end.
- **Annotations**: `// [!annotate] Explanation` at the end of a line adds a numbered popover note.
- **Code tabs**: `:::code-tabs` around two or more code blocks, each with `label="..."`. For tabs that also hold prose, use `<Tabs>` and `<TabItem label="...">`.
- **Inline code**: `` `channel.of(1, 2){:nextflow}` `` highlights inline code. The suffix goes inside the backticks.
- **Files from the repository**: import the file with `?raw` after the frontmatter and render it with `<Code>`. The import path is relative to the page:

  ```mdx
  import custom_config from "../../../../../nfcore-use/custom.config?raw";

  <Code
    code={custom_config}
    lang="nextflow"
    title="custom.config"
    meta="showLineNumbers"
  />
  ```

A Before/After example:

````mdx
<CodeWalkthrough>

```nextflow title="hello.nf" showLineNumbers step="Before"
process SAY_HELLO {
}
```

```nextflow {2} title="hello.nf" showLineNumbers step="After"
process SAY_HELLO {
    publishDir 'results'
}
```

</CodeWalkthrough>
````

A command with output:

````mdx
```bash
nextflow run hello.nf
# [!output]
 N E X T F L O W   ~  version 26.04.4

[65/7be2fa] sayHello | 1 of 1 ✔ # [!code highlight]
```
````

### Index page template

Course and module index pages render a block of collapsible information cards from structured frontmatter.
This keeps the index pages consistent and makes them easier to maintain.

To use it, add the fields below to the frontmatter and put `<CourseInfo frontmatter={frontmatter} />` where the cards should appear.

#### Basic structure

```mdx
---
title: Course Title
description: One-sentence summary of the course.
tableOfContents: false
page_type: index_page
index_type: course
additional_information:
  technical_requirements: true
  learning_objectives:
    - First objective
    - Second objective
  audience_prerequisites:
    - "**Audience:** Description of target audience"
    - "**Skills:** Required skills"
  videos_playlist: https://www.youtube.com/playlist?list=...
---

Summary paragraph describing the course.

<CourseInfo frontmatter={frontmatter} />

## Rest of page

Content after the component appears below the cards.
```

#### Frontmatter fields

| Field                    | Type                     | Required | Description                                                      |
| ------------------------ | ------------------------ | -------- | ---------------------------------------------------------------- |
| `page_type`              | `"index_page"`           | No       | Marks a course or module index page                              |
| `index_type`             | `"course"` or `"module"` | No       | Badge label displayed next to the "Additional information" title |
| `additional_information` | object                   | No       | Container for the collapsible cards                              |

#### Additional information fields

All fields within `additional_information` are optional:

| Field                    | Type             | Description                                                                                       |
| ------------------------ | ---------------- | ------------------------------------------------------------------------------------------------- |
| `technical_requirements` | `true` or string | If `true`, uses default text about GitHub/local installation. If a string, uses that custom text. |
| `learning_objectives`    | list of strings  | Rendered as a bulleted list. Must be a list, not `true`.                                          |
| `audience_prerequisites` | list of strings  | Rendered as a bulleted list. Supports markdown formatting. Must be a list, not `true`.            |
| `videos_playlist`        | URL string       | Uses default video description text plus a link to the playlist.                                  |
| `videos`                 | string           | Custom video description text (no link). Mutually exclusive with `videos_playlist`.               |

Text fields support Markdown.
Links in them resolve from the site root, so write `[Environment options](/envsetup/)`.

#### Default content

The default texts for `technical_requirements: true` and `videos_playlist` are UI strings in `docs/src/content/i18n/<lang>.json` (`courseInfo.defaultTechnicalRequirements` and `courseInfo.defaultVideos`), so they are translated with the rest of the interface.

#### Requirements

- The page must include `<CourseInfo frontmatter={frontmatter} />`
- `index_type` must be `course` or `module`
- `learning_objectives` and `audience_prerequisites` must be lists (not `true`)
- `videos_playlist` must be a URL
- `videos` and `videos_playlist` are mutually exclusive

The schema is in `docs/src/content.config.ts`.
The build fails with a descriptive error if these requirements are not met.

## TODO / FIXME

Please use `TODO` comments when you see something that needs coming back to.
I recommend the [Todo Tree VSCode extension](https://marketplace.visualstudio.com/items?itemName=Gruntfuggly.todo-tree) to find these comments easily.
In `.mdx` pages, write them as MDX comments: `{/* TODO: ... */}`.
