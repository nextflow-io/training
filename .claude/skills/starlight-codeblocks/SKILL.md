---
name: starlight-codeblocks
description: Write, choose and configure code block features on Astro Starlight sites that use the starlight-codeblocks plugin. Covers focus, line states, squiggles, comment notation directives, annotations, footnotes, inline callouts, side annotations, scrollycoding, code mentions, word-level diff, code walkthrough, hidden lines, expandable blocks, visible whitespace, colourised brackets, colour swatches, file icons, inline code highlighting, code links, API auto-linking, line permalinks, code tabs, fill-in placeholders, smart shell copy, open in playground and run code. Use when you install or configure starlight-codeblocks, when you write or edit code blocks in .md, .mdx or .mdoc pages of a site that has it, when you add these features to the existing pages of a site, or when you must choose how to explain, highlight, compare, shorten, link or run code in Starlight docs.
license: MIT
---

# starlight-codeblocks

starlight-codeblocks is a Starlight plugin. It adds 27 features to the code blocks that Expressive Code renders. A feature starts only when a code block uses its attribute or its directive, so the plugin does not change existing pages. A few features apply to every block of a type: they are in [Features that start on their own](#features-that-start-on-their-own).

Full docs: https://ewels.github.io/starlight-codeblocks/. Each page has a Markdown version at the same address with `.md` at the end. https://ewels.github.io/starlight-codeblocks/llms-full.txt has every page in one file.

## Set-up

1. Install the package: `npm install starlight-codeblocks` (or the pnpm or Yarn equivalent).
2. Add `codeblocks()` to the Starlight `plugins` list in `astro.config.mjs`:

   ```js
   import starlight from '@astrojs/starlight';
   import { defineConfig } from 'astro/config';
   import codeblocks from 'starlight-codeblocks';

   export default defineConfig({
     integrations: [starlight({ title: 'My docs', plugins: [codeblocks()] })],
   });
   ```

3. If the site has an `ec.config.mjs` file with a `plugins` list, add `pluginCodeblocks()` from `starlight-codeblocks/expressive-code` to that list. Keep `codeblocks()` in `astro.config.mjs`, with all the options. The build stops with a message if you forget this step.
4. Put `codeblocks()` before any plugin or theme that sets the Starlight `expressiveCode` option.
5. On an Astro site without Starlight, add `codeblocks()` from `starlight-codeblocks/astro` to `integrations` in place of `expressiveCode()`. [references/configuration.md](references/configuration.md) has the set-up.
6. If the site uses `starlight-links-validator`, give it `exclude: linksValidatorExclude`, imported from `starlight-codeblocks`.

Requirements: Astro 7 or later, Starlight 0.42 or later, Node.js 22.12 or later. The plugin works with Astro's default Markdown processor (Sätteri) and with `unified()`. [references/configuration.md](references/configuration.md) has every option and the set-up for other plugins, Markdoc and sites without Starlight.

## Add features to an existing site

To add features to pages that exist, follow [references/retrofit.md](references/retrofit.md). It goes through the pages one at a time, and through each block in the context of its page.

## Choose a feature

Find the goal, then use the feature in the same row. Each reference file has the full syntax, the options and the limits of its features.

### Explain code

[references/explain-code.md](references/explain-code.md)

| Goal | Use | Why |
|---|---|---|
| One short sentence about one name on one line | Inline callout: `[!callout /name/]` on or above the line | Always visible, with an arrow at the name. Use one or two in a block. |
| A long note, or a note for only some readers | Annotation: `[!annotate]` at the end of the line | A numbered marker. The note opens in a popover only when the reader asks. |
| Every reader needs every note, and the block is short | Footnotes: `[!ref]` on or above the line | A numbered list under the block, the same on phones and desktops. |
| Every reader needs every note, and the block is long | Side annotations: `[!annotate]` plus `annotations="side"` | Notes in a column beside the code on wide screens, a list on narrow screens. |
| A walkthrough of one block in prose steps (MDX only) | Scrollycoding: `<Scrollycoding>` with `<Step focus="...">` | The block stays in view and focuses the lines of each step as it scrolls past. A block between steps changes the code. |
| A file that grows or changes over the steps of a tutorial (MDX only) | Code walkthrough: `<CodeWalkthrough>` with `step="label"` on each block | One block with numbered steps. Code that stays moves, new code fades in. |
| A paragraph that names lines of the block below it | Code mentions: `[!mention name]` plus `[text](#mention:name)` | The link highlights the tagged lines. The page stays plain Markdown. |

Use one note style in a block. Annotations and side annotations use the same directive, so a change between them is only the fence line. The fence line is the first line of the code block, with the language.

Use scrollycoding when prose explains the code step by step as the reader scrolls. Use code walkthrough when readers step through the versions of a file with buttons. For versions that are alternatives, not steps in an order, use code tabs.

### Draw attention

[references/draw-attention.md](references/draw-attention.md)

| Goal | Use | Why |
|---|---|---|
| The text is about a few lines, but readers need the rest for context | Focus: `focus={4-7}` or `[!code focus]` | Blurs the other lines. Hover or keyboard focus makes every line sharp. |
| A line is wrong, needs care, or needs a note, as in an editor | Line states: `[!code error] message`, `[!code warning]`, `[!code info]`, `[!code success]` | A tint, a bar and an optional message after the code. |
| One name or expression has a problem, as a linter or compiler reports it | Squiggles: `[!error /name/] message`, `[!warning /name/]`, `[!info /name/]` | A wavy line under the text, as in VS Code. The message shows on hover and keyboard focus. |
| A neutral highlight with no meaning | `{3}` on the fence line or `[!code highlight]` | The `mark` of Expressive Code. |
| A state that is not error, warning, info or success, such as "To do" | A custom state in `lineStates.states` | Its name becomes an attribute and a directive. |

### Make code easier to read

[references/readability.md](references/readability.md)

| Goal | Use | Why |
|---|---|---|
| Imports and set-up that readers need to run the code but not to understand it | Hidden lines: `hidden={1-3}` or `[!code hide]` | Named lines go behind a marker. The copy button still copies them. |
| A long file that readers scan | Expandable block: `expandable` or `expandable={N}` | Shows the first lines, with a button that shows the rest. |
| Indentation or tabs change the meaning (Make, Python, YAML) | Visible whitespace: `whitespace` or `whitespace="all"` | Shows spaces and tabs as faint glyphs. |
| Deep nesting on dense lines | Colourised brackets: `brackets` | Colours bracket pairs by depth. |
| CSS colours in code | Colour swatches: on by default, `swatches=false` | A swatch of each colour before its value. |
| Readers must see the file type of a titled block | File icons: on by default, `icon="<name>"`, `icon=false` | The icon of the file type before the title: vscode-icons by default, or another set, such as the `<FileTree>` icons, with `fileIcons.set`. |
| Code inside a sentence of prose | Inline code highlighting: `` `code{:lang}` `` | Syntax colours for inline code. |
| An edit to a few lines, before and after in one block | A `diff` block, `ins={}` and `del={}`, or `[!code ++]` and `[!code --]` | Word-level diff then highlights the changed words. It starts on its own. |

Hidden lines remove chosen lines. Expandable blocks cut a block at a line count. Focus keeps every line in view.

### Link code

[references/link-code.md](references/link-code.md)

| Goal | Use | Why |
|---|---|---|
| One piece of text in the code links to a URL that you choose | Code links: `[!link /text/ url] description` on or above the line | For one-off links and names that no adapter knows. The description shows in a card. |
| Every library name in Python or Nextflow code links to its reference | API auto-linking | Starts on its own for `py`, `python`, `pycon`, `nextflow` and `nf` blocks. |
| Readers must link to one line or a range of lines | Line permalinks: `id="name"` | Line numbers that are links, as `#name-L2` or `#name-L2-L4`. |

### Adapt to the reader

[references/adapt-to-the-reader.md](references/adapt-to-the-reader.md)

| Goal | Use | Why |
|---|---|---|
| The same task in versions, and each version is one code block (npm, pnpm, Yarn) | Code tabs: `:::code-tabs{sync="pm"}` with `label="..."` on each block | Tabs in the title bar, or a menu with `control="menu"`. Works in `.md` and `.mdx`. |
| Each version needs prose, a list or more than one code block | Starlight `<Tabs>` with `syncKey` (MDX only) | Tabs can hold any content. |
| A value that each reader must change, such as `YOUR_TOKEN` | Fill-in placeholders: `placeholder="YOUR_TOKEN"` | Fields in the code. The value fills every block and the copied code. |

A code tabs block `sync` key and a `<Tabs>` `syncKey` do not switch together. Use one of the two for each kind of choice on a site.

### Copy and run

[references/copy-and-run.md](references/copy-and-run.md)

| Goal | Use | Why |
|---|---|---|
| A terminal session with prompts and output | Smart shell copy: start commands with `$ ` | Starts on its own. **Copy commands** copies the commands only. |
| Readers try the code in an online playground | Open in playground: `playground="typescript"` or `playground="rust"` | A title bar button that opens the playground with the code. |
| Readers run Python, JavaScript or TypeScript in the page | Run code: `runnable` | A **Run code** button and an output panel. Other languages need a runtime, or scripted output after an `[!output]` comment. |

## Syntax rules

Docs: https://ewels.github.io/starlight-codeblocks/comment-notation/, with every attribute in https://ewels.github.io/starlight-codeblocks/reference/attributes/ and every directive in https://ewels.github.io/starlight-codeblocks/reference/directives/

Attributes go on the fence line, after the language:

- A flag is a name on its own: `brackets`, `runnable`.
- A value is in quotes: `playground="rust"`, `id="cfg"`.
- A range or a number is in braces: `focus={2-4}`, `hidden={1, 4-6}`, `expandable={8}`.
- The attributes of Expressive Code, such as `title`, `frame`, `{2}`, `ins`, `del`, `"text"`, `wrap` and `showLineNumbers`, work as before.

Directives go in a comment, in the comment syntax of the language of the block. For example, JavaScript uses `//`, Python and shell use `#`, and HTML and Markdown use `<!-- -->`.

| Form | Meaning | Example |
|---|---|---|
| `[!code <name>]` | Applies to the line it is on. | `port: 3000, // [!code focus]` |
| `[!code <name>:N]` | Applies to its line and the next N-1 lines. | `// [!code ++:3]` |
| `[!code <name> <name>]` | Applies each name to the line. Text after it goes to the first name that takes text. | `// [!code focus ++] Loads the plugin` |
| `[!<name>] <text>` | Takes the text after it, to the end of the comment. | `// [!annotate] Runs once` |
| `[!<name> /<text>/]` | Points at the first whole-word match of the literal text on its target line, else the first match. | `// [!callout /port/] From the environment` |
| `[\!code <name>]` | Renders as the literal text `[!code <name>]`. | `// [\!code focus]` |

Rules that apply to every directive:

- A directive at the end of a line of code applies to that line. On a line with no code, it applies to the line below, and that line goes.
- The plugin removes directives from the rendered code and from the copied text. A comment that holds only directives goes completely.
- Ranges on the fence line count the lines that readers see. Lines that hold only directives do not count.
- `/text/` is literal text, not a regular expression.
- Features combine. One comment can hold several end-of-line directives, as in `// [!code ++] New [!annotate] Why`. Fence line attributes add to them.
- Inline code, links and bold in the text of a directive render as HTML. Other Markdown stays as text.
- One comment on each line can hold directives.
- JSON has no comments. Use `jsonc` for a block with directives.

## Features that start on their own

These features need no attribute. Know them, because they can change a block that you did not mean to change:

- Word-level diff applies to every removed line that has an added line below it. Turn it off for a block with `wordDiff=false`.
- API auto-linking links names in every `py`, `python`, `pycon`, `nextflow` and `nf` block. Turn it off for a block with `apiLinks=false`.
- Smart shell copy applies to every block with a terminal frame and a line that starts with a prompt (`$ ` or `> ` by default). It also applies to every `pycon` block with a line that starts with `>>> `. It applies to a `python` or `py` block if the first line starts with `>>> `.
- Colourised brackets apply to every block in the languages in `brackets.languages`, if the site sets that option.
- Colour swatches apply to every block with a CSS colour in it. Turn them off for a block with `swatches=false`.
- File icons apply to every block with a title in an editor frame. Turn them off for a block with `icon=false`, or for a file type or a language with `false` in `fileIcons.files` or `fileIcons.languages`.
- Expandable blocks apply to every block longer than `expandable.auto` lines, if the site sets that option. Turn it off for a block with `expandable=false`.
- Inline code highlighting applies to all inline code, if the site sets `inlineHighlighting.defaultLanguage`. Keep one piece plain with `{:txt}`.

If a page shows one feature, keep the others out of its examples with these attributes.

## Gotchas

- `<CodeWalkthrough>` and `<Scrollycoding>` work in MDX files only. Import them with `import { CodeWalkthrough, Scrollycoding, Step } from 'starlight-codeblocks/components';`. Put an empty line after the opening tag and before the closing tag.
- A `:::code-tabs` directive can contain only code blocks. A paragraph inside it fails the build.
- For inline code highlighting, put the suffix inside the backticks: `` `res.ok{:js}` ``. MDX reads `` `res.ok`{:js} `` as a JavaScript expression and the build fails.
- A directive with an unknown name, or of a feature that is off, stays in the code and logs a build warning. So does a directive outside a comment that the block's language reads. Read the build warnings after every change: each gives the file, the block and the line.
- Each `id` for line permalinks must be unique on the page, and must not match a heading id.
- Hidden lines still run in the copied code, in playgrounds and with `runnable`. Code for a playground or the **Run code** button must be complete.
- Fill-in placeholder values stay in the browser's `localStorage` by default. For secrets such as API tokens, set `placeholders: { storage: 'session' }`.
- In Markdoc (`.mdoc`) files, put attributes in a `meta` attribute: `` ```js {% meta="focus={2}" %} ``. Code tabs and inline code highlighting do not work there.
- The `<Code>` component of Starlight gets every feature inside a block. Give the attributes in its `meta` prop.
- In blocks from VitePress, put a space between the language and a range: `js {1,3}`. Use `"word"` on the fence line in place of `[!code word:...]`. Use `:::code-tabs` with `label="..."` in place of `::: code-group`.

## Check the result

1. Build the site and read every warning from `starlight-codeblocks`.
2. Open the page. Check that the copy button copies the code without directives.
3. For an interactive feature, check it with the keyboard as well as the mouse cursor.
