# Make code easier to read

Seven features make a block, or the code in the prose, easier to read. None of them changes the copied text.

- If lines are needed to run the code but not to understand it, hide them with hidden lines.
- If the whole file matters but readers scan it, cap it with an expandable block.
- If every line matters but a few matter more, use focus, in [draw-attention.md](draw-attention.md).
- If the exact whitespace changes the meaning, show it with visible whitespace.
- If brackets nest deeply on dense lines, colour them with colourised brackets.
- If code has CSS colours, colour swatches show them. They start on their own.
- If a block has a file name as its title, file icons show its file type. They start on their own.
- If code is inside a sentence, colour it with inline code highlighting.
- If the reader compares one before and one after of a few lines, use a diff in one block. Word-level diff starts on its own.

## Hidden lines

Docs: https://ewels.github.io/starlight-codeblocks/features/hidden-lines/

Use for imports, settings and boilerplate that readers need to run an example but not to understand it. Each run of hidden lines goes behind a dashed marker that the reader can open. Do not hide a line that the reader must change, such as an API key. Keep it visible, or use a fill-in placeholder.

| Syntax | Where |
|---|---|
| `hidden={range}` | Code block fence line, such as `hidden={1-3,6-7}` |
| `[!code hide]` | Comment at the end of the line |
| `[!code hide:N]` | Comment: its line and the next N-1 lines |

```ts
import { createServer } from 'node:http'; // [!code hide]
const port = Number(process.env.PORT ?? 3000);
createServer(handler).listen(port);
```

- The copy button always copies hidden lines, so the code still runs after a paste. Playgrounds and the **Run code** button get them too.
- The title bar gets a "Show N hidden lines" button that opens every run.
- A line permalink to a hidden line opens its marker.
- Hidden lines that show are dimmed to 75% opacity, which puts some syntax colours under 4.5:1 contrast. Set the `codeblocksHiddenLines.openOpacity` style setting to `1` if the site needs full contrast.
- Without JavaScript, hidden lines stay hidden.
- Option: `hiddenLines: false` ignores `hidden` and leaves `[!code hide]` in the code.
- Expressive Code's `collapse={range}` comes from `@expressive-code/plugin-collapsible-sections`. It keeps a summary line for each section and works without JavaScript. Hidden lines take only a thin marker, have a directive and a show-all button, and need JavaScript to open. The docs compare them: https://ewels.github.io/starlight-codeblocks/features/hidden-lines/#hidden-lines-or-collapsible-sections

## Expandable blocks

Docs: https://ewels.github.io/starlight-codeblocks/features/expandable-blocks/

Use when a full file is the clearest example, but a reader who scans the page does not need every line at once. The block shows its first lines with a fade, and a button that shows the rest. Do not use to hide chosen lines: use hidden lines.

| Syntax | Where |
|---|---|
| `expandable` | Code block fence line: caps the block at the site default, 12 lines |
| `expandable={N}` | Code block fence line: caps the block at `N` lines |
| `expandable=false` | Code block fence line: turns off the `expandable.auto` option for the block |

- A block collapses only if the collapse hides three lines or more.
- Lines that hidden lines remove do not count towards `N` or the total on the button.
- A block with Expressive Code's `collapse` attribute does not collapse, even with `expandable`.
- Find in the page still finds text in the collapsed lines, in browsers that support `hidden="until-found"`.
- Without JavaScript, and when the page prints, the block shows in full.
- Option: `expandable.lines` (default `12`) sets the count for the bare `expandable` attribute.
- Option: `expandable.auto` (default `false`) makes every block with more lines than this number expandable, at `expandable.lines`. It skips the variants of a code tabs block, runnable blocks, blocks with Expressive Code's `collapse` and the blocks in `<CodeWalkthrough>` and `<Scrollycoding>`.

## Visible whitespace

Docs: https://ewels.github.io/starlight-codeblocks/features/visible-whitespace/

Use when the exact whitespace character changes the meaning: a Makefile recipe needs a tab, Python needs consistent indentation, a patch has trailing spaces. Spaces show as middle dots and tabs as arrows. Do not use on blocks where whitespace does not matter, because the glyphs add noise.

| Syntax | Where |
|---|---|
| `whitespace` | Code block fence line: shows the leading whitespace of each line |
| `whitespace="all"` | Code block fence line: shows every space and tab, also between words and at the end of a line |

- Trailing whitespace shows only in fenced code blocks in Markdown and MDX. Expressive Code removes it from a block that the `<Code>` component renders.
- The glyphs are CSS. The copied text has the real spaces and tabs. Screen readers do not hear the glyphs.
- The plugin sets the Expressive Code `tabWidth` to `0` unless the site sets its own value, so tabs reach the block unchanged.
- Option: `whitespace: false` turns the feature off.

## Colourised brackets

Docs: https://ewels.github.io/starlight-codeblocks/features/colourised-brackets/

Use for lines with several levels of nested brackets. `()`, `[]` and `{}` cycle through three colours by depth. Hover over a bracket to outline it and its partner. Do not use on a short line with one level of nesting, where the colours add little.

| Syntax | Where |
|---|---|
| `brackets` | Code block fence line |
| `brackets=false` | Code block fence line. Turns off the colours for a language in `brackets.languages`. |

- Brackets in strings and comments, and brackets with no partner, keep their normal colour.
- The colours need no JavaScript. The outline on hover needs JavaScript.
- Option: `brackets.languages` (default `[]`) turns the colours on for every block in those languages, such as `['js', 'json']`. A name also covers the other names of its language (`js` covers `javascript`).
- Limits: the plugin finds strings and comments with the comment syntax of the language, not a full parser.

## Colour swatches

Docs: https://ewels.github.io/starlight-codeblocks/features/colour-swatches/

Starts on its own in every block. A small swatch in the colour goes before each CSS colour: hex, `rgb()`, `hsl()`, `hwb()`, `lab()`, `lch()`, `oklab()`, `oklch()`, `color()` and colour names. Hovering shows the colour text on a chip in its own colour. Clicking a colour copies it.

| Syntax | Where |
|---|---|
| `swatches=false` | Code block fence line. Turns off the swatches for one block. |
| `swatches` | Code block fence line. Turns on the swatches for a block outside `swatches.languages`. |
| `swatches.shape="square"`, `swatches.shape="rounded"`, `swatches.shape="circle"` | Code block fence line. Sets the shape for one block. |
| `swatches.match="all"`, `swatches.match="value"` | Code block fence line. `"all"` shows a swatch for every hex colour and colour function in one block. |

- In stylesheets, every colour in a declaration value gets a swatch. ID selectors, classes, variables and `url(#id)` do not.
- In other languages, a colour gets a swatch only in quotes or after `:`, `=` or `,`. Colour names need quotes.
- For a language that puts colours after other text, such as `| #4caf50 |`: set `swatches.delimiters` (`{ before: ['|'], after: ['|'] }`), or `swatches.match: 'all'`. `swatches.byLanguage` sets both per language, such as `{ mermaid: { match: 'all' } }`. In `'all'` mode, `#123` and `#add` still get no swatch.
- A colour function with `var()` or `calc()` in it has no swatch.
- The swatches copy no text. Copy needs JavaScript; the swatches do not.
- Options: `swatches.languages` (default `'all'`) and `swatches.formats` (default every format) choose what gets a swatch.
- Options: `swatches.shape` (`'square'`, `'rounded'` or `'circle'`, default `'rounded'`) and `swatches.size` (a CSS length, default `'0.8em'`).
- Options: `swatches.hover` and `swatches.copy` (default `true`), and `swatches.prose` (default `false`). With `swatches.prose`, colours in the text of a page, and inline code that is one colour, also get a swatch.
- A docs example about another feature with colours in it: add `swatches=false` so that it shows only its own feature.

## File icons

Docs: https://ewels.github.io/starlight-codeblocks/features/file-icons/

Starts on its own in every block with a title in an editor frame. An icon of the file type goes before the title. The plugin finds it from the title path, with the rules of the icon set, then from the language of the block. The Seti set uses the rules of Starlight's `<FileTree>`: the full file name, then the extension. Terminal frames and blocks with no title get no icon.

| Syntax | Where |
|---|---|
| `icon="<name>"` | Code block fence line. Sets the icon by name, such as `react`, `seti:vue` or a name in `fileIcons.icons`. |
| `icon=false`, `no-icon` | Code block fence line. Removes the icon. |
| `fileIcons.set="vscode-icons"`, `fileIcons.set="material"`, `fileIcons.set="catppuccin"`, `fileIcons.set="seti"` | Code block fence line. A coloured set, or the one-colour Seti icons. |
| `fileIcons.style="plain"`, `fileIcons.style="tile"` | Code block fence line. The icon alone, or on a square with rounded corners. |
| `fileIcons.colour="<colour>"` | Code block fence line. A CSS colour for the icon, or for the square of a tile. |

- The default icons are the coloured vscode-icons. Each coloured set has its own names. The Seti icons of `<FileTree>` have Starlight's names, with or without `seti:`. A Seti name works in every set.
- Option `fileIcons.set` picks the icons: the coloured `'vscode-icons'` (default), `'material'` or `'catppuccin'`, or `'seti'` (one colour). Material and Catppuccin need their package: `@iconify-json/material-icon-theme` or `@iconify-json/catppuccin`.
- Files in `.github/` and `.gitattributes` get the GitHub icon, from the whole path in the title.
- A coloured icon keeps its colours, on a neutral square in a tile. A Seti icon has the colour of the title, and its tile has the accent colour. A tile in a custom colour gets a black or white icon, whichever reads on it.
- Options: `fileIcons.style` (default `'plain'`). `fileIcons.languages` sets `icon`, `colour` and `style` for each language, such as `{ python: { colour: '#3776ab' } }`.
- Options: `fileIcons.languages.<lang>.icon` is an icon name or SVG markup, such as `siNextflow.svg` from `simple-icons`.
- Options: `false` in `fileIcons.files` or as `fileIcons.languages.<lang>.icon` gives no icon, such as `{ '.mmd': false }`. `icon="<name>"` on the fence line still sets one.
- Options: `fileIcons.icons` adds icons by name, as SVG markup or 24 by 24 path data. `fileIcons.files` maps a file name, an extension such as `.nf`, or a path pattern such as `docs/**/*.md` to an icon name. When two keys match, the later key wins.
- Code tabs show the icon of each file on its tab. A tab with only a `label` gets an icon only from `icon="..."`. The code tabs menu uses the same icon for each language, custom icons too.
- A docs example about another feature with a title: add `icon=false` if the icon distracts from the feature.

## Inline code highlighting

Docs: https://ewels.github.io/starlight-codeblocks/features/inline-code-highlighting/

Use for a short expression, command or declaration in a sentence. A `{:lang}` suffix gives inline code the syntax colours of the site's code blocks. The plugin removes the suffix.

| Syntax | Where |
|---|---|
| `` `<code>{:<lang>}` `` | Inside the backticks, at the end of the code. Use this form. |
| `` `<code>`{:<lang>} `` | Directly after the closing backtick. Fails in MDX. |

```md
Install the package with `pnpm add starlight-codeblocks{:sh}`, then restart the dev server.
```

- `<lang>` is any language name or alias that Expressive Code knows, such as `js`, `py`, `sh`, `css` or `html`.
- The suffix must follow the code with no space.
- MDX reads a `{:js}` after the closing backtick as a JavaScript expression, and the build fails. Put the suffix inside the backticks. `starlight-versions` and `starlight-md-txt` read `.md` files as MDX too.
- An unknown language shows as normal inline code, with a build warning.
- Inline code that contains a backtick keeps its suffix as text.
- Option: `inlineHighlighting: false` leaves the suffix in the prose as text.
- Option: `inlineHighlighting.defaultLanguage`, such as `'py'`, highlights inline code with no suffix in that language. A suffix wins over it. `{:txt}` keeps one piece of inline code plain.
- Limits: one line of code only. No attributes or directives. No token form, such as `{:.entity.name.function}`. Needs Starlight, and does not work in Markdoc files.

## Word-level diff

Docs: https://ewels.github.io/starlight-codeblocks/features/word-level-diff/

Word-level diff highlights the words that changed inside each pair of removed and added lines, so readers find a one-word edit in a long line. It has no attribute to turn it on. It applies to every removed line that has an added line below it, from any of these sources:

| Source | Where |
|---|---|
| A `diff` block, such as `diff lang="js"` | Code block fence line. `lang` keeps the syntax colours. |
| `ins={range}`, `del={range}` | Code block fence line, in a block of any language |
| `[!code ++]`, `[!code --]` | Comment at the end of the line |
| `wordDiff=false` | Code block fence line: turns the feature off for one block |

````md
```diff lang="js"
-const timeout = 5000;
+const timeout = options.timeout ?? 5000;
```
````

- Prefer `[!code ++]` and `[!code --]` when the block can change, because they move with their line.
- Changed words get a stronger tint and a bar under them. Removed words also have a line through them.
- A pair that is less than 40% similar keeps the whole-line tints only.
- Needs no JavaScript.
- Option: `wordDiff.minSimilarity` (default `0.4`, from 0 to 1). Raise it to show the highlight only for close edits. `wordDiff: false` turns it off for the site.
- Limits: pairing is by position, so a run of removed lines pairs one by one with the run of added lines below it. It compares visible text and does not parse the language. It skips long pairs, such as minified code.
- Add `wordDiff=false` to a block with `[!code ++]` and `[!code --]` that shows a different feature.
