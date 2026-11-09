# Explain code

Five features attach an explanation to lines of code: inline callouts, annotations, footnotes, side annotations and scrollycoding. Code walkthrough steps through versions of a block that changes. Code mentions link the prose around a block to its lines. All of them keep the explanation out of the copied code.

| Style | Reader sees the note | Note goes |
|---|---|---|
| Inline callout | At once | Above the line, with an arrow at one name |
| Annotation | After a click | In a popover that opens from a marker |
| Footnote | At once | In a list under the block |
| Side annotation | At once | In a column beside the code, on a wide screen |
| Scrollycoding | As each step scrolls past | In prose steps next to the block |

Rules to choose:

- If the note is about one name and fits in one short sentence, use a callout.
- If most readers can skip the note, use an annotation.
- If readers need every note, use footnotes for short blocks and side annotations for long blocks.
- If the explanation is several paragraphs that walk through one block, use scrollycoding.
- If the prose already explains the lines, use code mentions to connect the two.
- Use one style in a block, so that readers find every note in the same way.

## Inline callouts

Docs: https://ewels.github.io/starlight-codeblocks/features/inline-callouts/

Use for a note of one short sentence about one word or name on a line. The bubble is always visible, so keep to one or two callouts in a block. Do not use for a note about several lines: use annotations or footnotes.

| Syntax | Where |
|---|---|
| `[!callout /text/] note` | Comment at the end of the target line, or on its own line directly above it |
| `[!callout] note` | The same, with the arrow at the first character that is not a space |

```js
const controller = new AbortController();
// [!callout /signal/] Lets `controller.abort()` cancel the request.
const res = await fetch(url, { signal: controller.signal });
```

- The arrow points at the first whole-word match of `/text/` on the target line, else the first match.
- Two callout lines above one line give two bubbles, in source order.
- Between two lines with the same highlight (`{}`, `ins`, `del` or a line state), the bubble's row has the highlight too.
- The bubble is at most 60 characters wide, or 90% of the block. Longer notes wrap.
- Works without JavaScript. Screen readers read the bubble before its line.
- Option: `callouts: false` turns the feature off. The directive line then stays in the code.
- Limits: tabs after other characters on the line can move the arrow if the site sets a `tab-size` other than 2. With `wrap`, the arrow points at the first row of the line.

## Annotations

Docs: https://ewels.github.io/starlight-codeblocks/features/annotations/

Use when the note is long, or when most readers can skip it. The code stays as compact as a plain block. Do not use when readers need every note to follow the code: use footnotes or side annotations.

| Syntax | Where |
|---|---|
| `[!annotate] note` | Comment at the end of the line it explains |
| `startNoteNumber={N}` | Code block fence line: numbers the notes from `N`, to continue an earlier block. Works for side annotations too. |

```yaml
matrix:
  python: ["3.12", "3.13"]  # [!annotate] One job per version, run in parallel.
```

- Markers are numbered from 1 in each block, or from `startNoteNumber`. Hovering over a marker with a mouse shows its note until the mouse cursor leaves the marker and the note. Clicking a marker keeps its note open, and several can be open. Escape, or a click elsewhere, closes them.
- The note opens out of the marker, to the right of the line, when it fits there without covering code or another marker. Otherwise it opens under the marker, 340 px wide at most. Short notes on short lines fit beside the line.
- The note can hold inline code, links and bold text.
- Works without JavaScript, through the browser's `popover` attribute.
- When the page prints, each marker prints as its number, and the notes print as a numbered list under the block.
- Option: `annotations.style` is `'filled'` (default), for markers filled in the accent colour, or `'outline'`, for outlined markers in the theme's magenta. A block sets its own with `annotations.style="outline"` or `annotations.style="filled"` on the fence line. Side annotations use the same setting. `annotations: false` turns off annotations and side annotations.
- Limits: an annotation applies to one line. For a range, put it on the first line. Long notes make long source lines.

## Footnotes

Docs: https://ewels.github.io/starlight-codeblocks/features/footnotes/

Use when every reader needs every note, and the block is short enough that the list stays close to the lines. Each line gets a numbered badge, and the notes are a numbered list under the code.

| Syntax | Where |
|---|---|
| `[!ref] note` | Comment at the end of the line it explains, or on its own line directly above it |
| `footnotes="sticky"` | Code block fence line: keeps the list at the bottom of the window while the block is on screen |
| `footnotes="static"` | Code block fence line: turns the sticky list off for one block |
| `startNoteNumber={N}` | Code block fence line: numbers the footnotes from `N` |

```py
# [!ref] Creates the application object.
app = Flask(__name__)
```

- Hovering over a badge or a note highlights its line and its note. Clicking one keeps the highlight, and a second click clears it. Several can stay highlighted. Clicking a note does the same from the list.
- Without JavaScript, the badges and the numbers are plain links to each other.
- Option: `footnotes.sticky` (default `false`) makes every list sticky. `footnotes.style` is `'outline'` (default) or `'filled'`, the same two styles as annotations. A block sets its own with `footnotes.style="filled"` or `footnotes.style="outline"` on the fence line. `footnotes: false` turns the feature off.
- Limits: a footnote applies to one line. A sticky list covers the bottom of the block, so keep the notes short.

## Side annotations

Docs: https://ewels.github.io/starlight-codeblocks/features/side-annotations/

Use for a longer block that you explain step by step, where readers need every note. Write the notes with `[!annotate]`, and add `annotations="side"` to the fence line (the first line of the code block, with the language).

````md
```py title="report.py" annotations="side"
with path.open() as fh:  # [!annotate] Opens the file and closes it when the block ends.
    return list(csv.DictReader(fh))  # [!annotate] Each row becomes a dict keyed by the header line.
```
````

- The notes are a column beside the code when the container has space for the longest line. From its longest line, each block gets a width of 600, 800 or 1000 px. These hold about 42, 66 or 90 characters. The column sticks below the site header.
- In Starlight's default content column (45rem), only blocks with lines of about 42 characters or fewer get columns. Keep the lines of a block with side annotations short, or move a note off the longest line.
- For a page with longer lines, set `tableOfContents: false` in its frontmatter. Then a block that needs 800 or 1000 px spreads over the free space on each side of the content column. The text stays 45rem wide. This works only for a block directly on the page, not in tabs, asides, lists or components.
- In a narrower container, such as on a phone, the notes are a numbered list under the block.
- Hovering over a note, or focusing it, highlights its line. Hovering over a line highlights its note. Clicking a note or its number keeps the highlight, as for footnotes, and a second click clears it. Several can stay highlighted, and a click elsewhere clears them. Enter or Space does the same for a focused note.
- `codeSide="right"` on the fence line puts the code in the right column and the notes on the left.
- No options of its own. `annotations: false` turns it off.
- Limits: each note is next to its number, not next to its line. Lines of more than about 90 characters scroll inside the code column.

## Scrollycoding

Docs: https://ewels.github.io/starlight-codeblocks/features/scrollycoding/

Use for a walkthrough of a code block in prose steps. On a wide screen, the steps scroll in a column and the block stays in view. The step at the middle of the block sets the focus of the block. A code block between two steps is a new version of the code from the next step on. Use a code walkthrough instead when readers step through versions with buttons and need no prose between them.

| Syntax | Where |
|---|---|
| `<Scrollycoding>` ... `</Scrollycoding>` | Around one code block and its steps, in an MDX file |
| `<Step>` ... `</Step>` | After the code block, one for each step |
| `focus="<range>"` | On `<Step>`: a range without braces, such as `6-8` |
| `mark="<range>"` | On `<Step>`: marks lines as well |
| A code block between two steps | A new version of the code. The next steps count its lines. |
| `codeSide="left"` | On `<Scrollycoding>`: the code in the left column. The default is `right`. |

````mdx
import { Scrollycoding, Step } from 'starlight-codeblocks/components';

<Scrollycoding>

```js title="server.js"
import express from 'express';

const app = express();
app.listen(3000);
```

<Step focus="1">Import Express.</Step>
<Step focus="3">Create the app object that holds routes and middleware.</Step>
<Step focus="4">Start listening on port 3000.</Step>

</Scrollycoding>
````

- Put one code block first, then the steps. Leave an empty line after `<Scrollycoding>` and before `</Scrollycoding>`.
- A step can contain any Markdown. Text or components between the steps, but outside a `<Step>`, fail the build.
- A new version animates in as in a code walkthrough: code in both versions moves, new lines fade in. Under reduced motion, the code changes at once.
- The two columns need space for the longest line beside a text column of 12rem. From its longest line, the block gets a width of 600, 800 or 1000 px, which hold about 41, 65 or 89 characters. On a page with `tableOfContents: false`, a block directly on the page that needs 800 or 1000 px spreads past the content column.
- On a narrow screen, and without JavaScript, each step shows its own copy of the block, in the version of that step.
- Option: `scrollycoding: false` shows the narrow layout at every width, with no script. `walkthrough: false` changes versions without the animation.
- Limits: MDX only. A `focus` attribute on the fence line has no effect. A mark on the fence line shows in every step. Many steps with a long block make a long page on phones.

## Code walkthrough

Docs: https://ewels.github.io/starlight-codeblocks/features/code-walkthrough/

Use for a tutorial that shows the same file several times, with changes each time. The versions become one block with numbered steps in the title bar, and **Previous** and **Next** buttons under it. Code that stays moves to its new place, and new code fades in. Do not use for alternatives with no order: use code tabs.

| Syntax | Where |
|---|---|
| `<CodeWalkthrough>` ... `</CodeWalkthrough>` | Around two or more code blocks, in an MDX file |
| `step="<text>"` | Code block fence line of each step: the label of the step. Optional. |
| `start="last"` or `start={N}` | On `<CodeWalkthrough>`: the step that it opens on, counted from 1. Optional. |

````mdx
import { CodeWalkthrough } from 'starlight-codeblocks/components';

<CodeWalkthrough>

```js title="server.js" step="Create the app"
const app = express();

app.listen(3000);
```

```js title="server.js" step="Parse JSON bodies"
const app = express();
app.use(express.json());

app.listen(3000);
```

</CodeWalkthrough>
````

- Leave an empty line after `<CodeWalkthrough>` and before `</CodeWalkthrough>`. Each code block is one step, in order. Each keeps its other attributes, such as `title`.
- The walkthrough opens on its first step. For a Before and After pair, put Before first and add `start="last"`: the reader sees and copies After, and **Next** animates the change. `start="2"` in quotes works too. Any other value, or a number outside the steps, logs a build warning and opens the first step.
- The copy button copies the code of the current step.
- Arrow keys move between steps when a numbered step has focus. Under reduced motion, steps change without the animation.
- A line that is new in a step flashes the theme's green, then fades out in 1 second. There is no flash under reduced motion.
- Without JavaScript, and when the page prints, each step shows as a separate block with its label after the title.
- The page loads the animation library, about 3 kB, only on pages with `<CodeWalkthrough>`.
- Option: `walkthrough: { start: 'last' }` opens every walkthrough on its last step. `start` on a `<CodeWalkthrough>` overrides it. `walkthrough: false` shows each step as a separate block, with no script. The `codeblocksWalkthrough` style settings change the colours and `duration`.
- Limits: MDX only. During the animation, only the syntax colours show. Line states, line numbers and other decorations come back at the end.

## Code mentions

Docs: https://ewels.github.io/starlight-codeblocks/features/code-mentions/

Use when the prose explains the block line by line, such as "the base case stops the recursion". A link in the prose highlights the lines that it names. Do not use for lines that every reader must see at once: use focus.

| Syntax | Where |
|---|---|
| `[!mention <name>]` | Comment at the end of each line that belongs to the name |
| `[text](#mention:<name>)` | A link in the prose |

````md
The [base case](#mention:base) stops the recursion.

```py
def factorial(n):
    if n == 0:  # [!mention base]
        return 1  # [!mention base]
    return n * factorial(n - 1)
```
````

- A name is one word, such as `base` or `parse-args`. One line can have more than one tag.
- A link pairs with the next block in the same section that has the name. A section ends at the next heading. If no block follows, the link pairs with the nearest block before it.
- In a code tabs block, only the variant that shows counts. Tag the same name in each variant, and the link follows the reader's choice.
- A link with no matching block shows as plain text, with a build warning.
- `starlight-links-validator` reports `#mention:` links as broken. Give it `exclude: linksValidatorExclude`.
- Option: `mentions: false` turns the feature off.
- Limits: the build check does not see tags in a block from the `<Code>` component.
