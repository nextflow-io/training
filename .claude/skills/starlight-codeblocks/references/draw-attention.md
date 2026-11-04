# Draw attention

Focus, line states and squiggles mark code for every reader. Code mentions, in [explain-code.md](explain-code.md), mark lines only when the reader hovers over a link in the prose.

- If the text is about a few lines and the rest is context, use focus.
- If a line has a problem or a status, and the reader needs its name or a message, use a line state.
- If one name or expression on a line has the problem, use a squiggle.
- If the highlight has no meaning of its own, use `{3}` or `[!code highlight]`, the `mark` of Expressive Code.

## Focus

Docs: https://ewels.github.io/starlight-codeblocks/features/focus/

Use when a block shows a whole file but the text is about a few lines. Focus blurs and fades the other lines, which stay in the block for context. Do not use to point at one word: use an inline callout.

| Syntax | Where |
|---|---|
| `focus={range}` | Code block fence line, such as `focus={4-7}` or `focus={1, 4-6}` |
| `[!code focus]` | Comment at the end of the line |
| `[!code focus:N]` | Comment: its line and the next N-1 lines |

````md
```js title="src/config.js" focus={4-7}
import { defineConfig } from './lib.js';

export default defineConfig({
  cache: {
    dir: '.cache',
    maxAge: 3600,
  },
  retries: 2,
});
```
````

- Prefer the directive when the block can change, because it moves with its line.
- Hover over the block, or keyboard focus in it, makes every line sharp in 250 ms.
- Screen readers and the copy button get every line. Needs no JavaScript.
- Option: `focus.style` is `'blur'` (default) or `'dim'`. Use `'dim'` if the site does not need the blur, which is hard to read for some readers. `focus: false` leaves `[!code focus]` in the code.
- Style settings: `codeblocksFocus.blur`, `codeblocksFocus.opacity` and `codeblocksFocus.transitionDuration`.
- Limits: whole lines only. On a touch screen, the reader sees every line after a touch.

## Line states

Docs: https://ewels.github.io/starlight-codeblocks/features/line-states/

Use to tint a line as an error, a warning, a note or a success. An optional message shows after the code, as in a code editor. Do not use for a plain highlight: use `[!code highlight]`.

| Syntax | Where |
|---|---|
| `error={range}`, `warning={range}`, `info={range}`, `success={range}` | Code block fence line. Tints the lines, with no message. The first line of each group shows the name of the state. |
| `[!code error] message` | Comment. The message is optional. |
| `[!code warning] message` | Comment |
| `[!code info] message` | Comment. The label reads **Note**. |
| `[!code success] message` | Comment |
| `note={range}`, `[!code note]`, `warn={range}`, `[!code warn]` | Other names for `info` and `warning` |
| `[!code ++] message`, `[!code --] message`, `[!code highlight] message` | Comment. The message shows in a label in the colour of the marker, with no name before it. |
| `[!code <state>:N] message` | Comment: tints N lines, with the message on the first |
| `<state>={range}`, `[!code <state>]` | A custom state from the options |

```py
for name in sys.argv[1:]  # [!code error] SyntaxError: expected ':'
    print(name)
count = len(sys.argv)  # [!code warning] Includes the script name
```

- Use the attribute when the text around the block explains the problem. Use the directive to show a message.
- The label starts with the name of the state in bold: **Error**, **Warning**, **Note** or **Success**. Screen readers hear the name before the line.
- `lineStates: { prefix: false }`, or `lineStates.prefix=false` on one block, hides the name of the state on screen. A line with no message then shows only the tint, so colour alone marks the state. Screen readers still hear the name.
- Use `[!code ++] message` to say why a line was added, for example in a step of a tutorial. Without line states, the message stays in the code as a comment.
- A state on every line of a block tints the whole block, for example a one-line error message.
- The copy button leaves out the messages. Needs no JavaScript.
- Add a custom state in `lineStates.states`:

  ```js
  codeblocks({
    lineStates: {
      states: {
        todo: { label: 'To do', colour: { dark: '#c792ea', light: '#7c3aed' } },
      },
    },
  });
  ```

  Then use `todo={9}` or `[!code todo]`. A name uses lower-case letters, digits and hyphens, and cannot be the name of another attribute, such as `title` or `focus`. Each colour needs a contrast of 3:1 or more on the code background. A state with the name `error`, `warning`, `info` or `success` changes the built-in state. `note` and `warn` cannot be custom names.
- Option: `lineStates: false` leaves the directives in the code.
- Limits: a message is one line of text. `[!code <state>]` on a line with only a comment marks the line below, and the comment line goes.

## Squiggles

Docs: https://ewels.github.io/starlight-codeblocks/features/squiggles/

Use to draw a wavy error, warning or info line under one name or expression, as a code editor does for a linter problem. The message shows when the reader hovers over the text or focuses it. Do not use to mark a whole line with a message that is always visible: use a line state.

| Syntax | Where |
|---|---|
| `[!error /text/] message` | Comment on or above the line. The message is optional. |
| `[!warning /text/] message` | Comment on or above the line |
| `[!info /text/] message` | Comment on or above the line |
| `[!warn /text/]`, `[!note /text/]` | Other names for `warning` and `info` |
| `[!error] message` | Without `/text/`: the whole line, from the first character that is not a space |

```ts
import { join } from 'node:path'; // [!warning /join/] `join` is declared but its value is never read.
// [!error /readFle/] Cannot find name `readFle`. Did you mean `readFile`?
const config = await readFle('config.json');
```

- These directives have no `code` before the name. `[!code error]` is a line state, and `[!error]` is a squiggle.
- `/text/` is literal text, and the line goes under its first whole-word match, else its first match. Text that is not on the line logs a build warning, and the line renders without the squiggle.
- Several squiggles can go on one line. Put each directive on its own comment line above the line.
- The message starts with an icon and the name of the severity: **Error**, **Warning** or **Info**. Each squiggle is in the tab order, and Escape closes the message. Screen readers hear the severity and the message after the text.
- The copy button and a manual selection leave out the messages. In print, the messages are a list under the block. Without JavaScript, the wavy lines show but the messages do not open.
- Option: `squiggles: false` leaves the directives in the code.
- Style settings: `codeblocksSquiggles.error`, `codeblocksSquiggles.warning` and `codeblocksSquiggles.info`.
- Limits: one line per squiggle. The message opens only on hover or focus.
