# Configuration

Docs: https://ewels.github.io/starlight-codeblocks/configuration/ and https://ewels.github.io/starlight-codeblocks/reference/options/

With no options, `codeblocks()` turns on every feature. Pages that use no interactive feature load no JavaScript from the plugin.

## Options

`codeblocks()` takes one object with a key for each feature. Set a key to `false` to turn the feature off, or to an object to change its settings. A feature with no settings takes only `false`: leave the key out to keep it on. An unknown key or a value of the wrong type stops the build with a message that names the option.

```js
codeblocks({
  focus: { style: 'dim' },
  expandable: { lines: 20 },
  runnable: false,
});
```

| Key | Settings | Feature |
|---|---|---|
| `focus` | `style`: `'blur'` or `'dim'`, default `'blur'` | Focus |
| `lineStates` | `states`: custom states by name, each `{ label, colour: { dark, light } }`. `prefix`: default `true`. `false` hides the state name on screen | Line states |
| `notation` | `comments`: comment syntax for each language, added to the built-in map | Comment notation |
| `callouts` | None | Inline callouts |
| `squiggles` | None | Squiggles |
| `annotations` | `style`: `'filled'` or `'outline'`, default `'filled'` | Annotations and side annotations |
| `footnotes` | `sticky`: default `false`. `style`: `'filled'` or `'outline'`, default `'outline'` | Footnotes |
| `hiddenLines` | None | Hidden lines |
| `shellCopy` | `prompts`: default `['$ ', '> ']`. Python `>>>` sessions need no prompt setting | Smart shell copy |
| `wordDiff` | `minSimilarity`: from 0 to 1, default `0.4` | Word-level diff |
| `whitespace` | None | Visible whitespace |
| `brackets` | `languages`: default `[]` | Colourised brackets |
| `swatches` | `languages`: default `'all'`. `formats`: default every format. `shape`: `'square'`, `'rounded'` or `'circle'`, default `'rounded'`. `size`: default `'0.8em'`. `hover`, `copy`: default `true`. `prose`: default `false`. `match`: `'value'` or `'all'`, default `'value'`. `delimiters`: `{ before?, after? }` lists of extra text around a value. `byLanguage`: `match` and `delimiters` per language | Colour swatches |
| `fileIcons` | `set`: `'vscode-icons'`, `'material'`, `'catppuccin'` or `'seti'`, default `'vscode-icons'`. Material and Catppuccin need their package, such as `@iconify-json/catppuccin`. `style`: `'plain'` or `'tile'`, default `'plain'`. `languages`: `icon`, `colour` and `style` for each language. `files`: icon names, or `false` for no icon, by file name, extension or path pattern. `icons`: custom icons by name | File icons |
| `codeLinks` | None | Code links |
| `apiLinks` | `adapters`: default `[python(), nextflow()]` | API auto-linking |
| `expandable` | `lines`: default `12`. `auto`: a line count, default `false` | Expandable blocks |
| `playgrounds` | Custom playgrounds by name | Open in playground |
| `mentions` | None | Code mentions |
| `permalinks` | None | Line permalinks |
| `placeholders` | `storage`: `'local'`, `'session'` or `'none'`, default `'local'` | Fill-in placeholders |
| `codeTabs` | `control`: `'tabs'` or `'menu'`, default `'tabs'` | Code tabs |
| `walkthrough` | `start`: `1` (default) or `'last'`, the step that each walkthrough opens on | Code walkthrough |
| `scrollycoding` | None | Scrollycoding |
| `inlineHighlighting` | `defaultLanguage`: default `false` | Inline code highlighting |
| `runnable` | `runtimes`: language to module path, added to the built-in Python, JavaScript and TypeScript runtimes, or replacing them. `timeout`: default `10000` ms, at most `2147483647`. `label`: default `'Run code'`, or its translation. `againLabel`: default `'Run again'`, or its translation. Both replace the translations. `button`: `'below'` (default), `'title'` or `'both'`. `outputDelay`: default `200` ms | Run code |

### Options for one block

Put a setting on the fence line, with its name from `codeblocks()`, to change it for that block only. The fence line is the first line of the code block, with the language. The forms:

- `focus.style="dim"`
- `focus.style="blur"`
- `lineStates.prefix=false`
- `lineStates.prefix=true`
- `footnotes.sticky=false`
- `footnotes.sticky=true`
- `footnotes.style="filled"`
- `footnotes.style="outline"`
- `annotations.style="filled"`
- `annotations.style="outline"`
- `shellCopy.prompts="<text>"` (repeat it for more prompts)
- `wordDiff.minSimilarity=<0-1>`
- `expandable.lines=<N>`
- `placeholders.storage="local"`
- `placeholders.storage="session"`
- `placeholders.storage="none"`
- `swatches.shape="square"`, `swatches.shape="rounded"`, `swatches.shape="circle"`
- `swatches.match="all"`, `swatches.match="value"`
- `fileIcons.set="<set>"`, `fileIcons.style="plain"`, `fileIcons.style="tile"`, `fileIcons.colour="<colour>"`
- `runnable.label="<text>"`
- `runnable.againLabel="<text>"`
- `runnable.timeout=<ms>`
- `runnable.button="below"`, `runnable.button="title"`, `runnable.button="both"`
- `runnable.output="<text>"`, `runnable.outputDelay=<ms>`

A bad value gives a build warning, and the block uses the site setting. Options for the whole site (states, comment syntaxes, adapters, playgrounds, runtimes) have no fence line form.

A directive of a feature that is off stays in the code, with a build warning. With `notation: false`, the plugin reads no directives, and every comment renders as written. Attributes still work.

`notation.comments` adds or replaces the comment syntax of a language. A syntax with a space in it is a block comment: the opener, then the closer. An empty list removes a language.

```js
codeblocks({
  notation: { comments: { cypher: ['//'], jinja: ['{# #}'] } },
});
```

## Colours and sizes

Every colour and size is an Expressive Code style setting. The default colours come from the site's Expressive Code theme (for example `terminal.ansiBlue` for `codeblocks.accent`), with contrast correction, so they follow any theme. Change them in `styleOverrides` of the Expressive Code options. A string applies to both themes. A pair gives the dark value, then the light value.

The shared settings are in the `codeblocks` group. Each feature has its own group, such as `codeblocksFocus`. The list: https://ewels.github.io/starlight-codeblocks/reference/style-settings/

```js
starlight({
  expressiveCode: {
    styleOverrides: {
      codeblocks: { accent: ['#c792ea', '#7c3aed'] },
    },
  },
  plugins: [codeblocks()],
});
```

The defaults meet WCAG 2.2 AA contrast. A new colour must keep the contrast that the style settings page gives for it.

## Translations

Docs: https://ewels.github.io/starlight-codeblocks/configuration/#translate-the-plugin-texts

The buttons, labels and screen reader texts of the plugin are Starlight UI strings, with keys that start with `starlightCodeblocks.`. The package has texts in English and 17 more languages. The `lang` of a locale must be one of these codes:

`ca`, `de`, `eo`, `es`, `fr`, `hi`, `id`, `it`, `ja`, `ko`, `nb`, `no`, `pl`, `pt`, `pt-BR`, `pt-PT`, `ru`, `sv`, `tr`, `zh`, `zh-CN`, `zh-TW`

The code must match exactly: `de-AT` gets none. Portuguese uses Brazilian words, and Chinese uses simplified characters.

To change a text, or to add a language, a site adds the key in `src/content/i18n/<lang>.json`, next to Starlight's own strings. The site needs the `i18n` collection in `src/content.config.ts`. The site's texts come first. A key that a language does not have takes the text of the default language, then English.

```json
{
  "starlightCodeblocks.runnable.run": "Code starten",
  "starlightCodeblocks.hiddenLines.marker_one": "{{count}} versteckte Zeile",
  "starlightCodeblocks.hiddenLines.marker_other": "{{count}} versteckte Zeilen"
}
```

- A text with a count has a key for each plural form of the language, with the i18next suffixes. English has `_one` and `_other`, and Polish also has `_few` and `_many`. To change such a text, give all its plural keys.
- Keep `{{count}}` and every other name in double braces. The plugin fills in their values.
- `runnable.label` and `runnable.againLabel` replace the translations. Custom line states and playgrounds show the label from their options.
- Sites without Starlight show the English texts.

## ec.config.mjs

If `ec.config.mjs` has no `plugins` list, the plugin needs no change. If it has a `plugins` list, Expressive Code uses only that list. Add `pluginCodeblocks()` to it:

```js
import { pluginCollapsibleSections } from '@expressive-code/plugin-collapsible-sections';
import { pluginCodeblocks } from 'starlight-codeblocks/expressive-code';

export default {
  plugins: [pluginCollapsibleSections(), pluginCodeblocks()],
};
```

Keep `codeblocks()` in `astro.config.mjs`, and give it all the options. `pluginCodeblocks()` with no argument reads the same options. `styleOverrides` can stay in `ec.config.mjs`. `pluginCodeblocks()` can go before or after most other plugins. `pluginCollapsibleSections()` must come before it, or the build stops.

## Other plugins

Docs: https://ewels.github.io/starlight-codeblocks/reference/plugin-compatibility/

- List `codeblocks()` before any plugin or theme that sets the Starlight `expressiveCode` option. If a theme must come first, set `expressiveCode: {}` in the Starlight config.
- `starlight-links-validator`: give it `exclude: linksValidatorExclude`, imported from `starlight-codeblocks`. It then skips `#mention:` links and links to a block `id`, such as `#cfg-L2`.
- `starlight-llms-txt`: add `customSelectors: { all: ['.scb-deco'] }`, so the labels, buttons and notes of the plugin stay out of the code in the text files.
- `starlight-versions` and `starlight-md-txt` read `.md` files as MDX. Put the `{:lang}` suffix of inline code highlighting inside the backticks.
- `expressive-code-twoslash` removes some lines. A directive on a removed line has no effect, with a warning.

## Markdown processors and file types

- Sätteri, the default Markdown processor of Astro 7, needs no set-up.
- With `markdown: { processor: unified() }`, every feature works. The plugin adds its own remark plugin. Keep Sätteri unless another plugin needs `unified()`.
- In Markdoc (`.mdoc`) files, give the fence line attributes in a `meta` attribute: `` ```js {% title="app.js" meta="focus={2}" %} ``. Directives work. Code tabs, inline code highlighting, the check of code mention links and the check for duplicate block ids do not run.
- The `<Code>` component of Starlight gets the features inside a block. Give attributes in its `meta` prop. Code tabs do not work around it.

## Sites without Starlight

Docs: https://ewels.github.io/starlight-codeblocks/install/astro/

An Astro site without Starlight adds `codeblocks()` from `starlight-codeblocks/astro` to `integrations`, before `mdx()`. It adds Expressive Code itself, so remove `expressiveCode()`. Install `astro-expressive-code` next to the plugin.

```js
import mdx from '@astrojs/mdx';
import { defineConfig } from 'astro/config';
import codeblocks from 'starlight-codeblocks/astro';

export default defineConfig({
  integrations: [codeblocks({ expressiveCode: { themes: ['github-dark', 'github-light'] } }), mdx()],
});
```

It takes the same options as on Starlight, plus `expressiveCode` for the options of `astro-expressive-code`. Code tabs and inline code highlighting are off until the site sets `codeTabs: {}` and `inlineHighlighting: {}`, because they change Markdown outside the code blocks. API cards show on links in code blocks only. Side annotations and scrollycoding stay inside the content column. With `unified()`, code tabs need `remark-directive` in its `remarkPlugins`.

A site that keeps its own `expressiveCode()` gives the options to `pluginCodeblocks()` in `ec.config.mjs`, with `tabWidth: 0`. If the site has a `base`, give it as the second argument, `pluginCodeblocks(options, { base: '/docs' })`. Code tabs and inline code highlighting then do not work.

Without `codeblocks()`, `runnable.runtimes` values are URLs that the browser imports as they are, and there are no built-in runtimes. Map each language to the URL of its runtime module, such as `https://cdn.jsdelivr.net/npm/starlight-codeblocks/dist/runtimes/pyodide.mjs`. The JavaScript and TypeScript runtimes are `javascript.mjs` and `typescript.mjs` in the same folder. The single plugins, such as `pluginFocus()`, are in https://ewels.github.io/starlight-codeblocks/reference/expressive-code-plugins/
