# Link code

Three features add links to code blocks.

- If one piece of text needs a link that you choose, use a code link.
- If every library name in Python or Nextflow code needs a link to its reference, use API auto-linking. It starts on its own.
- If readers must share a link to one line or a range, use line permalinks.

## Code links

Docs: https://ewels.github.io/starlight-codeblocks/features/code-links/

Use for a one-off link from text in the code, or for a name that no API link adapter knows. Give each link a short description: it shows in the same card as API auto-linking. The linked text keeps its syntax colours, with an underline in the accent colour.

| Syntax | Where |
|---|---|
| `[!link /<text>/ <url>] <description>` | Comment at the end of the target line, or on its own line directly above it. The description shows in a card on hover and focus. |
| `[!link /<text>/ <url>]` | Same, as a plain link with no card |

```py
# [!link /linspace/ https://numpy.org/doc/stable/reference/generated/numpy.linspace.html] Returns evenly spaced numbers over an interval.
x = np.linspace(0, 1, 50)
```

- The plugin links the first whole-word match of the literal text on the target line, else the first match. Make the text longer to match a later place.
- Text after the directive shows in the same card as API auto-linking: the linked text, the description and the URL's domain. Use it to say what the link goes to. The card shows plain text only.
- Stack several `[!link]` lines above one line to link several texts on it.
- A URL that starts with `/` is a link inside the site. The plugin adds Astro's `base`.
- The URL is relative, `http` or `https`. If the text has no match, or the URL is missing or has another scheme, the build logs a warning.
- Names on the target line of a `[!link]` do not get API links, so a link is never inside another link.
- Option: `codeLinks: false` turns the feature off. The directive then stays in the code, with a warning.
- Limits: a link cannot span two lines. JSON has no comments: use `jsonc`.

## API auto-linking

Docs: https://ewels.github.io/starlight-codeblocks/features/api-auto-linking/

API auto-linking finds library names in a block and links each one to its reference page. A card on hover or focus shows the signature or the kind, a summary, and the source with the icon of its project. You write nothing in the block. An adapter links a name only when it is certain about it.

| Syntax | Where |
|---|---|
| `apiLinks=false` | Code block fence line: turns the feature off for one block |
| `pydocsBase="<base>"` | Code block fence line of a Python block: prefers the starlight-pydocs package at this base |
| `nextflow="config"` | Code block fence line of a Nextflow block: links it as a config file |

- The Python adapter runs on `py`, `python` and `pycon` blocks. It reads `import x`, `import x as y` and `from a import b as c`, then links the names that they bind and attribute chains on them, such as `os.path.join`.
- The Nextflow adapter runs on `nextflow` and `nf` blocks. It links channel factories, such as `channel.of`, and operators after them. With the `modules` option, it links processes and workflows from `include` statements.
- In a `process` body, the Nextflow adapter links the directives, such as `cpus`, and the `input:` and `output:` labels. It also links global names, such as `params` and `file()`, and the members of `workflow`, `log` and `nextflow`, such as `workflow.projectDir`.
- A Nextflow block is a config file if its title ends in `.config`. For a config block with another title or none, add `nextflow="config"` to the fence line. In a config file, the adapter links the config scopes, such as `docker` in `docker.enabled = true`, `params` and `profiles`, and the directives in the `process` scope.
- The build caches fetched indexes in `node_modules/.cache/starlight-codeblocks/`. A failed fetch logs a warning and does not fail the build.
- Option: `apiLinks.adapters` (default `[python(), nextflow()]`). A list replaces the default list, so give every adapter that you want:

  ```js
  import codeblocks from 'starlight-codeblocks';
  import { nextflow } from 'starlight-codeblocks/adapters/nextflow';
  import { python } from 'starlight-codeblocks/adapters/python';

  codeblocks({
    apiLinks: {
      adapters: [python({ inventories: ['https://numpy.org/doc/stable/objects.inv'] }), nextflow()],
    },
  });
  ```

- Python adapter options: `stdlib` (default `true`), `inventories` (Sphinx `objects.inv` URLs) and `summaries` (default `true`: fetch each documentation page for a short summary on the card).
- starlight-pydocs: install both plugins. The Python adapter links every package that starlight-pydocs documents, with signatures and summaries, and needs no configuration.
- `pydocsBase="<base>"` on a Python block links its names to the starlight-pydocs package at that base first, for example an older version at `1x/api/myproject`. Names that the package does not have link as usual. starlight-pydocs sets it on its docstring examples.
- Nextflow adapter option: `modules`, a function of `{ name, path }` that returns a URL, an object with `href`, or `undefined`.
- Another language needs an adapter of your own: https://ewels.github.io/starlight-codeblocks/extend/write-an-api-link-adapter/
- Another plugin can give its own links the card. Put `data-scb-api-links` on an element that holds the links. Put `data-scb-api-head` on each `a`, with optional `data-scb-api-summary`, `data-scb-api-source`, `data-scb-api-icon` and `data-scb-api-action`. `data-scb-api-icons` on the element that holds the links maps each icon key to an SVG path. This works on pages with no code blocks, and needs `codeblocks()`. starlight-pydocs does this for the types in its signatures.
- Limits: the Python adapter does not follow assignments and does not link built-ins such as `print`. A name that the block binds again does not link.

## Line permalinks

Docs: https://ewels.github.io/starlight-codeblocks/features/line-permalinks/

Use when readers must link to exact lines, for example to ask about one setting. The `id` gives the block line numbers, and each number is a link to its line. Do not use to highlight lines for the reader: use focus, line states or code mentions.

| Syntax | Where |
|---|---|
| `id="<id>"` | Code block fence line |

- Each number links to `#<id>-L<n>`. Shift and a second number link a range, `#<id>-L<a>-L<b>`. `#<id>` links to the whole block.
- Pick an `id` from the content, such as the file name. It must be unique on the page, and must not match a heading id such as `options`. A duplicate gives a build warning.
- `startLineNumber` of Expressive Code sets the first number.
- A link to a hidden line opens its marker. A link to a collapsed line opens an expandable block.
- `starlight-links-validator` reports these links as broken. Give it `exclude: linksValidatorExclude`.
- Option: `permalinks: false` ignores `id`.
- Limits: the line numbers take space, so long lines scroll sooner on phones. The page highlights one range at a time.
