# Adapt to the reader

Two features let the reader change what a block shows.

- If the same task comes in versions and each version is one code block, use code tabs.
- If a version needs prose, a list or more than one code block, use the `<Tabs>` component of Starlight in an MDX file.
- If the page is a `.md` file, use code tabs, or change the file to `.mdx` for `<Tabs>`.
- If only a value differs for each reader, such as a token or an ID, use fill-in placeholders in one block.

## Code tabs

Docs: https://ewels.github.io/starlight-codeblocks/features/code-tabs/

Use for the files of one example, install commands for each package manager, or one example in several languages. The variants become one block with a tab for each variant in the title bar, as in a code editor. The block takes no more space than a plain block.

| Syntax | Where |
|---|---|
| `:::code-tabs` ... `:::` | Around two or more code blocks. No space after the colons. |
| `sync="<key>"` | After `:::code-tabs`, in braces |
| `control="tabs"` or `control="menu"` | After `:::code-tabs`, in braces |
| `label="<text>"` | Code block fence line of each variant |

````md
:::code-tabs{sync="pm"}
```sh label="npm"
npm install starlight-codeblocks
```
```sh label="pnpm"
pnpm add starlight-codeblocks
```
:::
````

- The directive can contain only code blocks. A paragraph inside it fails the build.
- A tab shows the `title` of its variant, with its file icon. A variant with no `title` shows its `label`, with no icon unless the block has `icon="..."`. A variant with no `label` shows the name of its language, such as "Python" for `py`. Two variants with the same label fail the build: `sh` and `bash` both show "Shell", so give each one a `label`.
- Tabs use the editor frame, also for shell commands.
- `control="menu"` shows the `title` of the variant and a menu on the right in place of tabs. Option `codeTabs.control` (default `'tabs'`) sets it for the whole site.
- Blocks with the same `sync` key switch together, matched by exact label text. The browser keeps the choice in `localStorage`, so it applies across the site. Give every package manager block the same key, such as `pm`.
- The copy button copies the variant that shows. Without JavaScript, the first variant shows.
- `sync` and the `syncKey` of `<Tabs>` are separate, and do not switch together.
- Option: with `codeTabs: false`, each block in the directive renders on its own, and the directive lines do not show.
- Limits: works in `.md` and `.mdx` files, not in Markdoc and not around a `<Code>` component. Needs Starlight.

## Fill-in placeholders

Docs: https://ewels.github.io/starlight-codeblocks/features/fill-in-placeholders/

Use for values that each reader must change, such as `YOUR_TOKEN` or `WORKSPACE_ID`. Each match becomes a field. What a reader types fills every field with the same text on the site, and the copied code. Do not use for alternatives that change more than a value: use code tabs.

| Syntax | Where |
|---|---|
| `placeholder="<A>,<B>"` | Code block fence line: the literal texts, separated by commas |

````md
```sh placeholder="YOUR_TOKEN,WORKSPACE_ID"
curl -H "Authorization: Bearer YOUR_TOKEN" \
  https://api.example.com/workspaces/WORKSPACE_ID/runs
```
````

- Every match of each text becomes a field. Choose texts that appear nowhere else in the block.
- The copy button, playgrounds and the **Run code** button get the reader's values. An empty field copies its placeholder text.
- A text that is not in the code gives a build warning.
- Option: `placeholders.storage` is `'local'` (default, kept across visits), `'session'` (until the tab closes) or `'none'`. Use `'session'` on a site where readers type secrets, and tell readers where the values go.
- Limits: literal text only. Without JavaScript, fields are plain inputs that do not update each other.
