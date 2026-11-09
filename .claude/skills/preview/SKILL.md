---
name: preview
description: Start the local Astro dev server to view training material changes. Use when editing documentation, checking formatting, or testing how content renders before committing.
---

Start the local Astro (Starlight) dev server to view changes to training materials.

See [../shared/repo-conventions.md](../shared/repo-conventions.md) for directory structure.

All commands run in `docs/`.

**First, check if a server is already running:**

```bash
cd docs && npx astro dev status
```

- If it reports a running server, the site is already at the URL it prints (usually http://localhost:4321/)
- If it reports no server, proceed to start one

**Install dependencies (first time, or after `package-lock.json` changes):**

```bash
cd docs && npm ci
```

Needs Node 22 or later.

**Start the server in the background:**

```bash
cd docs && npm run dev -- --background
```

**After starting the server:**

- Open http://localhost:4321/ (a page such as `hello_nextflow/01_hello_world.mdx` is at http://localhost:4321/hello_nextflow/01_hello_world/)
- Translations are at http://localhost:4321/<lang>/...
- Pages reload when you save changes
- Read errors with `cd docs && npx astro dev logs`
- Use `/stop-preview` to stop the server

**Common issues:**

- If port 4321 is busy, check `npx astro dev status` and use `/stop-preview`, or start on another port with `npm run dev -- --background --port 4322`
- If a page shows an error overlay, the message names the file and line. Common causes are MDX syntax (a stray `{`, `<` or HTML comment in prose) or invalid frontmatter
- If a new page is missing from the navigation, add it to `docs/src/sidebar.mjs`
- The dev server does not run the links validator. Run a full build (`/validate`) to check links
