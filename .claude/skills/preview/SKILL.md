---
name: preview
description: Start local MkDocs preview server to view training material changes. Use when editing documentation, checking formatting, or testing how content renders before committing.
---

Start the local MkDocs preview server to view changes to training materials.

See [../shared/repo-conventions.md](../shared/repo-conventions.md) for directory structure.

**First, check if a server is already running:**

```bash
docker ps --filter "ancestor=ghcr.io/nextflow-io/training-mkdocs:latest" --format "{{.ID}}"
```

- If a container ID is returned, the server is already running at http://127.0.0.1:8000/
- If no output, proceed to start a new server

**Start the server using Docker (recommended):**

```bash
docker run --rm -d -p 8000:8000 -v ${PWD}:/docs -w /docs/docs/en ghcr.io/nextflow-io/training-mkdocs:latest
```

**Alternative - using uv:**

```bash
uv run _scripts/docs.py serve
```

**After starting the server:**

- Open http://127.0.0.1:8000/ (or http://0.0.0.0:8000/ for Docker)
- Pages auto-refresh when you save changes
- Check console for any build errors
- Use `/stop-preview` command or Ctrl+C to stop the server

**Common issues:**

- If port 8000 is busy, use `/stop-preview` to stop other servers
- If changes don't appear, check for markdown syntax errors in the console
