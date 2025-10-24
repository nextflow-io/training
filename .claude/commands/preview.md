---
description: Start local preview server to view training materials
---

Start the local MkDocs preview server to view changes to training materials.

**First, check if a server is already running:**

```bash
docker ps --filter "ancestor=ghcr.io/nextflow-io/training-mkdocs:latest" --format "{{.ID}}"
```

- If a container ID is returned, the server is already running at http://127.0.0.1:8000/
- If no output, proceed to start a new server

**Start the server using Docker (recommended):**

```bash
docker run --rm -d -p 8000:8000 -v ${PWD}:/docs ghcr.io/nextflow-io/training-mkdocs:latest
```

**Alternative options:**

If you encounter issues with social cards:

```bash
docker run --rm -p 8000:8000 -e 'CARDS=false' -v ${PWD}:/docs ghcr.io/nextflow-io/training-mkdocs:latest
```

If you have Python environment set up:

```bash
mkdocs serve
```

Or without social cards:

```bash
CARDS=false mkdocs serve
```

**After starting the server:**

- Open http://127.0.0.1:8000/ (or http://0.0.0.0:8000/ for Docker)
- Pages auto-refresh when you save changes
- Check console for any build errors
- Use `/stop-preview` command or Ctrl+C to stop the server

**Note:** The initial build may take a few minutes as it builds documentation in multiple languages (en, pt, es, fr, it, ko). Wait for the "Serving on" message before accessing the site.

**Common issues:**

- If social cards fail, use CARDS=false option
- If port 8000 is busy, use `/stop-preview` to stop other servers
- If changes don't appear, check for markdown syntax errors in the console
