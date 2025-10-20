---
description: Start local preview server to view training materials
---

Start the local MkDocs preview server to view changes to training materials.

Choose one of these methods:

**Option 1: Docker (Recommended)**
```bash
docker run --rm -it -p 8000:8000 -v ${PWD}:/docs ghcr.io/nextflow-io/training-mkdocs:latest
```

**Option 2: Docker without social cards** (if having issues)
```bash
docker run --rm -it -p 8000:8000 -e 'CARDS=false' -v ${PWD}:/docs ghcr.io/nextflow-io/training-mkdocs:latest
```

**Option 3: Python** (if Python environment is set up)
```bash
mkdocs serve
```

**Option 4: Python without social cards**
```bash
CARDS=false mkdocs serve
```

After starting the server:
- Open http://127.0.0.1:8000/ (or http://0.0.0.0:8000/ for Docker)
- Pages auto-refresh when you save changes
- Check console for any build errors
- Use Ctrl+C to stop the server

Common issues:
- If social cards fail, use CARDS=false option
- If port 8000 is busy, stop other preview servers
- If changes don't appear, check for markdown syntax errors in the console
