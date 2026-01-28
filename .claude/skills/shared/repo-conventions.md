# Repository Conventions

Shared reference for all training repository skills.

## Working Directory

All skill commands execute from the **repository root** - the directory containing:
- `mkdocs.yml`
- `docs/`
- `.github/`

Verify location before running commands. Use paths relative to repository root.

## Directory Mapping

| Content Type | Documentation | Working Directory | Solutions |
|--------------|---------------|-------------------|-----------|
| hello_nextflow | `docs/hello_nextflow/` | `hello-nextflow/` | `hello-nextflow/solutions/` |
| hello_nf-core | `docs/hello_nf-core/` | `hello-nf-core/` | `hello-nf-core/solutions/` |
| nf4_science/genomics | `docs/nf4_science/genomics/` | `nf4-science/genomics/` | `nf4-science/genomics/solutions/` |
| nf4_science/rnaseq | `docs/nf4_science/rnaseq/` | `nf4-science/rnaseq/` | `nf4-science/rnaseq/solutions/` |
| side_quests/* | `docs/side_quests/<name>.md` | `side-quests/<name>/` | `side-quests/solutions/<name>/` |

**Note**: Documentation uses underscores (`hello_nextflow`), working directories use hyphens (`hello-nextflow`).

## Common Paths

| Purpose | Path |
|---------|------|
| Heading checker | `.github/check_headings.py` |
| MkDocs config | `mkdocs.yml` |
| Site navigation | `mkdocs.yml` (nav section) |
| Devcontainer config | `.devcontainer/devcontainer.json` |
| Contributing guide | `CONTRIBUTING.md` |

## File Conventions

### Nextflow Scripts (.nf)

- Must start with `#!/usr/bin/env nextflow`
- DSL2 syntax only
- Process names in UPPERCASE
- Located in working directories, not docs/

### Markdown Files (.md)

- Heading numbering: `## 1.`, `### 1.1.` with trailing periods
- One sentence per line
- Code blocks with titles, linenums, and hl_lines where appropriate
- Admonitions indented 4 spaces

## Environment

### Nextflow Version

Read from devcontainer.json:
```bash
grep -o '"NXF_VER":\s*"[^"]*"' .devcontainer/devcontainer.json | cut -d'"' -f4
```

### Strict Syntax Parser

Always use `NXF_SYNTAX_PARSER=v2` for testing.
