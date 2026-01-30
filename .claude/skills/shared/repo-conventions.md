# Repository Conventions

Shared reference for all training repository skills.

## Working Directory

All skill commands execute from the **repository root** - the directory containing:

- `docs/en/mkdocs.yml`
- `docs/`
- `.github/`

Verify location before running commands. Use paths relative to repository root.

## Directory Mapping

| Content Type         | Documentation                        | Working Directory       | Solutions                         |
| -------------------- | ------------------------------------ | ----------------------- | --------------------------------- |
| hello_nextflow       | `docs/en/docs/hello_nextflow/`       | `hello-nextflow/`       | `hello-nextflow/solutions/`       |
| hello_nf-core        | `docs/en/docs/hello_nf-core/`        | `hello-nf-core/`        | `hello-nf-core/solutions/`        |
| nf4_science/genomics | `docs/en/docs/nf4_science/genomics/` | `nf4-science/genomics/` | `nf4-science/genomics/solutions/` |
| nf4_science/rnaseq   | `docs/en/docs/nf4_science/rnaseq/`   | `nf4-science/rnaseq/`   | `nf4-science/rnaseq/solutions/`   |
| side_quests/\*       | `docs/en/docs/side_quests/<name>.md` | `side-quests/<name>/`   | `side-quests/solutions/<name>/`   |

**Note**: Documentation uses underscores (`hello_nextflow`), working directories use hyphens (`hello-nextflow`).

## Common Paths

| Purpose                 | Path                               |
| ----------------------- | ---------------------------------- |
| Heading checker         | `.github/check_headings.py`        |
| MkDocs config (English) | `docs/en/mkdocs.yml`               |
| Site navigation         | `docs/en/mkdocs.yml` (nav section) |
| Translation configs     | `docs/{lang}/mkdocs.yml`           |
| Devcontainer config     | `.devcontainer/devcontainer.json`  |
| Contributing guide      | `CONTRIBUTING.md`                  |
| Translation guide       | `TRANSLATING.md`                   |

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
