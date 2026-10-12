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
| nfcore_build         | `docs/en/docs/nfcore_build/`         | `nfcore-build/`         | `nfcore-build/solutions/`         |
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

### Console Output Mode (CRITICAL)

Nextflow 26.04+ detects AI agent environments and switches to a machine-readable "agent mode" console format (`[PIPELINE]`, `[PROCESS]`, `[SUCCESS]`/`[FAILED]` tags) instead of the classic human-facing output (the `N E X T F L O W` banner, `executor >` lines, `[hash] process | N of M ✔` summaries) that learners actually see in their terminal.

Detection includes the `CLAUDECODE` environment variable, which Claude Code sets in every shell it runs. This means **any direct `nextflow run` invocation from a skill will silently render in agent mode** unless explicitly overridden, producing console output that does not match what a learner sees.

**Whenever a command's console output will be copied into learner-facing documentation (or compared against it), force human mode:**

```bash
env -u CLAUDECODE NXF_AGENT_MODE=false nextflow run ...
```

Running Nextflow inside a Docker container (`docker exec ...`) is not affected today, since a container does not inherit the host shell's environment unless explicitly passed with `-e`. Do not rely on that as the safeguard, though: apply the override explicitly any time Nextflow runs directly on the host, so the behavior doesn't silently break if the Docker invocation ever changes to forward host environment variables.
