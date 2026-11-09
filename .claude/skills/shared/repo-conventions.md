# Repository Conventions

Shared reference for all training repository skills.

## Working Directory

All skill commands execute from the **repository root** - the directory containing:

- `docs/astro.config.mjs`
- `docs/src/content/docs/`
- `.github/`

Verify location before running commands. Use paths relative to repository root.

## Directory Mapping

Pages live under `docs/src/content/docs/` (written `<pages>/` below).

| Content Type         | Documentation                      | Working Directory       | Solutions                                                          |
| -------------------- | ---------------------------------- | ----------------------- | ------------------------------------------------------------------ |
| hello_nextflow       | `<pages>/hello_nextflow/`          | `hello-nextflow/`       | `hello-nextflow/solutions/`                                        |
| nextflow_run         | `<pages>/nextflow_run/`            | `nextflow-run/`         | -                                                                  |
| nfcore_use           | `<pages>/nfcore_use/`              | `nfcore-use/`           | `nfcore-use/solutions/`                                            |
| nfcore_build         | `<pages>/nfcore_build/`            | `nfcore-build/`         | `nfcore-build/solutions/`                                          |
| config_exec          | `<pages>/config_exec/`             | `config-exec/`          | -                                                                  |
| nf4_science/genomics | `<pages>/nf4_science/genomics/`    | `nf4-science/genomics/` | `nf4-science/genomics/solutions/`                                  |
| nf4_science/rnaseq   | `<pages>/nf4_science/rnaseq/`      | `nf4-science/rnaseq/`   | `nf4-science/rnaseq/solutions/`                                    |
| side_quests/\*       | `<pages>/side_quests/<name>/*.mdx` | `side-quests/<name>/`   | `side-quests/solutions/<name>/` or `side-quests/<name>/solutions/` |

**Note**: Documentation uses underscores (`hello_nextflow`), working directories use hyphens (`hello-nextflow`).

Translations sit at the same relative paths under `<pages>/<lang>/` (ca de es fr hi it ko pl pt tr). They are AI-generated: never edit or review them, and leave them out of searches and checks.

## Common Paths

| Purpose              | Path                                   |
| -------------------- | -------------------------------------- |
| Heading checker      | `.github/check_headings.py`            |
| Code block audit     | `docs/scripts/check-fences.mjs`        |
| Site config          | `docs/astro.config.mjs`                |
| Site navigation      | `docs/src/sidebar.mjs`                 |
| Frontmatter schema   | `docs/src/content.config.ts`           |
| UI strings           | `docs/src/content/i18n/<lang>.json`    |
| Translation prompts  | `_scripts/translate/languages/<lang>/` |
| Devcontainer config  | `.devcontainer/devcontainer.json`      |
| Contributing guide   | `CONTRIBUTING.md`                      |
| Translation guide    | `TRANSLATING.md`                       |
| Code block reference | `.claude/skills/starlight-codeblocks/` |

## File Conventions

### Nextflow Scripts (.nf)

- Must start with `#!/usr/bin/env nextflow` (except module files, `**/modules/*.nf`)
- DSL2 syntax only
- Process names in UPPERCASE
- Located in working directories, not docs/

### Pages (.mdx)

- Frontmatter `title:` and `description:`; no H1 in the body
- Heading numbering: `## 1.`, `### 1.1.` with trailing periods
- One sentence per line
- Code blocks: ` ```nextflow ` for Nextflow, with `title="..."`, `showLineNumbers`, `startLineNumber=N` and `{3,5-7}` highlights where appropriate
- Asides `:::note` / `:::tip` / `:::caution` / `:::danger` ... `:::`; collapsibles `<details class="TYPE">` + `<summary>`, with blank lines around the body
- Before/After: `<CodeWalkthrough>` with `step="Before"` and `step="After"` blocks
- Command output: inside the `bash` block, after a `# [!output]` line
- Components are auto-imported: no import lines except `?raw` code includes
- New pages are listed in `docs/src/sidebar.mjs`

## Environment

### Nextflow Version

Read from devcontainer.json:

```bash
grep -o '"NXF_VER":\s*"[^"]*"' .devcontainer/devcontainer.json | cut -d'"' -f4
```

### Strict Syntax Parser

Always use `NXF_SYNTAX_PARSER=v2` for testing.

### Site Preview

The site needs Node 22+. From `docs/`: `npm ci` once, then `npm run dev` (http://localhost:4321/). In an agent session it starts in the background; `npx astro dev status`, `npx astro dev logs` and `npx astro dev stop` manage it.

### Console Output Mode (CRITICAL)

Nextflow 26.04+ detects AI agent environments and switches to a machine-readable "agent mode" console format (`[PIPELINE]`, `[PROCESS]`, `[SUCCESS]`/`[FAILED]` tags) instead of the classic human-facing output (the `N E X T F L O W` banner, `executor >` lines, `[hash] process | N of M ✔` summaries) that learners actually see in their terminal.

Detection includes the `CLAUDECODE` environment variable, which Claude Code sets in every shell it runs. This means **any direct `nextflow run` invocation from a skill will silently render in agent mode** unless explicitly overridden, producing console output that does not match what a learner sees.

**Whenever a command's console output will be copied into learner-facing documentation (or compared against it), force human mode:**

```bash
env -u CLAUDECODE NXF_AGENT_MODE=false nextflow run ...
```

Running Nextflow inside a Docker container (`docker exec ...`) is not affected today, since a container does not inherit the host shell's environment unless explicitly passed with `-e`. Do not rely on that as the safeguard, though: apply the override explicitly any time Nextflow runs directly on the host, so the behavior doesn't silently break if the Docker invocation ever changes to forward host environment variables.
