# Translation Guide

This document describes how translations work for the Nextflow training materials.

## How Translations Work

The Nextflow training uses an **AI-assisted translation system** with human review:

1. **LLM Translation**: Claude translates English content using language-specific prompts
2. **Automatic Post-processing**: Scripts fix common LLM mistakes (code blocks, links, etc.)
3. **Human Review**: Native speakers review and improve translations via PRs
4. **English Fallback**: Untranslated pages automatically display English content with a warning banner

## Directory Structure

```
docs/
├── en/                     # English (source language)
│   ├── mkdocs.yml          # Main config (inherited by other languages)
│   ├── overrides/          # Theme overrides
│   ├── hooks/              # MkDocs hooks
│   └── docs/               # English content
├── pt/                     # Portuguese
│   ├── mkdocs.yml          # Inherits from ../en/mkdocs.yml
│   ├── llm-prompt.md       # Translation rules & glossary
│   └── docs/               # Portuguese translations
├── es/                     # Spanish
│   └── ...
└── ...
```

## Languages

| Code | Language   | Status |
| ---- | ---------- | ------ |
| en   | English    | Source |
| pt   | Portuguese | Active |
| es   | Spanish    | Active |
| fr   | French     | Active |
| it   | Italian    | Active |
| ko   | Korean     | Active |
| pl   | Polish     | Active |
| tr   | Turkish    | Active |

## Contributing to Translations

### Option 1: Improve Translation Quality (Recommended)

The best way to improve translations is to **edit the language-specific prompt file** (`docs/{lang}/llm-prompt.md`). This ensures:

- Future translations follow better rules
- Consistency across all pages
- Less manual editing needed

**To improve the prompt:**

1. Fork the repository
2. Edit `docs/{lang}/llm-prompt.md`
3. Add or update translation rules, glossary terms, or grammar preferences
4. Submit a PR with your changes

### Option 2: Direct Translation Edits

You can also directly edit translated files in `docs/{lang}/docs/`:

1. Fork the repository
2. Edit the markdown file in your language directory
3. Submit a PR with your changes

**Note**: Direct edits may be overwritten if the English source changes and translations are regenerated. Prompt improvements are more durable.

### Option 3: Add a New Language

To add support for a new language:

1. Run `uv run python _scripts/docs.py new-lang <lang-code>`
2. Edit `docs/<lang>/llm-prompt.md` with language-specific rules
3. Add the language to `SUPPORTED_LANGS` in `_scripts/docs.py`
4. Add the language to `extra.alternate` in `docs/en/mkdocs.yml`
5. Submit a PR

## Translation Glossary Guidelines

Each language has a `llm-prompt.md` file containing:

### Terms to Keep in English

Technical terms that should NOT be translated:

- Nextflow syntax: `channel`, `process`, `workflow`, `emit`, `take`, etc.
- Directives: `publishDir`, `container`, `conda`, `memory`, etc.
- Operators: `map`, `filter`, `collect`, `join`, etc.
- Tools: `Nextflow`, `nf-core`, `Docker`, `Singularity`, etc.

### Terms to Translate

Common terms that SHOULD be translated:

| English     | Example Translation |
| ----------- | ------------------- |
| directory   | diretório (pt)      |
| environment | ambiente (pt)       |
| input       | entrada (pt)        |
| output      | saída (pt)          |
| task        | tarefa (pt)         |

### Admonition Titles

Standard translations for admonition types (Note, Tip, Warning, Exercise, Solution).

## Running Translations Locally

### Prerequisites

- Python 3.11+
- [uv](https://docs.astral.sh/uv/) (recommended) or pip
- `ANTHROPIC_API_KEY` environment variable

### Commands

```bash
cd _scripts

# List files missing translation
uv run python translate.py list-missing pt

# List outdated translations (English changed after translation)
uv run python translate.py list-outdated pt

# Translate a single page
uv run python translate.py translate-page --language pt --en-path docs/en/docs/index.md

# Translate all missing pages for a language
uv run python translate.py translate-lang --language pt

# Update outdated translations
uv run python translate.py update-outdated --language pt

# Add missing translations
uv run python translate.py add-missing --language pt

# Fix common issues in translations
uv run python translation_fixer.py fix-all pt
```

## Preview Translations

### Docker

```bash
# Preview Portuguese (or any language)
docker run --rm -it -p 8000:8000 -v ${PWD}:/docs \
  -w /docs/docs/pt \
  ghcr.io/nextflow-io/training-mkdocs:latest
```

### Python

```bash
# Preview specific language
cd _scripts
uv run python docs.py live pt

# Or directly with mkdocs
cd docs/pt
mkdocs serve --dev-addr 127.0.0.1:8008
```

## CI/CD Workflows

### Automatic Builds

When changes are pushed to `master`:

1. All languages are built in parallel
2. Built sites are combined (English at root, others in subdirectories)
3. Combined site is deployed with version tracking

### Translation Workflow

Use the **Translate** GitHub Action (Actions → Translate → Run workflow):

1. Select a language
2. Select a command:
   - `list-outdated`: Show translations older than English source
   - `list-missing`: Show untranslated files
   - `update-outdated`: Update stale translations
   - `add-missing`: Add new translations
   - `update-and-add`: Both update and add
3. The workflow creates a PR with the changes

## Notes

### Code Blocks

Only **comments** in code blocks are translated. The code itself remains unchanged:

```groovy
// This comment is translated
Channel.fromPath('data/*.fastq')  // Code stays the same
    .set { reads_ch }
```

### Admonitions

Admonition titles are translated, but keywords remain in English:

```markdown
!!! note "Nota"

    Conteúdo traduzido aqui.
```

### Headings

Heading text is translated, but anchors (`{ #anchor }`) are preserved for link compatibility:

```markdown
## 1. Primeiros Passos { #getting-started }
```

## References

- [FastAPI Translation System](https://github.com/fastapi/fastapi/tree/master/scripts) - Inspiration for this implementation
- [Portuguese Glossary (Google Sheets)](https://docs.google.com/spreadsheets/d/1HUa3BO2kwukhX4EXQ-1blXeP5iueUdM23OwDRpfarDg/edit)
