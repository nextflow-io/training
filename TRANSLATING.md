# Translation Guide

This document describes how translations work for the Nextflow training materials.

## Quick Reference

| Task                       | Command                                                         |
| -------------------------- | --------------------------------------------------------------- |
| Preview translations       | `uv run python docs.py serve <lang>`                            |
| List outdated translations | `uv run python translate.py list-outdated <lang>`               |
| List missing translations  | `uv run python translate.py list-missing <lang>`                |
| Update a single file       | `uv run python translate.py translate-page -l <lang> -p <path>` |
| Update all outdated        | `uv run python translate.py update-outdated -l <lang>`          |

---

## How Automatic Translation Updates Work

When English source files are modified, translations are automatically updated via GitHub Actions:

```
English file changed
        │
        ▼
┌───────────────────────────────────────┐
│  GitHub Actions: translate-check.yml  │
│  (Triggers on push to master)         │
└───────────────────────────────────────┘
        │
        ▼
┌───────────────────────────────────────┐
│  Detect outdated translations         │
│  (Compare git commit timestamps)      │
└───────────────────────────────────────┘
        │
        ▼
┌───────────────────────────────────────┐
│  AI updates only changed sections     │
│  (Minimal diffs for easy review)      │
└───────────────────────────────────────┘
        │
        ▼
┌───────────────────────────────────────┐
│  Create PR for each language          │
│  (One PR per language with changes)   │
└───────────────────────────────────────┘
        │
        ▼
┌───────────────────────────────────────┐
│  Human review + approve + merge       │
└───────────────────────────────────────┘
```

### Reviewing Automatic Translation PRs

When the bot creates a translation update PR:

1. **Check the diff** - The AI makes minimal changes, updating only sections that changed in English
2. **Verify accuracy** - Spot-check that translations are correct, especially for:
   - Technical terminology
   - Code block comments (should be translated)
   - Admonition titles (should match glossary)
3. **Test the build** - The CI will build the docs; check for errors
4. **Approve and merge** - If everything looks good, approve and merge

### Manual Trigger

You can also manually trigger translation updates:

1. Go to **Actions** → **Translate** → **Run workflow**
2. Select the language and command
3. The workflow will create a PR with changes

---

## How to Review and Fix Existing Translations

### Option 1: Improve the Translation Prompt (Recommended)

The best way to improve translations is to edit `docs/<lang>/llm-prompt.md`. This ensures:

- Future translations follow better rules
- Consistency across all pages
- Changes persist through re-translation

**Steps:**

1. Fork the repository
2. Edit `docs/<lang>/llm-prompt.md`
3. Add or update: translation rules, glossary terms, or grammar preferences
4. Submit a PR

### Option 2: Direct File Edits

For specific fixes, edit files directly in `docs/<lang>/docs/`:

1. Fork the repository
2. Edit the markdown file
3. Submit a PR

**Warning**: Direct edits may be overwritten when English source changes trigger re-translation. For permanent fixes, update the prompt file instead.

### Option 3: Report Issues

If you find translation errors but can't fix them:

1. Open an issue
2. Specify: language, file path, current text, suggested fix
3. A maintainer will update the translation or prompt

---

## How to Add a Missing Course for a Language

If a language exists but is missing a course (e.g., Portuguese exists but `nf4_science/` isn't translated):

### Using the Translation Script

```bash
cd _scripts

# Check what's missing
uv run python translate.py list-missing pt

# Translate specific files (one at a time recommended)
uv run python translate.py translate-page -l pt -p nf4_science/index.md
uv run python translate.py translate-page -l pt -p nf4_science/01_rnaseq.md
# ... continue for each file

# Or translate all missing at once (may timeout for large batches)
uv run python translate.py add-missing -l pt --include nf4_science
```

### Using GitHub Actions

1. Go to **Actions** → **Translate** → **Run workflow**
2. Select language: `pt`
3. Select command: `add-missing`
4. The workflow will translate all missing files and create a PR

### Manual Translation

For high-quality translations or languages where AI struggles:

1. Copy English source: `cp docs/en/docs/course/file.md docs/<lang>/docs/course/file.md`
2. Translate the content manually
3. Follow the glossary in `docs/<lang>/llm-prompt.md`
4. Submit a PR

---

## How to Add a New Language

### Step 1: Create Language Structure

```bash
cd _scripts
uv run python docs.py new-lang <lang-code>
```

This creates:

- `docs/<lang>/mkdocs.yml` - MkDocs config (inherits from English)
- `docs/<lang>/llm-prompt.md` - Translation prompt template
- `docs/<lang>/docs/.gitkeep` - Placeholder for content

### Step 2: Customize the Translation Prompt

Edit `docs/<lang>/llm-prompt.md` to add:

1. **Grammar preferences**: Formal/informal tone, regional variants
2. **Translation context rules**: When to translate vs. keep English
3. **Glossary**: Language-specific term translations
4. **Admonition titles**: Translations for Note, Tip, Warning, etc.

Example structure:

```markdown
# Translation Rules for <Language>

## Grammar Preferences

- Use formal/informal tone
- Regional spelling conventions

## Translation Context Rules

- In code blocks: Keep ALL Nextflow syntax in English
- In prose: Follow glossary for translations

## Glossary

### Terms to Keep in English

- Nextflow, Docker, GitHub, etc.

### Terms to Translate

| English | <Language>    |
| ------- | ------------- |
| channel | <translation> |
| process | <translation> |

### Admonition Titles

| English | <Language>    |
| ------- | ------------- |
| Note    | <translation> |
```

### Step 3: Register the Language

Add the language to:

1. `_scripts/docs.py` - Add to `SUPPORTED_LANGS` list
2. `docs/en/mkdocs.yml` - Add to `extra.alternate` for language switcher
3. `.github/workflows/translate.yml` - Add to language options

### Step 4: Translate Initial Content

```bash
cd _scripts

# Translate hello_nextflow (recommended starting point)
uv run python translate.py add-missing -l <lang> --include hello_nextflow

# Also translate supporting files
uv run python translate.py translate-page -l <lang> -p index.md
uv run python translate.py translate-page -l <lang> -p help.md
```

### Step 5: Submit PR

1. Create a branch: `lang-<code>` (e.g., `lang-de` for German)
2. Commit all changes
3. Create PR targeting the `lang` branch (if it exists) or `master`
4. Request review from a native speaker if possible

---

## Directory Structure

```
docs/
├── en/                     # English (source language)
│   ├── mkdocs.yml          # Main config (inherited by other languages)
│   ├── overrides/          # Theme overrides
│   ├── hooks/              # MkDocs hooks
│   └── docs/               # English content
│       ├── hello_nextflow/ # Course content
│       ├── envsetup/       # Setup guides
│       └── ...
├── pt/                     # Portuguese
│   ├── mkdocs.yml          # Inherits from ../en/mkdocs.yml
│   ├── llm-prompt.md       # Translation rules & glossary
│   └── docs/               # Portuguese translations
├── es/                     # Spanish
│   └── ...
└── ...
```

## Supported Languages

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

---

## Translation Script Reference

All commands run from `_scripts/` directory:

```bash
cd _scripts
```

### List Commands (Read-only)

```bash
# List files that need translation
uv run python translate.py list-missing <lang>

# List translations older than English source
uv run python translate.py list-outdated <lang>

# List translated files with no English source (orphans)
uv run python translate.py list-removable <lang>
```

### Translation Commands (Require ANTHROPIC_API_KEY)

```bash
# Translate a single file
uv run python translate.py translate-page -l <lang> -p <path-relative-to-en-docs>

# Translate all missing files
uv run python translate.py add-missing -l <lang>

# Translate missing files matching pattern
uv run python translate.py add-missing -l <lang> --include hello_nextflow

# Update outdated translations (smart diff)
uv run python translate.py update-outdated -l <lang>

# Remove orphaned translations
uv run python translate.py remove-removable -l <lang>
```

### Preview Commands

```bash
# Serve docs locally
uv run python docs.py serve <lang>

# Build docs
uv run python docs.py build-lang <lang>
```

---

## Translation Guidelines

### Code Blocks

Only **comments** are translated. Code stays in English:

```groovy
// Este comentário é traduzido
Channel.fromPath('data/*.fastq')  // Código permanece igual
    .set { reads_ch }
```

### Admonitions

Titles are translated, keywords stay in English:

```markdown
!!! note "Nota"

    Conteúdo traduzido aqui.
```

### Headings

Text is translated, anchors are preserved:

```markdown
## 1. Primeiros Passos { #getting-started }
```

### Links

Link text is translated, URLs and anchors stay unchanged:

```markdown
[Texto traduzido](../unchanged/path.md#unchanged-anchor)
```

---

## References

- [FastAPI Translation System](https://github.com/fastapi/fastapi/tree/master/scripts) - Inspiration for this implementation
- [Portuguese Glossary (Google Sheets)](https://docs.google.com/spreadsheets/d/1HUa3BO2kwukhX4EXQ-1blXeP5iueUdM23OwDRpfarDg/edit)
