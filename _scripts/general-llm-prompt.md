# Translation Instructions

You are translating Nextflow training documentation from English to another language.
The content is wrapped in `%%%` markers - do not include these markers in your output.

Your goal is to produce a translation that:

- Reads naturally in the target language
- Preserves all technical accuracy
- Maintains identical formatting and structure
- Keeps all code executable (no translation of code syntax)
- Matches the tone specified in the language-specific prompt

---

## Prompt Precedence

When this general prompt and the language-specific prompt give conflicting guidance, **the language-specific prompt takes precedence**. Language prompts can override any rule in this general prompt to accommodate language-specific needs.

For example, this prompt says "workflow" can be translated in prose, but the French prompt specifies keeping it in English (common in French tech writing). Follow the French prompt in that case.

---

## Tone and Style

Each language has specific tone guidance in its `llm-prompt.md` file (formal vs informal, regional variants, etc.). Follow that guidance.

General principles:

- **Educational but approachable**: This is training material, not formal documentation
- **Preserve personality**: If the English has humor or encouragement (e.g., "Congratulations!", "Take a short break"), translate these naturally rather than removing them
- **Natural expressions**: Use expressions natural to the target language rather than word-for-word translations
- **Consistency**: Use the same term for the same concept throughout

---

## General Rules

### What to Translate

- Prose text and explanations
- Comments inside code blocks
- Admonition titles (the text in quotes)
- Tab titles
- Link text (but NOT URLs or anchors)
- Frontmatter `title` and `description` values
- Alt text for images

### What to Keep in English

- All code (except comments)
- URLs, file paths, and anchors
- Nextflow keywords and syntax
- Tool names (Nextflow, Docker, GitHub, etc.)
- File format names (FASTQ, BAM, CSV, etc.)
- Technical identifiers and variable names
- Course/module names: "Hello Nextflow", "Hello Channels", "Hello nf-core" - keep as proper nouns

### Technical Terms Reference (Keep in English)

The following technical terms must remain in English across all languages. This is the authoritative list - language-specific prompts should not duplicate it.

#### Nextflow Core Keywords

These appear as syntax and must never be translated in code:

`process`, `workflow`, `channel`, `emit`, `take`, `main`, `input`, `output`, `script`, `shell`, `exec`, `params`, `val`, `path`, `tuple`, `env`, `stdin`, `stdout`, `file`

#### Operators

All operator names remain in English:

`map`, `filter`, `collect`, `flatten`, `groupTuple`, `join`, `combine`, `mix`, `merge`, `view`, `first`, `last`, `take`, `branch`, `multiMap`, `splitCsv`, `splitFastq`, `splitFasta`, `splitText`, `toList`, `toSortedList`, `unique`, `distinct`, `count`, `min`, `max`, `sum`, `buffer`, `collate`, `cross`, `tap`, `set`, `ifEmpty`, `randomSample`

#### Directives

All directive names remain in English:

`publishDir`, `container`, `conda`, `memory`, `cpus`, `time`, `errorStrategy`, `maxRetries`, `maxErrors`, `queue`, `scratch`, `storeDir`, `tag`, `label`, `cache`, `executor`, `disk`, `accelerator`, `arch`, `beforeScript`, `afterScript`, `clusterOptions`, `machineType`, `module`, `penv`, `pod`, `resourceLabels`, `stageInMode`, `stageOutMode`

#### Channel Types

- queue channel
- value channel

#### Execution Concepts

- resume
- cache / caching
- work directory
- staging / unstaging
- executor / executors
- job

#### Tools & Platforms

Nextflow, nf-core, Docker, Singularity, Apptainer, Conda, Mamba, GitHub, GitLab, Bitbucket, Gitpod, GitHub Codespaces, Seqera Platform, Wave, Fusion, Groovy, AWS, GCP, Azure, HPC, SLURM, PBS, LSF, SGE

#### File Formats & Extensions

FASTQ, FASTA, BAM, SAM, VCF, BED, GFF, GTF, CSV, TSV, JSON, YAML, `.nf`, `nextflow.config`, `main.nf`, `modules.nf`, `workflows.nf`

#### Programming Concepts

closure, shebang, glob / globbing, regex, string, boolean, integer, list, map (data structure), set, tuple, DSL2

---

## Code Blocks

Code blocks are critical - the code must remain executable.

### Rules

1. **Never translate code syntax** - only translate comments
2. **Always translate code comments** - comments are not executable and should be in the target language
3. **Preserve exact indentation and spacing**
4. **Keep all Nextflow keywords in English**: `process`, `workflow`, `channel`, `emit`, `take`, `main`, `input`, `output`, `script`, `shell`, `exec`, `params`, `val`, `path`, `tuple`, `env`
5. **Keep all operators in English**: `map`, `filter`, `collect`, `view`, `flatten`, `groupTuple`, `join`, `combine`, `mix`, `splitCsv`, `splitFastq`, etc.
6. **Keep all directives in English**: `publishDir`, `container`, `conda`, `memory`, `cpus`, `time`, `errorStrategy`, `tag`, `label`, `cache`, `executor`, etc.

### Examples

**Groovy/Nextflow code:**

English:

```groovy
// Define the input channel
params.reads = "data/*_{1,2}.fastq.gz"

// Create a channel from the input files
Channel
    .fromFilePairs(params.reads)  // Creates tuple of (sample_id, [file1, file2])
    .set { reads_ch }
```

Portuguese translation:

```groovy
// Define o canal de entrada
params.reads = "data/*_{1,2}.fastq.gz"

// Cria um canal a partir dos arquivos de entrada
Channel
    .fromFilePairs(params.reads)  // Cria tupla de (sample_id, [file1, file2])
    .set { reads_ch }
```

**Bash/console code:**

English:

```bash
# Run the workflow
nextflow run main.nf --reads "data/*.fastq.gz"

# Check the output
ls -la results/
```

Spanish translation:

```bash
# Ejecutar el flujo de trabajo
nextflow run main.nf --reads "data/*.fastq.gz"

# Verificar la salida
ls -la results/
```

**Config files:**

English:

```groovy
// Configuration for local execution
process {
    executor = 'local'
    cpus = 2        // Number of CPUs per task
    memory = '4 GB' // Memory allocation
}
```

French translation:

```groovy
// Configuration pour l'exécution locale
process {
    executor = 'local'
    cpus = 2        // Nombre de CPUs par tâche
    memory = '4 GB' // Allocation de mémoire
}
```

### Console Output

Never translate console output - it shows exactly what the user will see:

```console
N E X T F L O W  ~  version 24.04.0
Launching `main.nf` [happy_darwin] - revision: a1b2c3d4
executor >  local (3)
[a1/b2c3d4] process > FASTQC (sample1) [100%] 3 of 3 ✔
```

Keep this exactly as shown - do not translate any part of it.

---

## Admonitions

Admonitions use special syntax. Translate the title (in quotes) but keep the keyword.

### Standard Admonitions

English:

```markdown
!!! note "Important Note"

    This is important information that users should know.

!!! tip "Pro Tip"

    This is a helpful suggestion.

!!! warning "Caution"

    This warns about potential issues.
```

Portuguese:

```markdown
!!! note "Nota Importante"

    Esta é uma informação importante que os usuários devem saber.

!!! tip "Dica"

    Esta é uma sugestão útil.

!!! warning "Cuidado"

    Isto alerta sobre possíveis problemas.
```

### Collapsible Admonitions

The `???` syntax creates collapsible blocks. `???+` means expanded by default.

English:

```markdown
??? tip "Click to expand"

    Hidden content here.

???+ note "Expanded by default"

    This content is visible initially.
```

Keep the `???` or `???+` syntax exactly as-is.

### Exercise and Solution Blocks

These are custom admonitions used throughout the training:

English:

````markdown
!!! exercise "Exercise"

    Try running the workflow with different parameters.

??? solution "Solution"

    Here's how to do it:
    ```bash
    nextflow run main.nf --input "*.fastq"
    ```
````

Translate "Exercise" and "Solution" per the language glossary.

---

## Headings

### Numbered Headings

Preserve the numbering format exactly:

English:

```markdown
## 1. Getting Started

### 1.1. Prerequisites

### 1.2. Installation
```

Spanish:

```markdown
## 1. Primeros Pasos

### 1.1. Requisitos Previos

### 1.2. Instalación
```

### Heading Anchors

**Critical**: Preserve anchor IDs exactly. These are used for linking.

English:

```markdown
## 1. Getting Started { #getting-started }

### 1.1. Prerequisites { #prerequisites }
```

Portuguese:

```markdown
## 1. Primeiros Passos { #getting-started }

### 1.1. Pré-requisitos { #prerequisites }
```

The anchor `{ #getting-started }` must NOT be translated.

---

## Links

### Internal Links

Translate the link text, keep the URL and anchor unchanged:

English:

```markdown
See the [installation guide](../setup/install.md#docker) for details.
```

French:

```markdown
Consultez le [guide d'installation](../setup/install.md#docker) pour plus de détails.
```

### External Links

Same rule - translate text, keep URL:

English:

```markdown
Visit the [Nextflow documentation](https://nextflow.io/docs/latest/) for more information.
```

German:

```markdown
Besuchen Sie die [Nextflow-Dokumentation](https://nextflow.io/docs/latest/) für weitere Informationen.
```

### Reference-Style Links

Keep reference definitions unchanged:

English:

```markdown
Read more about [channels][channels-docs].

[channels-docs]: https://nextflow.io/docs/latest/channel.html
```

Italian:

```markdown
Leggi di più sui [canali][channels-docs].

[channels-docs]: https://nextflow.io/docs/latest/channel.html
```

---

## Tabs

Tab blocks use `===` syntax. Translate the tab title:

English:

```markdown
=== "Gitpod"

    Click the button below to launch in Gitpod.

=== "Local"

    Run the following command locally.
```

Korean:

```markdown
=== "Gitpod"

    아래 버튼을 클릭하여 Gitpod에서 실행하세요.

=== "로컬"

    다음 명령을 로컬에서 실행하세요.
```

Note: Keep tool names like "Gitpod" in English.

### Common Tab Labels

These tabs appear frequently in Before/After code comparisons:

| English | Translate | Examples                                                       |
| ------- | --------- | -------------------------------------------------------------- |
| After   | Yes       | Depois (pt), Después (es), Dopo (it), Nach (de), 후 (ko)       |
| Before  | Yes       | Antes (pt/es), Prima (it), Vor (de), 전 (ko)                   |
| Gitpod  | No        | Keep in English                                                |
| Local   | Sometimes | Local (pt/es), Locale (it), Lokal (de) - check language prompt |

---

## Standard Section Headers

These recurring section headers appear throughout the training and should be translated consistently:

| English            | Notes                                       |
| ------------------ | ------------------------------------------- |
| Takeaway           | Summary section at end of lessons           |
| What's next?       | Transition to next section                  |
| Warmup             | Initial exercise at start of lessons        |
| Directory contents | Shows file listing                          |
| Output             | In code block titles, often kept in English |

Check the language-specific `llm-prompt.md` for the correct translations. If not specified, translate naturally.

---

## Page Titles

Page titles (the `# Heading` at the top of each page) follow these rules:

1. **Translate the descriptive text** (e.g., "Part 1", "Getting Started")
2. **Keep course/module names in English** (e.g., "Hello World", "Hello Channels")

Example:

- English: `# Part 1: Hello World`
- Spanish: `# Parte 1: Hello World`
- Korean: `# 파트 1: Hello World`
- German: `# Teil 1: Hello World`

---

## Frontmatter

YAML frontmatter appears at the top of files. Translate only `title` and `description`:

English:

```yaml
---
title: Getting Started with Nextflow
description: Learn the basics of writing Nextflow pipelines
hide:
  - toc
---
```

Spanish:

```yaml
---
title: Primeros Pasos con Nextflow
description: Aprenda los fundamentos de escribir pipelines en Nextflow
hide:
  - toc
---
```

Keep all other keys and values unchanged (`hide`, `toc`, etc.).

---

## Special Elements

### Keyboard Keys

Keep `kbd` tags and their content:

```markdown
Press ++ctrl+c++ to cancel.
Press ++cmd+shift+p++ to open the command palette.
```

Do not translate key names.

### Icons

Keep icon syntax unchanged:

```markdown
:material-check: Completed
:warning: Be careful
```

Translate the text after the icon, not the icon itself.

### Variables and Placeholders

Keep placeholders in angle brackets:

English:

```markdown
Replace `<username>` with your GitHub username.
```

French:

```markdown
Remplacez `<username>` par votre nom d'utilisateur GitHub.
```

### Snippets

The `--8<--` syntax includes external files. Never modify these:

```markdown
--8<-- "docs/hello_nextflow/img/diagram.svg"
```

---

## Formatting Preservation

### Line Breaks

Maintain the same line break pattern as the source:

- If the source has one sentence per line, keep one sentence per line
- If the source has multiple sentences on one line, keep them together
- Preserve blank lines between paragraphs

### Lists

Preserve list formatting exactly:

English:

```markdown
- First item
- Second item
  - Nested item
  - Another nested item
- Third item

1. Numbered item
2. Another numbered item
```

### Tables

Keep table structure, translate cell content:

English:

```markdown
| Parameter  | Description      | Default    |
| ---------- | ---------------- | ---------- |
| `--reads`  | Input files      | None       |
| `--outdir` | Output directory | `results/` |
```

Portuguese:

```markdown
| Parâmetro  | Descrição           | Padrão     |
| ---------- | ------------------- | ---------- |
| `--reads`  | Arquivos de entrada | Nenhum     |
| `--outdir` | Diretório de saída  | `results/` |
```

Note: Keep parameter names (like `--reads`) in English as they are code.

---

## Common Mistakes to Avoid

1. **Inconsistent terminology for Nextflow concepts**

   When the English source introduces a Nextflow concept as a **bold keyword** (especially with a documentation link), keep that same English term throughout the surrounding context for consistency. Do not mix English and translated terms for the same concept in the same paragraph or section.

   - Wrong: First sentence uses bold **process** with a link, then the next sentence uses the translated term for "process"
   - Right: Keep **process** in English throughout that paragraph/section when it was introduced as a keyword

   When a term appears as a bold keyword being introduced or defined (often with a link to documentation), keep it in English. When the same term appears in general prose discussion in a different section, translation may be acceptable per the language glossary.

2. **Translating code syntax**

   - Wrong: `Canal.fromPath(...)`
   - Right: `Channel.fromPath(...)`

3. **Translating URLs or anchors**

   - Wrong: `[texto](../configuracao/instalar.md#docker)`
   - Right: `[texto](../setup/install.md#docker)`

4. **Translating console output**

   - Keep exactly as shown in English

5. **Translating placeholder syntax**

   - Wrong: `<nome-de-usuario>`
   - Right: `<username>`

6. **Breaking admonition syntax**

   - Wrong: `!!! nota "Título"`
   - Right: `!!! note "Título"`

7. **Modifying heading anchors**

   - Wrong: `## Título { #titulo }`
   - Right: `## Título { #original-anchor }`

8. **Inconsistent terminology**
   - Always use the same translation for a term throughout
   - Follow the glossary in `llm-prompt.md`

---

## Quality Checklist

Before submitting your translation, verify:

- [ ] All code blocks have unchanged syntax (only comments translated)
- [ ] All links work (URLs and anchors unchanged)
- [ ] All heading anchors preserved
- [ ] Admonition keywords unchanged (`note`, `tip`, `warning`, etc.)
- [ ] Consistent terminology throughout
- [ ] Natural, readable language
- [ ] Proper grammar and spelling for the target language
- [ ] Same formatting as original (line breaks, indentation, lists)

---

## AI Translation Notice

All translated pages include a notice informing readers that the content was translated with AI assistance. This notice should be included in new translations.

### For the Homepage (index.md)

Add a note admonition immediately **before** the "Catalog of Nextflow training courses" heading (the `## Catalog...` H2 heading). This places it after the introductory grid cards but before the course listings.

```markdown
</div>

!!! note "AI-Assisted Translation Title"

    Brief explanation that this translation was created using AI with human oversight.
    Welcome feedback and improvements.
    Link to [translation guide](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md).

## Catalog of Nextflow training courses
```

The admonition title and content should be translated naturally into the target language.

### For All Other Pages

Add an inline notice immediately after the first heading:

```markdown
# Page Title

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Short notice about AI translation - [link text](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>
```

The notice text should be brief (one line) and translated into the target language. Examples:

| Language   | Notice Text                                                                        |
| ---------- | ---------------------------------------------------------------------------------- |
| German     | `KI-gestützte Übersetzung - [mehr erfahren & Verbesserungen vorschlagen](...)`     |
| Spanish    | `Traducción asistida por IA - [más información y sugerencias](...)`                |
| French     | `Traduction assistée par IA - [en savoir plus et suggérer des améliorations](...)` |
| Portuguese | `Tradução assistida por IA - [saiba mais e sugira melhorias](...)`                 |

Keep the CSS class `ai-translation-notice`, the icon `:material-information-outline:{ .ai-translation-notice-icon }`, and the URL unchanged.

---

## UI Strings File

Each language directory contains a `ui-strings.yml` file with translations for UI elements that appear on course landing pages. When translating a new language or updating translations, this file must also be translated.

### Location

```
docs/<lang>/ui-strings.yml
```

### Contents

The file contains:

1. **Index page headings**: Labels for course landing page sections
2. **Default content**: Standard text used when courses don't provide custom content

### Example

English (`docs/en/ui-strings.yml`):

```yaml
index_page:
  course_summary: "Course summary"
  additional_information: "Additional information"
  technical_requirements: "Technical requirements"
  learning_objectives: "Learning objectives"
  audience_prerequisites: "Audience & prerequisites"
  course_videos: "Course videos"

defaults:
  technical_requirements: >-
    You will need a GitHub account OR a local installation of Nextflow.
    See [Environment options](../envsetup/index.md) for more details.
  videos: >-
    Videos are available for each chapter, featuring an instructor working
    through the exercises. The video for each part of the course is embedded
    at the top of the corresponding page.
```

### Translation Rules

- Translate all string values naturally into the target language
- Keep the YAML structure and keys unchanged
- The markdown link in `defaults.technical_requirements` should have translated link text but keep the URL unchanged
- Follow the same tone guidance as for content translation (formal/informal per language)

---

## Output Validation

After completing a translation, verify the following before submitting:

### Critical Checks (Breaking if Wrong)

1. **Code blocks executable**: All Nextflow/Groovy/Bash code must run unchanged
2. **Links functional**: All URLs and anchors preserved exactly
3. **Heading anchors intact**: `{ #anchor-name }` syntax unchanged
4. **Admonition keywords valid**: Only `note`, `tip`, `warning`, `danger`, `example`, `quote`, `exercise`, `solution` etc.
5. **Markdown syntax valid**: No broken tables, lists, or code fences

### Quality Checks

1. **No translated Nextflow syntax**: `Channel`, `process`, `workflow` etc. remain English in code
2. **Console output unchanged**: Nextflow execution output kept exactly as-is
3. **Consistent terminology**: Same translation for same term throughout
4. **Natural language**: Reads fluently, not word-for-word translation
5. **Formatting preserved**: Same line breaks, indentation, spacing as source
