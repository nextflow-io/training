# Translation Instructions

Translate English Markdown to the target language. Content is wrapped in `%%%` markers (don't include these in output).

## Rules

### Context-Dependent Translation

**Important**: Some technical terms have different translation rules depending on context:

1. **In code blocks**: Keep ALL Nextflow syntax in English (code must be executable)
2. **In prose text**: Follow the language-specific `llm-prompt.md` glossary for translations

For example, terms like "channel", "process", "workflow" may be translated in prose for some languages (e.g., Portuguese: "canal", "processo", "fluxo de trabalho") but must always remain in English in code blocks.

### Always Keep in English

- **Nextflow syntax** (in code): `channel`, `process`, `workflow`, `emit`, `take`, `main`, `output`, `input`, `script`, `shell`, `exec`, `params`
- **Directives**: `publishDir`, `container`, `conda`, `memory`, `cpus`, `time`, `errorStrategy`, `tag`, `label`, `cache`, `executor`
- **Operators**: `map`, `filter`, `collect`, `flatten`, `groupTuple`, `join`, `combine`, `mix`, `view`, `splitCsv`, `splitFastq`
- **Types** (in code): `val`, `path`, `env`, `stdin`, `stdout`, `tuple`
- **Concepts**: `DSL2`, `resume`, `cache`, `work directory`, `staging`
- **Tools**: `Nextflow`, `nf-core`, `Docker`, `Singularity`, `Conda`, `GitHub`, `Gitpod`, `Seqera Platform`
- **Formats**: `FASTQ`, `FASTA`, `BAM`, `VCF`, `CSV`, `JSON`, `YAML`

### Code blocks

Translate **only comments**, not code:

```groovy
// Translate this comment
Channel.fromPath('data/*.fastq')  // Keep code as-is
```

### Admonitions

Translate the title in quotes, keep the keyword:

```
!!! note "Translated Title"
    Translated content.
```

Same for collapsible: `??? tip "Title"` and `???+ warning "Title"`

### Headings

Translate text, preserve anchors and numbers:

```
## 1. Translated Title { #original-anchor }
```

### Links

Translate link text only. Never translate URLs or anchors:

```
[Translated text](../unchanged/path.md#unchanged-anchor)
```

### Tabs

Translate tab titles:

```
=== "Translated Tab"
    Translated content.
```

### Frontmatter

Translate `title` and `description` values. Keep keys and technical values unchanged.

### Formatting

- Preserve whitespace, line breaks, and indentation exactly
- Keep one sentence per line if the original does
