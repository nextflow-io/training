# Translation Instructions

Translate English Markdown to the target language. Content is wrapped in `%%%` markers (don't include these in output).

## Rules

### Keep in English (never translate)

- **Nextflow syntax**: `channel`, `process`, `workflow`, `emit`, `take`, `main`, `output`, `input`, `script`, `shell`, `exec`, `params`
- **Directives**: `publishDir`, `container`, `conda`, `memory`, `cpus`, `time`, `errorStrategy`, `tag`, `label`, `cache`, `executor`
- **Operators**: `map`, `filter`, `collect`, `flatten`, `groupTuple`, `join`, `combine`, `mix`, `view`, `splitCsv`, `splitFastq`
- **Types**: `val`, `path`, `env`, `stdin`, `stdout`, `tuple`, `queue channel`, `value channel`
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
