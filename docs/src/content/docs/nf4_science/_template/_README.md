# nf4_science Course Template

This is a generic template for creating a new domain-specific course under `nf4_science/`.
It is based on the common patterns found in the Genomics and RNAseq courses.

## How to use

1. Copy this directory and the corresponding script directory to create your course:

   ```bash
   # Docs
   cp -r docs/src/content/docs/nf4_science/_template docs/src/content/docs/nf4_science/{your_domain}
   rm docs/src/content/docs/nf4_science/{your_domain}/_README.md

   # Scripts
   cp -r nf4-science/_template nf4-science/{your_domain}
   ```

2. Rename the workflow script:

   ```bash
   mv nf4-science/{your_domain}/_template.nf nf4-science/{your_domain}/{your_domain}.nf
   ```

3. Rename the module files to match your tools:

   ```bash
   mv nf4-science/{your_domain}/modules/_tool_a.nf nf4-science/{your_domain}/modules/{tool_a}.nf
   mv nf4-science/{your_domain}/modules/_tool_b.nf nf4-science/{your_domain}/modules/{tool_b}.nf
   ```

4. Search for all `{PLACEHOLDER}` values in the docs and scripts and replace them with your content.
   In page prose, MDX needs the braces escaped, so placeholders there are written `\{PLACEHOLDER\}`: remove the backslashes along with the placeholder.
   In frontmatter, code spans and code blocks they are written `{PLACEHOLDER}`.

5. Add the course to `docs/src/sidebar.mjs` (see snippet below).

6. Populate the `data/` directory with test datasets.

7. Build out the `solutions/` directory with working code for each part.

8. Run `cd docs && npm run dev` and check every page renders, then run a full build (`npm run build`) to check the links.

## sidebar.mjs snippet

Add a group to the `items` of the "Nextflow for Science" group, next to Genomics and RNAseq.
Entries are page slugs without the `.mdx` extension; their labels come from the page titles.

```js
{
  "label": "{Your Domain}",
  "collapsed": true,
  "items": [
    { "label": "Overview", "slug": "nf4_science/{your_domain}" },
    "nf4_science/{your_domain}/00_orientation",
    "nf4_science/{your_domain}/01_method",
    "nf4_science/{your_domain}/02_single_sample",
    "nf4_science/{your_domain}/03_multi_sample",
    "nf4_science/{your_domain}/next_steps",
    "nf4_science/{your_domain}/survey"
  ]
}
```

Copy the `translations` map of the "Overview" label from the Genomics group.

## Placeholder reference

### Course-level placeholders

| Placeholder                  | Description                   | Example (Genomics)                         |
| ---------------------------- | ----------------------------- | ------------------------------------------ |
| `{DOMAIN}`                   | Domain name (title case)      | Genomics                                   |
| `{DOMAIN_DIR}`               | Directory name (lowercase)    | genomics                                   |
| `{METHOD}`                   | Analysis method name          | variant calling                            |
| `{METHOD_SHORT_DESCRIPTION}` | One-line method description   | variant calling with GATK                  |
| `{ACCESSORY_FILES}`          | Types of accessory files used | index files and reference genome resources |

### Tool placeholders

| Placeholder                           | Description                    | Example (Genomics)                                  |
| ------------------------------------- | ------------------------------ | --------------------------------------------------- |
| `{TOOL_A}` / `{TOOL_B}`               | Tool display names             | Samtools / GATK                                     |
| `{TOOL_A_URL}` / `{TOOL_B_URL}`       | Tool documentation URLs        | https://www.htslib.org/                             |
| `{TOOL_A_CONTAINER_URI}`              | Full container URI             | community.wave.seqera.io/library/samtools:1.20--... |
| `{TOOL_A_MODULE}` / `{TOOL_B_MODULE}` | Module filenames (without .nf) | samtools_index / gatk_haplotypecaller               |
| `{TOOL_A_PROCESS_NAME}`               | UPPERCASE process name         | SAMTOOLS_INDEX                                      |
| `{TOOL_A_COMMAND}`                    | The shell command to run       | samtools index '<input_bam>'                        |

### Input/Output placeholders

| Placeholder            | Description                 | Example (Genomics)   |
| ---------------------- | --------------------------- | -------------------- |
| `{PRIMARY_INPUT_TYPE}` | Primary input file type     | BAM file             |
| `{PRIMARY_PARAM_NAME}` | Nextflow param name         | reads_bam            |
| `{TEST_INPUT_PATH}`    | Path to test input in data/ | bam/reads_mother.bam |
| `{FINAL_OUTPUT_TYPE}`  | Final pipeline output type  | joint-called VCFs    |

### Content placeholders

| Placeholder             | Description                              |
| ----------------------- | ---------------------------------------- |
| `{PART2_SUMMARY}`       | One-line summary of Part 2 scope         |
| `{AGGREGATION_SUMMARY}` | One-line summary of the aggregation step |
| `{PROCESS_LIST}`        | Comma-separated list of process names    |
| `{TYPEFORM_ID}`         | Typeform embed ID for the survey         |

## Pedagogical structure

The template follows a three-part arc:

1. **Part 1 (01_method.mdx)**: Manual testing in Docker containers to understand the methodology
2. **Part 2 (02_single_sample.mdx)**: Wrap commands in Nextflow; single-sample, then batch
3. **Part 3 (03_multi_sample.mdx)**: Add multi-sample aggregation using channel operators

### Key conventions

- Use **`<CodeWalkthrough>`** with a `step="Before"` and a `step="After"` code block to show code changes, with `{...}` line highlights on the changed lines
- Put expected command output in the command's `bash` block, after a `# [!output]` line
- Use `<details class="abstract">` with `<summary>Directory contents</summary>` for collapsible directory trees
- End each major section with **Takeaway** and **What's next?** subsections
- Use `---` horizontal rules to separate major numbered sections
- Provide skeleton files (workflow + modules) for learners to fill in
- Organize solutions by part (`solutions/part2/`, `solutions/part3/`)

See `AGENTS.md` at the repository root for the full page syntax.

## Directory structure

```
docs/src/content/docs/nf4_science/{domain}/
├── index.mdx                   # Course overview with frontmatter
├── 00_orientation.mdx          # Environment setup
├── 01_method.mdx               # Manual testing in containers
├── 02_single_sample.mdx        # Single-sample Nextflow implementation
├── 03_multi_sample.mdx         # Multi-sample aggregation
├── survey.mdx                  # Typeform feedback survey
├── next_steps.mdx              # Course summary and suggestions
└── img/                        # Diagrams (.excalidraw.svg, .png)

nf4-science/{domain}/
├── {domain}.nf                 # Skeleton workflow file
├── nextflow.config             # Minimal config (docker.enabled = true)
├── data/                       # Test datasets and resources
│   └── samplesheet.csv         # Sample metadata
├── modules/                    # Skeleton module files
│   ├── {tool_a}.nf
│   └── {tool_b}.nf
└── solutions/
    ├── part2/                  # Complete Part 2 solution
    │   ├── {domain}-2.nf
    │   ├── nextflow.config
    │   └── modules/
    └── part3/                  # Complete Part 3 solution
        ├── {domain}-3.nf
        ├── nextflow.config
        └── modules/
```
