---
name: Check Inline Code Formatting
description: Review inline code in .mdx pages for proper formatting and syntax highlighting. Check for missing backticks on CLI flags, filenames, and code elements. Verify that inline syntax highlighting ({:nextflow} suffix) is used appropriately. Use when reviewing documentation or quiz questions.
---

# Check Inline Code Formatting

Review inline code in `.mdx` pages for proper formatting and appropriate use of syntax highlighting.

See [../shared/repo-conventions.md](../shared/repo-conventions.md) for file conventions.

## What to Check

### 1. Missing Backticks

The following should always be wrapped in backticks:

**CLI flags and options:**

- `-resume`, `-profile`, `-params-file`, `--input`, etc.

**Filenames and paths:**

- `nextflow.config`, `.command.sh`, `work/`, `modules/`, etc.

**Commands:**

- `nextflow run`, `nextflow log`, `docker run`, etc.

**Code elements:**

- Variable names: `params.input`, `outputDir`
- Process names: `SAYHELLO`
- Operators and methods: `collect()`, `map()`, `view()`
- Directives: `container`, `publishDir`, `memory`

### 2. Inline Syntax Highlighting

Add a `{:lang}` suffix inside the backticks to give inline code syntax colours: `` `code{:nextflow}` ``.
Use `{:nextflow}` for Nextflow code, `{:groovy}` only for real Groovy (plugin classes, `build.gradle`), and `{:bash}` for shell.

The suffix must sit inside the backticks. `` `code`{:nextflow} `` is a JavaScript expression in MDX and fails the build.

**DO use `{:nextflow}` for:**

- Variable interpolation: `` `${variable}{:nextflow}` ``
- Closures: `` `{ meta.id }{:nextflow}` `` or `` `.map { it * 2 }{:nextflow}` ``
- Method chains with closures: `` `channel.of(1,2,3).map { it * 2 }{:nextflow}` ``
- String interpolation: `` `"${greeting}-output.txt"{:nextflow}` ``
- Type declarations: `` `param: String = 'value'{:nextflow}` ``

**DO NOT use `{:nextflow}` for:**

- Single words or identifiers: use `` `outputDir` `` not `` `outputDir{:nextflow}` ``
- Simple directives: use `` `container 'uri'` `` not `` `container 'uri'{:nextflow}` ``
- Plain dot notation: use `` `processName.out` `` not `` `processName.out{:nextflow}` ``
- Keywords alone: use `` `include` `` not `` `include{:nextflow}` ``

**The rule:** Only use `{:nextflow}` when there's meaningful syntax structure (braces, interpolation, operators) that highlighting will visually enhance. A single word or simple expression gains nothing from syntax highlighting.

### 3. Consistency

- Similar code references should be formatted the same way throughout a document
- Quiz answer options should have consistent formatting

### 4. Leftover MkDocs syntax

`` `#!groovy code` `` is the old MkDocs form and renders literally. Replace it with `` `code{:nextflow}` `` (or `{:groovy}` for real Groovy), following the rules above.

## Examples

### Good inline code formatting

```mdx
Run the workflow with `nextflow run main.nf -resume`.

The `${greeting}{:nextflow}` variable is interpolated into the string.

Set memory with `process { withName: 'FOO' { memory = '4.GB' } }{:nextflow}`.

The `container` directive specifies the Docker image.

Access outputs with `processName.out.outputName`.
```

### Bad inline code formatting

```mdx
{/* Missing backticks */}

Run the workflow with nextflow run main.nf -resume.

{/* Unnecessary {:nextflow} on single word */}

The `container{:nextflow}` directive specifies the image.

{/* Missing {:nextflow} on closure */}

Transform with `.map { it.toUpperCase() }`.

{/* Suffix outside the backticks: fails the build */}

Transform with `.map { it.toUpperCase() }`{:nextflow}.
```
