# Translation Rules for Italian

The target language for this translation is **Italian** (`it`).

## Grammar Preferences

- Use first person plural inclusive (noi/voi) rather than formal Lei - this creates a more collaborative, educational tone
- Example: "ci addentriamo", "apriamo", "eseguiamo", "vediamo insieme"
- Follow standard Italian spelling conventions
- Prefer active voice when possible
- Use natural Italian expressions (e.g., "Voilà!" is acceptable as it's commonly used in Italian)

## Translation Context Rules

**Important distinction**: Some technical terms have different translation rules depending on context:

1. **In code blocks**: Keep ALL Nextflow syntax in English (the code must run)
2. **In code comments**: TRANSLATE comments to Italian (they are not executable)
3. **In prose/explanatory text**: Follow the glossary below for translations

For example:

- In prose: "Il canale di input riceve i file..." (translate "channel" to "canale")
- In code: `channel.fromPath('*.fastq')` (keep "channel" in English)
- In comments: `// emit a greeting` → `// emette un saluto`

## Code Comments

**Always translate code comments to Italian.** Comments are not executable code and should be in the target language for better comprehension.

```groovy
// English original
params.greeting = "Hello" // set default greeting

// Italian translation
params.greeting = "Hello" // imposta il saluto predefinito
```

## Glossary

### Terms to Keep in English (DO NOT TRANSLATE)

Note: Italian tech writing commonly uses English terms. Many of these terms are kept in English even in prose because they are standard in Italian technical documentation.

#### Nextflow Core Concepts

- Nextflow
- DSL2
- pipeline / pipelines
- script
- shell
- exec
- emit
- take
- main
- params
- output (commonly kept in English in Italian tech writing)
- input (can also use "ingresso", but English is common)
- file (kept in English, not "archivio")
- directory (kept in English, not "cartella")

#### Channel Types

- queue channel
- value channel

#### Data Types

- val
- path
- env
- stdin
- stdout
- tuple
- file (when referring to Nextflow type)

#### Operators

- map
- filter
- collect
- flatten
- groupTuple
- join
- combine
- mix
- merge
- view
- first
- last
- take (operator)
- branch
- multiMap
- splitCsv
- splitFastq
- splitFasta
- splitText

#### Directives

- publishDir
- container
- conda
- memory
- cpus
- time
- errorStrategy
- maxRetries
- maxErrors
- queue
- scratch
- storeDir
- tag
- label
- cache
- executor

#### Execution Concepts

- resume
- cache / caching
- work directory
- staging
- unstaging
- executor / executors
- job

#### Tools & Platforms

- Nextflow
- nf-core
- Docker
- Singularity
- Apptainer
- Conda
- Mamba
- GitHub
- GitLab
- Bitbucket
- Gitpod
- GitHub Codespaces
- Seqera Platform
- Wave
- Fusion

#### File Formats & Extensions

- `.nf`
- `nextflow.config`
- `main.nf`
- `modules.nf`
- `workflows.nf`
- FASTQ
- FASTA
- BAM
- SAM
- VCF
- BED
- GFF
- GTF
- CSV
- TSV
- JSON
- YAML

#### Programming Concepts

- closure
- Groovy
- shebang
- glob / globbing
- regex
- string
- boolean
- integer
- list
- map (data structure)
- set

### Terms to Translate

| English         | Italian                          | Notes                                 |
| --------------- | -------------------------------- | ------------------------------------- |
| workflow        | flusso di lavoro                 | Always translate in prose             |
| process         | processo                         | Always translate in prose             |
| channel         | canale / canali                  | Translate in prose                    |
| channel factory | fabbrica di canali               |                                       |
| queue channel   | canale di coda                   |                                       |
| work directory  | directory di lavoro              | Note: "directory" stays English       |
| alignment       | allineamento                     |                                       |
| command         | comando                          |                                       |
| container       | container                        | Keep in English (standard in Italian) |
| directive       | direttiva                        |                                       |
| directory       | directory                        | Keep in English (standard in Italian) |
| environment     | ambiente                         |                                       |
| file (general)  | file                             | Keep in English (standard in Italian) |
| index           | indice                           |                                       |
| input           | input / ingresso                 | English is common, both acceptable    |
| module          | modulo                           |                                       |
| operator        | operatore                        |                                       |
| output          | output                           | Keep in English (standard in Italian) |
| parameter       | parametro                        |                                       |
| qualifier       | qualificatore                    |                                       |
| reference       | riferimento                      |                                       |
| run             | eseguire / esecuzione / lanciare |                                       |
| sample          | campione                         |                                       |
| task            | attività                         |                                       |
| training        | formazione                       |                                       |
| tuple           | tupla                            |                                       |
| greeting        | saluto                           | In context of Hello World examples    |
| log             | registri / log                   | Both acceptable                       |

### Admonition Titles

| English  | Italian             |
| -------- | ------------------- |
| Note     | Nota                |
| Tip      | Suggerimento        |
| Warning  | Avviso / Attenzione |
| Exercise | Esercizio           |
| Solution | Soluzione           |
| Example  | Esempio             |

### Section Headers

These recurring section headers should be translated consistently:

| English            | Italian            |
| ------------------ | ------------------ | --------------------------------------------------------- |
| Takeaway           | Takeaway           | (Keep in English - common in Italian educational content) |
| What's next?       | Cosa c'è dopo?     |
| Warmup             | Riscaldamento      |
| Directory contents | Directory contents | (Keep in code block titles)                               |
| Output             | Output             | (Keep in code block titles)                               |

### Tab Labels

| English | Italian |
| ------- | ------- |
| After   | Dopo    |
| Before  | Prima   |
| Gitpod  | Gitpod  |
| Local   | Locale  |

### Common Expressions

| English            | Italian                        |
| ------------------ | ------------------------------ |
| Congratulations!   | Congratulazioni!               |
| Good news          | Buone notizie                  |
| Voilà!             | Voilà! (acceptable in Italian) |
| Take a short break | Prendetevi una piccola pausa   |
| you've earned it   | ve la siete meritata           |
