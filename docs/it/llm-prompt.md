# Translation Rules for Italian

The target language for this translation is **Italian** (`it`).

## Grammar Preferences

- Use formal tone (Lei instead of tu)
- Follow standard Italian spelling conventions
- Prefer active voice when possible

## Translation Context Rules

**Important distinction**: Some technical terms have different translation rules depending on context:

1. **In code blocks**: Keep ALL Nextflow syntax in English (the code must run)
2. **In prose/explanatory text**: Follow the glossary below for translations

For example:

- In prose: "Il canale di input riceve i file..." (translate "channel" to "canale")
- In code: `channel.fromPath('*.fastq')` (keep "channel" in English)

## Glossary

### Terms to Keep in English (DO NOT TRANSLATE)

#### Nextflow Core Concepts

- Nextflow
- DSL2
- channel / channels
- process / processes
- workflow / workflows
- pipeline / pipelines
- script
- shell
- exec
- emit
- take
- main
- params

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

| English        | Italian               |
| -------------- | --------------------- |
| alignment      | allineamento          |
| command        | comando               |
| container      | container             |
| directive      | direttiva             |
| directory      | directory             |
| environment    | ambiente              |
| file (general) | file                  |
| index          | indice                |
| input          | input                 |
| module         | modulo                |
| operator       | operatore             |
| output         | output                |
| parameter      | parametro             |
| reference      | riferimento           |
| run            | eseguire / esecuzione |
| sample         | campione              |
| task           | attività              |
| training       | formazione            |
| tuple          | tupla                 |

### Admonition Titles

| English  | Italian      |
| -------- | ------------ |
| Note     | Nota         |
| Tip      | Suggerimento |
| Warning  | Avviso       |
| Exercise | Esercizio    |
| Solution | Soluzione    |
| Example  | Esempio      |
