# Translation Rules for Polish

The target language for this translation is **Polish** (`pl`).

## Grammar Preferences

- Use informal tone (Ty for one person, Wy for multiple, not Pan/Pani); always capitalize second-person pronouns, including possessives (Twój, Wasz)
- Follow standard Polish spelling conventions
- Prefer active voice when possible
- Use natural Polish syntax and sentence structure, feel free to rephrase the English text wherever it improves clarity and naturality
- Always use correct declension, even with English terms. Use apostrophe endings (e.g. "workflow'u") where necessary, and in particular for terms in the list below.
- Use inanimate genders when referring to programming terms. Exception: use virile (animate masculine) to refer to Nextflow (i.e. Nextflow'a).

### Terms with apostrophe endings

- pipeline (pipeline'u, pipeline'owi, etc.)
- workflow (workflow'u, workflow'owi, etc.)

## Translation Context Rules

**Important distinction**: Some technical terms have different translation rules depending on context:

1. **In code blocks**: Keep ALL Nextflow syntax in English (the code must run)
2. **In code comments**: TRANSLATE comments to Polish (they are not executable)
3. **In prose/explanatory text**: Follow the glossary below for translations; use Polish words where available unless referring to code keywords verbatim; always use Markdown inline code for keywords

For example:

- In prose: "Kanał wejściowy otrzymuje pliki..." (translate "channel" to "kanał"); BUT: "kanały można tworzyć przy pomocy przestrzeni nazw `channel`" (verbatim syntax reference, do not translate)
- In code: `channel.fromPath('*.fastq')` (keep "channel" in English)
- In comments: `// emit a greeting` → `// Wyemituj powitanie`

## Glossary

### Technical terms (translate in prose if a translation is provided, keep as-is otherwise)

If a translation is available, it is written in the same line after a dash. Assume no translation otherwise.

#### Nextflow Core Concepts

- Nextflow
- DSL2
- channel / channels – kanał / kanały
- process / processes – proces / procesy
- workflow / workflows – workflow / workflow'y
- pipeline / pipelines – pipeline / pipeline'y
- script – skrypt
- shell – powłoka
- exec
- emit
- take
- main
- params – parametry

#### Channel Types

- queue channel – kanał kolejki
- value channel – kanał wartości

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

| English         | Polish                   |
| --------------- | ------------------------ |
| alignment       | dopasowanie              |
| command         | polecenie                |
| container       | kontener                 |
| directive       | dyrektywa                |
| directory       | katalog                  |
| environment     | środowisko               |
| file (general)  | plik                     |
| index           | indeks                   |
| input           | wejście                  |
| module          | moduł                    |
| operator        | operator                 |
| output          | wyjście                  |
| parameter       | parametr                 |
| reference       | referencja               |
| run             | uruchomić / uruchomienie |
| sample          | próbka                   |
| task            | zadanie                  |
| training        | szkolenie                |
| tuple           | krotka                   |
| parallelization | paralelizacja            |
| parallelize     | paralelizować            |

### Admonition Titles

| English  | Polish      |
| -------- | ----------- |
| Note     | Uwaga       |
| Tip      | Wskazówka   |
| Warning  | Ostrzeżenie |
| Exercise | Ćwiczenie   |
| Solution | Rozwiązanie |
| Example  | Przykład    |
