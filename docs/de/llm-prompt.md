# Translation Rules for German

The target language for this translation is **German** (`de`).

## Grammar Preferences

- Use informal tone (du instead of Sie)
- Use standard German (Hochdeutsch) spelling conventions
- Prefer active voice when possible
- Use gender-neutral language where appropriate

## Translation Context Rules

**Important distinction**: Some technical terms have different translation rules depending on context:

1. **In code blocks**: Keep ALL Nextflow syntax in English (the code must run)
2. **In code comments**: TRANSLATE comments to German (they are not executable)
3. **In prose/explanatory text**: Follow the glossary below for translations

For example:

- In prose: "Der Eingabekanal empfängt die Dateien..." (translate "channel" to "Kanal")
- In code: `channel.fromPath('*.fastq')` (keep "channel" in English)
- In comments: `// emit a greeting` → `// Eine Begrüßung ausgeben`

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

| English        | German                 |
| -------------- | ---------------------- |
| alignment      | Alignment              |
| command        | Befehl                 |
| container      | Container              |
| directive      | Direktive              |
| directory      | Verzeichnis            |
| environment    | Umgebung               |
| file (general) | Datei                  |
| index          | Index                  |
| input          | Eingabe                |
| module         | Modul                  |
| operator       | Operator               |
| output         | Ausgabe                |
| parameter      | Parameter              |
| reference      | Referenz               |
| run            | ausführen / Ausführung |
| sample         | Probe                  |
| task           | Aufgabe                |
| training       | Training               |
| tuple          | Tupel                  |

### Admonition Titles

| English  | German   |
| -------- | -------- |
| Note     | Hinweis  |
| Tip      | Tipp     |
| Warning  | Warnung  |
| Exercise | Übung    |
| Solution | Lösung   |
| Example  | Beispiel |
