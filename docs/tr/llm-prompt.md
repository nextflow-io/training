# Translation Rules for Turkish

The target language for this translation is **Turkish** (`tr`).

## Grammar Preferences

- Use formal tone (siz instead of sen)
- Follow Turkish Language Association (TDK) spelling conventions
- Prefer active voice when possible
- Pay attention to vowel harmony in suffixes

## Translation Context Rules

**Important distinction**: Some technical terms have different translation rules depending on context:

1. **In code blocks**: Keep ALL Nextflow syntax in English (the code must run)
2. **In prose/explanatory text**: Follow the glossary below for translations

For example:

- In prose: "Giriş kanalı dosyaları alır..." (translate "channel" to "kanal")
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

| English        | Turkish                  |
| -------------- | ------------------------ |
| alignment      | hizalama                 |
| command        | komut                    |
| container      | konteyner                |
| directive      | yönerge                  |
| directory      | dizin                    |
| environment    | ortam                    |
| file (general) | dosya                    |
| index          | dizin                    |
| input          | girdi                    |
| module         | modül                    |
| operator       | operatör                 |
| output         | çıktı                    |
| parameter      | parametre                |
| reference      | referans                 |
| run            | çalıştırmak / çalıştırma |
| sample         | örnek                    |
| task           | görev                    |
| training       | eğitim                   |
| tuple          | demet                    |

### Admonition Titles

| English  | Turkish   |
| -------- | --------- |
| Note     | Not       |
| Tip      | İpucu     |
| Warning  | Uyarı     |
| Exercise | Alıştırma |
| Solution | Çözüm     |
| Example  | Örnek     |
