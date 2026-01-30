# Translation Rules for Korean

The target language for this translation is **Korean** (`ko`).

## Grammar Preferences

- Use formal polite speech level (합쇼체/하십시오체)
- Maintain consistent honorific usage throughout
- Follow standard Korean spelling conventions (한글 맞춤법)

## Translation Context Rules

**Important distinction**: Some technical terms have different translation rules depending on context:

1. **In code blocks**: Keep ALL Nextflow syntax in English (the code must run)
2. **In prose/explanatory text**: Follow the glossary below for translations

For example:

- In prose: "입력 채널이 파일을 받습니다..." (translate "channel" to "채널")
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

| English        | Korean   |
| -------------- | -------- |
| alignment      | 정렬     |
| command        | 명령     |
| container      | 컨테이너 |
| directive      | 지시문   |
| directory      | 디렉토리 |
| environment    | 환경     |
| file (general) | 파일     |
| index          | 인덱스   |
| input          | 입력     |
| module         | 모듈     |
| operator       | 연산자   |
| output         | 출력     |
| parameter      | 매개변수 |
| reference      | 참조     |
| run            | 실행     |
| sample         | 샘플     |
| task           | 작업     |
| training       | 교육     |
| tuple          | 튜플     |

### Admonition Titles

| English  | Korean |
| -------- | ------ |
| Note     | 참고   |
| Tip      | 팁     |
| Warning  | 경고   |
| Exercise | 연습   |
| Solution | 해결책 |
| Example  | 예제   |
