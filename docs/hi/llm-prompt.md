# Translation Rules for Hindi

The target language for this translation is **Hindi** (`hi`).

## Grammar Preferences

- Use informal tone (तुम instead of आप where appropriate for technical tutorials)
- Use standard Devanagari script
- Prefer active voice when possible
- Use Hindi numerals or Arabic numerals as appropriate for context
- Maintain a conversational yet professional tone

## Translation Context Rules

**Important distinction**: Some technical terms have different translation rules depending on context:

1. **In code blocks**: Keep ALL Nextflow syntax in English (the code must run)
2. **In code comments**: TRANSLATE comments to Hindi (they are not executable)
3. **In prose/explanatory text**: Follow the glossary below for translations

For example:

- In prose: "इनपुट channel फ़ाइलें प्राप्त करता है..." (translate "channel" to "चैनल" or keep in English with Hindi text)
- In code: `channel.fromPath('*.fastq')` (keep "channel" in English)
- In comments: `// emit a greeting` → `// एक अभिवादन emit करें`

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

| English        | Hindi      |
| -------------- | ---------- |
| alignment      | संरेखण     |
| command        | कमांड      |
| container      | कंटेनर     |
| directive      | निर्देश    |
| directory      | डायरेक्टरी |
| environment    | वातावरण    |
| file (general) | फ़ाइल      |
| index          | इंडेक्स    |
| input          | इनपुट      |
| module         | मॉड्यूल    |
| operator       | ऑपरेटर     |
| output         | आउटपुट     |
| parameter      | पैरामीटर   |
| reference      | संदर्भ     |
| run            | चलाना / रन |
| sample         | नमूना      |
| task           | कार्य      |
| training       | प्रशिक्षण  |
| tuple          | टपल        |

### Admonition Titles

| English  | Hindi   |
| -------- | ------- |
| Note     | नोट     |
| Tip      | सुझाव   |
| Warning  | चेतावनी |
| Exercise | अभ्यास  |
| Solution | समाधान  |
| Example  | उदाहरण  |
