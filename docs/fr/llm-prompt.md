# Translation Rules for French

The target language for this translation is **French** (`fr`).

## Grammar Preferences

- Use formal tone (vous instead of tu)
- Follow standard French spelling conventions
- Prefer active voice when possible
- Standard quotation marks ("") are acceptable; French quotation marks (« ») are optional

### Inclusive Writing (Écriture Inclusive)

Consider using inclusive writing with middle dots for gender-neutral terms:

- "débutant·es" (beginners, both genders)
- "prêt·e" (ready, both genders)
- "chercheur·euses" (researchers, both genders)
- "formateur·trice" (trainer, both genders)

This is optional but aligns with modern French writing practices.

## Translation Context Rules

**Important distinction**: Some technical terms have different translation rules depending on context:

1. **In code blocks**: Keep ALL Nextflow syntax in English (the code must run)
2. **In code comments**: TRANSLATE comments to French (they are not executable)
3. **In prose/explanatory text**: Follow the glossary below for translations

For example:

- In prose: "Le canal d'entrée reçoit les fichiers..." (translate "channel" to "canal")
- In code: `channel.fromPath('*.fastq')` (keep "channel" in English)
- In comments: `// emit a greeting` → `// émet un message de bienvenue`

## Code Comments

**Always translate code comments to French.** Comments are not executable code and should be in the target language for better comprehension.

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

### Terms to Keep in English

Note: French technical writing commonly keeps certain English terms, especially in software contexts:

- **workflow** - commonly kept in English in French tech docs (not "flux de travail")
- **pipeline** - commonly kept in English (not "pipeline de traitement")
- Course/module names: "Hello Nextflow", "Hello nf-core" - keep as proper nouns

### Terms to Translate

| English         | French               | Notes                            |
| --------------- | -------------------- | -------------------------------- |
| training        | formation            |                                  |
| workshop        | atelier              |                                  |
| Side Quests     | Quêtes secondaires   | Creative translation, keep theme |
| quick reference | Référence rapide     |                                  |
| alignment       | alignement           |                                  |
| command         | commande             |                                  |
| container       | conteneur            |                                  |
| directive       | directive            |                                  |
| directory       | répertoire           |                                  |
| environment     | environnement        |                                  |
| file (general)  | fichier              |                                  |
| index           | index                |                                  |
| input           | entrée               |                                  |
| module          | module               |                                  |
| operator        | opérateur            |                                  |
| output          | sortie               |                                  |
| parameter       | paramètre            |                                  |
| reference       | référence            |                                  |
| run             | exécuter / exécution |                                  |
| sample          | échantillon          |                                  |
| task            | tâche                |                                  |
| tuple           | tuple                |                                  |
| channel         | canal                | In prose; keep English in code   |
| process         | processus            | In prose; keep English in code   |

### Admonition Titles

| English  | French        |
| -------- | ------------- |
| Note     | Note          |
| Tip      | Astuce        |
| Warning  | Avertissement |
| Exercise | Exercice      |
| Solution | Solution      |
| Example  | Exemple       |

### Section Headers

| English           | French                            |
| ----------------- | --------------------------------- |
| Environment Setup | Configuration de l'environnement  |
| Getting Started   | Premiers pas                      |
| (course titles)   | (keep in English as proper nouns) |
