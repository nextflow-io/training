# Translation Rules for Korean

The target language for this translation is **Korean** (`ko`).

## Grammar Preferences

- Use formal polite speech level (합쇼체/하십시오체)
- Sentence endings: "~입니다", "~습니다", "~하세요", "~수 있습니다"
- Maintain consistent honorific usage throughout
- Follow standard Korean spelling conventions (한글 맞춤법)
- Preserve personality and humor where present (e.g., "그동안 차 한 잔을 준비하거나" for suggesting a tea break)

## Translation Context Rules

**Important distinction**: Some technical terms have different translation rules depending on context:

1. **In code blocks**: Keep ALL Nextflow syntax in English (the code must run)
2. **In code comments**: TRANSLATE comments to Korean (they are not executable)
3. **In prose/explanatory text**: Follow the glossary below for translations

For example:

- In prose: "입력 채널이 파일을 받습니다..." (translate "channel" to "채널")
- In code: `channel.fromPath('*.fastq')` (keep "channel" in English)
- In comments: `// emit a greeting` → `// 인사말을 내보냅니다`

## Code Comments

**Always translate code comments to Korean.** Comments are not executable code and should be in the target language for better comprehension.

## English Terms with Korean Particles

Korean uses particles attached to words. When English terms appear in Korean text, attach Korean particles directly:

- "README.md 파일" (README.md file)
- "VSCode IDE에서" (in VSCode IDE)
- "Docker를" (Docker + object marker)

## UI Element Pattern

For UI elements, consider including English in parentheses to help users match translations to what they see in the (often English) interface:

- "**사이드바(Sidebar)**"
- "**파일 탐색기(File Explorer)**"

## Image Alt Text

Translate image alt text to Korean for accessibility:

- "GitHub Codespace 세션 리스트"
- "GitHub Codespaces 환영 화면"

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

Note: Some technical terms use transliteration (음차) while others use actual Korean translations. Follow these patterns:

#### Transliterations (음차)

| English   | Korean     | Notes                  |
| --------- | ---------- | ---------------------- |
| workflow  | 워크플로우 | Common transliteration |
| pipeline  | 파이프라인 | Common transliteration |
| terminal  | 터미널     | Common transliteration |
| container | 컨테이너   | Common transliteration |
| sidebar   | 사이드바   | Common transliteration |
| module    | 모듈       | Common transliteration |
| sample    | 샘플       | Common transliteration |
| tuple     | 튜플       | Common transliteration |
| index     | 인덱스     | Common transliteration |
| core      | 코어       | Common transliteration |

#### Actual Translations

| English        | Korean      | Notes                                         |
| -------------- | ----------- | --------------------------------------------- |
| standalone     | 단독        | Not 독립 (political independence)             |
| mini-course    | 단기 과정   | Not 미니 과정 (too casual)                    |
| wrap           | 적용하다    | Not 래핑하다 (sounds like physical packaging) |
| training       | 교육        | Actual translation                            |
| environment    | 환경        | Actual translation                            |
| repository     | 저장소      | Actual translation                            |
| file explorer  | 파일 탐색기 | Mixed                                         |
| main editor    | 메인 편집기 | Mixed                                         |
| directory      | 디렉토리    | Transliteration acceptable                    |
| file (general) | 파일        | Transliteration                               |
| alignment      | 정렬        | Actual translation                            |
| command        | 명령        | Actual translation                            |
| directive      | 지시문      | Actual translation                            |
| input          | 입력        | Actual translation                            |
| output         | 출력        | Actual translation                            |
| operator       | 연산자      | Actual translation                            |
| parameter      | 매개변수    | Actual translation                            |
| reference      | 참조        | Actual translation                            |
| run            | 실행        | Actual translation                            |
| task           | 작업        | Actual translation                            |
| lowercase      | 소문자      | Actual translation                            |
| materials      | 자료        | Actual translation                            |
| course         | 과정        | Actual translation                            |

### Admonition Titles

| English  | Korean |
| -------- | ------ |
| Note     | 참고   |
| Tip      | 팁     |
| Warning  | 경고   |
| Exercise | 연습   |
| Solution | 해결책 |
| Example  | 예제   |

### Section Headers

| English           | Korean                            |
| ----------------- | --------------------------------- |
| Environment Setup | 환경 설정                         |
| Getting Started   | 시작하기                          |
| (course titles)   | (keep in English as proper nouns) |
