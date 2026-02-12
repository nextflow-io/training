# nf4_science 과정 템플릿

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } AI 지원 번역 - [자세히 알아보기 및 개선 제안](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

이것은 `nf4_science/` 하위에 새로운 도메인별 과정을 생성하기 위한 범용 템플릿입니다.
Genomics 및 RNAseq 과정에서 발견되는 공통 패턴을 기반으로 합니다.

## 사용 방법

1. 이 디렉토리와 해당 스크립트 디렉토리를 복사하여 과정을 생성하세요:

   ```bash
   # 문서
   cp -r docs/en/docs/nf4_science/_template docs/en/docs/nf4_science/{your_domain}

   # 스크립트
   cp -r nf4-science/_template nf4-science/{your_domain}
   ```

2. 워크플로우 스크립트 이름을 변경하세요:

   ```bash
   mv nf4-science/{your_domain}/_template.nf nf4-science/{your_domain}/{your_domain}.nf
   ```

3. 모듈 파일 이름을 도구에 맞게 변경하세요:

   ```bash
   mv nf4-science/{your_domain}/modules/_tool_a.nf nf4-science/{your_domain}/modules/{tool_a}.nf
   mv nf4-science/{your_domain}/modules/_tool_b.nf nf4-science/{your_domain}/modules/{tool_b}.nf
   ```

4. 문서와 스크립트에서 모든 `{PLACEHOLDER}` 값을 검색하여 내용으로 교체하세요.

5. `mkdocs.yml` nav에 추가하세요 (아래 코드 조각 참조).

6. `data/` 디렉토리를 테스트 데이터셋으로 채우세요.

7. 각 파트에 대한 작동하는 코드로 `solutions/` 디렉토리를 구성하세요.

## mkdocs.yml nav 코드 조각

```yaml
- { Your Domain }:
    - nf4_science/{your_domain}/index.md
    - nf4_science/{your_domain}/00_orientation.md
    - nf4_science/{your_domain}/01_method.md
    - nf4_science/{your_domain}/02_single_sample.md
    - nf4_science/{your_domain}/03_multi_sample.md
    - nf4_science/{your_domain}/survey.md
    - nf4_science/{your_domain}/next_steps.md
```

## 플레이스홀더 참조

### 과정 수준 플레이스홀더

| 플레이스홀더                 | 설명                       | 예시 (Genomics)                            |
| ---------------------------- | -------------------------- | ------------------------------------------ |
| `{DOMAIN}`                   | 도메인 이름 (제목 형식)    | Genomics                                   |
| `{DOMAIN_DIR}`               | 디렉토리 이름 (소문자)     | genomics                                   |
| `{METHOD}`                   | 분석 방법 이름             | variant calling                            |
| `{METHOD_SHORT_DESCRIPTION}` | 한 줄 방법 설명            | variant calling with GATK                  |
| `{ACCESSORY_FILES}`          | 사용되는 보조 파일 유형    | index files and reference genome resources |

### 도구 플레이스홀더

| 플레이스홀더                          | 설명                       | 예시 (Genomics)                                     |
| ------------------------------------- | -------------------------- | --------------------------------------------------- |
| `{TOOL_A}` / `{TOOL_B}`               | 도구 표시 이름             | Samtools / GATK                                     |
| `{TOOL_A_URL}` / `{TOOL_B_URL}`       | 도구 문서 URL              | https://www.htslib.org/                             |
| `{TOOL_A_CONTAINER_URI}`              | 전체 컨테이너 URI          | community.wave.seqera.io/library/samtools:1.20--... |
| `{TOOL_A_MODULE}` / `{TOOL_B_MODULE}` | 모듈 파일 이름 (.nf 제외)  | samtools_index / gatk_haplotypecaller               |
| `{TOOL_A_PROCESS_NAME}`               | 대문자 프로세스 이름       | SAMTOOLS_INDEX                                      |
| `{TOOL_A_COMMAND}`                    | 실행할 셸 명령             | samtools index '<input_bam>'                        |

### 입력/출력 플레이스홀더

| 플레이스홀더           | 설명                        | 예시 (Genomics)      |
| ---------------------- | --------------------------- | -------------------- |
| `{PRIMARY_INPUT_TYPE}` | 주요 입력 파일 유형         | BAM file             |
| `{PRIMARY_PARAM_NAME}` | Nextflow 매개변수 이름      | reads_bam            |
| `{TEST_INPUT_PATH}`    | data/ 내 테스트 입력 경로   | bam/reads_mother.bam |
| `{FINAL_OUTPUT_TYPE}`  | 최종 파이프라인 출력 유형   | joint-called VCFs    |

### 콘텐츠 플레이스홀더

| 플레이스홀더            | 설명                                 |
| ----------------------- | ------------------------------------ |
| `{PART2_SUMMARY}`       | Part 2 범위에 대한 한 줄 요약        |
| `{AGGREGATION_SUMMARY}` | 집계 단계에 대한 한 줄 요약          |
| `{PROCESS_LIST}`        | 쉼표로 구분된 프로세스 이름 목록     |
| `{TYPEFORM_ID}`         | 설문조사용 Typeform 임베드 ID        |

## 교육 구조

템플릿은 세 부분으로 구성된 구조를 따릅니다:

1. **Part 1 (01_method.md)**: Docker 컨테이너에서 수동 테스트를 통해 방법론 이해
2. **Part 2 (02_single_sample.md)**: 명령을 Nextflow로 래핑; 단일 샘플, 그 다음 배치
3. **Part 3 (03_multi_sample.md)**: 채널 연산자를 사용하여 다중 샘플 집계 추가

### 주요 규칙

- **전/후 탭 코드 블록**을 `hl_lines`와 함께 사용하여 코드 변경 사항 표시
- 접을 수 있는 예상 출력에 `??? success "Command output"` 사용
- 접을 수 있는 디렉토리 트리에 `??? abstract "Directory contents"` 사용
- 각 주요 섹션을 **핵심 정리** 및 **다음 단계** 하위 섹션으로 마무리
- 주요 번호가 매겨진 섹션을 구분하기 위해 `---` 수평선 사용
- 학습자가 채울 수 있도록 스켈레톤 파일(워크플로우 + 모듈) 제공
- 파트별로 해결책 구성 (`solutions/part2/`, `solutions/part3/`)

## 디렉토리 구조

```
docs/en/docs/nf4_science/{domain}/
├── index.md                    # Course overview with frontmatter
├── 00_orientation.md           # Environment setup
├── 01_method.md                # Manual testing in containers
├── 02_single_sample.md         # Single-sample Nextflow implementation
├── 03_multi_sample.md          # Multi-sample aggregation
├── survey.md                   # Typeform feedback survey
├── next_steps.md               # Course summary and suggestions
└── img/                        # Diagrams (.excalidraw.svg, .png)

nf4-science/{domain}/
├── {domain}.nf                 # Skeleton workflow file
├── nextflow.config             # Minimal config (docker.enabled = true)
├── data/                       # Test datasets and resources
│   └── samplesheet.csv         # Sample metadata
├── modules/                    # Skeleton module files
│   ├── {tool_a}.nf
│   └── {tool_b}.nf
└── solutions/
    ├── part2/                  # Complete Part 2 solution
    │   ├── {domain}-2.nf
    │   ├── nextflow.config
    │   └── modules/
    └── part3/                  # Complete Part 3 solution
        ├── {domain}-3.nf
        ├── nextflow.config
        └── modules/
```
