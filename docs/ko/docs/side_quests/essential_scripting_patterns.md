# 필수 Nextflow 스크립팅 패턴

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } AI 지원 번역 - [자세히 알아보기 및 개선 제안](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Nextflow는 Java Virtual Machine에서 실행되는 프로그래밍 언어입니다. Nextflow는 [Groovy](http://groovy-lang.org/)를 기반으로 구축되어 많은 구문을 공유하지만, Nextflow는 단순히 "확장 기능이 있는 Groovy"가 아닙니다. 완전히 명세화된 [구문](https://nextflow.io/docs/latest/reference/syntax.html)과 [표준 라이브러리](https://nextflow.io/docs/latest/reference/stdlib.html)를 갖춘 독립적인 언어입니다.

변수, 맵, 리스트에 대한 기본 구문을 넘어서지 않고도 많은 Nextflow 코드를 작성할 수 있습니다. 대부분의 Nextflow 튜토리얼은 워크플로 오케스트레이션(채널, 프로세스, 데이터 흐름)에 초점을 맞추며, 이것만으로도 놀라울 정도로 많은 것을 할 수 있습니다.

그러나 데이터를 조작하거나, 복잡한 파일 이름을 파싱하거나, 조건부 로직을 구현하거나, 견고한 프로덕션 워크플로를 구축해야 할 때는 코드의 두 가지 측면을 구분하여 생각하는 것이 도움이 됩니다: **데이터플로**(채널, 연산자, 프로세스, 워크플로)와 **스크립팅**(클로저, 함수, 프로세스 스크립트 내부의 코드). 이 구분은 다소 임의적이지만(모두 Nextflow 코드입니다), 파이프라인을 오케스트레이션할 때와 데이터를 조작할 때를 이해하는 데 유용한 멘탈 모델을 제공합니다. 두 가지를 모두 마스터하면 명확하고 유지보수 가능한 워크플로를 작성하는 능력이 크게 향상됩니다.

### 학습 목표

이 사이드 퀘스트는 기본 개념부터 프로덕션 준비 패턴까지 실습을 통한 여정을 안내합니다.
간단한 CSV 읽기 워크플로를 정교한 생물정보학 파이프라인으로 변환하며, 현실적인 과제를 통해 단계별로 발전시킵니다:

- **경계 이해하기:** 데이터플로 작업과 스크립팅을 구분하고 이들이 함께 작동하는 방식 이해하기
- **데이터 조작:** 강력한 연산자를 사용하여 맵과 컬렉션 추출, 변환, 부분집합 만들기
- **문자열 처리:** 정규식 패턴으로 복잡한 파일 명명 규칙 파싱 및 변수 보간 마스터하기
- **재사용 가능한 함수:** 복잡한 로직을 명명된 함수로 추출하여 더 깔끔하고 유지보수 가능한 워크플로 만들기
- **동적 로직:** 다양한 입력 타입에 적응하는 프로세스 구축 및 동적 리소스 할당을 위한 클로저 사용
- **조건부 라우팅:** 메타데이터 특성에 따라 샘플을 다양한 프로세스로 지능적으로 라우팅
- **안전한 작업:** null-safe 연산자로 누락된 데이터를 우아하게 처리하고 명확한 오류 메시지로 입력 검증
- **설정 기반 핸들러:** 로깅, 알림, 라이프사이클 관리를 위한 워크플로 이벤트 핸들러 사용

### 사전 요구사항

이 사이드 퀘스트를 시작하기 전에:

- [Hello Nextflow](../hello_nextflow/README.md) 튜토리얼 또는 이에 상응하는 초급 과정을 완료해야 합니다.
- 기본 Nextflow 개념과 메커니즘(프로세스, 채널, 연산자, 파일 작업, 메타데이터)을 편안하게 사용할 수 있어야 합니다.
- 일반적인 프로그래밍 구조(변수, 맵, 리스트)에 대한 기본적인 친숙도가 필요합니다.

이 튜토리얼은 프로그래밍 개념을 만나면서 설명하므로 광범위한 프로그래밍 경험이 필요하지 않습니다.
기본 개념부터 시작하여 고급 패턴까지 단계적으로 구축합니다.

---

## 0. 시작하기

#### 교육 코드스페이스 열기

아직 완료하지 않았다면, [환경 설정](../envsetup/index.md)에 설명된 대로 교육 환경을 엽니다.

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

#### 프로젝트 디렉토리로 이동

이 튜토리얼의 파일이 있는 디렉토리로 이동합니다.

```bash
cd side-quests/essential_scripting_patterns
```

#### 자료 검토

메인 워크플로 파일과 예제 데이터 파일이 포함된 `data` 디렉토리를 찾을 수 있습니다.

```console title="디렉토리 내용"
.
├── collect.nf
├── data
│   ├── samples.csv
│   └── sequences
│       ├── SAMPLE_001_S1_L001_R1_001.fastq
│       ├── SAMPLE_002_S2_L001_R1_001.fastq
│       └── SAMPLE_003_S3_L001_R1_001.fastq
├── main.nf
├── modules
│   ├── fastp.nf
│   ├── generate_report.nf
│   └── trimgalore.nf
└── nextflow.config
```

샘플 CSV에는 특성에 따라 다른 처리가 필요한 생물학적 샘플에 대한 정보가 포함되어 있습니다:

```console title="samples.csv"
sample_id,organism,tissue_type,sequencing_depth,file_path,quality_score
SAMPLE_001,human,liver,30000000,data/sequences/SAMPLE_001_S1_L001_R1_001.fastq,38.5
SAMPLE_002,mouse,brain,25000000,data/sequences/SAMPLE_002_S2_L001_R1_001.fastq,35.2
SAMPLE_003,human,kidney,45000000,data/sequences/SAMPLE_003_S3_L001_R1_001.fastq,42.1
```

이 현실적인 데이터셋을 사용하여 실제 생물정보학 워크플로에서 마주치게 될 실용적인 프로그래밍 기법을 탐색합니다.

<!-- TODO: Can we make this more domain-agnostic? -->

<!-- TODO: add an assignment statement? #### Review the assignment -->

#### 준비 체크리스트

시작할 준비가 되었나요?

- [ ] 이 과정의 목표와 사전 요구사항을 이해했습니다
- [ ] 코드스페이스가 실행 중입니다
- [ ] 작업 디렉토리를 적절하게 설정했습니다
<!-- - [ ] I understand the assignment -->

모든 항목을 체크할 수 있다면, 시작할 준비가 되었습니다.

---

## 1. 데이터플로 vs 스크립팅: 경계 이해하기

### 1.1. 무엇이 무엇인지 식별하기

Nextflow 워크플로를 작성할 때 **데이터플로**(채널과 프로세스를 통한 데이터 이동 방식)와 **스크립팅**(데이터를 조작하고 결정을 내리는 코드)을 구분하는 것이 중요합니다. 이들이 함께 작동하는 방식을 보여주는 워크플로를 구축해 봅시다.

#### 1.1.1. 기본 Nextflow 워크플로

CSV 파일을 읽는 간단한 워크플로부터 시작합니다(`main.nf`에 이미 준비되어 있습니다):

```groovy title="main.nf" linenums="1"
workflow {
    ch_samples = channel.fromPath("./data/samples.csv")
        .splitCsv(header: true)
        .view()
}
```

`workflow` 블록은 파이프라인 구조를 정의하고, `channel.fromPath()`는 파일 경로에서 채널을 생성합니다. `.splitCsv()` 연산자는 CSV 파일을 처리하고 각 행을 맵 데이터 구조로 변환합니다.

이 워크플로를 실행하여 원시 CSV 데이터를 확인합니다:

```bash
nextflow run main.nf
```

??? success "명령 출력"

    ```console
    Launching `main.nf` [marvelous_tuckerman] DSL2 - revision: 6113e05c17

    [sample_id:SAMPLE_001, organism:human, tissue_type:liver, sequencing_depth:30000000, file_path:data/sequences/SAMPLE_001_S1_L001_R1_001.fastq, quality_score:38.5]
    [sample_id:SAMPLE_002, organism:mouse, tissue_type:brain, sequencing_depth:25000000, file_path:data/sequences/SAMPLE_002_S2_L001_R1_001.fastq, quality_score:35.2]
    [sample_id:SAMPLE_003, organism:human, tissue_type:kidney, sequencing_depth:45000000, file_path:data/sequences/SAMPLE_003_S3_L001_R1_001.fastq, quality_score:42.1]
    ```

#### 1.1.2. Map 연산자 추가하기

이제 이미 익숙할 `.map()` 연산자를 사용하여 데이터를 변환하는 스크립팅을 추가합니다. 이 연산자는 각 항목을 변환할 수 있는 코드를 작성할 수 있는 '클로저'를 받습니다.

!!! note

    **클로저**는 전달되고 나중에 실행될 수 있는 코드 블록입니다. 인라인으로 정의하는 함수로 생각하면 됩니다. 클로저는 중괄호 `{ }`로 작성되며 매개변수를 받을 수 있습니다. Nextflow 연산자 작동 방식의 기본이며, Nextflow를 작성한 지 오래되었다면 이미 알지 못한 채 사용하고 있었을 수 있습니다!

다음은 map 연산의 모습입니다:

=== "변경 후"

    ```groovy title="main.nf" linenums="2" hl_lines="3-6"
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map { row ->
                return row
            }
            .view()
    ```

=== "변경 전"

    ```groovy title="main.nf" linenums="2" hl_lines="3"
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .view()
    ```

이것이 첫 번째 **클로저**입니다 - 인수로 전달할 수 있는 익명 함수입니다(Python의 람다나 JavaScript의 화살표 함수와 유사합니다). 클로저는 Nextflow 연산자 작업에 필수적입니다.

클로저 `{ row -> return row }`는 매개변수 `row`를 받습니다(어떤 이름이든 가능: `item`, `sample` 등).

`.map()` 연산자가 각 채널 항목을 처리할 때 해당 항목을 클로저에 전달합니다. 여기서 `row`는 한 번에 하나의 CSV 행을 보유합니다.

이 변경 사항을 적용하고 워크플로를 실행합니다:

```bash
nextflow run main.nf
```

map 연산자가 입력을 변경 없이 반환하므로 이전과 동일한 출력이 표시됩니다. 이는 map 연산자가 올바르게 작동하는 것을 확인합니다. 이제 데이터 변환을 시작해 봅시다.

#### 1.1.3. Map 데이터 구조 생성하기

이제 클로저 내부에 **스크립팅** 로직을 작성하여 각 데이터 행을 변환합니다. 여기서 데이터 흐름을 오케스트레이션하는 것이 아니라 개별 데이터 항목을 처리합니다.

=== "변경 후"

    ```groovy title="main.nf" linenums="2" hl_lines="4-12"
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map { row ->
                // 데이터 변환을 위한 스크립팅
                def sample_meta = [
                    id: row.sample_id.toLowerCase(),
                    organism: row.organism,
                    tissue: row.tissue_type.replaceAll('_', ' ').toLowerCase(),
                    depth: row.sequencing_depth.toInteger(),
                    quality: row.quality_score.toDouble()
                ]
                return sample_meta
            }
            .view()
    ```

=== "변경 전"

    ```groovy title="main.nf" linenums="2" hl_lines="4"
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map { row ->
                return row
            }
            .view()
    ```

`sample_meta` 맵은 관련 정보를 저장하는 키-값 데이터 구조입니다(Python의 딕셔너리, JavaScript의 객체, Ruby의 해시와 유사): 샘플 ID, 유기체, 조직 유형, 시퀀싱 깊이, 품질 점수.

`.toLowerCase()`와 `.replaceAll()` 같은 문자열 조작 메서드를 사용하여 데이터를 정리하고, `.toInteger()`와 `.toDouble()` 같은 타입 변환 메서드를 사용하여 CSV의 문자열 데이터를 적절한 숫자 타입으로 변환합니다.

이 변경 사항을 적용하고 워크플로를 실행합니다:

```bash
nextflow run main.nf
```

??? success "명령 출력"

    ```console
    [id:sample_001, organism:human, tissue:liver, depth:30000000, quality:38.5]
    [id:sample_002, organism:mouse, tissue:brain, depth:25000000, quality:35.2]
    [id:sample_003, organism:human, tissue:kidney, depth:45000000, quality:42.1]
    ```

#### 1.1.4. 조건부 로직 추가하기

이제 더 많은 스크립팅을 추가합니다 - 이번에는 삼항 연산자를 사용하여 데이터 값을 기반으로 결정을 내립니다.

다음 변경을 수행합니다:

=== "변경 후"

    ```groovy title="main.nf" linenums="2" hl_lines="11-12"
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map { row ->
                def sample_meta = [
                    id: row.sample_id.toLowerCase(),
                    organism: row.organism,
                    tissue: row.tissue_type.replaceAll('_', ' ').toLowerCase(),
                    depth: row.sequencing_depth.toInteger(),
                    quality: row.quality_score.toDouble()
                ]
                def priority = sample_meta.quality > 40 ? 'high' : 'normal'
                return sample_meta + [priority: priority]
            }
            .view()
    ```

=== "변경 전"

    ```groovy title="main.nf" linenums="2" hl_lines="11"
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map { row ->
                def sample_meta = [
                    id: row.sample_id.toLowerCase(),
                    organism: row.organism,
                    tissue: row.tissue_type.replaceAll('_', ' ').toLowerCase(),
                    depth: row.sequencing_depth.toInteger(),
                    quality: row.quality_score.toDouble()
                ]
                return sample_meta
            }
            .view()
    ```

삼항 연산자는 `condition ? value_if_true : value_if_false` 패턴을 따르는 if/else 문의 축약형입니다. 이 줄은 "품질이 40보다 크면 'high'를 사용하고, 그렇지 않으면 'normal'을 사용한다"는 의미입니다. 이와 유사한 **Elvis 연산자**(`?:`)는 무언가가 null이거나 비어있을 때 기본값을 제공합니다 - 이 패턴은 이 튜토리얼의 뒷부분에서 살펴보겠습니다.

맵 덧셈 연산자 `+`는 기존 맵을 수정하는 것이 아니라 **새로운 맵**을 생성합니다. 이 줄은 `sample_meta`의 모든 키-값 쌍과 새로운 `priority` 키를 포함하는 새로운 맵을 생성합니다.

!!! Note

    클로저에 전달된 맵을 절대 수정하지 마세요 - 항상 `+`(예를 들어)를 사용하여 새로운 맵을 생성하세요. Nextflow에서는 동일한 데이터가 종종 여러 작업을 통해 동시에 흐릅니다. 맵을 제자리에서 수정하면 다른 작업이 동일한 객체를 참조할 때 예측할 수 없는 부작용이 발생할 수 있습니다. 새로운 맵을 생성하면 각 작업이 자체적인 깨끗한 복사본을 갖도록 보장합니다.

수정된 워크플로를 실행합니다:

```bash
nextflow run main.nf
```

??? success "명령 출력"

    ```console
    [id:sample_001, organism:human, tissue:liver, depth:30000000, quality:38.5, priority:normal]
    [id:sample_002, organism:mouse, tissue:brain, depth:25000000, quality:35.2, priority:normal]
    [id:sample_003, organism:human, tissue:kidney, depth:45000000, quality:42.1, priority:high]
    ```

품질 점수를 기반으로 우선순위 수준으로 메타데이터를 풍부하게 하는 조건부 로직을 성공적으로 추가했습니다.

#### 1.1.5. `.subMap()`으로 맵 부분집합 만들기

`+` 연산자가 맵에 키를 추가하는 반면, 때로는 그 반대가 필요합니다 - 특정 키만 추출하기. `.subMap()` 메서드가 이에 완벽합니다.

식별 필드만 포함하는 메타데이터의 단순화된 버전을 생성하는 줄을 추가해 봅시다:

=== "변경 후"

    ```groovy title="main.nf" linenums="2" hl_lines="12-15"
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map { row ->
                // 데이터 변환을 위한 스크립팅
                def sample_meta = [
                    id: row.sample_id.toLowerCase(),
                    organism: row.organism,
                    tissue: row.tissue_type.replaceAll('_', ' ').toLowerCase(),
                    depth: row.sequencing_depth.toInteger(),
                    quality: row.quality_score.toDouble()
                ]
                def id_only = sample_meta.subMap(['id', 'organism', 'tissue'])
                println "ID fields only: ${id_only}"

                def priority = sample_meta.quality > 40 ? 'high' : 'normal'
                return sample_meta + [priority: priority]
            }
            .view()
    ```

=== "변경 전"

    ```groovy title="main.nf" linenums="2" hl_lines="12"
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map { row ->
                // 데이터 변환을 위한 스크립팅
                def sample_meta = [
                    id: row.sample_id.toLowerCase(),
                    organism: row.organism,
                    tissue: row.tissue_type.replaceAll('_', ' ').toLowerCase(),
                    depth: row.sequencing_depth.toInteger(),
                    quality: row.quality_score.toDouble()
                ]
                def priority = sample_meta.quality > 40 ? 'high' : 'normal'
                return sample_meta + [priority: priority]
            }
            .view()
    ```

수정된 워크플로를 실행합니다:

```bash
nextflow run main.nf
```

??? success "명령 출력"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [peaceful_cori] DSL2 - revision: 4cc4a8340f

    ID fields only: [id:sample_001, organism:human, tissue:liver]
    ID fields only: [id:sample_002, organism:mouse, tissue:brain]
    ID fields only: [id:sample_003, organism:human, tissue:kidney]
    [id:sample_001, organism:human, tissue:liver, depth:30000000, quality:38.5, priority:normal]
    [id:sample_002, organism:mouse, tissue:brain, depth:25000000, quality:35.2, priority:normal]
    [id:sample_003, organism:human, tissue:kidney, depth:45000000, quality:42.1, priority:high]
    ```

이는 `view()` 작업으로 표시된 전체 메타데이터와 `println`으로 출력한 추출된 부분집합을 모두 보여줍니다.

`.subMap()` 메서드는 키 목록을 받아 해당 키만 포함하는 새로운 맵을 반환합니다. 키가 원본 맵에 존재하지 않으면 결과에 포함되지 않습니다.

이는 다양한 프로세스에 대해 다양한 메타데이터 버전을 생성해야 할 때 특히 유용합니다 - 일부는 전체 메타데이터가 필요하고 다른 일부는 최소한의 식별 필드만 필요할 수 있습니다.

이제 앞으로 필요하지 않으므로 워크플로를 이전 상태로 복원하기 위해 println 문을 제거합니다.

!!! tip "맵 작업 요약"

    - **키 추가**: `map1 + [new_key: value]` - 추가 키가 있는 새로운 맵 생성
    - **키 추출**: `map1.subMap(['key1', 'key2'])` - 지정된 키만 있는 새로운 맵 생성
    - **두 작업 모두 새로운 맵을 생성** - 원본 맵은 변경되지 않음

#### 1.1.6. 맵 결합 및 결과 반환

지금까지는 Nextflow 커뮤니티가 '메타 맵'이라고 부르는 것만 반환했고, 해당 메타데이터와 관련된 파일은 무시했습니다. 하지만 Nextflow 워크플로를 작성한다면 아마도 그 파일로 무언가를 하고 싶을 것입니다.

2개 요소로 구성된 튜플을 포함하는 채널 구조를 출력해 봅시다: 풍부한 메타데이터 맵과 해당 파일 경로. 이는 프로세스에 데이터를 전달하는 Nextflow의 일반적인 패턴입니다.

=== "변경 후"

    ```groovy title="main.nf" linenums="2" hl_lines="12"
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map { row ->
                def sample_meta = [
                    id: row.sample_id.toLowerCase(),
                    organism: row.organism,
                    tissue: row.tissue_type.replaceAll('_', ' ').toLowerCase(),
                    depth: row.sequencing_depth.toInteger(),
                    quality: row.quality_score.toDouble()
                ]
                def priority = sample_meta.quality > 40 ? 'high' : 'normal'
                return tuple( sample_meta + [priority: priority], file(row.file_path) )
            }
            .view()
    ```

=== "변경 전"

    ```groovy title="main.nf" linenums="2" hl_lines="12"
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map { row ->
                def sample_meta = [
                    id: row.sample_id.toLowerCase(),
                    organism: row.organism,
                    tissue: row.tissue_type.replaceAll('_', ' ').toLowerCase(),
                    depth: row.sequencing_depth.toInteger(),
                    quality: row.quality_score.toDouble()
                ]
                def priority = sample_meta.quality > 40 ? 'high' : 'normal'
                return sample_meta + [priority: priority]
            }
            .view()
    ```

이 변경 사항을 적용하고 워크플로를 실행합니다:

```bash
nextflow run main.nf
```

??? success "명령 출력"

    ```console
    [[id:sample_001, organism:human, tissue:liver, depth:30000000, quality:38.5, priority:normal], /workspaces/training/side-quests/essential_scripting_patterns/data/sequences/SAMPLE_001_S1_L001_R1_001.fastq]
    [[id:sample_002, organism:mouse, tissue:brain, depth:25000000, quality:35.2, priority:normal], /workspaces/training/side-quests/essential_scripting_patterns/data/sequences/SAMPLE_002_S2_L001_R1_001.fastq]
    [[id:sample_003, organism:human, tissue:kidney, depth:45000000, quality:42.1, priority:high], /workspaces/training/side-quests/essential_scripting_patterns/data/sequences/SAMPLE_003_S3_L001_R1_001.fastq]
    ```

이 `[meta, file]` 튜플 구조는 메타데이터와 관련 파일을 프로세스에 전달하는 Nextflow의 일반적인 패턴입니다.

!!! note

    **맵과 메타데이터**: 맵은 Nextflow에서 메타데이터 작업의 기본입니다. 메타데이터 맵 작업에 대한 자세한 설명은 [메타데이터 작업](./metadata.md) 사이드 퀘스트를 참조하세요.

우리 워크플로는 핵심 패턴을 보여줍니다: **데이터플로 작업**(`workflow`, `channel.fromPath()`, `.splitCsv()`, `.map()`, `.view()`)은 파이프라인을 통한 데이터 이동 방식을 오케스트레이션하고, **스크립팅**(맵 `[key: value]`, 문자열 메서드, 타입 변환, 삼항 연산자)은 `.map()` 클로저 내부에서 개별 데이터 항목의 변환을 처리합니다.

### 1.2. 다른 타입 이해하기: Channel vs List

지금까지는 데이터플로 작업과 스크립팅을 구분할 수 있었습니다. 하지만 두 컨텍스트 모두에 동일한 메서드 이름이 존재하는 경우는 어떻게 될까요?

완벽한 예는 Nextflow 표준 라이브러리의 채널 타입과 List 타입 모두에 존재하는 `collect` 메서드입니다. List의 `collect()` 메서드는 각 요소를 변환하는 반면, 채널의 `collect()` 연산자는 모든 채널 방출을 단일 항목 채널로 수집합니다.

채널 `collect()` 연산자가 무엇을 하는지 복습하면서 샘플 데이터로 이를 시연해 봅시다. `collect.nf`를 확인합니다:

```groovy title="collect.nf" linenums="1"
def sample_ids = ['sample_001', 'sample_002', 'sample_003']

// channel.collect() - 여러 채널 방출을 하나로 그룹화
ch_input = channel.fromList(sample_ids)
ch_input.view { sample -> "Individual channel item: ${sample}" }
ch_collected = ch_input.collect()
ch_collected.view { list -> "channel.collect() result: ${list} (${list.size()} items grouped into 1)" }
```

단계:

- 샘플 ID 목록 정의
- 각 샘플 ID를 별도로 방출하는 `fromList()`로 채널 생성
- 흐르는 각 항목을 `view()`로 출력
- 채널의 `collect()` 연산자로 모든 항목을 단일 리스트로 수집
- 두 번째 `view()`로 수집된 결과(모든 샘플 ID를 포함하는 단일 항목) 출력

채널의 구조를 변경했지만 데이터 자체는 변경하지 않았습니다.

이를 확인하기 위해 워크플로를 실행합니다:

```bash
nextflow run collect.nf
```

??? success "명령 출력"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `collect.nf` [loving_mendel] DSL2 - revision: e8d054a46e

    Individual channel item: sample_001
    Individual channel item: sample_002
    Individual channel item: sample_003
    channel.collect() result: [sample_001, sample_002, sample_003] (3 items grouped into 1)
    ```

`view()`는 모든 채널 방출에 대해 출력을 반환하므로, 이 단일 출력에는 하나의 리스트로 그룹화된 원래 3개 항목이 모두 포함되어 있음을 알 수 있습니다.

이제 List의 `collect` 메서드가 작동하는 것을 봅시다. `collect.nf`를 수정하여 원래 샘플 ID 리스트에 List의 `collect` 메서드를 적용합니다:

=== "변경 후"

    ```groovy title="main.nf" linenums="1" hl_lines="9-13"
    def sample_ids = ['sample_001', 'sample_002', 'sample_003']

    // channel.collect() - 여러 채널 방출을 하나로 그룹화
    ch_input = channel.fromList(sample_ids)
    ch_input.view { sample -> "Individual channel item: ${sample}" }
    ch_collected = ch_input.collect()
    ch_collected.view { list -> "channel.collect() result: ${list} (${list.size()} items grouped into 1)" }

    // List.collect() - 각 요소 변환, 구조 유지
    def formatted_ids = sample_ids.collect { id ->
        id.toUpperCase().replace('SAMPLE_', 'SPECIMEN_')
    }
    println "List.collect() result: ${formatted_ids} (${sample_ids.size()} items transformed into ${formatted_ids.size()})"
    ```

=== "변경 전"

    ```groovy title="main.nf" linenums="1"
    def sample_ids = ['sample_001', 'sample_002', 'sample_003']

    // channel.collect() - 여러 채널 방출을 하나로 그룹화
    ch_input = channel.fromList(sample_ids)
    ch_input.view { sample -> "Individual channel item: ${sample}" }
    ch_collected = ch_input.collect()
    ch_collected.view { list -> "channel.collect() result: ${list} (${list.size()} items grouped into 1)" }
    ```

이 새로운 스니펫에서 우리는:

- List의 `collect` 메서드를 사용하여 원래 리스트의 각 샘플 ID를 변환하는 새로운 변수 `formatted_ids`를 정의
- `println`을 사용하여 결과 출력

수정된 워크플로를 실행합니다:

```bash
nextflow run collect.nf
```

??? success "명령 출력"

    ```console hl_lines="5"
    N E X T F L O W   ~  version 25.10.2

    Launching `collect.nf` [cheeky_stonebraker] DSL2 - revision: 2d5039fb47

    List.collect() result: [SPECIMEN_001, SPECIMEN_002, SPECIMEN_003] (3 items transformed into 3)
    Individual channel item: sample_001
    Individual channel item: sample_002
    Individual channel item: sample_003
    channel.collect() result: [sample_001, sample_002, sample_003] (3 items grouped into 1)
    ```

이번에는 일부 샘플을 제외하는 필터를 선택했기 때문에 데이터 구조를 변경하지 않았고, 리스트에 여전히 3개 항목이 있지만, List의 `collect` 메서드를 사용하여 각 항목을 변환하여 수정된 값이 있는 새로운 리스트를 생성했습니다. 이는 채널에서 `map` 연산자를 사용하는 것과 유사하지만 채널이 아닌 List 데이터 구조에서 작동합니다.

`collect`는 여기서 포인트를 만들기 위해 사용하는 극단적인 경우입니다. 핵심 교훈은 워크플로를 작성할 때 항상 **데이터 구조**(Lists, Maps 등)와 **채널**(데이터플로 구조)을 구분해야 한다는 것입니다. 작업은 이름을 공유할 수 있지만 호출되는 타입에 따라 완전히 다르게 작동합니다.

### 1.3. 스프레드 연산자(`*.`) - 속성 추출의 축약형

List의 `collect` 메서드와 관련된 것은 스프레드 연산자(`*.`)로, 컬렉션에서 속성을 추출하는 간결한 방법을 제공합니다. 이는 기본적으로 일반적인 `collect` 패턴의 문법적 설탕입니다.

`collect.nf` 파일에 데모를 추가해 봅시다:

=== "변경 후"

    ```groovy title="collect.nf" linenums="1" hl_lines="15-18"
    def sample_ids = ['sample_001', 'sample_002', 'sample_003']

    // channel.collect() - 여러 채널 방출을 하나로 그룹화
    ch_input = channel.fromList(sample_ids)
    ch_input.view { sample -> "Individual channel item: ${sample}" }
    ch_collected = ch_input.collect()
    ch_collected.view { list -> "channel.collect() result: ${list} (${list.size()} items grouped into 1)" }

    // List.collect() - 각 요소 변환, 구조 유지
    def formatted_ids = sample_ids.collect { id ->
        id.toUpperCase().replace('SAMPLE_', 'SPECIMEN_')
    }
    println "List.collect() result: ${formatted_ids} (${sample_ids.size()} items transformed into ${formatted_ids.size()})"

    // 스프레드 연산자 - 간결한 속성 접근
    def sample_data = [[id: 's1', quality: 38.5], [id: 's2', quality: 42.1], [id: 's3', quality: 35.2]]
    def all_ids = sample_data*.id
    println "Spread operator result: ${all_ids}"
    ```

=== "변경 전"

    ```groovy title="collect.nf" linenums="1"
    def sample_ids = ['sample_001', 'sample_002', 'sample_003']

    // channel.collect() - 여러 채널 방출을 하나로 그룹화
    ch_input = channel.fromList(sample_ids)
    ch_input.view { sample -> "Individual channel item: ${sample}" }
    ch_collected = ch_input.collect()
    ch_collected.view { list -> "channel.collect() result: ${list} (${list.size()} items grouped into 1)" }

    // List.collect() - 각 요소 변환, 구조 유지
    def formatted_ids = sample_ids.collect { id ->
        id.toUpperCase().replace('SAMPLE_', 'SPECIMEN_')
    }
    println "List.collect() result: ${formatted_ids} (${sample_ids.size()} items transformed into ${formatted_ids.size()})"
    ```

업데이트된 워크플로를 실행합니다:

```bash title="스프레드 연산자 테스트"
nextflow run collect.nf
```

??? success "명령 출력"

    ```console hl_lines="6"
    N E X T F L O W   ~  version 25.10.2

    Launching `collect.nf` [cranky_galileo] DSL2 - revision: 5f3c8b2a91

    List.collect() result: [SPECIMEN_001, SPECIMEN_002, SPECIMEN_003] (3 items transformed into 3)
    Spread operator result: [s1, s2, s3]
    Individual channel item: sample_001
    Individual channel item: sample_002
    Individual channel item: sample_003
    channel.collect() result: [sample_001, sample_002, sample_003] (3 items grouped into 1)
    ```

스프레드 연산자 `*.`는 일반적인 collect 패턴의 축약형입니다:

```groovy
// 다음은 동등합니다:
def ids = samples*.id
def ids = samples.collect { it.id }

// 메서드 호출과도 작동합니다:
def names = files*.getName()
def names = files.collect { it.getName() }
```

스프레드 연산자는 객체 리스트에서 단일 속성을 추출해야 할 때 특히 유용합니다 - 전체 `collect` 클로저를 작성하는 것보다 더 읽기 쉽습니다.

!!! tip "스프레드 vs Collect 사용 시기"

    - **스프레드(`*.`) 사용**: 간단한 속성 접근: `samples*.id`, `files*.name`
    - **collect 사용**: 변환 또는 복잡한 로직: `samples.collect { it.id.toUpperCase() }`, `samples.collect { [it.id, it.quality > 40] }`

### 요약

이 섹션에서 배운 내용:

- **데이터플로 vs 스크립팅**: 채널 연산자는 파이프라인을 통한 데이터 흐름 방식을 오케스트레이션하고, 스크립팅은 개별 데이터 항목을 변환합니다
- **타입 이해**: 동일한 메서드 이름(`collect` 같은)이 호출되는 타입(Channel vs List)에 따라 다르게 작동할 수 있습니다
- **컨텍스트의 중요성**: 채널(데이터플로) 또는 데이터 구조(스크립팅) 중 무엇을 다루고 있는지 항상 인식해야 합니다

이러한 경계를 이해하는 것은 디버깅, 문서화, 유지보수 가능한 워크플로 작성에 필수적입니다.

다음으로는 실제 데이터 처리에 필수적인 문자열 처리 기능을 더 깊이 살펴보겠습니다.

---

## 2. 문자열 처리 및 동적 스크립트 생성

문자열 처리를 마스터하는 것은 취약한 워크플로와 견고한 파이프라인을 구분합니다. 이 섹션에서는 복잡한 파일 이름 파싱, 동적 스크립트 생성, 변수 보간을 다룹니다.

### 2.1. 패턴 매칭과 정규 표현식

생물정보학 파일은 종종 메타데이터를 인코딩하는 복잡한 명명 규칙을 가지고 있습니다. 정규 표현식으로 패턴 매칭을 사용하여 이를 자동으로 추출해 봅시다.

`main.nf` 워크플로로 돌아가서 파일 이름에서 추가 샘플 정보를 추출하는 패턴 매칭 로직을 추가합니다. 데이터셋의 FASTQ 파일은 `SAMPLE_001_S1_L001_R1_001.fastq.gz`와 같은 이름으로 Illumina 스타일 명명 규칙을 따릅니다. 이러한 이름은 암호처럼 보일 수 있지만, 실제로는 샘플 ID, 레인 번호, 읽기 방향과 같은 유용한 메타데이터를 인코딩합니다. 정규식 기능을 사용하여 이러한 이름을 파싱합니다.

기존 `main.nf` 워크플로를 다음과 같이 변경합니다:

=== "변경 후"

    ```groovy title="main.nf" linenums="4" hl_lines="10-21"
            .map { row ->
                // 데이터 변환을 위한 스크립팅
                def sample_meta = [
                    id: row.sample_id.toLowerCase(),
                    organism: row.organism,
                    tissue: row.tissue_type.replaceAll('_', ' ').toLowerCase(),
                    depth: row.sequencing_depth.toInteger(),
                    quality: row.quality_score.toDouble()
                ]
                def fastq_path = file(row.file_path)

                def m = (fastq_path.name =~ /^(.+)_S(\d+)_L(\d{3})_(R[12])_(\d{3})\.fastq(?:\.gz)?$/)
                def file_meta = m ? [
                    sample_num: m[0][2].toInteger(),
                    lane: m[0][3],
                    read: m[0][4],
                    chunk: m[0][5]
                ] : [:]

                def priority = sample_meta.quality > 40 ? 'high' : 'normal'
                return tuple(sample_meta + file_meta + [priority: priority], fastq_path)
            }
    ```

=== "변경 전"

    ```groovy title="main.nf" linenums="4" hl_lines="10-11"
            .map { row ->
                // 데이터 변환을 위한 스크립팅
                def sample_meta = [
                    id: row.sample_id.toLowerCase(),
                    organism: row.organism,
                    tissue: row.tissue_type.replaceAll('_', ' ').toLowerCase(),
                    depth: row.sequencing_depth.toInteger(),
                    quality: row.quality_score.toDouble()
                ]
                def priority = sample_meta.quality > 40 ? 'high' : 'normal'
                return tuple(sample_meta + [priority: priority], file(row.file_path))
            }
    ```

이는 핵심 **문자열 처리 개념**을 보여줍니다:

1. **정규 표현식 리터럴** `~/pattern/` 구문 사용 - 백슬래시를 이스케이프할 필요 없이 정규식 패턴 생성
2. **패턴 매칭** `=~` 연산자 사용 - 문자열을 정규식 패턴과 매칭 시도
3. **Matcher 객체** `[0][1]`, `[0][2]` 등으로 캡처 그룹 - `[0]`은 전체 매칭을 참조하고, `[1]`, `[2]` 등은 괄호 안의 캡처 그룹을 참조

정규식 패턴 `^(.+)_S(\d+)_L(\d{3})_(R[12])_(\d{3})\.fastq(?:\.gz)?$`을 분석해 봅시다:

| 패턴                | 매칭                                  | 캡처                         |
| ------------------- | ------------------------------------- | ---------------------------- |
| `^(.+)`             | 시작부터 샘플 이름                    | 그룹 1: 샘플 이름            |
| `_S(\d+)`           | 샘플 번호 `_S1`, `_S2` 등             | 그룹 2: 샘플 번호            |
| `_L(\d{3})`         | 레인 번호 `_L001`                     | 그룹 3: 레인 (3자리)         |
| `_(R[12])`          | 읽기 방향 `_R1` 또는 `_R2`            | 그룹 4: 읽기 방향            |
| `_(\d{3})`          | 청크 번호 `_001`                      | 그룹 5: 청크 (3자리)         |
| `\.fastq(?:\.gz)?$` | 파일 확장자 `.fastq` 또는 `.fastq.gz` | 캡처되지 않음 (?: 는 비캡처) |

이는 메타데이터를 자동으로 추출하기 위해 Illumina 스타일 명명 규칙을 파싱합니다.

수정된 워크플로를 실행합니다:

```bash title="패턴 매칭 테스트"
nextflow run main.nf
```

??? success "명령 출력"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [clever_pauling] DSL2 - revision: 605d2058b4

    [[id:sample_001, organism:human, tissue:liver, depth:30000000, quality:38.5, sample_num:1, lane:001, read:R1, chunk:001, priority:normal], /workspaces/training/side-quests/essential_scripting_patterns/data/sequences/SAMPLE_001_S1_L001_R1_001.fastq]
    [[id:sample_002, organism:mouse, tissue:brain, depth:25000000, quality:35.2, sample_num:2, lane:001, read:R1, chunk:001, priority:normal], /workspaces/training/side-quests/essential_scripting_patterns/data/sequences/SAMPLE_002_S2_L001_R1_001.fastq]
    [[id:sample_003, organism:human, tissue:kidney, depth:45000000, quality:42.1, sample_num:3, lane:001, read:R1, chunk:001, priority:high], /workspaces/training/side-quests/essential_scripting_patterns/data/sequences/SAMPLE_003_S3_L001_R1_001.fastq]
    ```

파일 이름에서 풍부해진 메타데이터가 표시됩니다.

### 2.2. 프로세스의 동적 스크립트 생성

프로세스 스크립트 블록은 본질적으로 셸에 전달되는 다중 라인 문자열입니다. **조건부 로직**(if/else, 삼항 연산자)을 사용하여 입력 특성에 따라 다양한 스크립트 문자열을 동적으로 생성할 수 있습니다. 이는 프로세스 정의를 복제하지 않고 단일 말단 대 쌍 말단 시퀀싱 읽기와 같은 다양한 입력 타입을 처리하는 데 필수적입니다.

이 패턴을 보여주는 프로세스를 워크플로에 추가해 봅시다. `modules/fastp.nf`를 열어 확인합니다:

```groovy title="modules/fastp.nf" linenums="1"
process FASTP {
    container 'community.wave.seqera.io/library/fastp:0.24.0--62c97b06e8447690'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*_trimmed*.fastq.gz"), emit: reads

    script:
    """
    fastp \\
        --in1 ${reads[0]} \\
        --in2 ${reads[1]} \\
        --out1 ${meta.id}_trimmed_R1.fastq.gz \\
        --out2 ${meta.id}_trimmed_R2.fastq.gz \\
        --json ${meta.id}.fastp.json \\
        --html ${meta.id}.fastp.html \\
        --thread $task.cpus
    """
}
```

이 프로세스는 FASTQ 파일을 입력으로 받아 `fastp` 도구를 실행하여 어댑터를 트리밍하고 저품질 읽기를 필터링합니다. 불행히도 이 프로세스를 작성한 사람은 예제 데이터셋에 있는 단일 말단 읽기를 고려하지 않았습니다. 워크플로에 추가하고 무슨 일이 일어나는지 봅시다:

먼저 `main.nf` 워크플로의 맨 첫 번째 줄에 모듈을 포함합니다:

```groovy title="main.nf" linenums="1"
include { FASTP } from './modules/fastp.nf'
```

그런 다음 `workflow` 블록을 수정하여 `ch_samples` 채널을 `FASTP` 프로세스에 연결합니다:

=== "변경 후"

    ```groovy title="main.nf" linenums="25" hl_lines="27"
    workflow {

        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map { row ->
                def sample_meta = [
                    id: row.sample_id.toLowerCase(),
                    organism: row.organism,
                    tissue: row.tissue_type.replaceAll('_', ' ').toLowerCase(),
                    depth: row.sequencing_depth.toInteger(),
                    quality: row.quality_score.toDouble()
                ]
                def fastq_path = file(row.file_path)

                def m = (fastq_path.name =~ /^(.+)_S(\d+)_L(\d{3})_(R[12])_(\d{3})\.fastq(?:\.gz)?$/)
                def file_meta = m ? [
                    sample_num: m[0][2].toInteger(),
                    lane: m[0][3],
                    read: m[0][4],
                    chunk: m[0][5]
                ] : [:]

                def priority = sample_meta.quality > 40 ? 'high' : 'normal'
                return tuple(sample_meta + file_meta + [priority: priority], fastq_path)
            }

        ch_fastp = FASTP(ch_samples)
    }
    ```

=== "변경 전"

    ```groovy title="main.nf" linenums="25" hl_lines="26"
    workflow {

        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map { row ->
                def sample_meta = [
                    id: row.sample_id.toLowerCase(),
                    organism: row.organism,
                    tissue: row.tissue_type.replaceAll('_', ' ').toLowerCase(),
                    depth: row.sequencing_depth.toInteger(),
                    quality: row.quality_score.toDouble()
                ]
                def fastq_path = file(row.file_path)

                def m = (fastq_path.name =~ /^(.+)_S(\d+)_L(\d{3})_(R[12])_(\d{3})\.fastq(?:\.gz)?$/)
                def file_meta = m ? [
                    sample_num: m[0][2].toInteger(),
                    lane: m[0][3],
                    read: m[0][4],
                    chunk: m[0][5]
                ] : [:]

                def priority = sample_meta.quality > 40 ? 'high' : 'normal'
                return [sample_meta + file_meta + [priority: priority], file(row.file_path)]
            }
            .view()
    }
    ```

이 수정된 워크플로를 실행합니다:

```bash
nextflow run main.nf
```

??? failure "명령 출력"

    ```console
    ERROR ~ Error executing process > 'FASTP (3)'

    Caused by:
      Process `FASTP (3)` terminated with an error exit status (255)


    Command executed:

      fastp \
          --in1 SAMPLE_003_S3_L001_R1_001.fastq \
          --in2 null \
          --out1 sample_003_trimmed_R1.fastq.gz \
          --out2 sample_003_trimmed_R2.fastq.gz \
          --json sample_003.fastp.json \
          --html sample_003.fastp.html \
          --thread 2

    Command exit status:
      255

    Command output:
      (empty)
    ```

프로세스가 두 번째 입력 파일에 대해 `null` 값으로 `fastp`를 실행하려고 시도하여 실패하는 것을 볼 수 있습니다. 이는 데이터셋에 단일 말단 읽기가 포함되어 있지만 프로세스가 쌍 말단 읽기(한 번에 두 개의 입력 파일)를 기대하도록 하드코딩되어 있기 때문입니다.

`FASTP` 프로세스 `script:` 블록에 조건부 로직을 추가하여 이를 수정합니다. if/else 문이 읽기 파일 개수를 확인하고 그에 따라 명령을 조정합니다.

=== "변경 후"

    ```groovy title="main.nf" linenums="10" hl_lines="3-27"
        script:
        // 간단한 단일 말단 vs 쌍 말단 감지
        def is_single = reads instanceof List ? reads.size() == 1 : true

        if (is_single) {
            def input_file = reads instanceof List ? reads[0] : reads
            """
            fastp \\
                --in1 ${input_file} \\
                --out1 ${meta.id}_trimmed.fastq.gz \\
                --json ${meta.id}.fastp.json \\
                --html ${meta.id}.fastp.html \\
                --thread $task.cpus
            """
