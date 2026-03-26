# 필수 Nextflow 스크립팅 패턴

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } AI 지원 번역 - [자세히 알아보기 및 개선 제안](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Nextflow는 Java Virtual Machine에서 실행되는 프로그래밍 언어입니다. Nextflow는 [Groovy](http://groovy-lang.org/)를 기반으로 하며 많은 문법을 공유하지만, 단순히 "확장된 Groovy"가 아닙니다. Nextflow는 완전히 명세된 [문법](https://nextflow.io/docs/latest/reference/syntax.html)과 [표준 라이브러리](https://nextflow.io/docs/latest/reference/stdlib.html)를 갖춘 독립적인 언어입니다.

변수, map, 리스트에 대한 기본 문법만으로도 Nextflow 코드를 많이 작성할 수 있습니다. 대부분의 Nextflow 튜토리얼은 워크플로우 오케스트레이션(채널, 프로세스, 데이터 흐름)에 초점을 맞추며, 그것만으로도 상당히 많은 것을 할 수 있습니다.

그러나 데이터를 조작하거나, 복잡한 파일명을 파싱하거나, 조건부 로직을 구현하거나, 견고한 프로덕션 워크플로우를 구축해야 할 때는 코드의 두 가지 측면을 구분하는 것이 도움이 됩니다. 바로 **데이터플로우**(채널, 연산자, 프로세스, 워크플로우)와 **스크립팅**(closure, 함수, 프로세스 스크립트 내부의 코드)입니다. 이 구분이 다소 임의적으로 보일 수 있지만—결국 모두 Nextflow 코드이므로—파이프라인을 오케스트레이션하는 시점과 데이터를 조작하는 시점을 이해하는 데 유용한 사고 모델을 제공합니다. 두 가지를 모두 숙달하면 명확하고 유지보수하기 쉬운 워크플로우를 작성하는 능력이 크게 향상됩니다.

### 학습 목표

이 사이드 퀘스트는 기본 개념부터 프로덕션 수준의 패턴까지 실습 중심으로 진행됩니다.
간단한 CSV 읽기 워크플로우를 정교한 바이오인포매틱스 파이프라인으로 발전시키며, 현실적인 과제를 단계적으로 해결합니다.

- **경계 이해하기:** 데이터플로우 연산과 스크립팅을 구분하고, 두 가지가 어떻게 함께 작동하는지 이해합니다
- **데이터 조작:** 강력한 연산자를 사용하여 map과 컬렉션을 추출, 변환, 부분 선택합니다
- **문자열 처리:** regex 패턴으로 복잡한 파일 명명 체계를 파싱하고 변수 보간을 마스터합니다
- **재사용 가능한 함수:** 복잡한 로직을 명명된 함수로 추출하여 더 깔끔하고 유지보수하기 쉬운 워크플로우를 만듭니다
- **동적 로직:** 다양한 입력 유형에 적응하는 프로세스를 구축하고 동적 리소스 할당을 위한 closure를 사용합니다
- **조건부 라우팅:** 메타데이터 특성에 따라 샘플을 다른 프로세스로 지능적으로 라우팅합니다
- **안전한 연산:** null-safe 연산자로 누락된 데이터를 우아하게 처리하고 명확한 오류 메시지로 입력을 검증합니다
- **설정 기반 핸들러:** 로깅, 알림, 라이프사이클 관리를 위한 워크플로우 이벤트 핸들러를 사용합니다

### 사전 요구사항

이 사이드 퀘스트를 시작하기 전에 다음을 갖추어야 합니다.

- [Hello Nextflow](../hello_nextflow/README.md) 튜토리얼 또는 동급의 입문 과정을 완료했어야 합니다.
- 기본적인 Nextflow 개념과 메커니즘(프로세스, 채널, 연산자, 파일 작업, 메타 데이터)에 익숙해야 합니다.
- 일반적인 프로그래밍 구조(변수, map, 리스트)에 대한 기본적인 이해가 있어야 합니다.

이 튜토리얼은 프로그래밍 개념을 접할 때마다 설명하므로, 광범위한 프로그래밍 경험이 없어도 됩니다.
기본 개념부터 시작하여 고급 패턴까지 단계적으로 학습합니다.

---

## 0. 시작하기

#### 교육 코드스페이스 열기

아직 열지 않았다면, [환경 설정](../envsetup/index.md)에 설명된 대로 교육 환경을 열어야 합니다.

[![GitHub Codespaces에서 열기](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

#### 프로젝트 디렉토리로 이동

이 튜토리얼의 파일이 있는 디렉토리로 이동합니다.

```bash
cd side-quests/essential_scripting_patterns
```

#### 자료 검토

메인 워크플로우 파일과 예제 데이터 파일이 포함된 `data` 디렉토리를 확인할 수 있습니다.

```console title="Directory contents"
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

샘플 CSV에는 특성에 따라 다른 처리가 필요한 생물학적 샘플에 대한 정보가 포함되어 있습니다.

```console title="samples.csv"
sample_id,organism,tissue_type,sequencing_depth,file_path,quality_score
SAMPLE_001,human,liver,30000000,data/sequences/SAMPLE_001_S1_L001_R1_001.fastq,38.5
SAMPLE_002,mouse,brain,25000000,data/sequences/SAMPLE_002_S2_L001_R1_001.fastq,35.2
SAMPLE_003,human,kidney,45000000,data/sequences/SAMPLE_003_S3_L001_R1_001.fastq,42.1
```

이 현실적인 데이터셋을 사용하여 실제 바이오인포매틱스 워크플로우에서 접하게 될 실용적인 프로그래밍 기법을 살펴봅니다.

<!-- TODO: Can we make this more domain-agnostic? -->

<!-- TODO: add an assignment statement? #### Review the assignment -->

#### 준비 체크리스트

시작할 준비가 되었나요?

- [ ] 이 과정의 목표와 사전 요구사항을 이해했습니다
- [ ] 코드스페이스가 실행 중입니다
- [ ] 작업 디렉토리를 적절히 설정했습니다
<!-- - [ ] I understand the assignment -->

모든 항목을 확인했다면 시작할 준비가 된 것입니다.

---

## 1. 데이터플로우 vs 스크립팅: 경계 이해하기

### 1.1. 무엇이 무엇인지 파악하기

Nextflow 워크플로우를 작성할 때, **데이터플로우**(데이터가 채널과 프로세스를 통해 이동하는 방식)와 **스크립팅**(데이터를 조작하고 결정을 내리는 코드)을 구분하는 것이 중요합니다. 두 가지가 어떻게 함께 작동하는지 보여주는 워크플로우를 구축합니다.

#### 1.1.1. 기본 Nextflow 워크플로우

CSV 파일을 읽는 간단한 워크플로우부터 시작합니다(`main.nf`에 이미 작성되어 있습니다).

```groovy title="main.nf" linenums="1"
workflow {
    ch_samples = channel.fromPath("./data/samples.csv")
        .splitCsv(header: true)
        .view()
}
```

`workflow` 블록은 파이프라인 구조를 정의하고, `channel.fromPath()`는 파일 경로에서 채널을 생성합니다. `.splitCsv()` 연산자는 CSV 파일을 처리하여 각 행을 map 데이터 구조로 변환합니다.

이 워크플로우를 실행하여 원시 CSV 데이터를 확인합니다.

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

이제 `.map()` 연산자를 사용하여 데이터를 변환하는 스크립팅을 추가합니다. 이 연산자는 이미 익숙할 것입니다. 이 연산자는 각 항목을 변환하는 코드를 작성할 수 있는 'closure'를 받습니다.

!!! note "참고"

    **closure**는 나중에 전달되고 실행될 수 있는 코드 블록입니다. 인라인으로 정의하는 함수라고 생각하면 됩니다. closure는 중괄호 `{ }`로 작성되며 매개변수를 받을 수 있습니다. Nextflow 연산자가 작동하는 방식의 핵심이며, Nextflow를 한동안 작성해왔다면 이미 사용하고 있었을 수도 있습니다!

map 연산이 어떻게 생겼는지 확인합니다.

=== "후"

    ```groovy title="main.nf" linenums="2" hl_lines="3-6"
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map { row ->
                return row
            }
            .view()
    ```

=== "전"

    ```groovy title="main.nf" linenums="2" hl_lines="3"
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .view()
    ```

이것이 첫 번째 **closure**입니다. 인자로 전달할 수 있는 익명 함수입니다(Python의 람다나 JavaScript의 화살표 함수와 유사합니다). closure는 Nextflow 연산자를 사용하는 데 필수적입니다.

closure `{ row -> return row }`는 매개변수 `row`를 받습니다(이름은 `item`, `sample` 등 무엇이든 될 수 있습니다).

`.map()` 연산자가 각 채널 항목을 처리할 때, 해당 항목을 closure에 전달합니다. 여기서 `row`는 한 번에 하나의 CSV 행을 보유합니다.

이 변경 사항을 적용하고 워크플로우를 실행합니다.

```bash
nextflow run main.nf
```

입력을 그대로 반환하고 있으므로 이전과 동일한 출력이 표시됩니다. 이는 map 연산자가 올바르게 작동하고 있음을 확인합니다. 이제 데이터 변환을 시작합니다.

#### 1.1.3. Map 데이터 구조 생성하기

이제 closure 내부에 **스크립팅** 로직을 작성하여 각 데이터 행을 변환합니다. 여기서는 데이터 흐름을 오케스트레이션하는 것이 아니라 개별 데이터 항목을 처리합니다.

=== "후"

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

=== "전"

    ```groovy title="main.nf" linenums="2" hl_lines="4"
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map { row ->
                return row
            }
            .view()
    ```

`sample_meta` map은 키-값 데이터 구조입니다(Python의 딕셔너리, JavaScript의 객체, Ruby의 해시와 유사합니다). 샘플 ID, 생물체, 조직 유형, 시퀀싱 깊이, 품질 점수 등 관련 정보를 저장합니다.

`.toLowerCase()`와 `.replaceAll()` 같은 문자열 조작 메서드를 사용하여 데이터를 정리하고, `.toInteger()`와 `.toDouble()` 같은 타입 변환 메서드를 사용하여 CSV의 문자열 데이터를 적절한 숫자 타입으로 변환합니다.

이 변경 사항을 적용하고 워크플로우를 실행합니다.

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

이제 삼항 연산자를 사용하여 데이터 값에 따라 결정을 내리는 스크립팅을 추가합니다.

다음과 같이 변경합니다.

=== "후"

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

=== "전"

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

삼항 연산자는 `조건 ? 참일_때_값 : 거짓일_때_값` 패턴을 따르는 if/else 문의 축약형입니다. 이 줄은 "품질이 40보다 크면 'high'를 사용하고, 그렇지 않으면 'normal'을 사용한다"는 의미입니다. 이와 관련된 **Elvis 연산자**(`?:`)는 값이 null이거나 비어 있을 때 기본값을 제공합니다. 이 패턴은 이 튜토리얼 후반부에서 살펴봅니다.

map 추가 연산자 `+`는 기존 map을 수정하는 것이 아니라 **새로운 map**을 생성합니다. 이 줄은 `sample_meta`의 모든 키-값 쌍과 새로운 `priority` 키를 포함하는 새 map을 생성합니다.

!!! Note "참고"

    closure에 전달된 map은 절대 수정하지 마세요. 항상 `+`를 사용하여 새로운 map을 생성하세요. Nextflow에서는 동일한 데이터가 여러 연산을 동시에 통과하는 경우가 많습니다. map을 직접 수정하면 다른 연산이 동일한 객체를 참조할 때 예측할 수 없는 부작용이 발생할 수 있습니다. 새로운 map을 생성하면 각 연산이 자체적인 깨끗한 복사본을 갖게 됩니다.

수정된 워크플로우를 실행합니다.

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

#### 1.1.5. `.subMap()`으로 Map 부분 선택하기

`+` 연산자가 map에 키를 추가하는 반면, 때로는 반대로 특정 키만 추출해야 할 때가 있습니다. `.subMap()` 메서드가 이에 적합합니다.

식별 필드만 포함하는 메타데이터의 단순화된 버전을 생성하는 줄을 추가합니다.

=== "후"

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

=== "전"

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

수정된 워크플로우를 실행합니다.

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

이는 `view()` 연산으로 표시된 전체 메타데이터와 `println`으로 출력한 추출된 부분 집합을 모두 보여줍니다.

`.subMap()` 메서드는 키 목록을 받아 해당 키만 포함하는 새 map을 반환합니다. 원본 map에 키가 없으면 결과에 포함되지 않습니다.

이는 다른 프로세스에 대해 다른 메타데이터 버전을 생성해야 할 때 특히 유용합니다. 일부 프로세스는 전체 메타데이터가 필요하고 다른 프로세스는 최소한의 식별 필드만 필요할 수 있습니다.

이제 println 문을 제거하여 워크플로우를 이전 상태로 복원합니다. 앞으로는 필요하지 않습니다.

!!! tip "팁: Map 연산 요약"

    - **키 추가**: `map1 + [new_key: value]` - 추가 키가 있는 새 map 생성
    - **키 추출**: `map1.subMap(['key1', 'key2'])` - 지정된 키만 있는 새 map 생성
    - **두 연산 모두 새 map을 생성합니다** - 원본 map은 변경되지 않습니다

#### 1.1.6. Map 결합 및 결과 반환하기

지금까지는 Nextflow 커뮤니티에서 'meta map'이라고 부르는 것만 반환했으며, 해당 메타데이터와 관련된 파일은 무시했습니다. 하지만 Nextflow 워크플로우를 작성한다면 그 파일로 무언가를 하고 싶을 것입니다.

풍부해진 메타데이터 map과 해당 파일 경로, 두 요소로 구성된 튜플을 포함하는 채널 구조를 출력합니다. 이는 프로세스에 데이터를 전달하는 Nextflow의 일반적인 패턴입니다.

=== "후"

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

=== "전"

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

이 변경 사항을 적용하고 워크플로우를 실행합니다.

```bash
nextflow run main.nf
```

??? success "명령 출력"

    ```console
    [[id:sample_001, organism:human, tissue:liver, depth:30000000, quality:38.5, priority:normal], /workspaces/training/side-quests/essential_scripting_patterns/data/sequences/SAMPLE_001_S1_L001_R1_001.fastq]
    [[id:sample_002, organism:mouse, tissue:brain, depth:25000000, quality:35.2, priority:normal], /workspaces/training/side-quests/essential_scripting_patterns/data/sequences/SAMPLE_002_S2_L001_R1_001.fastq]
    [[id:sample_003, organism:human, tissue:kidney, depth:45000000, quality:42.1, priority:high], /workspaces/training/side-quests/essential_scripting_patterns/data/sequences/SAMPLE_003_S3_L001_R1_001.fastq]
    ```

이 `[meta, file]` 튜플 구조는 메타데이터와 관련 파일을 모두 프로세스에 전달하는 Nextflow의 일반적인 패턴입니다.

!!! note "참고"

    **Map과 메타데이터**: Map은 Nextflow에서 메타데이터를 다루는 데 기본적입니다. 메타데이터 map 작업에 대한 자세한 설명은 [메타데이터 작업](../metadata/) 사이드 퀘스트를 참조하세요.

우리의 워크플로우는 핵심 패턴을 보여줍니다. **데이터플로우 연산**(`workflow`, `channel.fromPath()`, `.splitCsv()`, `.map()`, `.view()`)은 데이터가 파이프라인을 통해 이동하는 방식을 오케스트레이션하고, `.map()` closure 내부의 **스크립팅**(map `[key: value]`, 문자열 메서드, 타입 변환, 삼항 연산자)은 개별 데이터 항목의 변환을 처리합니다.

### 1.2. 다양한 타입 이해하기: 채널 vs 리스트

지금까지는 데이터플로우 연산과 스크립팅을 구분할 수 있었습니다. 그런데 동일한 메서드 이름이 두 컨텍스트 모두에 존재하는 경우는 어떨까요?

완벽한 예시는 `collect` 메서드입니다. Nextflow 표준 라이브러리에서 채널 타입과 List 타입 모두에 존재합니다. List의 `collect()` 메서드는 각 요소를 변환하는 반면, 채널의 `collect()` 연산자는 모든 채널 방출을 단일 항목 채널로 수집합니다.

샘플 데이터로 이를 시연합니다. 먼저 채널 `collect()` 연산자가 무엇을 하는지 복습합니다. `collect.nf`를 확인합니다.

```groovy title="collect.nf" linenums="1"
def sample_ids = ['sample_001', 'sample_002', 'sample_003']

// channel.collect() - 여러 채널 방출을 하나로 그룹화
ch_input = channel.fromList(sample_ids)
ch_input.view { sample -> "Individual channel item: ${sample}" }
ch_collected = ch_input.collect()
ch_collected.view { list -> "channel.collect() result: ${list} (${list.size()} items grouped into 1)" }
```

단계:

- 샘플 ID의 List 정의
- 각 샘플 ID를 별도로 방출하는 `fromList()`로 채널 생성
- 각 항목이 통과할 때 `view()`로 출력
- 채널의 `collect()` 연산자로 모든 항목을 단일 리스트로 수집
- 두 번째 `view()`로 수집된 결과(모든 샘플 ID를 포함하는 단일 항목) 출력

채널의 구조는 변경했지만 데이터 자체는 변경하지 않았습니다.

워크플로우를 실행하여 확인합니다.

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

`view()`는 모든 채널 방출에 대해 출력을 반환하므로, 이 단일 출력에 원래 3개 항목이 하나의 리스트로 그룹화되어 있음을 알 수 있습니다.

이제 List의 `collect` 메서드를 실제로 확인합니다. `collect.nf`를 수정하여 원래 샘플 ID 리스트에 List의 `collect` 메서드를 적용합니다.

=== "후"

    ```groovy title="main.nf" linenums="1" hl_lines="9-13"
    def sample_ids = ['sample_001', 'sample_002', 'sample_003']

    // channel.collect() - 여러 채널 방출을 하나로 그룹화
    ch_input = channel.fromList(sample_ids)
    ch_input.view { sample -> "Individual channel item: ${sample}" }
    ch_collected = ch_input.collect()
    ch_collected.view { list -> "channel.collect() result: ${list} (${list.size()} items grouped into 1)" }

    // List.collect() - 각 요소를 변환하고 구조를 유지
    def formatted_ids = sample_ids.collect { id ->
        id.toUpperCase().replace('SAMPLE_', 'SPECIMEN_')
    }
    println "List.collect() result: ${formatted_ids} (${sample_ids.size()} items transformed into ${formatted_ids.size()})"
    ```

=== "전"

    ```groovy title="main.nf" linenums="1"
    def sample_ids = ['sample_001', 'sample_002', 'sample_003']

    // channel.collect() - 여러 채널 방출을 하나로 그룹화
    ch_input = channel.fromList(sample_ids)
    ch_input.view { sample -> "Individual channel item: ${sample}" }
    ch_collected = ch_input.collect()
    ch_collected.view { list -> "channel.collect() result: ${list} (${list.size()} items grouped into 1)" }
    ```

이 새로운 코드에서:

- List의 `collect` 메서드를 사용하여 원래 리스트의 각 샘플 ID를 변환하는 새 변수 `formatted_ids`를 정의합니다
- `println`을 사용하여 결과를 출력합니다

수정된 워크플로우를 실행합니다.

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

이번에는 데이터의 구조를 변경하지 않았습니다. 리스트에 여전히 3개의 항목이 있지만, List의 `collect` 메서드를 사용하여 각 항목을 변환하여 수정된 값이 있는 새 리스트를 생성했습니다. 이는 채널에서 `map` 연산자를 사용하는 것과 유사하지만, 채널이 아닌 List 데이터 구조에서 작동합니다.

`collect`는 요점을 설명하기 위해 사용한 극단적인 예시입니다. 핵심 교훈은 워크플로우를 작성할 때 항상 **데이터 구조**(List, Map 등)와 **채널**(데이터플로우 구조)을 구분해야 한다는 것입니다. 연산은 이름을 공유할 수 있지만 호출되는 타입에 따라 완전히 다르게 동작합니다.

### 1.3. 스프레드 연산자(`*.`) - 속성 추출의 축약형

List의 `collect` 메서드와 관련된 것이 스프레드 연산자(`*.`)입니다. 이는 컬렉션에서 속성을 추출하는 간결한 방법을 제공합니다. 본질적으로 일반적인 `collect` 패턴의 문법적 설탕입니다.

`collect.nf` 파일에 시연을 추가합니다.

=== "후"

    ```groovy title="collect.nf" linenums="1" hl_lines="15-18"
    def sample_ids = ['sample_001', 'sample_002', 'sample_003']

    // channel.collect() - 여러 채널 방출을 하나로 그룹화
    ch_input = channel.fromList(sample_ids)
    ch_input.view { sample -> "Individual channel item: ${sample}" }
    ch_collected = ch_input.collect()
    ch_collected.view { list -> "channel.collect() result: ${list} (${list.size()} items grouped into 1)" }

    // List.collect() - 각 요소를 변환하고 구조를 유지
    def formatted_ids = sample_ids.collect { id ->
        id.toUpperCase().replace('SAMPLE_', 'SPECIMEN_')
    }
    println "List.collect() result: ${formatted_ids} (${sample_ids.size()} items transformed into ${formatted_ids.size()})"

    // 스프레드 연산자 - 간결한 속성 접근
    def sample_data = [[id: 's1', quality: 38.5], [id: 's2', quality: 42.1], [id: 's3', quality: 35.2]]
    def all_ids = sample_data*.id
    println "Spread operator result: ${all_ids}"
    ```

=== "전"

    ```groovy title="collect.nf" linenums="1"
    def sample_ids = ['sample_001', 'sample_002', 'sample_003']

    // channel.collect() - 여러 채널 방출을 하나로 그룹화
    ch_input = channel.fromList(sample_ids)
    ch_input.view { sample -> "Individual channel item: ${sample}" }
    ch_collected = ch_input.collect()
    ch_collected.view { list -> "channel.collect() result: ${list} (${list.size()} items grouped into 1)" }

    // List.collect() - 각 요소를 변환하고 구조를 유지
    def formatted_ids = sample_ids.collect { id ->
        id.toUpperCase().replace('SAMPLE_', 'SPECIMEN_')
    }
    println "List.collect() result: ${formatted_ids} (${sample_ids.size()} items transformed into ${formatted_ids.size()})"
    ```

업데이트된 워크플로우를 실행합니다.

```bash title="Test spread operator"
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

스프레드 연산자 `*.`는 일반적인 collect 패턴의 축약형입니다.

```groovy
// 다음 두 표현은 동일합니다:
def ids = samples*.id
def ids = samples.collect { it.id }

// 메서드 호출에도 작동합니다:
def names = files*.getName()
def names = files.collect { it.getName() }
```

스프레드 연산자는 객체 리스트에서 단일 속성을 추출해야 할 때 특히 유용합니다. 전체 `collect` closure를 작성하는 것보다 가독성이 높습니다.

!!! tip "팁: 스프레드 vs Collect 사용 시기"

    - **스프레드(`*.`) 사용**: 단순 속성 접근: `samples*.id`, `files*.name`
    - **collect 사용**: 변환이나 복잡한 로직: `samples.collect { it.id.toUpperCase() }`, `samples.collect { [it.id, it.quality > 40] }`

### 핵심 정리

이 섹션에서 다음을 학습했습니다.

- **데이터플로우 vs 스크립팅**: 채널 연산자는 데이터가 파이프라인을 통해 흐르는 방식을 오케스트레이션하고, 스크립팅은 개별 데이터 항목을 변환합니다
- **타입 이해하기**: 동일한 메서드 이름(예: `collect`)이 호출되는 타입(채널 vs List)에 따라 다르게 동작할 수 있습니다
- **컨텍스트가 중요합니다**: 채널(데이터플로우)로 작업하는지 데이터 구조(스크립팅)로 작업하는지 항상 인식하세요

이러한 경계를 이해하는 것은 디버깅, 문서화, 유지보수 가능한 워크플로우 작성에 필수적입니다.

다음으로 실제 데이터를 처리하는 데 필수적인 문자열 처리 기능을 더 깊이 살펴봅니다.

---

## 2. 문자열 처리와 동적 스크립트 생성

문자열 처리를 마스터하면 취약한 워크플로우와 견고한 파이프라인을 구분할 수 있습니다. 이 섹션에서는 복잡한 파일명 파싱, 동적 스크립트 생성, 변수 보간을 다룹니다.

### 2.1. 패턴 매칭과 정규 표현식

바이오인포매틱스 파일은 종종 메타데이터를 인코딩하는 복잡한 명명 규칙을 가집니다. 정규 표현식을 사용한 패턴 매칭으로 이를 자동으로 추출합니다.

`main.nf` 워크플로우로 돌아가서 파일명에서 추가 샘플 정보를 추출하는 패턴 매칭 로직을 추가합니다. 데이터셋의 FASTQ 파일은 `SAMPLE_001_S1_L001_R1_001.fastq.gz`와 같은 Illumina 스타일 명명 규칙을 따릅니다. 이 이름들이 암호처럼 보일 수 있지만, 실제로는 샘플 ID, 레인 번호, 읽기 방향과 같은 유용한 메타데이터를 인코딩합니다. regex 기능을 사용하여 이 이름들을 파싱합니다.

기존 `main.nf` 워크플로우에 다음과 같이 변경합니다.

=== "후"

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

=== "전"

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

이는 핵심 **문자열 처리 개념**을 보여줍니다.

1. **정규 표현식 리터럴** `~/pattern/` 문법 사용 - 백슬래시를 이스케이프할 필요 없이 regex 패턴을 생성합니다
2. **패턴 매칭** `=~` 연산자 사용 - 문자열을 regex 패턴과 매칭하려고 시도합니다
3. **매처 객체**가 `[0][1]`, `[0][2]` 등으로 그룹을 캡처합니다 - `[0]`은 전체 매치를 참조하고, `[1]`, `[2]` 등은 괄호 안의 캡처된 그룹을 참조합니다

regex 패턴 `^(.+)_S(\d+)_L(\d{3})_(R[12])_(\d{3})\.fastq(?:\.gz)?$`를 분석합니다.

| 패턴                | 매칭                                  | 캡처                             |
| ------------------- | ------------------------------------- | -------------------------------- |
| `^(.+)`             | 시작부터 샘플 이름                    | 그룹 1: 샘플 이름                |
| `_S(\d+)`           | 샘플 번호 `_S1`, `_S2` 등             | 그룹 2: 샘플 번호                |
| `_L(\d{3})`         | 레인 번호 `_L001`                     | 그룹 3: 레인 (3자리)             |
| `_(R[12])`          | 읽기 방향 `_R1` 또는 `_R2`            | 그룹 4: 읽기 방향                |
| `_(\d{3})`          | 청크 번호 `_001`                      | 그룹 5: 청크 (3자리)             |
| `\.fastq(?:\.gz)?$` | 파일 확장자 `.fastq` 또는 `.fastq.gz` | 캡처되지 않음 (?:는 비캡처 그룹) |

이는 Illumina 스타일 명명 규칙을 파싱하여 메타데이터를 자동으로 추출합니다.

수정된 워크플로우를 실행합니다.

```bash title="Test pattern matching"
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

파일명에서 풍부해진 메타데이터를 확인할 수 있습니다.

### 2.2. 프로세스에서의 동적 스크립트 생성

프로세스 script 블록은 본질적으로 셸에 전달되는 여러 줄 문자열입니다. **조건부 로직**(if/else, 삼항 연산자)을 사용하여 입력 특성에 따라 다른 스크립트 문자열을 동적으로 생성할 수 있습니다. 이는 프로세스 정의를 중복하지 않고 단일 엔드 vs 페어드 엔드 시퀀싱 읽기와 같은 다양한 입력 유형을 처리하는 데 필수적입니다.

이 패턴을 보여주는 프로세스를 워크플로우에 추가합니다. `modules/fastp.nf`를 열어 확인합니다.

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

이 프로세스는 FASTQ 파일을 입력으로 받아 `fastp` 도구를 실행하여 어댑터를 트리밍하고 저품질 읽기를 필터링합니다. 안타깝게도 이 프로세스를 작성한 사람은 예제 데이터셋에 있는 단일 엔드 읽기를 허용하지 않았습니다. 워크플로우에 추가하고 어떻게 되는지 확인합니다.

먼저 `main.nf` 워크플로우의 첫 번째 줄에 모듈을 포함합니다.

```groovy title="main.nf" linenums="1"
include { FASTP } from './modules/fastp.nf'
```

그런 다음 `workflow` 블록을 수정하여 `ch_samples` 채널을 `FASTP` 프로세스에 연결합니다.

=== "후"

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

=== "전"

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

수정된 워크플로우를 실행합니다.

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

프로세스가 두 번째 입력 파일에 `null` 값으로 `fastp`를 실행하려고 하여 실패하는 것을 볼 수 있습니다. 데이터셋에 단일 엔드 읽기가 포함되어 있지만, 프로세스는 페어드 엔드 읽기(한 번에 두 개의 입력 파일)를 기대하도록 하드코딩되어 있기 때문입니다.

`FASTP` 프로세스 `script:` 블록에 조건부 로직을 추가하여 이를 수정합니다. if/else 문이 읽기 파일 수를 확인하고 그에 따라 명령을 조정합니다.

=== "후"

    ```groovy title="main.nf" linenums="10" hl_lines="3-27"
        script:
        // 단일 엔드 vs 페어드 엔드 감지
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
        } else {
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

=== "전"

    ```groovy title="main.nf" linenums="10" hl_lines="2-11"
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

이제 워크플로우가 단일 엔드와 페어드 엔드 읽기를 모두 우아하게 처리할 수 있습니다. 조건부 로직이 입력 파일 수를 확인하고 `fastp`에 적합한 명령을 구성합니다. 작동하는지 확인합니다.

```bash
nextflow run main.nf
```

??? success "명령 출력"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [adoring_rosalind] DSL2 - revision: 04b1cd93e9

    executor >  local (3)
    [31/a8ad4d] process > FASTP (3) [100%] 3 of 3 ✔
    ```

잘 작동합니다! 실행된 실제 명령을 확인합니다(작업 해시에 맞게 조정하세요).

```console title="Check commands executed"
cat work/31/a8ad4d95749e685a6d842d3007957f/.command.sh
```

Nextflow가 단일 엔드 읽기에 올바른 명령을 선택했음을 확인할 수 있습니다.

```bash title=".command.sh"
#!/bin/bash -ue
fastp \
    --in1 SAMPLE_003_S3_L001_R1_001.fastq \
    --out1 sample_003_trimmed.fastq.gz \
    --json sample_003.fastp.json \
    --html sample_003.fastp.html \
    --thread 2
```

동적 스크립트 로직의 또 다른 일반적인 사용 예는 [Nextflow for Science Genomics 모듈](../../nf4science/genomics/02_joint_calling)에서 볼 수 있습니다. 해당 모듈에서 호출되는 GATK 프로세스는 여러 입력 파일을 받을 수 있지만, 올바른 명령줄을 구성하려면 각 파일 앞에 `-V`를 붙여야 합니다. 프로세스는 스크립팅을 사용하여 입력 파일 컬렉션(`all_gvcfs`)을 올바른 명령 인자로 변환합니다.

```groovy title="command line manipulation for GATK" linenums="1" hl_lines="2 5"
    script:
    def gvcfs_line = all_gvcfs.collect { gvcf -> "-V ${gvcf}" }.join(' ')
    """
    gatk GenomicsDBImport \
        ${gvcfs_line} \
        -L ${interval_list} \
        --genomicsdb-workspace-path ${cohort_name}_gdb
    """
```

프로세스 script 블록에서 스크립팅을 사용하는 이러한 패턴은 매우 강력하며 다양한 시나리오에 적용할 수 있습니다. 가변 입력 유형 처리부터 파일 컬렉션에서 복잡한 명령줄 인자 구축까지, 프로세스를 실제 데이터의 다양한 요구사항에 진정으로 적응 가능하게 만듭니다.

### 2.3. 변수 보간: Nextflow와 셸 변수

프로세스 스크립트는 Nextflow 변수, 셸 변수, 명령 치환을 혼합하며, 각각 다른 보간 문법을 사용합니다. 잘못된 문법을 사용하면 오류가 발생합니다. 처리 보고서를 생성하는 프로세스로 이를 살펴봅니다.

모듈 파일 `modules/generate_report.nf`를 확인합니다.

```groovy title="modules/generate_report.nf" linenums="1"
process GENERATE_REPORT {

    publishDir 'results/reports', mode: 'copy'

    input:
    tuple val(meta), path(reads)

    output:
    path "${meta.id}_report.txt"

    script:
    """
    echo "Processing ${reads}" > ${meta.id}_report.txt
    echo "Sample: ${meta.id}" >> ${meta.id}_report.txt
    """
}
```

이 프로세스는 샘플 ID와 파일명이 포함된 간단한 보고서를 작성합니다. 이제 다양한 유형의 변수를 혼합해야 할 때 어떻게 되는지 실행하여 확인합니다.

`main.nf`에 프로세스를 포함하고 워크플로우에 추가합니다.

=== "후"

    ```groovy title="main.nf" linenums="1" hl_lines="2 30"
    include { FASTP } from './modules/fastp.nf'
    include { GENERATE_REPORT } from './modules/generate_report.nf'

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
        GENERATE_REPORT(ch_samples)
    }
    ```

=== "전"

    ```groovy title="main.nf" linenums="1" hl_lines="1 10-29"
    include { FASTP } from './modules/fastp.nf'

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

워크플로우를 실행하고 `results/reports/`에 생성된 보고서를 확인합니다. 각 샘플에 대한 기본 정보가 포함되어 있어야 합니다.

<!-- TODO: add the run command -->

??? success "명령 출력"

    ```console
    <!-- TODO: output -->
    ```

처리가 언제 어디서 이루어졌는지에 대한 정보를 추가하려면 어떻게 해야 할까요? 프로세스를 수정하여 **셸** 변수와 명령 치환을 사용하여 현재 사용자, 호스트명, 날짜를 보고서에 포함합니다.

=== "후"

    ```groovy title="modules/generate_report.nf" linenums="10" hl_lines="5-7"
        script:
        """
        echo "Processing ${reads}" > ${meta.id}_report.txt
        echo "Sample: ${meta.id}" >> ${meta.id}_report.txt
        echo "Processed by: ${USER}" >> ${meta.id}_report.txt
        echo "Hostname: $(hostname)" >> ${meta.id}_report.txt
        echo "Date: $(date)" >> ${meta.id}_report.txt
        """
    ```

=== "전"

    ```groovy title="modules/generate_report.nf" linenums="10"
        script:
        """
        echo "Processing ${reads}" > ${meta.id}_report.txt
        echo "Sample: ${meta.id}" >> ${meta.id}_report.txt
        """
    ```

실행하면 오류가 발생합니다. Nextflow가 `${USER}`를 존재하지 않는 Nextflow 변수로 해석하려고 하기 때문입니다.

??? failure "명령 출력"

    ```console
    Error modules/generate_report.nf:15:27: `USER` is not defined
    │  15 |     echo "Processed by: ${USER}" >> ${meta.id}_report.txt
    ╰     |                           ^^^^

    ERROR ~ Script compilation failed
    ```

Bash가 처리할 수 있도록 이스케이프해야 합니다.

백슬래시(`\`)로 셸 변수와 명령 치환을 이스케이프하여 수정합니다.

=== "후"

    ```groovy title="modules/generate_report.nf" linenums="10" hl_lines="5-7"
        script:
        """
        echo "Processing ${reads}" > ${meta.id}_report.txt
        echo "Sample: ${meta.id}" >> ${meta.id}_report.txt
        echo "Processed by: \${USER}" >> ${meta.id}_report.txt
        echo "Hostname: \$(hostname)" >> ${meta.id}_report.txt
        echo "Date: \$(date)" >> ${meta.id}_report.txt
        """
    ```

=== "전"

    ```groovy title="modules/generate_report.nf" linenums="10"
        script:
        """
        echo "Processing ${reads}" > ${meta.id}_report.txt
        echo "Sample: ${meta.id}" >> ${meta.id}_report.txt
        echo "Processed by: ${USER}" >> ${meta.id}_report.txt
        echo "Hostname: $(hostname)" >> ${meta.id}_report.txt
        echo "Date: $(date)" >> ${meta.id}_report.txt
        """
    ```

이제 작동합니다! 백슬래시(`\`)는 Nextflow에게 "이것을 해석하지 말고 Bash에 그대로 전달하라"고 알립니다.

### 핵심 정리

이 섹션에서 **문자열 처리** 기법을 학습했습니다.

- **파일 파싱을 위한 정규 표현식**: `=~` 연산자와 regex 패턴(`~/pattern/`)을 사용하여 복잡한 파일 명명 규칙에서 메타데이터를 추출합니다
- **동적 스크립트 생성**: 조건부 로직(if/else, 삼항 연산자)을 사용하여 입력 특성에 따라 다른 스크립트 문자열을 생성합니다
- **변수 보간**: Nextflow가 문자열을 해석하는 시점과 셸이 해석하는 시점을 이해합니다
  - `${var}` - Nextflow 변수 (워크플로우 컴파일 시 Nextflow가 보간)
  - `\${var}` - 셸 환경 변수 (이스케이프되어 런타임에 bash에 전달)
  - `\$(cmd)` - 셸 명령 치환 (이스케이프되어 런타임에 bash가 실행)

이러한 문자열 처리 및 생성 패턴은 실제 바이오인포매틱스 워크플로우에서 접하게 될 다양한 파일 형식과 명명 규칙을 처리하는 데 필수적입니다.

---

## 3. 재사용 가능한 함수 만들기

채널 연산자나 프로세스 정의에 인라인으로 복잡한 워크플로우 로직을 작성하면 가독성과 유지보수성이 떨어집니다. **함수**를 사용하면 이 로직을 명명된 재사용 가능한 컴포넌트로 추출할 수 있습니다.

map 연산이 길고 복잡해졌습니다. `def` 키워드를 사용하여 재사용 가능한 함수로 추출합니다.

기존 워크플로우에서 어떻게 보이는지 설명하기 위해, `def`를 사용하여 `separateMetadata`라는 재사용 가능한 함수를 정의하는 아래 수정 사항을 적용합니다.

=== "후"

    ```groovy title="main.nf" linenums="1" hl_lines="4-24 29"
    include { FASTP } from './modules/fastp.nf'
    include { GENERATE_REPORT } from './modules/generate_report.nf'

    def separateMetadata(row) {
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

    workflow {
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map{ row -> separateMetadata(row) }

        ch_fastp = FASTP(ch_samples)
        GENERATE_REPORT(ch_samples)
    }
    ```

=== "전"

    ```groovy title="main.nf" linenums="1" hl_lines="7-27"
    include { FASTP } from './modules/fastp.nf'
    include { GENERATE_REPORT } from './modules/generate_report.nf'

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
        GENERATE_REPORT(ch_samples)
    }
    ```

이 로직을 함수로 추출함으로써 실제 워크플로우 로직이 훨씬 깔끔해졌습니다.

```groovy title="minimal workflow"
    ch_samples = channel.fromPath("./data/samples.csv")
        .splitCsv(header: true)
        .map{ row -> separateMetadata(row) }

    ch_fastp = FASTP(ch_samples)
    GENERATE_REPORT(ch_samples)
```

이렇게 하면 워크플로우 로직을 한눈에 훨씬 쉽게 읽고 이해할 수 있습니다. `separateMetadata` 함수는 메타데이터를 파싱하고 풍부하게 만드는 모든 복잡한 로직을 캡슐화하여 재사용 가능하고 테스트 가능하게 만듭니다.

워크플로우를 실행하여 여전히 작동하는지 확인합니다.

```bash
nextflow run main.nf
```

??? success "명령 출력"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [admiring_panini] DSL2 - revision: 8cc832e32f

    executor >  local (6)
    [8c/2e3f91] process > FASTP (3)           [100%] 3 of 3 ✔
    [7a/1b4c92] process > GENERATE_REPORT (3) [100%] 3 of 3 ✔
    ```

출력에 두 프로세스가 성공적으로 완료된 것이 표시되어야 합니다. 워크플로우가 훨씬 깔끔해졌으며, 복잡한 메타데이터 처리 로직이 모두 `separateMetadata` 함수에 캡슐화되어 유지보수하기 쉬워졌습니다.

### 핵심 정리

이 섹션에서 **함수 생성**을 학습했습니다.

- **`def`로 함수 정의하기**: 명명된 함수를 생성하는 키워드입니다(Python의 `def`나 JavaScript의 `function`과 유사합니다)
- **함수 범위**: 스크립트 수준에서 정의된 함수는 Nextflow 워크플로우 전체에서 접근 가능합니다
- **반환 값**: 함수는 마지막 표현식을 자동으로 반환하거나 명시적 `return`을 사용합니다
- **더 깔끔한 코드**: 복잡한 로직을 함수로 추출하는 것은 모든 언어에서 기본적인 소프트웨어 엔지니어링 관행입니다

다음으로 동적 리소스 할당을 위한 프로세스 지시문에서 closure를 사용하는 방법을 살펴봅니다.

---

## 4. closure를 사용한 동적 리소스 지시문

지금까지는 프로세스의 `script` 블록에서 스크립팅을 사용했습니다. 하지만 **closure**(섹션 1.1에서 소개됨)는 프로세스 지시문에서도 매우 유용합니다. 특히 동적 리소스 할당에 유용합니다. 샘플 특성에 따라 적응하는 리소스 지시문을 FASTP 프로세스에 추가합니다.

### 4.1. 샘플별 리소스 할당

현재 FASTP 프로세스는 기본 리소스를 사용합니다. 고심도 샘플에 더 많은 CPU를 할당하여 더 스마트하게 만듭니다. `modules/fastp.nf`를 편집하여 동적 `cpus` 지시문과 정적 `memory` 지시문을 포함합니다.

=== "후"

    ```groovy title="modules/fastp.nf" linenums="1" hl_lines="4-5"
    process FASTP {
        container 'community.wave.seqera.io/library/fastp:0.24.0--62c97b06e8447690'

        cpus { meta.depth > 40000000 ? 2 : 1 }
        memory 2.GB

        input:
        tuple val(meta), path(reads)
    ```

=== "전"

    ```groovy title="modules/fastp.nf" linenums="1"
    process FASTP {
        container 'community.wave.seqera.io/library/fastp:0.24.0--62c97b06e8447690'

        input:
        tuple val(meta), path(reads)
    ```

closure `{ meta.depth > 40000000 ? 2 : 1 }`는 **삼항 연산자**(섹션 1.1에서 다룸)를 사용하며 각 작업에 대해 평가되어 샘플별 리소스 할당이 가능합니다. 고심도 샘플(>4천만 읽기)은 CPU 2개를 받고, 나머지는 CPU 1개를 받습니다.

!!! note "참고: 지시문에서 입력 변수 접근하기"

    Nextflow가 각 작업 실행 컨텍스트에서 이 closure를 평가하기 때문에, closure는 모든 입력 변수(여기서는 `meta`)에 접근할 수 있습니다.

작업 해시를 더 쉽게 볼 수 있도록 `-ansi-log false` 옵션으로 워크플로우를 다시 실행합니다.

```bash
nextflow run main.nf -ansi-log false
```

??? success "명령 출력"

    ```console
    N E X T F L O W  ~  version 25.10.2
    Launching `main.nf` [fervent_albattani] DSL2 - revision: fa8f249759
    [bd/ff3d41] Submitted process > FASTP (2)
    [a4/a3aab2] Submitted process > FASTP (1)
    [48/6db0c9] Submitted process > FASTP (3)
    [ec/83439d] Submitted process > GENERATE_REPORT (3)
    [bd/15d7cc] Submitted process > GENERATE_REPORT (2)
    [42/699357] Submitted process > GENERATE_REPORT (1)
    ```

특정 작업에 대해 실행된 실제 `docker` 명령을 확인하여 CPU 할당을 볼 수 있습니다.

```console title="Check docker command"
cat work/48/6db0c9e9d8aa65e4bb4936cd3bd59e/.command.run | grep "docker run"
```

다음과 같은 내용이 표시됩니다.

```bash title="docker command"
    docker run -i --cpu-shares 4096 --memory 2048m -e "NXF_TASK_WORKDIR" -v /workspaces/training/side-quests/essential_scripting_patterns:/workspaces/training/side-quests/essential_scripting_patterns -w "$NXF_TASK_WORKDIR" --name $NXF_BOXID community.wave.seqera.io/library/fastp:0.24.0--62c97b06e8447690 /bin/bash -ue /workspaces/training/side-quests/essential_scripting_patterns/work/48/6db0c9e9d8aa65e4bb4936cd3bd59e/.command.sh
```

이 예시에서는 고심도 샘플이었기 때문에 CPU 2개를 요청했습니다(`--cpu-shares 2048`). 샘플 심도에 따라 다른 CPU 할당이 표시될 것입니다. 다른 작업에 대해서도 확인해 보세요.

### 4.2. 재시도 전략

또 다른 강력한 패턴은 재시도 전략을 위한 `task.attempt` 사용입니다. 이것이 왜 유용한지 보여주기 위해, 먼저 FASTP에 필요한 것보다 적은 메모리를 할당합니다. `modules/fastp.nf`의 `memory` 지시문을 `1.GB`로 변경합니다.

=== "후"

    ```groovy title="modules/fastp.nf" linenums="1" hl_lines="5"
    process FASTP {
        container 'community.wave.seqera.io/library/fastp:0.24.0--62c97b06e8447690'

        cpus { meta.depth > 40000000 ? 4 : 2 }
        memory 1.GB

        input:
        tuple val(meta), path(reads)
    ```

=== "전"

    ```groovy title="modules/fastp.nf" linenums="1" hl_lines="5"
    process FASTP {
        container 'community.wave.seqera.io/library/fastp:0.24.0--62c97b06e8447690'

        cpus { meta.depth > 40000000 ? 4 : 2 }
        memory 2.GB

        input:
        tuple val(meta), path(reads)
    ```

... 그리고 워크플로우를 다시 실행합니다.

```bash
nextflow run main.nf
```

??? failure "명령 출력"

    ```console hl_lines="2 11"
    Command exit status:
      137

    Command output:
      (empty)

    Command error:
      Detecting adapter sequence for read1...
      No adapter detected for read1

      .command.sh: line 7:   101 Killed                  fastp --in1 SAMPLE_002_S2_L001_R1_001.fastq --out1 sample_002_trimmed.fastq.gz --json sample_002.fastp.json --html sample_002.fastp.html --thread 2
    ```

이는 프로세스가 메모리 한도를 초과하여 종료되었음을 나타냅니다.

이는 실제 워크플로우에서 매우 일반적인 시나리오입니다. 실행하기 전까지는 작업에 얼마나 많은 메모리가 필요한지 알 수 없는 경우가 있습니다.

워크플로우를 더 견고하게 만들기 위해, Groovy closure를 사용하여 각 시도마다 메모리 할당을 늘리는 재시도 전략을 구현할 수 있습니다. `memory` 지시문을 수정하여 기본 메모리에 `task.attempt`를 곱하고, `errorStrategy 'retry'`와 `maxRetries 2` 지시문을 추가합니다.

=== "후"

    ```groovy title="modules/fastp.nf" linenums="1" hl_lines="5-7"
    process FASTP {
        container 'community.wave.seqera.io/library/fastp:0.24.0--62c97b06e8447690'

        cpus { meta.depth > 40000000 ? 4 : 2 }
        memory { 1.GB * task.attempt }
        errorStrategy 'retry'
        maxRetries 2

        input:
        tuple val(meta), path(reads)
    ```

=== "전"

    ```groovy title="modules/fastp.nf" linenums="1" hl_lines="5"
    process FASTP {
        container 'community.wave.seqera.io/library/fastp:0.24.0--62c97b06e8447690'

        cpus { meta.depth > 40000000 ? 4 : 2 }
        memory 2.GB

        input:
        tuple val(meta), path(reads)
    ```

이제 메모리 부족으로 프로세스가 실패하면 Nextflow가 더 많은 메모리로 재시도합니다.

- 첫 번째 시도: 1 GB (task.attempt = 1)
- 두 번째 시도: 2 GB (task.attempt = 2)

... `maxRetries` 한도까지 계속됩니다.

### 핵심 정리

closure를 사용한 동적 지시문으로 다음을 할 수 있습니다.

- 입력 특성에 따라 리소스 할당
- 리소스를 늘리는 자동 재시도 전략 구현
- 여러 요소(메타데이터, 시도 횟수, 우선순위) 결합
- 복잡한 리소스 계산을 위한 조건부 로직 사용

이렇게 하면 워크플로우가 더 효율적(과도한 할당 방지)이고 더 견고(더 많은 리소스로 자동 재시도)해집니다.

---

## 5. 조건부 로직과 프로세스 제어

이전에는 `.map()`과 스크립팅을 사용하여 채널 데이터를 변환했습니다. 이제 조건부 로직을 사용하여 데이터에 따라 어떤 프로세스가 실행될지 제어합니다. 이는 다양한 샘플 유형에 적응하는 유연한 워크플로우에 필수적입니다.

Nextflow의 [데이터플로우 연산자](https://www.nextflow.io/docs/latest/reference/operator.html)는 런타임에 평가되는 closure를 받아, 채널 내용에 따라 워크플로우 결정을 내리는 조건부 로직을 가능하게 합니다.

### 5.1. `.branch()`로 라우팅하기

예를 들어, 시퀀싱 샘플이 특정 임계값 이상의 커버리지를 가진 인간 샘플인 경우에만 FASTP로 트리밍해야 한다고 가정합니다. 마우스 샘플이나 저커버리지 샘플은 대신 Trimgalore로 실행해야 합니다(이는 가상의 예시이지만 요점을 잘 설명합니다).

`modules/trimgalore.nf`에 간단한 Trimgalore 프로세스를 제공했습니다. 원하면 확인해 보세요. 세부 사항은 이 실습에서 중요하지 않습니다. 핵심은 메타데이터에 따라 샘플을 라우팅하고 싶다는 것입니다.

`modules/trimgalore.nf`에서 새 모듈을 포함합니다.

=== "후"

    ```groovy title="main.nf" linenums="1" hl_lines="2"
    include { FASTP } from './modules/fastp.nf'
    include { TRIMGALORE } from './modules/trimgalore.nf'
    ```

=== "전"

    ```groovy title="main.nf" linenums="1"
    include { FASTP } from './modules/fastp.nf'
    ```

... 그런 다음 `main.nf` 워크플로우를 수정하여 메타데이터에 따라 샘플을 분기하고 적절한 트리밍 프로세스로 라우팅합니다.

=== "후"

    ```groovy title="main.nf" linenums="28" hl_lines="5-12"
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map { row -> separateMetadata(row) }

        trim_branches = ch_samples
            .branch { meta, reads ->
                fastp: meta.organism == 'human' && meta.depth >= 30000000
                trimgalore: true
            }

        ch_fastp = FASTP(trim_branches.fastp)
        ch_trimgalore = TRIMGALORE(trim_branches.trimgalore)
        GENERATE_REPORT(ch_samples)
    ```

=== "전"

    ```groovy title="main.nf" linenums="28" hl_lines="5"
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map { row -> separateMetadata(row) }

        ch_fastp = FASTP(ch_samples)
        GENERATE_REPORT(ch_samples)
    ```

수정된 워크플로우를 실행합니다.

```bash
nextflow run main.nf
```

??? success "명령 출력"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [adoring_galileo] DSL2 - revision: c9e83aaef1

    executor >  local (6)
    [1d/0747ac] process > FASTP (2)           [100%] 2 of 2 ✔
    [cc/c44caf] process > TRIMGALORE (1)      [100%] 1 of 1 ✔
    [34/bd5a9f] process > GENERATE_REPORT (1) [100%] 3 of 3 ✔
    ```

여기서는 `.branch{}` 연산자 내부의 작지만 강력한 조건부 표현식을 사용하여 메타데이터에 따라 샘플을 라우팅했습니다. 고커버리지 인간 샘플은 `FASTP`를 통과하고, 다른 모든 샘플은 `TRIMGALORE`를 통과합니다.

### 5.2. 진리값(Truthiness)을 사용한 `.filter()`

워크플로우 실행을 제어하는 또 다른 강력한 패턴은 `.filter()` 연산자입니다. 이 연산자는 closure를 사용하여 파이프라인을 계속 진행할 항목을 결정합니다. filter closure 내부에서 어떤 항목이 통과할지 결정하는 **boolean 표현식**을 작성합니다.

Nextflow(많은 동적 언어와 마찬가지로)는 boolean 컨텍스트에서 어떤 값이 `true` 또는 `false`로 평가되는지를 결정하는 **"진리값(truthiness)"** 개념을 가지고 있습니다.

- **참(Truthy)**: null이 아닌 값, 비어 있지 않은 문자열, 0이 아닌 숫자, 비어 있지 않은 컬렉션
- **거짓(Falsy)**: `null`, 빈 문자열 `""`, 0, 빈 컬렉션 `[]` 또는 `[:]`, `false`

이는 `meta.id` 단독으로(명시적 `!= null` 없이) ID가 존재하고 비어 있지 않은지 확인한다는 것을 의미합니다. 이를 사용하여 품질 요구사항을 충족하지 않는 샘플을 필터링합니다.

분기 연산 전에 다음을 추가합니다.

=== "후"

    ```groovy title="main.nf" linenums="28" hl_lines="5-11"
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map { row -> separateMetadata(row) }

        // 유효하지 않거나 저품질 샘플 필터링
        ch_valid_samples = ch_samples
            .filter { meta, reads ->
                meta.id && meta.organism && meta.depth >= 25000000
            }

        trim_branches = ch_valid_samples
            .branch { meta, reads ->
                fastp: meta.organism == 'human' && meta.depth >= 30000000
                trimgalore: true
            }
    ```

=== "전"

    ```groovy title="main.nf" linenums="28" hl_lines="5"
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map { row -> separateMetadata(row) }

        trim_branches = ch_samples
            .branch { meta, reads ->
                fastp: meta.organism == 'human' && meta.depth >= 30000000
                trimgalore: true
            }
    ```

워크플로우를 다시 실행합니다.

```bash
nextflow run main.nf
```

??? success "명령 출력"

    ```console
    N E X T F L O W  ~  version 25.10.2
    Launching `main.nf` [lonely_williams] DSL2 - revision: d0b3f121ec
    [94/b48eac] Submitted process > FASTP (2)
    [2c/d2b28f] Submitted process > GENERATE_REPORT (2)
    [65/2e3be4] Submitted process > GENERATE_REPORT (1)
    [94/b48eac] NOTE: Process `FASTP (2)` terminated with an error exit status (137) -- Execution is retried (1)
    [3e/0d8664] Submitted process > TRIMGALORE (1)
    [6a/9137b0] Submitted process > FASTP (1)
    [6a/9137b0] NOTE: Process `FASTP (1)` terminated with an error exit status (137) -- Execution is retried (1)
    [83/577ac0] Submitted process > GENERATE_REPORT (3)
    [a2/5117de] Re-submitted process > FASTP (1)
    [1f/a1a4ca] Re-submitted process > FASTP (2)
    ```

일부 샘플을 제외하는 필터를 선택했기 때문에 더 적은 작업이 실행되었습니다.

필터 표현식 `meta.id && meta.organism && meta.depth >= 25000000`은 진리값과 명시적 비교를 결합합니다.

- `meta.id && meta.organism`은 두 필드가 존재하고 비어 있지 않은지 확인합니다(진리값 사용)
- `meta.depth >= 25000000`은 명시적 비교로 충분한 시퀀싱 심도를 보장합니다

!!! note "참고: 실제 진리값 활용"

    `meta.id && meta.organism` 표현식은 다음보다 더 간결합니다.
    ```groovy
    meta.id != null && meta.id != '' && meta.organism != null && meta.organism != ''
    ```

    이렇게 하면 필터링 로직이 훨씬 깔끔하고 읽기 쉬워집니다.

### 핵심 정리

이 섹션에서는 `.branch{}`와 `.filter{}`와 같은 Nextflow 연산자의 closure 인터페이스를 사용하여 워크플로우 실행을 제어하는 조건부 로직을 사용하는 방법을 학습했습니다. 진리값을 활용하여 간결한 조건부 표현식을 작성했습니다.

파이프라인이 이제 샘플을 적절한 프로세스로 지능적으로 라우팅하지만, 프로덕션 워크플로우는 유효하지 않은 데이터를 우아하게 처리해야 합니다. 누락되거나 null인 값에 대해 워크플로우를 견고하게 만듭니다.

---

## 6. 안전한 탐색과 Elvis 연산자

`separateMetadata` 함수는 현재 모든 CSV 필드가 존재하고 유효하다고 가정합니다. 하지만 불완전한 데이터가 있으면 어떻게 될까요?

### 6.1. 문제: 존재하지 않는 속성 접근하기

선택적 시퀀싱 실행 정보에 대한 지원을 추가하려고 한다고 가정합니다. 일부 실험실에서는 샘플에 시퀀싱 실행 ID나 배치 번호에 대한 추가 필드가 있을 수 있지만, 현재 CSV에는 이 열이 없습니다. 어쨌든 접근해 봅니다.

`separateMetadata` 함수를 수정하여 run_id 필드를 포함합니다.

=== "후"

    ```groovy title="main.nf" linenums="5" hl_lines="9"
    def separateMetadata(row) {
        def sample_meta = [
            id: row.sample_id.toLowerCase(),
            organism: row.organism,
            tissue: row.tissue_type.replaceAll('_', ' ').toLowerCase(),
            depth: row.sequencing_depth.toInteger(),
            quality: row.quality_score.toDouble()
        ]
        def run_id = row.run_id.toUpperCase()
    ```

=== "전"

    ```groovy title="main.nf" linenums="5"
    def separateMetadata(row) {
        def sample_meta = [
            id: row.sample_id.toLowerCase(),
            organism: row.organism,
            tissue: row.tissue_type.replaceAll('_', ' ').toLowerCase(),
            depth: row.sequencing_depth.toInteger(),
            quality: row.quality_score.toDouble()
        ]
    ```

워크플로우를 실행합니다.

```bash
nextflow run main.nf
```

??? failure "명령 출력"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [trusting_torvalds] DSL2 - revision: b56fbfbce2

    ERROR ~ Cannot invoke method toUpperCase() on null object

    -- Check script 'main.nf' at line: 13 or see '.nextflow.log' file for more details
    ```

NullPointerException으로 충돌합니다.

문제는 `run_id` 열이 CSV에 없기 때문에 `row.run_id`가 `null`을 반환한다는 것입니다. `null`에 `.toUpperCase()`를 호출하려고 하면 충돌합니다. 여기서 안전한 탐색 연산자가 도움이 됩니다.

### 6.2. 안전한 탐색 연산자(`?.`)

안전한 탐색 연산자(`?.`)는 `null` 값에 호출될 때 예외를 던지는 대신 `null`을 반환합니다. `?.` 앞의 객체가 `null`이면 메서드를 실행하지 않고 전체 표현식이 `null`로 평가됩니다.

안전한 탐색을 사용하도록 함수를 업데이트합니다.

=== "후"

    ```groovy title="main.nf" linenums="4" hl_lines="9"
    def separateMetadata(row) {
        def sample_meta = [
            id: row.sample_id.toLowerCase(),
            organism: row.organism,
            tissue: row.tissue_type.replaceAll('_', ' ').toLowerCase(),
            depth: row.sequencing_depth.toInteger(),
            quality: row.quality_score.toDouble()
        ]
        def run_id = row.run_id?.toUpperCase()
    ```

=== "전"

    ```groovy title="main.nf" linenums="4" hl_lines="9"
    def separateMetadata(row) {
        def sample_meta = [
            id: row.sample_id.toLowerCase(),
            organism: row.organism,
            tissue: row.tissue_type.replaceAll('_', ' ').toLowerCase(),
            depth: row.sequencing_depth.toInteger(),
            quality: row.quality_score.toDouble()
        ]
        def run_id = row.run_id.toUpperCase()
    ```

다시 실행합니다.

```bash
nextflow run main.nf
```

??? success "명령 출력"

    ```console
    <!-- TODO: output -->
    ```

충돌이 없습니다! 워크플로우가 이제 누락된 필드를 우아하게 처리합니다. `row.run_id`가 `null`이면 `?.` 연산자가 `.toUpperCase()` 호출을 방지하고, `run_id`는 예외를 발생시키는 대신 `null`이 됩니다.

### 6.3. 기본값을 위한 Elvis 연산자(`?:`)

Elvis 연산자(`?:`)는 왼쪽이 "거짓(falsy)"일 때(앞서 설명한 대로) 기본값을 제공합니다. 옆으로 보면 `?:`가 Elvis Presley의 유명한 헤어스타일과 눈처럼 보이기 때문에 이런 이름이 붙었습니다!

이제 안전한 탐색을 사용하고 있으므로, 해당 필드가 없는 샘플의 경우 `run_id`가 `null`이 됩니다. Elvis 연산자를 사용하여 기본값을 제공하고 `sample_meta` map에 추가합니다.

=== "후"

    ```groovy title="main.nf" linenums="5" hl_lines="9-10"
    def separateMetadata(row) {
        def sample_meta = [
            id: row.sample_id.toLowerCase(),
            organism: row.organism,
            tissue: row.tissue_type.replaceAll('_', ' ').toLowerCase(),
            depth: row.sequencing_depth.toInteger(),
            quality: row.quality_score.toDouble()
        ]
        def run_id = row.run_id?.toUpperCase() ?: 'UNSPECIFIED'
        sample_meta.run = run_id
    ```

=== "전"

    ```groovy title="main.nf" linenums="5" hl_lines="9"
    def separateMetadata(row) {
        def sample_meta = [
            id: row.sample_id.toLowerCase(),
            organism: row.organism,
            tissue: row.tissue_type.replaceAll('_', ' ').toLowerCase(),
            depth: row.sequencing_depth.toInteger(),
            quality: row.quality_score.toDouble()
        ]
        def run_id = row.run_id?.toUpperCase()
    ```

결과를 확인하기 위해 워크플로우에 `view()` 연산자도 추가합니다.

=== "후"

    ```groovy title="main.nf" linenums="30" hl_lines="4"
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map{ row -> separateMetadata(row) }
            .view()
    ```

=== "전"

    ```groovy title="main.nf" linenums="30"
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map{ row -> separateMetadata(row) }
    ```

워크플로우를 실행합니다.

```bash
nextflow run main.nf
```

??? success "명령 출력"

    ```console
    [[id:sample_001, organism:human, tissue:liver, depth:30000000, quality:38.5, run:UNSPECIFIED, sample_num:1, lane:001, read:R1, chunk:001, priority:normal], /workspaces/training/side-quests/essential_scripting_patterns/data/sequences/SAMPLE_001_S1_L001_R1_001.fastq]
    [[id:sample_002, organism:mouse, tissue:brain, depth:25000000, quality:35.2, run:UNSPECIFIED, sample_num:2, lane:001, read:R1, chunk:001, priority:normal], /workspaces/training/side-quests/essential_scripting_patterns/data/sequences/SAMPLE_002_S2_L001_R1_001.fastq]
    [[id:sample_003, organism:human, tissue:kidney, depth:45000000, quality:42.1, run:UNSPECIFIED, sample_num:3, lane:001, read:R1, chunk:001, priority:high], /workspaces/training/side-quests/essential_scripting_patterns/data/sequences/SAMPLE_003_S3_L001_R1_001.fastq]
    ```

완벽합니다! 이제 모든 샘플에 실제 실행 ID(대문자) 또는 기본값 'UNSPECIFIED'가 있는 `run` 필드가 있습니다. `?.`와 `?:`의 조합은 안전성(충돌 없음)과 합리적인 기본값을 모두 제공합니다.

작동이 확인되었으므로 `.view()` 연산자를 제거합니다.

!!! tip "팁: 안전한 탐색과 Elvis 결합하기"

    `value?.method() ?: 'default'` 패턴은 프로덕션 워크플로우에서 일반적입니다.

    - `value?.method()` - 메서드를 안전하게 호출하고, `value`가 `null`이면 `null`을 반환합니다
    - `?: 'default'` - 결과가 `null`이면 대체값을 제공합니다

    이 패턴은 누락되거나 불완전한 데이터를 우아하게 처리합니다.

함수, 연산자 closure(`.map{}`, `.filter{}`), 프로세스 스크립트, 설정 파일에서 이 연산자들을 일관되게 사용하세요. 실제 데이터를 처리할 때 충돌을 방지합니다.

### 핵심 정리

- **안전한 탐색(`?.`)**: null 값에서의 충돌 방지 - 예외를 던지는 대신 null을 반환합니다
- **Elvis 연산자(`?:`)**: 기본값 제공 - `value ?: 'default'`
- **결합**: `value?.method() ?: 'default'`가 일반적인 패턴입니다

이 연산자들은 불완전한 데이터에 대해 워크플로우를 탄력적으로 만듭니다. 실제 작업에 필수적입니다.

---

## 7. `error()`와 `log.warn`을 사용한 검증

때로는 입력 매개변수가 유효하지 않은 경우 워크플로우를 즉시 중지해야 합니다. Nextflow에서는 `error()`와 `log.warn` 같은 내장 함수와 `if` 문 및 boolean 로직 같은 표준 프로그래밍 구조를 사용하여 검증 로직을 구현할 수 있습니다. 워크플로우에 검증을 추가합니다.

워크플로우 블록 앞에 검증 함수를 생성하고, 워크플로우에서 호출하고, CSV 파일 경로에 매개변수를 사용하도록 채널 생성을 변경합니다. 매개변수가 없거나 파일이 존재하지 않으면 `error()`를 호출하여 명확한 메시지와 함께 실행을 중지합니다.

=== "후"

    ```groovy title="main.nf" linenums="1" hl_lines="5-20 23-24"
    include { FASTP } from './modules/fastp.nf'
    include { TRIMGALORE } from './modules/trimgalore.nf'
    include { GENERATE_REPORT } from './modules/generate_report.nf'

    def validateInputs() {
        // 입력 매개변수가 제공되었는지 확인
        if (!params.input) {
            error("Input CSV file path not provided. Please specify --input <file.csv>")
        }

        // CSV 파일이 존재하는지 확인
        if (!file(params.input).exists()) {
            error("Input CSV file not found: ${params.input}")
        }
    }
    ...
    workflow {
        validateInputs()
        ch_samples = channel.fromPath(params.input)
    ```

=== "전"

    ```groovy title="main.nf" linenums="1"
    include { FASTP } from './modules/fastp.nf'
    include { TRIMGALORE } from './modules/trimgalore.nf'
    include { GENERATE_REPORT } from './modules/generate_report.nf'

    ...
    workflow {
        ch_samples = channel.fromPath("./data/samples.csv")
    ```

CSV 파일 없이 실행해 봅니다.

```bash
nextflow run main.nf
```

??? failure "명령 출력"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [confident_coulomb] DSL2 - revision: 07059399ed

    WARN: Access to undefined parameter `input` -- Initialise it to a default value eg. `params.input = some_value`
    Input CSV file path not provided. Please specify --input <file.csv>
    ```

나중에 불가사의하게 실패하는 대신 명확한 오류 메시지와 함께 워크플로우가 즉시 중지됩니다.

존재하지 않는 파일로 실행합니다.

```bash
nextflow run main.nf --input ./data/nonexistent.csv
```

??? failure "명령 출력"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [cranky_gates] DSL2 - revision: 26839ae3eb

    Input CSV file not found: ./data/nonexistent.csv
    ```

마지막으로 올바른 파일로 실행합니다.

```bash
nextflow run main.nf --input ./data/samples.csv
```

??? success "명령 출력"

    ```console
    <!-- TODO: output -->
    ```

이번에는 성공적으로 실행됩니다.

`separateMetadata` 함수 내에도 검증을 추가할 수 있습니다. 치명적이지 않은 `log.warn`을 사용하여 낮은 시퀀싱 심도를 가진 샘플에 대한 경고를 발행하되, 워크플로우는 계속 진행합니다.

=== "후"

    ```groovy title="main.nf" linenums="1" hl_lines="3-6"
        def priority = sample_meta.quality > 40 ? 'high' : 'normal'

        // 데이터가 합리적인지 검증
        if (sample_meta.depth < 30000000) {
            log.warn "Low sequencing depth for ${sample_meta.id}: ${sample_meta.depth}"
        }

        return tuple(sample_meta + file_meta + [priority: priority], fastq_path)
    }
    ```

=== "전"

    ```groovy title="main.nf" linenums="1"
        def priority = sample_meta.quality > 40 ? 'high' : 'normal'

        return tuple(sample_meta + file_meta + [priority: priority], fastq_path)
    }
    ```

원래 CSV로 워크플로우를 다시 실행합니다.

```bash
nextflow run main.nf --input ./data/samples.csv
```

??? warning "명령 출력"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [awesome_goldwasser] DSL2 - revision: a31662a7c1

    executor >  local (5)
    [ce/df5eeb] process > FASTP (2)           [100%] 2 of 2 ✔
    [-        ] process > TRIMGALORE          -
    [d1/7d2b4b] process > GENERATE_REPORT (3) [100%] 3 of 3 ✔
    WARN: Low sequencing depth for sample_002: 25000000
    ```

샘플 중 하나의 낮은 시퀀싱 심도에 대한 경고가 표시됩니다.

### 핵심 정리

- **`error()`**: 명확한 메시지와 함께 워크플로우를 즉시 중지합니다
- **`log.warn`**: 워크플로우를 중지하지 않고 경고를 발행합니다
- **조기 검증**: 처리 전에 입력을 확인하여 명확한 오류 메시지와 함께 빠르게 실패합니다
- **검증 함수**: 워크플로우 시작 시 호출할 수 있는 재사용 가능한 검증 로직을 생성합니다

적절한 검증은 명확한 오류 메시지로 문제를 조기에 발견하여 워크플로우를 더 견고하고 사용자 친화적으로 만듭니다.

---

## 8. 워크플로우 이벤트 핸들러

지금까지는 워크플로우 스크립트와 프로세스 정의에 코드를 작성했습니다. 하지만 알아야 할 중요한 기능이 하나 더 있습니다. 바로 워크플로우 이벤트 핸들러입니다.

이벤트 핸들러는 워크플로우 라이프사이클의 특정 시점에 실행되는 closure입니다. 로깅, 알림, 정리 작업을 추가하는 데 적합합니다. 이 핸들러는 워크플로우 정의와 함께 워크플로우 스크립트에 정의해야 합니다.

### 8.1. `onComplete` 핸들러

가장 일반적으로 사용되는 이벤트 핸들러는 `onComplete`입니다. 워크플로우가 완료될 때(성공 또는 실패 여부에 관계없이) 실행됩니다. 파이프라인 결과를 요약하는 핸들러를 추가합니다.

`main.nf` 파일의 워크플로우 정의 내부에 이벤트 핸들러를 추가합니다.

=== "후"

    ```groovy title="main.nf" linenums="66" hl_lines="5-16"
        ch_fastp = FASTP(trim_branches.fastp)
        ch_trimgalore = TRIMGALORE(trim_branches.trimgalore)
        GENERATE_REPORT(ch_samples)

        workflow.onComplete = {
            println ""
            println "Pipeline execution summary:"
            println "=========================="
            println "Completed at: ${workflow.complete}"
            println "Duration    : ${workflow.duration}"
            println "Success     : ${workflow.success}"
            println "workDir     : ${workflow.workDir}"
            println "exit status : ${workflow.exitStatus}"
            println ""
        }
    }
    ```

=== "전"

    ```groovy title="main.nf" linenums="66" hl_lines="4"
        ch_fastp = FASTP(trim_branches.fastp)
        ch_trimgalore = TRIMGALORE(trim_branches.trimgalore)
        GENERATE_REPORT(ch_samples)
    }
    ```

이 closure는 워크플로우가 완료될 때 실행됩니다. 내부에서 실행에 대한 유용한 속성을 제공하는 `workflow` 객체에 접근할 수 있습니다.

워크플로우를 실행하면 마지막에 이 요약이 표시됩니다!

```bash
nextflow run main.nf --input ./data/samples.csv -ansi-log false
```

??? success "명령 출력"

    ```console
    N E X T F L O W  ~  version 25.10.2
    Launching `main.nf` [marvelous_boltzmann] DSL2 - revision: a31662a7c1
    WARN: Low sequencing depth for sample_002: 25000000
    [9b/d48e40] Submitted process > FASTP (2)
    [6a/73867a] Submitted process > GENERATE_REPORT (2)
    [79/ad0ac5] Submitted process > GENERATE_REPORT (1)
    [f3/bda6cb] Submitted process > FASTP (1)
    [34/d5b52f] Submitted process > GENERATE_REPORT (3)

    Pipeline execution summary:
    ==========================
    Completed at: 2025-10-10T12:14:24.885384+01:00
    Duration    : 2.9s
    Success     : true
    workDir     : /workspaces/training/side-quests/essential_scripting_patterns/work
    exit status : 0
    ```

조건부 로직을 추가하여 더 유용하게 만듭니다.

=== "후"

    ```groovy title="main.nf" linenums="66" hl_lines="5-22"
        ch_fastp = FASTP(trim_branches.fastp)
        ch_trimgalore = TRIMGALORE(trim_branches.trimgalore)
        GENERATE_REPORT(ch_samples)

        workflow.onComplete = {
            println ""
            println "Pipeline execution summary:"
            println "=========================="
            println "Completed at: ${workflow.complete}"
            println "Duration    : ${workflow.duration}"
            println "Success     : ${workflow.success}"
            println "workDir     : ${workflow.workDir}"
            println "exit status : ${workflow.exitStatus}"
            println ""

            if (workflow.success) {
                println "✅ Pipeline completed successfully!"
            } else {
                println "❌ Pipeline failed!"
                println "Error: ${workflow.errorMessage}"
            }
        }
    }
    ```

=== "전"

    ```groovy title="main.nf" linenums="66" hl_lines="5-16"
        ch_fastp = FASTP(trim_branches.fastp)
        ch_trimgalore = TRIMGALORE(trim_branches.trimgalore)
        GENERATE_REPORT(ch_samples)

        workflow.onComplete = {
            println ""
            println "Pipeline execution summary:"
            println "=========================="
            println "Completed at: ${workflow.complete}"
            println "Duration    : ${workflow.duration}"
            println "Success     : ${workflow.success}"
            println "workDir     : ${workflow.workDir}"
            println "exit status : ${workflow.exitStatus}"
            println ""
        }
    }
    ```

이제 성공/실패 메시지를 포함하는 더 유익한 요약을 얻을 수 있습니다.

<!-- TODO: add run command -->

??? success "명령 출력"

    ```console
    N E X T F L O W  ~  version 25.10.2
    Launching `main.nf` [boring_linnaeus] DSL2 - revision: a31662a7c1
    WARN: Low sequencing depth for sample_002: 25000000
    [e5/242efc] Submitted process > FASTP (2)
    [3b/74047c] Submitted process > GENERATE_REPORT (3)
    [8a/7a57e6] Submitted process > GENERATE_REPORT (1)
    [a8/b1a31f] Submitted process > GENERATE_REPORT (2)
    [40/648429] Submitted process > FASTP (1)

    Pipeline execution summary:
    ==========================
    Completed at: 2025-10-10T12:16:00.522569+01:00
    Duration    : 3.6s
    Success     : true
    workDir     : /workspaces/training/side-quests/essential_scripting_patterns/work
    exit status : 0

    ✅ Pipeline completed successfully!
    ```

파일 작업을 사용하여 요약을 파일에 작성할 수도 있습니다.

```groovy title="main.nf - Writing summary to file"
workflow {
    // ... 워크플로우 코드 ...

    workflow.onComplete = {
        def summary = """
        Pipeline Execution Summary
        ===========================
        Completed: ${workflow.complete}
        Duration : ${workflow.duration}
        Success  : ${workflow.success}
        Command  : ${workflow.commandLine}
        """

        println summary

        // 로그 파일에 작성
        def log_file = file("${workflow.launchDir}/pipeline_summary.txt")
        log_file.text = summary
    }
}
```

### 8.2. `onError` 핸들러

`onComplete` 외에도 사용할 수 있는 이벤트 핸들러가 하나 더 있습니다. `onError`는 워크플로우가 실패할 때만 실행됩니다.

```groovy title="main.nf - onError handler"
workflow {
    // ... 워크플로우 코드 ...

    workflow.onError = {
        println "="* 50
        println "Pipeline execution failed!"
        println "Error message: ${workflow.errorMessage}"
        println "="* 50

        // 상세 오류 로그 작성
        def error_file = file("${workflow.launchDir}/error.log")
        error_file.text = """
        Workflow Error Report
        =====================
        Time: ${new Date()}
        Error: ${workflow.errorMessage}
        Error report: ${workflow.errorReport ?: 'No detailed report available'}
        """

        println "Error details written to: ${error_file}"
    }
}
```

워크플로우 스크립트에서 여러 핸들러를 함께 사용할 수 있습니다.

```groovy title="main.nf - Combined handlers"
workflow {
    // ... 워크플로우 코드 ...

    workflow.onError = {
        println "Workflow failed: ${workflow.errorMessage}"
    }

    workflow.onComplete = {
        def duration_mins = workflow.duration.toMinutes().round(2)
        def status = workflow.success ? "SUCCESS ✅" : "FAILED ❌"

        println """
        Pipeline finished: ${status}
        Duration: ${duration_mins} minutes
        """
    }
}
```

### 핵심 정리

이 섹션에서 다음을 학습했습니다.

- **이벤트 핸들러 closure**: 워크플로우 스크립트에서 다양한 라이프사이클 시점에 실행되는 closure입니다
- **`onComplete` 핸들러**: 실행 요약 및 결과 보고에 사용합니다
- **`onError` 핸들러**: 오류 처리 및 실패 로깅에 사용합니다
- **워크플로우 객체 속성**: `workflow.success`, `workflow.duration`, `workflow.errorMessage` 등에 접근합니다

이벤트 핸들러는 워크플로우 스크립트 내에서 Nextflow 언어의 모든 기능을 사용하여 정교한 로깅 및 알림 기능을 추가하는 방법을 보여줍니다.

---

## 요약

축하합니다, 완료했습니다!

이 사이드 퀘스트 전반에 걸쳐 기본적인 메타데이터 처리에서 정교한 프로덕션 수준의 워크플로우로 발전하는 포괄적인 샘플 처리 파이프라인을 구축했습니다.
각 섹션은 이전 섹션을 기반으로 하여, 프로그래밍 구조가 단순한 워크플로우를 강력한 데이터 처리 시스템으로 어떻게 변환하는지 보여주었습니다. 다음과 같은 이점이 있습니다.

- **더 명확한 코드**: 데이터플로우와 스크립팅을 이해하면 더 체계적인 워크플로우를 작성할 수 있습니다
- **견고한 처리**: 안전한 탐색과 Elvis 연산자로 누락된 데이터에 대해 워크플로우를 탄력적으로 만듭니다
- **유연한 처리**: 조건부 로직으로 워크플로우가 다양한 샘플 유형을 적절히 처리할 수 있습니다
- **적응형 리소스**: 동적 지시문으로 입력 특성에 따라 리소스 사용을 최적화합니다

이 발전 과정은 소수의 샘플을 처리하는 연구 프로토타입에서 실험실과 기관 전반에 걸쳐 수천 개의 샘플을 처리하는 프로덕션 시스템으로의 실제 바이오인포매틱스 파이프라인 진화를 반영합니다.
해결한 모든 과제와 학습한 패턴은 Nextflow 워크플로우를 확장할 때 개발자들이 직면하는 실제 문제를 반영합니다.

이러한 패턴을 자신의 작업에 적용하면 견고하고 프로덕션 수준의 워크플로우를 구축할 수 있습니다.

### 핵심 패턴

1.  **데이터플로우 vs 스크립팅:** 데이터플로우 연산(채널 오케스트레이션)과 스크립팅(데이터를 조작하는 코드)을 구분하는 방법을 학습했습니다. 채널 vs List에서의 `collect`와 같이 다른 타입에서의 연산 간의 중요한 차이점도 포함합니다.

    - 데이터플로우: 채널 오케스트레이션

    ```groovy
    channel.fromPath('*.fastq').splitCsv(header: true)
    ```

    - 스크립팅: 컬렉션에 대한 데이터 처리

    ```groovy
    sample_data.collect { it.toUpperCase() }
    ```

2.  **고급 문자열 처리**: 파일명 파싱을 위한 정규 표현식, 프로세스에서의 동적 스크립트 생성, 변수 보간(Nextflow vs Bash vs 셸)을 마스터했습니다.

    - 패턴 매칭

    ```groovy
    filename =~ ~/^(\w+)_(\w+)_(\d+)\.fastq$/
    ```

    - 조건부 반환이 있는 함수

    ```groovy
    def parseSample(filename) {
        def matcher = filename =~ pattern
        return matcher ? [valid: true, data: matcher[0]] : [valid: false]
    }
    ```

    - 파일 컬렉션을 명령 인자로 변환 (프로세스 script 블록에서)

    ```groovy
    script:
    def file_args = input_files.collect { file -> "--input ${file}" }.join(' ')
    """
    analysis_tool ${file_args} --output results.txt
    """
    ```

3.  **재사용 가능한 함수 만들기**: 복잡한 로직을 채널 연산자에서 호출할 수 있는 명명된 함수로 추출하여 워크플로우를 더 읽기 쉽고 유지보수하기 쉽게 만드는 방법을 학습했습니다.

    - 명명된 함수 정의

    ```groovy
    def separateMetadata(row) {
        def sample_meta = [ /* code hidden for brevity */ ]
        def fastq_path = file(row.file_path)
        def m = (fastq_path.name =~ /^(.+)_S(\d+)_L(\d{3})_(R[12])_(\d{3})\.fastq(?:\.gz)?$/)
        def file_meta = m ? [ /* code hidden for brevity */ ] : [:]
        def priority = sample_meta.quality > 40 ? 'high' : 'normal'

        return tuple(sample_meta + file_meta + [priority: priority], fastq_path)
    }
    ```

    - 워크플로우에서 명명된 함수 호출

    ```groovy
    workflow {
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map{ row -> separateMetadata(row) }

        ch_fastp = FASTP(ch_samples)
    }
    ```

4.  **closure를 사용한 동적 리소스 지시문**: 입력 특성에 따른 적응형 리소스 할당을 위해 프로세스 지시문에서 closure를 사용하는 방법을 살펴봤습니다.

    - 명명된 closure와 합성

    ```groovy
    def enrichData = normalizeId >> addQualityCategory >> addFlags
    def processor = generalFunction.curry(fixedParam)
    ```

    - 범위 접근이 있는 closure

    ```groovy
    def collectStats = { data -> stats.count++; return data }
    ```

5.  **조건부 로직과 프로세스 제어**: `.branch()`와 `.filter()` 연산자를 사용한 지능적인 라우팅을 추가하고, 간결한 조건부 표현식을 위한 진리값을 활용했습니다.

    - `.branch()`를 사용하여 데이터를 다른 워크플로우 분기로 라우팅

    ```groovy
    trim_branches = ch_samples
    .branch { meta, reads ->
        fastp: meta.organism == 'human' && meta.depth >= 30000000
        trimgalore: true
    }

    ch_fastp = FASTP(trim_branches.fastp)
    ch_trimgalore = TRIMGALORE(trim_branches.trimgalore)
    ```

    - Groovy 진리값을 사용한 boolean 평가

    ```groovy
    if (sample.files) println "Has files"
    ```

    - '진리값'으로 데이터를 부분 선택하는 `filter()` 사용

    ```groovy
    ch_valid_samples = ch_samples
        .filter { meta, reads ->
            meta.id && meta.organism && meta.depth >= 25000000
        }
    ```

6.  **안전한 탐색과 Elvis 연산자**: null-safe 속성 접근을 위한 `?.`와 기본값 제공을 위한 `?:`를 사용하여 파이프라인을 누락된 데이터에 대해 견고하게 만들었습니다.

    ```groovy
    def id = data?.sample?.id ?: 'unknown'
    ```

7.  **error()와 log.warn을 사용한 검증**: 입력을 조기에 검증하고 명확한 오류 메시지와 함께 빠르게 실패하는 방법을 학습했습니다.

    ```groovy
    try {
        def errors = validateSample(sample)
        if (errors) throw new RuntimeException("Invalid: ${errors.join(', ')}")
    } catch (Exception e) {
        println "Error: ${e.message}"
    }
    ```

8.  **설정 이벤트 핸들러**: 로깅, 알림, 라이프사이클 관리를 위한 워크플로우 이벤트 핸들러(`onComplete`와 `onError`)를 사용하는 방법을 학습했습니다.

    - `onComplete`를 사용하여 로그 및 알림

    ```groovy
    workflow.onComplete = {
        println "Success     : ${workflow.success}"
        println "exit status : ${workflow.exitStatus}"

        if (workflow.success) {
            println "✅ Pipeline completed successfully!"
        } else {
            println "❌ Pipeline failed!"
            println "Error: ${workflow.errorMessage}"
        }
    }
    ```

    - `onError`를 사용하여 실패 시 특정 작업 수행

    ```groovy
    workflow.onError = {
        // 상세 오류 로그 작성
        def error_file = file("${workflow.launchDir}/error.log")
        error_file.text = """
        Time: ${new Date()}
        Error: ${workflow.errorMessage}
        Error report: ${workflow.errorReport ?: 'No detailed report available'}
        """

        println "Error details written to: ${error_file}"
    }
    ```

### 추가 자료

- [Nextflow 언어 참조](https://nextflow.io/docs/latest/reference/syntax.html)
- [Nextflow 연산자](https://www.nextflow.io/docs/latest/operator.html)
- [Nextflow 스크립트 문법](https://www.nextflow.io/docs/latest/script.html)
- [Nextflow 표준 라이브러리](https://nextflow.io/docs/latest/reference/stdlib.html)

더 고급 기능을 살펴볼 때 이 자료들을 참고하세요.

다음을 위해 기술을 연습하고 확장하면 도움이 됩니다.

- 데이터플로우와 스크립팅을 적절히 분리하여 더 깔끔한 워크플로우 작성
- Nextflow, Bash, 셸 변수와 관련된 일반적인 함정을 피하기 위한 변수 보간 마스터
- 효율적이고 적응형 워크플로우를 위한 동적 리소스 지시문 사용
- 파일 컬렉션을 올바르게 형식화된 명령줄 인자로 변환
- regex와 문자열 처리를 사용하여 다양한 파일 명명 규칙과 입력 형식을 우아하게 처리
- 고급 closure 패턴과 함수형 프로그래밍을 사용하여 재사용 가능하고 유지보수하기 쉬운 코드 구축
- 컬렉션 연산을 사용하여 복잡한 데이터셋 처리 및 구성
- 워크플로우를 프로덕션 수준으로 만들기 위한 검증, 오류 처리, 로깅 추가
- 이벤트 핸들러를 사용한 워크플로우 라이프사이클 관리 구현

---

## 다음 단계

[사이드 퀘스트 메뉴](../)로 돌아가거나 페이지 오른쪽 하단의 버튼을 클릭하여 목록의 다음 주제로 이동하세요.
