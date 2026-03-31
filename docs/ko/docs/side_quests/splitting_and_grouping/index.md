# 분할 및 그룹화

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } AI 지원 번역 - [자세히 알아보기 및 개선 제안](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Nextflow는 데이터를 유연하게 처리할 수 있는 강력한 도구를 제공합니다. 핵심 기능 중 하나는 데이터를 여러 스트림으로 분할한 다음 관련 항목을 다시 그룹화하는 것입니다. 이는 분석을 위해 결과를 결합하기 전에 서로 다른 유형의 샘플을 독립적으로 처리해야 하는 생물정보학 워크플로우에서 특히 유용합니다.

우편물 분류에 비유하면 이해하기 쉽습니다. 목적지별로 편지를 분류하고, 각 묶음을 다르게 처리한 다음, 같은 사람에게 가는 항목을 다시 합치는 것과 같습니다. Nextflow는 특수 연산자를 사용하여 과학 데이터에 이 작업을 수행합니다. 이 접근 방식은 분산 컴퓨팅 및 생물정보학 워크플로우에서 **scatter/gather** 패턴으로도 널리 알려져 있습니다.

Nextflow의 채널 시스템은 이러한 유연성의 핵심입니다. 채널은 워크플로우의 여러 부분을 연결하여 데이터가 분석 과정을 통해 흐를 수 있도록 합니다. 단일 데이터 소스에서 여러 채널을 생성하고, 각 채널을 다르게 처리한 다음, 필요할 때 채널을 다시 병합할 수 있습니다. 이 접근 방식을 통해 복잡한 생물정보학 분석의 분기 및 수렴 경로를 자연스럽게 반영하는 워크플로우를 설계할 수 있습니다.

### 학습 목표

이 사이드 퀘스트에서는 Nextflow의 채널 연산자를 사용하여 데이터를 분할하고 그룹화하는 방법을 학습합니다.
샘플 정보와 관련 데이터 파일이 포함된 CSV 파일로 시작하여 이 데이터를 조작하고 재구성합니다.

이 사이드 퀘스트를 마치면 다음 기법을 사용하여 데이터 스트림을 효과적으로 분리하고 결합할 수 있습니다.

- `splitCsv`를 사용하여 파일에서 데이터 읽기
- `filter`와 `map`으로 데이터 필터링 및 변환
- `join`과 `groupTuple`을 사용하여 관련 데이터 결합
- 병렬 처리를 위한 `combine`으로 데이터 조합 생성
- `subMap`과 중복 제거 전략을 사용하여 데이터 구조 최적화
- 채널 구조 조작을 위한 이름 있는 closure로 재사용 가능한 함수 작성

이러한 기술은 깔끔하고 유지 관리하기 쉬운 코드 구조를 유지하면서 여러 입력 파일과 다양한 유형의 데이터를 효율적으로 처리하는 워크플로우를 구축하는 데 도움이 됩니다.

### 사전 요구 사항

이 사이드 퀘스트를 시작하기 전에 다음을 준비하세요.

- [Hello Nextflow](../hello_nextflow/README.md) 튜토리얼 또는 동급의 입문 과정을 완료해야 합니다.
- 기본 Nextflow 개념과 메커니즘(프로세스, 채널, 연산자, 파일 작업, 메타 데이터)에 익숙해야 합니다.

**선택 사항:** [워크플로우의 메타데이터](../metadata/) 사이드 퀘스트를 먼저 완료하는 것을 권장합니다.
해당 퀘스트에서는 `splitCsv`로 CSV 파일을 읽고 meta map을 생성하는 기본 사항을 다루며, 여기서도 이를 많이 활용합니다.

---

## 0. 시작하기

#### 교육 코드스페이스 열기

아직 열지 않았다면 [환경 설정](../envsetup/index.md)에 설명된 대로 교육 환경을 열어 주세요.

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

#### 프로젝트 디렉토리로 이동

이 튜토리얼의 파일이 있는 디렉토리로 이동합니다.

```bash
cd side-quests/splitting_and_grouping
```

VSCode에서 이 디렉토리에 포커스를 설정할 수 있습니다.

```bash
code .
```

#### 자료 검토

메인 워크플로우 파일과 `samplesheet.csv`라는 샘플시트가 포함된 `data` 디렉토리를 확인할 수 있습니다.

```console title="Directory contents"
.
├── data
│   └── samplesheet.csv
└── main.nf
```

샘플시트에는 환자 ID, 샘플 반복 번호, 유형(정상 또는 종양), 가상 데이터 파일 경로(실제로는 존재하지 않지만 존재한다고 가정합니다)를 포함한 여러 환자의 샘플 정보가 포함되어 있습니다.

```console title="samplesheet.csv"
id,repeat,type,bam
patientA,1,normal,patientA_rep1_normal.bam
patientA,1,tumor,patientA_rep1_tumor.bam
patientA,2,normal,patientA_rep2_normal.bam
patientA,2,tumor,patientA_rep2_tumor.bam
patientB,1,normal,patientB_rep1_normal.bam
patientB,1,tumor,patientB_rep1_tumor.bam
patientC,1,normal,patientC_rep1_normal.bam
patientC,1,tumor,patientC_rep1_tumor.bam
```

이 샘플시트에는 세 명의 환자(A, B, C)로부터 채취한 8개의 샘플이 나열되어 있습니다.

각 환자에 대해 `tumor` 유형(일반적으로 종양 생검에서 유래)과 `normal` 유형(건강한 조직이나 혈액에서 채취)의 샘플이 있습니다.
암 분석에 익숙하지 않더라도, 이것이 대조 분석을 수행하기 위해 쌍을 이룬 종양/정상 샘플을 사용하는 실험 모델에 해당한다는 것만 알면 됩니다.

특히 환자 A의 경우 두 세트의 기술적 반복(repeat)이 있습니다.

!!! note "참고"

    이 실험 설계에 익숙하지 않더라도 걱정하지 마세요. 이 튜토리얼을 이해하는 데 필수적인 내용은 아닙니다.

#### 과제 검토

여러분의 과제는 다음을 수행하는 Nextflow 워크플로우를 작성하는 것입니다.

1. **읽기**: CSV 파일에서 샘플 데이터를 읽고 meta map으로 구조화
2. **분리**: 유형(정상 vs 종양)에 따라 샘플을 다른 채널로 분리
3. **결합**: 환자 ID와 반복 번호로 매칭된 종양/정상 쌍을 join
4. **분배**: 병렬 처리를 위해 유전체 구간에 샘플 분배
5. **그룹화**: 다운스트림 분석을 위해 관련 샘플을 다시 그룹화

이는 독립적인 처리를 위해 데이터를 분할한 다음 비교 분석을 위해 관련 항목을 재결합해야 하는 일반적인 생물정보학 패턴입니다.

#### 준비 체크리스트

시작할 준비가 되었나요?

- [ ] 이 과정의 목표와 사전 요구 사항을 이해했습니다
- [ ] 코드스페이스가 실행 중입니다
- [ ] 작업 디렉토리를 적절히 설정했습니다
- [ ] 과제를 이해했습니다

모든 항목을 체크할 수 있다면 시작할 준비가 된 것입니다.

---

## 1. 샘플 데이터 읽기

### 1.1. `splitCsv`로 샘플 데이터를 읽고 meta map 생성

`splitCsv`로 샘플 데이터를 읽고 meta map 패턴으로 구성하는 것부터 시작합니다. `main.nf`에서 이미 워크플로우가 시작된 것을 확인할 수 있습니다.

```groovy title="main.nf" linenums="1" hl_lines="2"
workflow {
    ch_samplesheet = channel.fromPath("./data/samplesheet.csv")
}
```

!!! note "참고"

    이 튜토리얼 전반에 걸쳐 모든 채널 변수에 `ch_` 접두사를 사용하여 Nextflow 채널임을 명확히 나타냅니다.

[워크플로우의 메타데이터](../metadata/) 사이드 퀘스트를 완료했다면 이 패턴을 알아볼 것입니다. `splitCsv`를 사용하여 CSV를 읽고 meta map으로 데이터를 구조화하여 메타데이터와 파일 경로를 분리합니다.

!!! info "정보"

    이 교육에서는 `map`이라는 두 가지 다른 개념을 접하게 됩니다.

    - **데이터 구조**: 키-값 쌍을 저장하는 Groovy map(다른 언어의 딕셔너리/해시에 해당)
    - **채널 연산자**: 채널의 항목을 변환하는 `.map()` 연산자

    문맥에 따라 어느 것을 의미하는지 명확히 하겠지만, Nextflow로 작업할 때 이 구분을 이해하는 것이 중요합니다.

`main.nf`에 다음 변경 사항을 적용합니다.

=== "후"

    ```groovy title="main.nf" linenums="2" hl_lines="2-6"
        ch_samples = channel.fromPath("./data/samplesheet.csv")
            .splitCsv(header: true)
            .map{ row ->
              [[id:row.id, repeat:row.repeat, type:row.type], row.bam]
            }
            .view()
    ```

=== "전"

    ```groovy title="main.nf" linenums="2" hl_lines="1"
        ch_samplesheet = channel.fromPath("./data/samplesheet.csv")
    ```

이는 `splitCsv` 작업(헤더가 있는 CSV 읽기)과 `map` 작업(데이터를 `[meta, file]` 튜플로 구조화)을 한 단계로 결합합니다. 변경 사항을 적용하고 파이프라인을 실행합니다.

```bash
nextflow run main.nf
```

??? success "명령 출력"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [deadly_mercator] DSL2 - revision: bd6b0224e9

    [[id:patientA, repeat:1, type:normal], patientA_rep1_normal.bam]
    [[id:patientA, repeat:1, type:tumor], patientA_rep1_tumor.bam]
    [[id:patientA, repeat:2, type:normal], patientA_rep2_normal.bam]
    [[id:patientA, repeat:2, type:tumor], patientA_rep2_tumor.bam]
    [[id:patientB, repeat:1, type:normal], patientB_rep1_normal.bam]
    [[id:patientB, repeat:1, type:tumor], patientB_rep1_tumor.bam]
    [[id:patientC, repeat:1, type:normal], patientC_rep1_normal.bam]
    [[id:patientC, repeat:1, type:tumor], patientC_rep1_tumor.bam]
    ```

이제 각 항목이 `[meta, file]` 튜플인 채널이 생성되었습니다. 메타데이터가 파일 경로와 분리되어 있습니다. 이 구조를 통해 메타데이터 필드를 기반으로 작업을 분할하고 그룹화할 수 있습니다.

---

## 2. 데이터 필터링 및 변환

### 2.1. `filter`로 데이터 필터링

[`filter` 연산자](https://www.nextflow.io/docs/latest/operator.html#filter)를 사용하여 조건에 따라 데이터를 필터링할 수 있습니다. 정상 샘플만 처리하려는 경우, `type` 필드를 기반으로 데이터를 필터링할 수 있습니다. `view` 연산자 앞에 이를 삽입합니다.

=== "후"

    ```groovy title="main.nf" linenums="2" hl_lines="6"
        ch_samples = channel.fromPath("./data/samplesheet.csv")
            .splitCsv(header: true)
            .map{ row ->
              [[id:row.id, repeat:row.repeat, type:row.type], row.bam]
            }
            .filter { meta, file -> meta.type == 'normal' }
            .view()
    ```

=== "전"

    ```groovy title="main.nf" linenums="2"
        ch_samples = channel.fromPath("./data/samplesheet.csv")
            .splitCsv(header: true)
            .map{ row ->
              [[id:row.id, repeat:row.repeat, type:row.type], row.bam]
            }
            .view()
    ```

워크플로우를 다시 실행하여 필터링된 결과를 확인합니다.

```bash
nextflow run main.nf
```

??? success "명령 출력"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [admiring_brown] DSL2 - revision: 194d61704d

    [[id:patientA, repeat:1, type:normal], patientA_rep1_normal.bam]
    [[id:patientA, repeat:2, type:normal], patientA_rep2_normal.bam]
    [[id:patientB, repeat:1, type:normal], patientB_rep1_normal.bam]
    [[id:patientC, repeat:1, type:normal], patientC_rep1_normal.bam]
    ```

정상 샘플만 포함하도록 데이터를 성공적으로 필터링했습니다. 작동 방식을 정리해 보겠습니다.

`filter` 연산자는 채널의 각 요소에 적용되는 closure를 받습니다. closure가 `true`를 반환하면 해당 요소가 포함되고, `false`를 반환하면 제외됩니다.

여기서는 `meta.type == 'normal'`인 샘플만 유지하려 합니다. closure는 각 샘플을 참조하기 위해 튜플 `meta,file`을 사용하고, `meta.type`으로 샘플 유형에 접근하여 `'normal'`과 같은지 확인합니다.

이는 위에서 소개한 단일 closure로 수행됩니다.

```groovy title="main.nf" linenums="7"
    .filter { meta, file -> meta.type == 'normal' }
```

### 2.2. 별도의 필터링된 채널 생성

현재 CSV에서 직접 생성된 채널에 filter를 적용하고 있지만, 여러 방식으로 필터링하려 하므로 정상 샘플을 위한 별도의 필터링된 채널을 생성하도록 로직을 재작성합니다.

=== "후"

    ```groovy title="main.nf" linenums="2" hl_lines="6 8"
        ch_samples = channel.fromPath("./data/samplesheet.csv")
            .splitCsv(header: true)
            .map{ row ->
                [[id:row.id, repeat:row.repeat, type:row.type], row.bam]
            }
        ch_normal_samples = ch_samples
            .filter { meta, file -> meta.type == 'normal' }
        ch_normal_samples
            .view()
    ```

=== "전"

    ```groovy title="main.nf" linenums="2"
        ch_samples = channel.fromPath("./data/samplesheet.csv")
            .splitCsv(header: true)
            .map{ row ->
              [[id:row.id, repeat:row.repeat, type:row.type], row.bam]
            }
            .filter { meta, file -> meta.type == 'normal' }
            .view()
    ```

파이프라인을 실행하여 결과를 확인합니다.

```bash
nextflow run main.nf
```

??? success "명령 출력"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [trusting_poisson] DSL2 - revision: 639186ee74

    [[id:patientA, repeat:1, type:normal], patientA_rep1_normal.bam]
    [[id:patientA, repeat:2, type:normal], patientA_rep2_normal.bam]
    [[id:patientB, repeat:1, type:normal], patientB_rep1_normal.bam]
    [[id:patientC, repeat:1, type:normal], patientC_rep1_normal.bam]
    ```

데이터를 성공적으로 필터링하고 정상 샘플을 위한 별도의 채널을 생성했습니다.

종양 샘플을 위한 필터링된 채널도 생성합니다.

=== "후"

    ```groovy title="main.nf" linenums="7" hl_lines="3-8"
        ch_normal_samples = ch_samples
            .filter { meta, file -> meta.type == 'normal' }
        ch_tumor_samples = ch_samples
            .filter { meta, file -> meta.type == 'tumor' }
        ch_normal_samples
            .view{'Normal sample: ' + it}
        ch_tumor_samples
            .view{'Tumor sample: ' + it}
    ```

=== "전"

    ```groovy title="main.nf" linenums="7" hl_lines="3 4"
        ch_normal_samples = ch_samples
            .filter { meta, file -> meta.type == 'normal' }
        ch_normal_samples
            .view()
    ```

```bash
nextflow run main.nf
```

??? success "명령 출력"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [maniac_boltzmann] DSL2 - revision: 3636b6576b

    Tumor sample: [[id:patientA, repeat:1, type:tumor], patientA_rep1_tumor.bam]
    Tumor sample: [[id:patientA, repeat:2, type:tumor], patientA_rep2_tumor.bam]
    Normal sample: [[id:patientA, repeat:1, type:normal], patientA_rep1_normal.bam]
    Normal sample: [[id:patientA, repeat:2, type:normal], patientA_rep2_normal.bam]
    Normal sample: [[id:patientB, repeat:1, type:normal], patientB_rep1_normal.bam]
    Normal sample: [[id:patientC, repeat:1, type:normal], patientC_rep1_normal.bam]
    Tumor sample: [[id:patientB, repeat:1, type:tumor], patientB_rep1_tumor.bam]
    Tumor sample: [[id:patientC, repeat:1, type:tumor], patientC_rep1_tumor.bam]
    ```

정상 샘플과 종양 샘플을 두 개의 다른 채널로 분리하고, `view()`에 제공된 closure를 사용하여 출력에서 다르게 레이블을 지정했습니다: `ch_tumor_samples.view{'Tumor sample: ' + it}`.

### 핵심 정리

이 섹션에서 학습한 내용:

- **데이터 필터링**: `filter`로 데이터를 필터링하는 방법
- **데이터 분할**: 조건에 따라 데이터를 다른 채널로 분할하는 방법
- **데이터 확인**: `view`를 사용하여 데이터를 출력하고 다른 채널의 출력에 레이블을 지정하는 방법

이제 정상 샘플과 종양 샘플을 두 개의 다른 채널로 분리했습니다. 다음으로 `id` 필드를 기준으로 정상 샘플과 종양 샘플을 join합니다.

---

## 3. 식별자로 채널 결합

이전 섹션에서 정상 샘플과 종양 샘플을 두 개의 다른 채널로 분리했습니다. 이들은 유형에 따라 특정 프로세스나 워크플로우를 사용하여 독립적으로 처리될 수 있습니다. 하지만 같은 환자의 정상 샘플과 종양 샘플을 비교하려면 어떻게 해야 할까요? 이 시점에서 `id` 필드를 기반으로 샘플을 매칭하여 다시 결합해야 합니다.

Nextflow에는 채널을 결합하는 다양한 방법이 있지만, 이 경우 가장 적합한 연산자는 [`join`](https://www.nextflow.io/docs/latest/operator.html#join)입니다. SQL에 익숙하다면 `JOIN` 작업처럼 결합할 키와 수행할 join 유형을 지정하는 방식으로 동작합니다.

### 3.1. `map`과 `join`을 사용하여 환자 ID 기반으로 결합

[`join`](https://www.nextflow.io/docs/latest/operator.html#join) 문서를 확인하면 기본적으로 각 튜플의 첫 번째 항목을 기반으로 두 채널을 결합한다는 것을 알 수 있습니다.

#### 3.1.1. 데이터 구조 확인

콘솔 출력이 더 이상 없다면 파이프라인을 실행하여 데이터 구조를 확인하고 `id` 필드로 join하기 위해 어떻게 수정해야 하는지 살펴봅니다.

```bash
nextflow run main.nf
```

??? success "명령 출력"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [maniac_boltzmann] DSL2 - revision: 3636b6576b

    Tumor sample: [[id:patientA, repeat:1, type:tumor], patientA_rep1_tumor.bam]
    Tumor sample: [[id:patientA, repeat:2, type:tumor], patientA_rep2_tumor.bam]
    Normal sample: [[id:patientA, repeat:1, type:normal], patientA_rep1_normal.bam]
    Normal sample: [[id:patientA, repeat:2, type:normal], patientA_rep2_normal.bam]
    Normal sample: [[id:patientB, repeat:1, type:normal], patientB_rep1_normal.bam]
    Normal sample: [[id:patientC, repeat:1, type:normal], patientC_rep1_normal.bam]
    Tumor sample: [[id:patientB, repeat:1, type:tumor], patientB_rep1_tumor.bam]
    Tumor sample: [[id:patientC, repeat:1, type:tumor], patientC_rep1_tumor.bam]
    ```

`id` 필드가 각 meta map의 첫 번째 요소임을 알 수 있습니다. `join`이 작동하려면 각 튜플에서 `id` 필드를 분리해야 합니다. 그런 다음 `join` 연산자를 사용하여 두 채널을 결합할 수 있습니다.

#### 3.1.2. `id` 필드 분리

`id` 필드를 분리하기 위해 [`map` 연산자](https://www.nextflow.io/docs/latest/operator.html#map)를 사용하여 `id` 필드를 첫 번째 요소로 하는 새 튜플을 생성할 수 있습니다.

=== "후"

    ```groovy title="main.nf" linenums="7" hl_lines="3 6"
        ch_normal_samples = ch_samples
            .filter { meta, file -> meta.type == 'normal' }
            .map { meta, file -> [meta.id, meta, file] }
        ch_tumor_samples = ch_samples
            .filter { meta, file -> meta.type == 'tumor' }
            .map { meta, file -> [meta.id, meta, file] }
        ch_normal_samples
            .view{'Normal sample: ' + it}
        ch_tumor_samples
            .view{'Tumor sample: ' + it}
    ```

=== "전"

    ```groovy title="main.nf" linenums="7"
        ch_normal_samples = ch_samples
            .filter { meta, file -> meta.type == 'normal' }
        ch_tumor_samples = ch_samples
            .filter { meta, file -> meta.type == 'tumor' }
        ch_normal_samples
            .view{'Normal sample: ' + it}
        ch_tumor_samples
            .view{'Tumor sample: ' + it}
    ```

```bash
nextflow run main.nf
```

??? success "명령 출력"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [mad_lagrange] DSL2 - revision: 9940b3f23d

    Tumor sample: [patientA, [id:patientA, repeat:1, type:tumor], patientA_rep1_tumor.bam]
    Tumor sample: [patientA, [id:patientA, repeat:2, type:tumor], patientA_rep2_tumor.bam]
    Normal sample: [patientA, [id:patientA, repeat:1, type:normal], patientA_rep1_normal.bam]
    Normal sample: [patientA, [id:patientA, repeat:2, type:normal], patientA_rep2_normal.bam]
    Tumor sample: [patientB, [id:patientB, repeat:1, type:tumor], patientB_rep1_tumor.bam]
    Tumor sample: [patientC, [id:patientC, repeat:1, type:tumor], patientC_rep1_tumor.bam]
    Normal sample: [patientB, [id:patientB, repeat:1, type:normal], patientB_rep1_normal.bam]
    Normal sample: [patientC, [id:patientC, repeat:1, type:normal], patientC_rep1_normal.bam]
    ```

미묘할 수 있지만, 각 튜플의 첫 번째 요소가 `id` 필드임을 확인할 수 있습니다.

#### 3.1.3. 두 채널 결합

이제 `join` 연산자를 사용하여 `id` 필드를 기반으로 두 채널을 결합할 수 있습니다.

`view`를 사용하여 결합된 출력을 출력합니다.

=== "후"

    ```groovy title="main.nf" linenums="7" hl_lines="7-9"
        ch_normal_samples = ch_samples
            .filter { meta, file -> meta.type == 'normal' }
            .map { meta, file -> [meta.id, meta, file] }
        ch_tumor_samples = ch_samples
            .filter { meta, file -> meta.type == 'tumor' }
            .map { meta, file -> [meta.id, meta, file] }
        ch_joined_samples = ch_normal_samples
            .join(ch_tumor_samples)
        ch_joined_samples.view()
    ```

=== "전"

    ```groovy title="main.nf" linenums="7" hl_lines="7-10"
        ch_normal_samples = ch_samples
            .filter { meta, file -> meta.type == 'normal' }
            .map { meta, file -> [meta.id, meta, file] }
        ch_tumor_samples = ch_samples
            .filter { meta, file -> meta.type == 'tumor' }
            .map { meta, file -> [meta.id, meta, file] }
        ch_normal_samples
            .view{'Normal sample: ' + it}
        ch_tumor_samples
            .view{'Tumor sample: ' + it}
    ```

```bash
nextflow run main.nf
```

??? success "명령 출력"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [soggy_wiles] DSL2 - revision: 3bc1979889

    [patientA, [id:patientA, repeat:1, type:normal], patientA_rep1_normal.bam, [id:patientA, repeat:1, type:tumor], patientA_rep1_tumor.bam]
    [patientA, [id:patientA, repeat:2, type:normal], patientA_rep2_normal.bam, [id:patientA, repeat:2, type:tumor], patientA_rep2_tumor.bam]
    [patientB, [id:patientB, repeat:1, type:normal], patientB_rep1_normal.bam, [id:patientB, repeat:1, type:tumor], patientB_rep1_tumor.bam]
    [patientC, [id:patientC, repeat:1, type:normal], patientC_rep1_normal.bam, [id:patientC, repeat:1, type:tumor], patientC_rep1_tumor.bam]
    ```

너무 넓어서 읽기 어려울 수 있지만, `id` 필드를 기준으로 샘플이 결합된 것을 확인할 수 있습니다. 각 튜플의 형식은 다음과 같습니다.

- `id`: 샘플 ID
- `normal_meta_map`: 유형, 반복 번호, BAM 파일 경로를 포함한 정상 샘플 메타 데이터
- `normal_sample_file`: 정상 샘플 파일
- `tumor_meta_map`: 유형, 반복 번호, BAM 파일 경로를 포함한 종양 샘플 메타 데이터
- `tumor_sample`: 유형, 반복 번호, BAM 파일 경로를 포함한 종양 샘플

!!! warning "경고"

    `join` 연산자는 매칭되지 않는 튜플을 버립니다. 이 예제에서는 모든 샘플이 종양과 정상으로 매칭되도록 했지만, 그렇지 않은 경우 매칭되지 않는 튜플을 유지하려면 `remainder: true` 매개변수를 사용해야 합니다. 자세한 내용은 [문서](https://www.nextflow.io/docs/latest/operator.html#join)를 확인하세요.

이제 `map`을 사용하여 튜플에서 필드를 분리하는 방법과 `join`을 사용하여 첫 번째 필드를 기반으로 튜플을 결합하는 방법을 알게 되었습니다.
이 지식을 바탕으로 공유 필드를 기반으로 채널을 성공적으로 결합할 수 있습니다.

다음으로 여러 필드로 join하려는 상황을 살펴보겠습니다.

### 3.2. 여러 필드로 Join

sampleA에는 2개의 반복이 있지만 sampleB와 sampleC에는 1개만 있습니다. 이 경우 `id` 필드를 사용하여 효과적으로 join할 수 있었지만, 동기화가 맞지 않으면 어떻게 될까요? 서로 다른 반복의 정상 샘플과 종양 샘플이 혼합될 수 있습니다!

이를 방지하기 위해 여러 필드로 join할 수 있습니다. 이를 달성하는 방법은 여러 가지가 있지만, 샘플 `id`와 `replicate` 번호를 모두 포함하는 새 joining 키를 생성하는 방법에 집중합니다.

새 joining 키를 생성하는 것부터 시작합니다. [`map` 연산자](https://www.nextflow.io/docs/latest/operator.html#map)를 사용하여 `id`와 `repeat` 필드를 첫 번째 요소로 하는 새 튜플을 생성하면 됩니다.

=== "후"

    ```groovy title="main.nf" linenums="7" hl_lines="3 6"
        ch_normal_samples = ch_samples
            .filter { meta, file -> meta.type == 'normal' }
            .map { meta, file -> [[meta.id, meta.repeat], meta, file] }
        ch_tumor_samples = ch_samples
            .filter { meta, file -> meta.type == 'tumor' }
            .map { meta, file -> [[meta.id, meta.repeat], meta, file] }
    ```

=== "전"

    ```groovy title="main.nf" linenums="7" hl_lines="3 6"
        ch_normal_samples = ch_samples
            .filter { meta, file -> meta.type == 'normal' }
            .map { meta, file -> [meta.id, meta, file] }
        ch_tumor_samples = ch_samples
            .filter { meta, file -> meta.type == 'tumor' }
            .map { meta, file -> [meta.id, meta, file] }
    ```

이제 `id`와 `repeat` 필드를 모두 사용하여 join이 수행되는 것을 확인할 수 있습니다. 워크플로우를 실행합니다.

```bash
nextflow run main.nf
```

??? success "명령 출력"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [prickly_wing] DSL2 - revision: 3bebf22dee

    [[patientA, 1], [id:patientA, repeat:1, type:normal], patientA_rep1_normal.bam, [id:patientA, repeat:1, type:tumor], patientA_rep1_tumor.bam]
    [[patientA, 2], [id:patientA, repeat:2, type:normal], patientA_rep2_normal.bam, [id:patientA, repeat:2, type:tumor], patientA_rep2_tumor.bam]
    [[patientB, 1], [id:patientB, repeat:1, type:normal], patientB_rep1_normal.bam, [id:patientB, repeat:1, type:tumor], patientB_rep1_tumor.bam]
    [[patientC, 1], [id:patientC, repeat:1, type:normal], patientC_rep1_normal.bam, [id:patientC, repeat:1, type:tumor], patientC_rep1_tumor.bam]
    ```

각 결합 결과의 첫 번째 요소로 두 요소(`id`와 `repeat` 필드)로 구성된 튜플이 있음을 확인합니다. 이는 복잡한 항목을 joining 키로 사용하여 동일한 조건의 샘플 간에 상당히 복잡한 매칭을 가능하게 하는 방법을 보여줍니다.

다른 키로 join하는 방법을 더 알아보려면 [join 연산자 문서](https://www.nextflow.io/docs/latest/operator.html#join)에서 추가 옵션과 예제를 확인하세요.

### 3.3. `subMap`을 사용하여 새 joining 키 생성

이전 접근 방식은 joining 키에서 필드 이름을 잃게 됩니다. `id`와 `repeat` 필드가 단순한 값 목록이 됩니다. 나중에 이름으로 접근할 수 있도록 필드 이름을 유지하려면 [`subMap` 메서드](<https://docs.groovy-lang.org/latest/html/groovy-jdk/java/util/Map.html#subMap(java.util.Collection)>)를 사용할 수 있습니다.

`subMap` 메서드는 map에서 지정된 키-값 쌍만 추출합니다. 여기서는 joining 키를 생성하기 위해 `id`와 `repeat` 필드만 추출합니다.

=== "후"

    ```groovy title="main.nf" linenums="7" hl_lines="3 6"
        ch_normal_samples = ch_samples
            .filter { meta, file -> meta.type == 'normal' }
            .map { meta, file -> [meta.subMap(['id', 'repeat']), meta, file] }
        ch_tumor_samples = ch_samples
            .filter { meta, file -> meta.type == 'tumor' }
            .map { meta, file -> [meta.subMap(['id', 'repeat']), meta, file] }
    ```

=== "전"

    ```groovy title="main.nf" linenums="7" hl_lines="3 6"
        ch_normal_samples = ch_samples
            .filter { meta, file -> meta.type == 'normal' }
            .map { meta, file -> [[meta.id, meta.repeat], meta, file] }
        ch_tumor_samples = ch_samples
            .filter { meta, file -> meta.type == 'tumor' }
            .map { meta, file -> [[meta.id, meta.repeat], meta, file] }
    ```

```bash
nextflow run main.nf
```

??? success "명령 출력"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [reverent_wing] DSL2 - revision: 847016c3b7

    [[id:patientA, repeat:1], [id:patientA, repeat:1, type:normal], patientA_rep1_normal.bam, [id:patientA, repeat:1, type:tumor], patientA_rep1_tumor.bam]
    [[id:patientA, repeat:2], [id:patientA, repeat:2, type:normal], patientA_rep2_normal.bam, [id:patientA, repeat:2, type:tumor], patientA_rep2_tumor.bam]
    [[id:patientB, repeat:1], [id:patientB, repeat:1, type:normal], patientB_rep1_normal.bam, [id:patientB, repeat:1, type:tumor], patientB_rep1_tumor.bam]
    [[id:patientC, repeat:1], [id:patientC, repeat:1, type:normal], patientC_rep1_normal.bam, [id:patientC, repeat:1, type:tumor], patientC_rep1_tumor.bam]
    ```

이제 `id`와 `repeat` 필드를 포함할 뿐만 아니라 나중에 `meta.id`, `meta.repeat`와 같이 이름으로 접근할 수 있도록 필드 이름도 유지하는 새 joining 키가 생성되었습니다.

### 3.4. map에서 이름 있는 closure 사용

중복을 방지하고 오류를 줄이기 위해 이름 있는 closure를 사용할 수 있습니다. 이름 있는 closure를 사용하면 여러 곳에서 호출할 수 있는 재사용 가능한 함수를 만들 수 있습니다.

먼저 closure를 새 변수로 정의합니다.

=== "후"

    ```groovy title="main.nf" linenums="2" hl_lines="7"
        ch_samples = channel.fromPath("./data/samplesheet.csv")
            .splitCsv(header: true)
            .map{ row ->
                [[id:row.id, repeat:row.repeat, type:row.type], row.bam]
            }

        getSampleIdAndReplicate = { meta, bam -> [ meta.subMap(['id', 'repeat']), meta, file(bam) ] }

        ch_normal_samples = ch_samples
            .filter { meta, file -> meta.type == 'normal' }
    ```

=== "전"

    ```groovy title="main.nf" linenums="2"
        ch_samples = channel.fromPath("./data/samplesheet.csv")
            .splitCsv(header: true)
            .map{ row ->
                [[id:row.id, repeat:row.repeat, type:row.type], row.bam]
            }
        ch_normal_samples = ch_samples
            .filter { meta, file -> meta.type == 'normal' }
    ```

map 변환을 재사용 가능한 이름 있는 변수로 정의했습니다.

또한 `file()`을 사용하여 파일 경로를 Path 객체로 변환하여 이 채널을 받는 모든 프로세스가 파일을 올바르게 처리할 수 있도록 합니다(자세한 내용은 [파일 작업](../working_with_files/) 참조).

워크플로우에 closure를 적용합니다.

=== "후"

    ```groovy title="main.nf" linenums="10" hl_lines="3 6"
        ch_normal_samples = ch_samples
            .filter { meta, file -> meta.type == 'normal' }
             .map ( getSampleIdAndReplicate )
        ch_tumor_samples = ch_samples
            .filter { meta, file -> meta.type == 'tumor' }
             .map ( getSampleIdAndReplicate )

    ```

=== "전"

    ```groovy title="main.nf" linenums="10" hl_lines="3 6"
        ch_normal_samples = ch_samples
            .filter { meta, file -> meta.type == 'normal' }
            .map { meta, file -> [meta.subMap(['id', 'repeat']), meta, file] }
        ch_tumor_samples = ch_samples
            .filter { meta, file -> meta.type == 'tumor' }
            .map { meta, file -> [meta.subMap(['id', 'repeat']), meta, file] }
    ```

!!! note "참고"

    `map` 연산자가 `{ }`에서 `( )`로 변경되어 closure를 인자로 전달합니다. `map` 연산자는 closure를 인자로 기대하며 `{ }`는 익명 closure를 정의하는 데 사용됩니다. 이름 있는 closure를 호출할 때는 `( )` 구문을 사용합니다.

워크플로우를 다시 실행하여 모든 것이 여전히 작동하는지 확인합니다.

```bash
nextflow run main.nf
```

??? success "명령 출력"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [angry_meninsky] DSL2 - revision: 2edc226b1d

    [[id:patientA, repeat:1], [id:patientA, repeat:1, type:normal], patientA_rep1_normal.bam, [id:patientA, repeat:1, type:tumor], patientA_rep1_tumor.bam]
    [[id:patientA, repeat:2], [id:patientA, repeat:2, type:normal], patientA_rep2_normal.bam, [id:patientA, repeat:2, type:tumor], patientA_rep2_tumor.bam]
    [[id:patientB, repeat:1], [id:patientB, repeat:1, type:normal], patientB_rep1_normal.bam, [id:patientB, repeat:1, type:tumor], patientB_rep1_tumor.bam]
    [[id:patientC, repeat:1], [id:patientC, repeat:1, type:normal], patientC_rep1_normal.bam, [id:patientC, repeat:1, type:tumor], patientC_rep1_tumor.bam]
    ```

이름 있는 closure를 사용하면 여러 곳에서 동일한 변환을 재사용할 수 있어 오류 위험을 줄이고 코드를 더 읽기 쉽고 유지 관리하기 쉽게 만들 수 있습니다.

### 3.5. 데이터 중복 줄이기

워크플로우에 중복된 데이터가 많습니다. 결합된 샘플의 각 항목은 `id`와 `repeat` 필드를 반복합니다. 이 정보는 이미 그룹화 키에서 사용 가능하므로 이 중복을 피할 수 있습니다. 현재 데이터 구조는 다음과 같습니다.

```groovy
[
  [
    "id": "sampleC",
    "repeat": "1",
  ],
  [
    "id": "sampleC",
    "repeat": "1",
    "type": "normal",
  ],
  "sampleC_rep1_normal.bam"
  [
    "id": "sampleC",
    "repeat": "1",
    "type": "tumor",
  ],
  "sampleC_rep1_tumor.bam"
]
```

`id`와 `repeat` 필드가 그룹화 키에서 사용 가능하므로, 중복을 피하기 위해 각 채널 항목의 나머지 부분에서 이를 제거합니다. `subMap` 메서드를 사용하여 `type` 필드만 포함하는 새 map을 생성할 수 있습니다. 이 접근 방식을 통해 데이터 구조의 중복을 제거하면서 필요한 모든 정보를 유지할 수 있습니다.

=== "후"

    ```groovy title="main.nf" linenums="8" hl_lines="1"
        getSampleIdAndReplicate = { meta, bam -> [ meta.subMap(['id', 'repeat']), meta.subMap(['type']), file(bam) ] }
    ```

=== "전"

    ```groovy title="main.nf" linenums="8" hl_lines="1"
        getSampleIdAndReplicate = { meta, bam -> [ meta.subMap(['id', 'repeat']), meta, file(bam) ] }
    ```

이제 closure는 첫 번째 요소에 `id`와 `repeat` 필드가 포함되고 두 번째 요소에 `type` 필드만 포함된 튜플을 반환합니다. 그룹화 키에 `id`와 `repeat` 정보를 한 번만 저장하여 중복을 제거하면서 필요한 모든 정보를 유지합니다.

워크플로우를 실행하여 결과를 확인합니다.

```bash
nextflow run main.nf
```

??? success "명령 출력"

    ```console
    [[id:patientA, repeat:1], [type:normal], /workspaces/training/side-quests/splitting_and_grouping/patientA_rep1_normal.bam, [type:tumor], /workspaces/training/side-quests/splitting_and_grouping/patientA_rep1_tumor.bam]
    [[id:patientA, repeat:2], [type:normal], /workspaces/training/side-quests/splitting_and_grouping/patientA_rep2_normal.bam, [type:tumor], /workspaces/training/side-quests/splitting_and_grouping/patientA_rep2_tumor.bam]
    [[id:patientB, repeat:1], [type:normal], /workspaces/training/side-quests/splitting_and_grouping/patientB_rep1_normal.bam, [type:tumor], /workspaces/training/side-quests/splitting_and_grouping/patientB_rep1_tumor.bam]
    [[id:patientC, repeat:1], [type:normal], /workspaces/training/side-quests/splitting_and_grouping/patientC_rep1_normal.bam, [type:tumor], /workspaces/training/side-quests/splitting_and_grouping/patientC_rep1_tumor.bam]
    ```

그룹화 키에 `id`와 `repeat` 필드를 한 번만 명시하고 샘플 데이터에 `type` 필드가 있음을 확인할 수 있습니다. 정보를 잃지 않으면서 채널 내용을 더 간결하게 만들었습니다.

### 3.6. 불필요한 정보 제거

위에서 중복된 정보를 제거했지만, 채널에 여전히 불필요한 정보가 있습니다.

처음에 `filter`를 사용하여 정상 샘플과 종양 샘플을 분리한 다음 `id`와 `repeat` 키를 기반으로 join했습니다. `join` 연산자는 튜플이 병합되는 순서를 유지하므로, 왼쪽에 정상 샘플, 오른쪽에 종양 샘플이 있는 경우 결과 채널은 `id, <정상 요소>, <종양 요소>` 구조를 유지합니다.

채널의 각 요소 위치를 알고 있으므로 `[type:normal]`과 `[type:tumor]` 메타데이터를 제거하여 구조를 더 단순화할 수 있습니다.

=== "후"

    ```groovy title="main.nf" linenums="8" hl_lines="1"
        getSampleIdAndReplicate = { meta, file -> [ meta.subMap(['id', 'repeat']), file ] }
    ```

=== "전"

    ```groovy title="main.nf" linenums="8" hl_lines="1"
        getSampleIdAndReplicate = { meta, file -> [ meta.subMap(['id', 'repeat']), meta.subMap(['type']), file ] }
    ```

다시 실행하여 결과를 확인합니다.

```bash
nextflow run main.nf
```

??? success "명령 출력"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [confident_leavitt] DSL2 - revision: a2303895bd

    [[id:patientA, repeat:1], patientA_rep1_normal.bam, patientA_rep1_tumor.bam]
    [[id:patientA, repeat:2], patientA_rep2_normal.bam, patientA_rep2_tumor.bam]
    [[id:patientB, repeat:1], patientB_rep1_normal.bam, patientB_rep1_tumor.bam]
    [[id:patientC, repeat:1], patientC_rep1_normal.bam, patientC_rep1_tumor.bam]
    ```

### 핵심 정리

이 섹션에서 학습한 내용:

- **튜플 조작**: `map`을 사용하여 튜플에서 필드를 분리하는 방법
- **튜플 결합**: `join`을 사용하여 첫 번째 필드를 기반으로 튜플을 결합하는 방법
- **Joining 키 생성**: `subMap`을 사용하여 새 joining 키를 생성하는 방법
- **이름 있는 Closure**: map에서 이름 있는 closure를 사용하는 방법
- **여러 필드 Join**: 더 정확한 매칭을 위해 여러 필드로 join하는 방법
- **데이터 구조 최적화**: 불필요한 정보를 제거하여 채널 구조를 간소화하는 방법

이제 샘플시트를 분할하고, 정상 샘플과 종양 샘플을 필터링하고, 샘플 ID와 반복 번호로 join한 다음 결과를 출력하는 워크플로우가 완성되었습니다.

이는 독립적으로 처리한 후 샘플이나 다른 유형의 데이터를 매칭해야 하는 생물정보학 워크플로우에서 일반적인 패턴이므로 유용한 기술입니다. 다음으로 샘플을 여러 번 반복하는 방법을 살펴보겠습니다.

## 4. 샘플을 구간에 분배

생물정보학 워크플로우의 핵심 패턴은 유전체 영역에 걸쳐 분석을 분배하는 것입니다. 예를 들어, 변이 호출은 게놈을 구간(염색체 또는 더 작은 영역)으로 나누어 병렬화할 수 있습니다. 이 병렬화 전략은 여러 코어 또는 노드에 계산 부하를 분산하여 전체 실행 시간을 줄임으로써 파이프라인 효율성을 크게 향상시킵니다.

다음 섹션에서는 샘플 데이터를 여러 유전체 구간에 분배하는 방법을 보여줍니다. 각 샘플을 모든 구간과 쌍으로 만들어 서로 다른 유전체 영역의 병렬 처리를 가능하게 합니다. 이렇게 하면 구간 수만큼 데이터셋 크기가 늘어나 나중에 다시 합칠 수 있는 여러 독립적인 분석 단위가 생성됩니다.

### 4.1. `combine`을 사용하여 구간에 샘플 분배

구간 채널을 생성하는 것부터 시작합니다. 간단하게 수동으로 정의한 3개의 구간을 사용합니다. 실제 워크플로우에서는 파일 입력에서 읽거나 많은 구간 파일이 있는 채널을 생성할 수도 있습니다.

=== "후"

    ```groovy title="main.nf" linenums="17" hl_lines="2"
            .join(ch_tumor_samples)
        ch_intervals = channel.of('chr1', 'chr2', 'chr3')
    ```

=== "전"

    ```groovy title="main.nf" linenums="17" hl_lines="2"
            .join(ch_tumor_samples)
        ch_joined_samples.view()
    ```

각 샘플을 각 구간에 대해 반복하려 합니다. 이를 샘플과 구간의 카르테시안 곱이라고도 합니다. [`combine` 연산자](https://www.nextflow.io/docs/latest/operator.html#combine)를 사용하여 이를 달성할 수 있습니다. 채널 1의 모든 항목을 가져와 채널 2의 각 항목에 대해 반복합니다. 워크플로우에 combine 연산자를 추가합니다.

=== "후"

    ```groovy title="main.nf" linenums="18" hl_lines="3-5"
        ch_intervals = channel.of('chr1', 'chr2', 'chr3')

        ch_combined_samples = ch_joined_samples
            .combine(ch_intervals)
            .view()
    ```

=== "전"

    ```groovy title="main.nf" linenums="18"
        ch_intervals = channel.of('chr1', 'chr2', 'chr3')
    ```

실행하여 결과를 확인합니다.

```bash
nextflow run main.nf
```

??? success "명령 출력"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [mighty_tesla] DSL2 - revision: ae013ab70b

    [[id:patientA, repeat:1], patientA_rep1_normal.bam, patientA_rep1_tumor.bam, chr1]
    [[id:patientA, repeat:1], patientA_rep1_normal.bam, patientA_rep1_tumor.bam, chr2]
    [[id:patientA, repeat:1], patientA_rep1_normal.bam, patientA_rep1_tumor.bam, chr3]
    [[id:patientA, repeat:2], patientA_rep2_normal.bam, patientA_rep2_tumor.bam, chr1]
    [[id:patientA, repeat:2], patientA_rep2_normal.bam, patientA_rep2_tumor.bam, chr2]
    [[id:patientA, repeat:2], patientA_rep2_normal.bam, patientA_rep2_tumor.bam, chr3]
    [[id:patientB, repeat:1], patientB_rep1_normal.bam, patientB_rep1_tumor.bam, chr1]
    [[id:patientB, repeat:1], patientB_rep1_normal.bam, patientB_rep1_tumor.bam, chr2]
    [[id:patientB, repeat:1], patientB_rep1_normal.bam, patientB_rep1_tumor.bam, chr3]
    [[id:patientC, repeat:1], patientC_rep1_normal.bam, patientC_rep1_tumor.bam, chr1]
    [[id:patientC, repeat:1], patientC_rep1_normal.bam, patientC_rep1_tumor.bam, chr2]
    [[id:patientC, repeat:1], patientC_rep1_normal.bam, patientC_rep1_tumor.bam, chr3]
    ```

성공입니다! 3개의 구간 목록에서 모든 구간에 대해 모든 샘플을 반복했습니다. 채널의 항목 수가 효과적으로 3배가 되었습니다.

읽기가 조금 어려우므로 다음 섹션에서 정리하겠습니다.

### 4.2. 채널 구성

`map` 연산자를 사용하여 샘플 데이터를 정리하고 리팩토링하여 이해하기 쉽게 만들 수 있습니다. 구간 string을 첫 번째 요소의 joining map으로 이동합니다.

=== "후"

    ```groovy title="main.nf" linenums="20" hl_lines="3-9"
        ch_combined_samples = ch_joined_samples
            .combine(ch_intervals)
            .map { grouping_key, normal, tumor, interval ->
                [
                    grouping_key + [interval: interval],
                    normal,
                    tumor
                ]
            }
            .view()
    ```

=== "전"

    ```groovy title="main.nf" linenums="20"
        ch_combined_samples = ch_joined_samples
            .combine(ch_intervals)
            .view()
    ```

이 map 작업이 단계별로 수행하는 작업을 분석해 보겠습니다.

먼저 이름 있는 매개변수를 사용하여 코드를 더 읽기 쉽게 만듭니다. `grouping_key`, `normal`, `tumor`, `interval`이라는 이름을 사용하면 인덱스 대신 이름으로 튜플의 요소를 참조할 수 있습니다.

```groovy
        .map { grouping_key, normal, tumor, interval ->
```

다음으로 `grouping_key`와 `interval` 필드를 결합합니다. `grouping_key`는 `id`와 `repeat` 필드를 포함하는 map입니다. Groovy의 map 덧셈(`+`)을 사용하여 `interval`이 있는 새 map을 생성하고 병합합니다.

```groovy
                grouping_key + [interval: interval],
```

마지막으로 결합된 메타데이터 map, 정상 샘플 파일, 종양 샘플 파일의 세 요소로 구성된 튜플로 반환합니다.

```groovy
            [
                grouping_key + [interval: interval],
                normal,
                tumor
            ]
```

다시 실행하여 채널 내용을 확인합니다.

```bash
nextflow run main.nf
```

??? success "명령 출력"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [sad_hawking] DSL2 - revision: 1f6f6250cd

    [[id:patientA, repeat:1, interval:chr1], patientA_rep1_normal.bam, patientA_rep1_tumor.bam]
    [[id:patientA, repeat:1, interval:chr2], patientA_rep1_normal.bam, patientA_rep1_tumor.bam]
    [[id:patientA, repeat:1, interval:chr3], patientA_rep1_normal.bam, patientA_rep1_tumor.bam]
    [[id:patientA, repeat:2, interval:chr1], patientA_rep2_normal.bam, patientA_rep2_tumor.bam]
    [[id:patientA, repeat:2, interval:chr2], patientA_rep2_normal.bam, patientA_rep2_tumor.bam]
    [[id:patientA, repeat:2, interval:chr3], patientA_rep2_normal.bam, patientA_rep2_tumor.bam]
    [[id:patientB, repeat:1, interval:chr1], patientB_rep1_normal.bam, patientB_rep1_tumor.bam]
    [[id:patientB, repeat:1, interval:chr2], patientB_rep1_normal.bam, patientB_rep1_tumor.bam]
    [[id:patientB, repeat:1, interval:chr3], patientB_rep1_normal.bam, patientB_rep1_tumor.bam]
    [[id:patientC, repeat:1, interval:chr1], patientC_rep1_normal.bam, patientC_rep1_tumor.bam]
    [[id:patientC, repeat:1, interval:chr2], patientC_rep1_normal.bam, patientC_rep1_tumor.bam]
    [[id:patientC, repeat:1, interval:chr3], patientC_rep1_normal.bam, patientC_rep1_tumor.bam]
    ```

`map`을 사용하여 데이터를 올바른 구조로 변환하는 것은 까다로울 수 있지만, 효과적인 데이터 조작에 매우 중요합니다.

이제 모든 샘플이 모든 유전체 구간에 걸쳐 반복되어 병렬로 처리할 수 있는 여러 독립적인 분석 단위가 생성되었습니다. 하지만 관련 샘플을 다시 합치려면 어떻게 해야 할까요? 다음 섹션에서는 공통 속성을 공유하는 샘플을 그룹화하는 방법을 학습합니다.

### 핵심 정리

이 섹션에서 학습한 내용:

- **구간에 샘플 분배**: `combine`을 사용하여 구간에 샘플을 반복하는 방법
- **카르테시안 곱 생성**: 샘플과 구간의 모든 조합을 생성하는 방법
- **채널 구조 구성**: `map`을 사용하여 가독성을 위해 데이터를 재구성하는 방법
- **병렬 처리 준비**: 분산 분석을 위한 데이터 설정 방법

## 5. `groupTuple`을 사용하여 샘플 집계

이전 섹션에서는 입력 파일에서 데이터를 분할하고 특정 필드(정상 및 종양 샘플)로 필터링하는 방법을 학습했습니다. 하지만 이는 단일 유형의 결합만 다룹니다. 특정 속성으로 샘플을 그룹화하려면 어떻게 해야 할까요? 예를 들어, 매칭된 정상-종양 쌍을 join하는 대신 유형에 관계없이 "sampleA"의 모든 샘플을 함께 처리하려 할 수 있습니다. 이 패턴은 마지막에 결과를 비교하거나 결합하기 전에 효율성을 위해 관련 샘플을 독립적으로 처리하려는 생물정보학 워크플로우에서 일반적입니다.

Nextflow에는 이를 위한 내장 메서드가 있으며, 주로 살펴볼 것은 `groupTuple`입니다.

동일한 `id`와 `interval` 필드를 가진 모든 샘플을 그룹화하는 것부터 시작합니다. 이는 기술적 반복을 그룹화하되 의미 있게 다른 샘플은 분리된 상태로 유지하려는 분석에서 일반적입니다.

이를 위해 그룹화 변수를 분리하여 독립적으로 사용할 수 있도록 해야 합니다.

첫 번째 단계는 이전 섹션에서 수행한 것과 유사합니다. 그룹화 변수를 튜플의 첫 번째 요소로 분리해야 합니다. 현재 첫 번째 요소는 `id`, `repeat`, `interval` 필드의 map임을 기억하세요.

```groovy title="main.nf" linenums="1"
{
  "id": "sampleA",
  "repeat": "1",
  "interval": "chr1"
}
```

이전의 `subMap` 메서드를 재사용하여 map에서 `id`와 `interval` 필드를 분리할 수 있습니다. 이전과 마찬가지로 `map` 연산자를 사용하여 각 샘플의 튜플 첫 번째 요소에 `subMap` 메서드를 적용합니다.

=== "후"

    ```groovy title="main.nf" linenums="20" hl_lines="11-19"
        ch_combined_samples = ch_joined_samples
            .combine(ch_intervals)
            .map { grouping_key, normal, tumor, interval ->
                [
                    grouping_key + [interval: interval],
                    normal,
                    tumor
                ]
            }

        ch_grouped_samples = ch_combined_samples
            .map { grouping_key, normal, tumor ->
                [
                    grouping_key.subMap('id', 'interval'),
                    normal,
                    tumor
                ]
              }
              .view()
    ```

=== "전"

    ```groovy title="main.nf" linenums="20" hl_lines="10"
        ch_combined_samples = ch_joined_samples
            .combine(ch_intervals)
            .map { grouping_key, normal, tumor, interval ->
                [
                    grouping_key + [interval: interval],
                    normal,
                    tumor
                ]
            }
            .view()
    ```

다시 실행하여 채널 내용을 확인합니다.

```bash
nextflow run main.nf
```

??? success "명령 출력"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [hopeful_brenner] DSL2 - revision: 7f4f7fea76

    [[id:patientA, interval:chr1], patientA_rep1_normal.bam, patientA_rep1_tumor.bam]
    [[id:patientA, interval:chr2], patientA_rep1_normal.bam, patientA_rep1_tumor.bam]
    [[id:patientA, interval:chr3], patientA_rep1_normal.bam, patientA_rep1_tumor.bam]
    [[id:patientA, interval:chr1], patientA_rep2_normal.bam, patientA_rep2_tumor.bam]
    [[id:patientA, interval:chr2], patientA_rep2_normal.bam, patientA_rep2_tumor.bam]
    [[id:patientA, interval:chr3], patientA_rep2_normal.bam, patientA_rep2_tumor.bam]
    [[id:patientB, interval:chr1], patientB_rep1_normal.bam, patientB_rep1_tumor.bam]
    [[id:patientB, interval:chr2], patientB_rep1_normal.bam, patientB_rep1_tumor.bam]
    [[id:patientB, interval:chr3], patientB_rep1_normal.bam, patientB_rep1_tumor.bam]
    [[id:patientC, interval:chr1], patientC_rep1_normal.bam, patientC_rep1_tumor.bam]
    [[id:patientC, interval:chr2], patientC_rep1_normal.bam, patientC_rep1_tumor.bam]
    [[id:patientC, interval:chr3], patientC_rep1_normal.bam, patientC_rep1_tumor.bam]
    ```

`id`와 `interval` 필드를 성공적으로 분리했지만 아직 샘플을 그룹화하지 않았습니다.

!!! note "참고"

    여기서 `replicate` 필드를 버립니다. 이후 다운스트림 처리에 필요하지 않기 때문입니다. 이 튜토리얼을 완료한 후 나중의 그룹화에 영향을 주지 않고 이를 포함할 수 있는지 확인해 보세요!

이제 [`groupTuple` 연산자](https://www.nextflow.io/docs/latest/operator.html#grouptuple)를 사용하여 이 새 그룹화 요소로 샘플을 그룹화합니다.

=== "후"

    ```groovy title="main.nf" linenums="30" hl_lines="9"
        ch_grouped_samples = ch_combined_samples
            .map { grouping_key, normal, tumor ->
                [
                    grouping_key.subMap('id', 'interval'),
                    normal,
                    tumor
                ]
              }
              .groupTuple()
              .view()
    ```

=== "전"

    ```groovy title="main.nf" linenums="30"
        ch_grouped_samples = ch_combined_samples
            .map { grouping_key, normal, tumor ->
                [
                    grouping_key.subMap('id', 'interval'),
                    normal,
                    tumor
                ]
              }
              .view()
    ```

이것이 전부입니다! 코드 한 줄만 추가했습니다. 실행하면 어떻게 되는지 확인합니다.

```bash
nextflow run main.nf
```

??? success "명령 출력"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [friendly_jang] DSL2 - revision: a1bee1c55d

    [[id:patientA, interval:chr1], [patientA_rep1_normal.bam, patientA_rep2_normal.bam], [patientA_rep1_tumor.bam, patientA_rep2_tumor.bam]]
    [[id:patientA, interval:chr2], [patientA_rep1_normal.bam, patientA_rep2_normal.bam], [patientA_rep1_tumor.bam, patientA_rep2_tumor.bam]]
    [[id:patientA, interval:chr3], [patientA_rep1_normal.bam, patientA_rep2_normal.bam], [patientA_rep1_tumor.bam, patientA_rep2_tumor.bam]]
    [[id:patientB, interval:chr1], [patientB_rep1_normal.bam], [patientB_rep1_tumor.bam]]
    [[id:patientB, interval:chr2], [patientB_rep1_normal.bam], [patientB_rep1_tumor.bam]]
    [[id:patientB, interval:chr3], [patientB_rep1_normal.bam], [patientB_rep1_tumor.bam]]
    [[id:patientC, interval:chr1], [patientC_rep1_normal.bam], [patientC_rep1_tumor.bam]]
    [[id:patientC, interval:chr2], [patientC_rep1_normal.bam], [patientC_rep1_tumor.bam]]
    [[id:patientC, interval:chr3], [patientC_rep1_normal.bam], [patientC_rep1_tumor.bam]]
    ```

데이터 구조가 변경되었으며 각 채널 요소 내의 파일이 이제 `[patientA_rep1_normal.bam, patientA_rep2_normal.bam]`과 같은 튜플에 포함되어 있습니다. `groupTuple`을 사용하면 Nextflow가 그룹의 각 샘플에 대한 단일 파일을 결합하기 때문입니다. 다운스트림에서 데이터를 처리할 때 이 점을 기억하는 것이 중요합니다.

!!! note "참고"

    [`transpose`](https://www.nextflow.io/docs/latest/reference/operator.html#transpose)는 groupTuple의 반대입니다. 채널의 항목을 풀고 평탄화합니다. `transpose`를 추가하여 위에서 수행한 그룹화를 되돌려 보세요!

### 핵심 정리

이 섹션에서 학습한 내용:

- **관련 샘플 그룹화**: `groupTuple`을 사용하여 공통 속성으로 샘플을 집계하는 방법
- **그룹화 키 분리**: `subMap`을 사용하여 그룹화를 위한 특정 필드를 추출하는 방법
- **그룹화된 데이터 구조 처리**: `groupTuple`로 생성된 내포된 구조로 작업하는 방법
- **기술적 반복 처리**: 동일한 실험 조건을 공유하는 샘플을 그룹화하는 방법

---

## 요약

이 사이드 퀘스트에서는 채널을 사용하여 데이터를 분할하고 그룹화하는 방법을 학습했습니다.

파이프라인을 통해 흐르는 데이터를 수정함으로써 루프나 while 문을 사용하지 않고 확장 가능한 파이프라인을 구축할 수 있으며, 이는 더 전통적인 접근 방식에 비해 여러 장점을 제공합니다.

- 추가 코드 없이 원하는 만큼 많거나 적은 입력으로 확장 가능
- 반복 대신 파이프라인을 통한 데이터 흐름 처리에 집중
- 필요에 따라 복잡하거나 단순하게 구성 가능
- 파이프라인이 더 선언적이 되어 어떻게 해야 하는지보다 무엇이 일어나야 하는지에 집중
- Nextflow가 독립적인 작업을 병렬로 실행하여 실행을 최적화

이러한 채널 작업을 마스터하면 루프나 반복 프로그래밍에 의존하지 않고 복잡한 데이터 관계를 처리하는 유연하고 확장 가능한 파이프라인을 구축할 수 있으며, Nextflow가 실행을 최적화하고 독립적인 작업을 자동으로 병렬화할 수 있습니다.

### 핵심 패턴

1.  **구조화된 입력 데이터 생성:** meta map이 있는 CSV 파일에서 시작([워크플로우의 메타데이터](../metadata/)의 패턴 기반)

    ```groovy
    ch_samples = channel.fromPath("./data/samplesheet.csv")
        .splitCsv(header: true)
        .map{ row ->
          [[id:row.id, repeat:row.repeat, type:row.type], row.bam]
        }
    ```

2.  **데이터를 별도 채널로 분할:** `filter`를 사용하여 `type` 필드를 기반으로 데이터를 독립적인 스트림으로 분리

    ```groovy
    channel.filter { it.type == 'tumor' }
    ```

3.  **매칭된 샘플 결합:** `join`을 사용하여 `id`와 `repeat` 필드를 기반으로 관련 샘플을 재결합

    - 키(튜플의 첫 번째 요소)로 두 채널 join

    ```groovy
    tumor_ch.join(normal_ch)
    ```

    - joining 키를 추출하고 이 값으로 join

    ```groovy
    tumor_ch.map { meta, file -> [meta.id, meta, file] }
        .join(
          normal_ch.map { meta, file -> [meta.id, meta, file] }
        )
    ```

    - subMap을 사용하여 여러 필드로 join

    ```groovy
    tumor_ch.map { meta, file -> [meta.subMap(['id', 'repeat']), meta, file] }
        .join(
          normal_ch.map { meta, file -> [meta.subMap(['id', 'repeat']), meta, file] }
        )
    ```

4.  **구간에 분배:** `combine`을 사용하여 병렬 처리를 위해 유전체 구간과 샘플의 카르테시안 곱을 생성

    ```groovy
    samples_ch.combine(intervals_ch)
    ```

5.  **그룹화 키로 집계:** `groupTuple`을 사용하여 각 튜플의 첫 번째 요소로 그룹화하여 `id`와 `interval` 필드를 공유하는 샘플을 수집하고 기술적 반복을 병합

    ```groovy
    channel.groupTuple()
    ```

6.  **데이터 구조 최적화:** `subMap`을 사용하여 특정 필드를 추출하고 변환을 재사용 가능하게 만들기 위한 이름 있는 closure를 생성

    - map에서 특정 필드 추출

    ```groovy
    meta.subMap(['id', 'repeat'])
    ```

    - 재사용 가능한 변환을 위한 이름 있는 closure 사용

    ```groovy
    getSampleIdAndReplicate = { meta, file -> [meta.subMap(['id', 'repeat']), file] }
    channel.map(getSampleIdAndReplicate)
    ```

### 추가 자료

- [filter](https://www.nextflow.io/docs/latest/operator.html#filter)
- [map](https://www.nextflow.io/docs/latest/operator.html#map)
- [join](https://www.nextflow.io/docs/latest/operator.html#join)
- [groupTuple](https://www.nextflow.io/docs/latest/operator.html#grouptuple)
- [combine](https://www.nextflow.io/docs/latest/operator.html#combine)

---

## 다음 단계

[사이드 퀘스트 메뉴](../)로 돌아가거나 페이지 오른쪽 하단의 버튼을 클릭하여 목록의 다음 주제로 이동하세요.
