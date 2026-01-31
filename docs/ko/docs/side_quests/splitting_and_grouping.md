# 데이터 분할 및 그룹화

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } AI 지원 번역 - [자세히 알아보기 및 개선 제안](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Nextflow는 데이터를 유연하게 다룰 수 있는 강력한 도구를 제공합니다. 주요 기능 중 하나는 데이터를 여러 스트림으로 분할한 다음 관련 항목을 다시 그룹화하는 것입니다. 이는 특히 생물정보학 워크플로에서 유용합니다. 이러한 워크플로에서는 서로 다른 유형의 샘플을 별도로 처리한 후 분석을 위해 결과를 결합해야 합니다.

우편물을 분류하는 것으로 생각해 보십시오. 목적지별로 편지를 분리하고, 각 묶음을 다르게 처리한 다음, 동일한 사람에게 보내지는 항목을 다시 결합합니다. Nextflow는 과학 데이터로 이를 수행하기 위해 특수 연산자를 사용합니다. 이 접근 방식은 분산 컴퓨팅 및 생물정보학 워크플로에서 일반적으로 **scatter/gather** 패턴으로 알려져 있습니다.

Nextflow의 채널 시스템은 이러한 유연성의 핵심입니다. 채널은 워크플로의 여러 부분을 연결하여 데이터가 분석을 통해 흐를 수 있도록 합니다. 단일 데이터 소스에서 여러 채널을 생성하고, 각 채널을 다르게 처리한 다음, 필요할 때 채널을 다시 병합할 수 있습니다. 이 접근 방식을 사용하면 복잡한 생물정보학 분석의 분기 및 수렴 경로를 자연스럽게 반영하는 워크플로를 설계할 수 있습니다.

### 학습 목표

이 사이드 퀘스트에서는 Nextflow의 채널 연산자를 사용하여 데이터를 분할하고 그룹화하는 방법을 배웁니다.
샘플 정보 및 관련 데이터 파일이 포함된 CSV 파일로 시작한 다음 이 데이터를 조작하고 재구성합니다.

이 사이드 퀘스트가 끝나면 다음 기술을 사용하여 데이터 스트림을 효과적으로 분리하고 결합할 수 있습니다:

- `splitCsv`를 사용하여 파일에서 데이터 읽기
- `filter` 및 `map`으로 데이터 필터링 및 변환
- `join` 및 `groupTuple`을 사용하여 관련 데이터 결합
- 병렬 처리를 위해 `combine`으로 데이터 조합 생성
- `subMap` 및 중복 제거 전략을 사용하여 데이터 구조 최적화
- 명명된 클로저를 사용하여 채널 구조를 조작하는 데 도움이 되는 재사용 가능한 함수 구축

이러한 기술은 깨끗하고 유지 관리 가능한 코드 구조를 유지하면서 여러 입력 파일과 다양한 유형의 데이터를 효율적으로 처리할 수 있는 워크플로를 구축하는 데 도움이 됩니다.

### 전제 조건

이 사이드 퀘스트를 시작하기 전에 다음을 수행해야 합니다:

- [Hello Nextflow](../hello_nextflow/README.md) 튜토리얼 또는 이에 상응하는 초급 과정을 완료했습니다.
- 기본적인 Nextflow 개념 및 메커니즘(프로세스, 채널, 연산자, 파일 작업, 메타데이터)을 사용하는 데 익숙합니다.

**선택 사항:** 먼저 [워크플로의 메타데이터](./metadata.md) 사이드 퀘스트를 완료하는 것이 좋습니다.
여기에서는 `splitCsv`로 CSV 파일을 읽고 메타 맵을 생성하는 기본 사항을 다룹니다. 이는 여기에서 많이 사용할 것입니다.

---

## 0. 시작하기

#### 교육 codespace 열기

아직 수행하지 않았다면 [환경 설정](../envsetup/index.md)에 설명된 대로 교육 환경을 열어야 합니다.

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

#### 프로젝트 디렉토리로 이동

이 튜토리얼의 파일이 있는 디렉토리로 이동하겠습니다.

```bash
cd side-quests/splitting_and_grouping
```

VSCode가 이 디렉토리에 집중하도록 설정할 수 있습니다:

```bash
code .
```

#### 자료 검토

메인 워크플로 파일과 `samplesheet.csv`라는 샘플시트가 포함된 `data` 디렉토리를 찾을 수 있습니다.

```console title="디렉토리 내용"
.
├── data
│   └── samplesheet.csv
└── main.nf
```

샘플시트에는 환자 ID, 샘플 반복 번호, 유형(정상 또는 종양) 및 가상 데이터 파일 경로(실제로 존재하지 않지만 존재한다고 가정함)를 포함하여 여러 환자의 샘플에 대한 정보가 포함되어 있습니다.

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

이 샘플시트에는 세 명의 환자(A, B, C)로부터 얻은 8개의 샘플이 나열되어 있습니다.

각 환자에 대해 `tumor`(일반적으로 종양 생검에서 유래) 또는 `normal`(건강한 조직 또는 혈액에서 채취) 유형의 샘플이 있습니다.
암 분석에 익숙하지 않다면 이것이 대조 분석을 수행하기 위해 쌍을 이루는 종양/정상 샘플을 사용하는 실험 모델에 해당한다는 것만 알아두십시오.

특히 환자 A의 경우 두 세트의 기술적 복제(반복)가 있습니다.

!!! note "참고"

    이 실험 설계에 익숙하지 않더라도 걱정하지 마십시오. 이 튜토리얼을 이해하는 데 중요하지 않습니다.

#### 과제 검토

귀하의 과제는 다음을 수행하는 Nextflow 워크플로를 작성하는 것입니다:

1. CSV 파일에서 샘플 데이터를 **읽고** 메타 맵으로 구조화합니다
2. 유형(정상 대 종양)에 따라 샘플을 여러 채널로 **분리**합니다
3. 환자 ID 및 복제 번호로 일치하는 종양/정상 쌍을 **조인**합니다
4. 병렬 처리를 위해 게놈 간격에 샘플을 **분산**합니다
5. 다운스트림 분석을 위해 관련 샘플을 다시 **그룹화**합니다

이것은 독립적인 처리를 위해 데이터를 분할한 다음 비교 분석을 위해 관련 항목을 다시 결합해야 하는 일반적인 생물정보학 패턴을 나타냅니다.

#### 준비 상태 체크리스트

준비가 되었다고 생각하십니까?

- [ ] 이 과정의 목표와 전제 조건을 이해합니다
- [ ] codespace가 작동 중입니다
- [ ] 작업 디렉토리를 적절하게 설정했습니다
- [ ] 과제를 이해합니다

모든 항목을 체크할 수 있다면 시작할 준비가 된 것입니다.

---

## 1. 샘플 데이터 읽기

### 1.1. `splitCsv`로 샘플 데이터를 읽고 메타 맵 생성

먼저 `splitCsv`로 샘플 데이터를 읽고 메타 맵 패턴으로 구성해 보겠습니다. `main.nf`에서 워크플로를 이미 시작했음을 알 수 있습니다.

```groovy title="main.nf" linenums="1" hl_lines="2"
workflow {
    ch_samplesheet = channel.fromPath("./data/samplesheet.csv")
}
```

!!! note "참고"

    이 튜토리얼 전체에서 모든 채널 변수에 `ch_` 접두사를 사용하여 Nextflow 채널임을 명확하게 표시합니다.

[워크플로의 메타데이터](./metadata.md) 사이드 퀘스트를 완료했다면 이 패턴을 알아볼 것입니다. `splitCsv`를 사용하여 CSV를 읽고 즉시 메타 맵으로 데이터를 구조화하여 메타데이터를 파일 경로와 분리합니다.

!!! info

    이 교육에서 `map`이라는 두 가지 개념을 접하게 됩니다:

    - **데이터 구조**: 키-값 쌍을 저장하는 Groovy 맵(다른 언어의 딕셔너리/해시에 해당)
    - **채널 연산자**: 채널의 항목을 변환하는 `.map()` 연산자

    컨텍스트에서 어떤 것을 의미하는지 명확히 하겠지만 이 구분은 Nextflow 작업 시 이해하는 것이 중요합니다.

`main.nf`에 다음 변경 사항을 적용하십시오:

=== "수정 후"

    ```groovy title="main.nf" linenums="2" hl_lines="2-6"
        ch_samples = channel.fromPath("./data/samplesheet.csv")
            .splitCsv(header: true)
            .map{ row ->
              [[id:row.id, repeat:row.repeat, type:row.type], row.bam]
            }
            .view()
    ```

=== "수정 전"

    ```groovy title="main.nf" linenums="2" hl_lines="1"
        ch_samplesheet = channel.fromPath("./data/samplesheet.csv")
    ```

이것은 `splitCsv` 작업(헤더가 있는 CSV 읽기)과 `map` 작업(데이터를 `[meta, file]` 튜플로 구조화)을 한 단계로 결합합니다. 해당 변경 사항을 적용하고 파이프라인을 실행하십시오:

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

이제 각 항목이 `[meta, file]` 튜플인 채널이 있습니다. 메타데이터가 파일 경로와 분리되었습니다. 이 구조를 통해 메타데이터 필드를 기반으로 작업을 분할하고 그룹화할 수 있습니다.

---

## 2. 데이터 필터링 및 변환

### 2.1. `filter`로 데이터 필터링

[`filter` 연산자](https://www.nextflow.io/docs/latest/operator.html#filter)를 사용하여 조건에 따라 데이터를 필터링할 수 있습니다. 정상 샘플만 처리하려고 한다고 가정해 보겠습니다. `type` 필드를 기반으로 데이터를 필터링하여 이를 수행할 수 있습니다. `view` 연산자 앞에 이것을 삽입해 보겠습니다.

=== "수정 후"

    ```groovy title="main.nf" linenums="2" hl_lines="6"
        ch_samples = channel.fromPath("./data/samplesheet.csv")
            .splitCsv(header: true)
            .map{ row ->
              [[id:row.id, repeat:row.repeat, type:row.type], row.bam]
            }
            .filter { meta, file -> meta.type == 'normal' }
            .view()
    ```

=== "수정 전"

    ```groovy title="main.nf" linenums="2"
        ch_samples = channel.fromPath("./data/samplesheet.csv")
            .splitCsv(header: true)
            .map{ row ->
              [[id:row.id, repeat:row.repeat, type:row.type], row.bam]
            }
            .view()
    ```

워크플로를 다시 실행하여 필터링된 결과를 확인하십시오:

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

정상 샘플만 포함하도록 데이터를 성공적으로 필터링했습니다. 이것이 어떻게 작동하는지 요약해 보겠습니다.

`filter` 연산자는 채널의 각 요소에 적용되는 클로저를 사용합니다. 클로저가 `true`를 반환하면 요소가 포함되고, `false`를 반환하면 요소가 제외됩니다.

우리의 경우 `meta.type == 'normal'`인 샘플만 유지하려고 합니다. 클로저는 튜플 `meta,file`을 사용하여 각 샘플을 참조하고, `meta.type`으로 샘플 유형에 액세스하고, 그것이 `'normal'`과 같은지 확인합니다.

이는 위에서 소개한 단일 클로저로 수행됩니다:

```groovy title="main.nf" linenums="7"
    .filter { meta, file -> meta.type == 'normal' }
```

### 2.2. 별도의 필터링된 채널 생성

현재 CSV에서 직접 생성된 채널에 필터를 적용하고 있지만 이것을 하나 이상의 방법으로 필터링하려고 하므로 정상 샘플에 대한 별도의 필터링된 채널을 생성하도록 로직을 다시 작성해 보겠습니다:

=== "수정 후"

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

=== "수정 전"

    ```groovy title="main.nf" linenums="2"
        ch_samples = channel.fromPath("./data/samplesheet.csv")
            .splitCsv(header: true)
            .map{ row ->
              [[id:row.id, repeat:row.repeat, type:row.type], row.bam]
            }
            .filter { meta, file -> meta.type == 'normal' }
            .view()
    ```

파이프라인을 실행하여 결과를 확인하십시오:

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

데이터를 성공적으로 필터링하고 정상 샘플에 대한 별도의 채널을 생성했습니다.

종양 샘플에 대한 필터링된 채널도 생성해 보겠습니다:

=== "수정 후"

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

=== "수정 전"

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

정상 및 종양 샘플을 두 개의 다른 채널로 분리하고 `view()`에 제공된 클로저를 사용하여 출력에서 다르게 레이블을 지정했습니다: `ch_tumor_samples.view{'Tumor sample: ' + it}`.

### 핵심 요약

이 섹션에서는 다음을 배웠습니다:

- **데이터 필터링**: `filter`로 데이터를 필터링하는 방법
- **데이터 분할**: 조건에 따라 데이터를 여러 채널로 분할하는 방법
- **데이터 보기**: `view`를 사용하여 데이터를 출력하고 여러 채널의 출력에 레이블을 지정하는 방법

이제 정상 및 종양 샘플을 두 개의 다른 채널로 분리했습니다. 다음으로 `id` 필드에서 정상 및 종양 샘플을 조인하겠습니다.

---

## 3. 식별자로 채널 조인

이전 섹션에서는 정상 및 종양 샘플을 두 개의 다른 채널로 분리했습니다. 이들은 유형에 따라 특정 프로세스 또는 워크플로를 사용하여 독립적으로 처리될 수 있습니다. 그러나 동일한 환자의 정상 및 종양 샘플을 비교하려면 어떻게 해야 합니까? 이 시점에서 `id` 필드를 기반으로 샘플을 일치시켜 다시 조인해야 합니다.

Nextflow에는 채널을 결합하는 여러 방법이 포함되어 있지만 이 경우 가장 적절한 연산자는 [`join`](https://www.nextflow.io/docs/latest/operator.html#join)입니다. SQL에 익숙하다면 `JOIN` 연산처럼 작동하며, 조인할 키와 수행할 조인 유형을 지정합니다.

### 3.1. `map` 및 `join`을 사용하여 환자 ID를 기반으로 결합

[`join`](https://www.nextflow.io/docs/latest/operator.html#join) 문서를 확인하면 기본적으로 각 튜플의 첫 번째 항목을 기반으로 두 채널을 조인한다는 것을 알 수 있습니다.

#### 3.1.1. 데이터 구조 확인

콘솔 출력을 계속 사용할 수 없는 경우 파이프라인을 실행하여 데이터 구조를 확인하고 `id` 필드에서 조인하도록 수정해야 하는 방법을 살펴보겠습니다.

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

`id` 필드가 각 메타 맵의 첫 번째 요소임을 알 수 있습니다. `join`이 작동하려면 각 튜플에서 `id` 필드를 분리해야 합니다. 그런 다음 `join` 연산자를 사용하여 두 채널을 결합하면 됩니다.

#### 3.1.2. `id` 필드 분리

`id` 필드를 분리하려면 [`map` 연산자](https://www.nextflow.io/docs/latest/operator.html#map)를 사용하여 `id` 필드를 첫 번째 요소로 하는 새 튜플을 생성할 수 있습니다.

=== "수정 후"

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

=== "수정 전"

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

미묘할 수 있지만 각 튜플의 첫 번째 요소가 `id` 필드임을 알 수 있어야 합니다.

#### 3.1.3. 두 채널 결합

이제 `join` 연산자를 사용하여 `id` 필드를 기반으로 두 채널을 결합할 수 있습니다.

다시 한 번 `view`를 사용하여 조인된 출력을 출력합니다.

=== "수정 후"

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

=== "수정 전"

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

너무 넓어서 알기 어렵지만 샘플이 `id` 필드로 조인되었음을 알 수 있어야 합니다. 이제 각 튜플은 다음 형식을 갖습니다:

- `id`: 샘플 ID
- `normal_meta_map`: 유형, 복제 및 bam 파일 경로를 포함한 정상 샘플 메타데이터
- `normal_sample_file`: 정상 샘플 파일
- `tumor_meta_map`: 유형, 복제 및 bam 파일 경로를 포함한 종양 샘플 메타데이터
- `tumor_sample`: 유형, 복제 및 bam 파일 경로를 포함한 종양 샘플

!!! warning "경고"

    `join` 연산자는 일치하지 않는 튜플을 삭제합니다. 이 예제에서는 모든 샘플이 종양 및 정상에 대해 일치하는지 확인했지만 그렇지 않은 경우 일치하지 않는 튜플을 유지하려면 매개변수 `remainder: true`를 사용해야 합니다. 자세한 내용은 [문서](https://www.nextflow.io/docs/latest/operator.html#join)를 확인하십시오.

이제 `map`을 사용하여 튜플의 필드를 분리하는 방법과 `join`을 사용하여 첫 번째 필드를 기반으로 튜플을 결합하는 방법을 알았습니다.
이 지식을 통해 공유 필드를 기반으로 채널을 성공적으로 결합할 수 있습니다.

다음으로 여러 필드에서 조인하려는 상황을 고려하겠습니다.

### 3.2. 여러 필드에서 조인

sampleA에는 2개의 복제가 있지만 sampleB와 sampleC에는 1개만 있습니다. 이 경우 `id` 필드를 사용하여 효과적으로 조인할 수 있었지만 동기화되지 않으면 어떻게 될까요? 서로 다른 복제의 정상 및 종양 샘플을 혼동할 수 있습니다!

이를 방지하려면 여러 필드에서 조인할 수 있습니다. 실제로 이를 달성하는 여러 방법이 있지만 샘플 `id`와 `replicate` 번호를 모두 포함하는 새 조인 키를 만드는 데 집중하겠습니다.

새 조인 키를 만드는 것부터 시작하겠습니다. 이전과 같은 방법으로 [`map` 연산자](https://www.nextflow.io/docs/latest/operator.html#map)를 사용하여 `id` 및 `repeat` 필드를 첫 번째 요소로 하는 새 튜플을 생성할 수 있습니다.

=== "수정 후"

    ```groovy title="main.nf" linenums="7" hl_lines="3 6"
        ch_normal_samples = ch_samples
            .filter { meta, file -> meta.type == 'normal' }
            .map { meta, file -> [[meta.id, meta.repeat], meta, file] }
        ch_tumor_samples = ch_samples
            .filter { meta, file -> meta.type == 'tumor' }
            .map { meta, file -> [[meta.id, meta.repeat], meta, file] }
    ```

=== "수정 전"

    ```groovy title="main.nf" linenums="7" hl_lines="3 6"
        ch_normal_samples = ch_samples
            .filter { meta, file -> meta.type == 'normal' }
            .map { meta, file -> [meta.id, meta, file] }
        ch_tumor_samples = ch_samples
            .filter { meta, file -> meta.type == 'tumor' }
            .map { meta, file -> [meta.id, meta, file] }
    ```

이제 `id` 및 `repeat` 필드를 모두 사용하여 조인이 발생하는 것을 볼 수 있어야 합니다. 워크플로를 실행하십시오:

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

각 조인 결과의 첫 번째 요소로 두 요소(`id` 및 `repeat` 필드)의 튜플이 어떻게 있는지 주목하십시오. 이것은 복잡한 항목을 조인 키로 사용할 수 있음을 보여주며 동일한 조건의 샘플 간에 상당히 복잡한 일치를 가능하게 합니다.

다른 키에서 조인하는 더 많은 방법을 탐색하려면 추가 옵션 및 예제에 대한 [join 연산자 문서](https://www.nextflow.io/docs/latest/operator.html#join)를 확인하십시오.

### 3.3. `subMap`을 사용하여 새 조인 키 생성

이전 접근 방식은 조인 키에서 필드 이름을 잃어버립니다. `id` 및 `repeat` 필드는 값 목록이 됩니다. 나중에 액세스할 수 있도록 필드 이름을 유지하려면 [`subMap` 메서드](<https://docs.groovy-lang.org/latest/html/groovy-jdk/java/util/Map.html#subMap(java.util.Collection)>)를 사용할 수 있습니다.

`subMap` 메서드는 맵에서 지정된 키-값 쌍만 추출합니다. 여기서는 `id` 및 `repeat` 필드만 추출하여 조인 키를 만듭니다.

=== "수정 후"

    ```groovy title="main.nf" linenums="7" hl_lines="3 6"
        ch_normal_samples = ch_samples
            .filter { meta, file -> meta.type == 'normal' }
            .map { meta, file -> [meta.subMap(['id', 'repeat']), meta, file] }
        ch_tumor_samples = ch_samples
            .filter { meta, file -> meta.type == 'tumor' }
            .map { meta, file -> [meta.subMap(['id', 'repeat']), meta, file] }
    ```

=== "수정 전"

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

이제 `id` 및 `repeat` 필드를 포함할 뿐만 아니라 필드 이름도 유지하여 나중에 이름으로 액세스할 수 있는 새 조인 키가 있습니다(예: `meta.id` 및 `meta.repeat`).

### 3.4. map에서 명명된 클로저 사용

중복을 피하고 오류를 줄이기 위해 명명된 클로저를 사용할 수 있습니다. 명명된 클로저를 사용하면 여러 곳에서 호출할 수 있는 재사용 가능한 함수를 만들 수 있습니다.

이를 위해 먼저 클로저를 새 변수로 정의합니다:

=== "수정 후"

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

=== "수정 전"

    ```groovy title="main.nf" linenums="2"
        ch_samples = channel.fromPath("./data/samplesheet.csv")
            .splitCsv(header: true)
            .map{ row ->
                [[id:row.id, repeat:row.repeat, type:row.type], row.bam]
            }
        ch_normal_samples = ch_samples
            .filter { meta, file -> meta.type == 'normal' }
    ```

재사용할 수 있는 명명된 변수로 맵 변환을 정의했습니다.

또한 `file()`을 사용하여 파일 경로를 Path 객체로 변환하여 이 채널을 수신하는 모든 프로세스가 파일을 올바르게 처리할 수 있도록 합니다(자세한 내용은 [파일 작업](./working_with_files.md) 참조).

워크플로에서 클로저를 구현해 보겠습니다:

=== "수정 후"

    ```groovy title="main.nf" linenums="10" hl_lines="3 6"
        ch_normal_samples = ch_samples
            .filter { meta, file -> meta.type == 'normal' }
             .map ( getSampleIdAndReplicate )
        ch_tumor_samples = ch_samples
            .filter { meta, file -> meta.type == 'tumor' }
             .map ( getSampleIdAndReplicate )

    ```

=== "수정 전"

    ```groovy title="main.nf" linenums="10" hl_lines="3 6"
        ch_normal_samples = ch_samples
            .filter { meta, file -> meta.type == 'normal' }
            .map { meta, file -> [meta.subMap(['id', 'repeat']), meta, file] }
        ch_tumor_samples = ch_samples
            .filter { meta, file -> meta.type == 'tumor' }
            .map { meta, file -> [meta.subMap(['id', 'repeat']), meta, file] }
    ```

!!! note "참고"

    `map` 연산자는 클로저를 인수로 전달하기 위해 `{ }`를 사용하는 것에서 `( )`를 사용하는 것으로 전환되었습니다. `map` 연산자는 클로저를 인수로 예상하고 `{ }`는 익명 클로저를 정의하는 데 사용되기 때문입니다. 명명된 클로저를 호출할 때는 `( )` 구문을 사용하십시오.

모든 것이 여전히 작동하는지 확인하기 위해 워크플로를 한 번 더 실행하십시오:

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

명명된 클로저를 사용하면 여러 곳에서 동일한 변환을 재사용할 수 있으므로 오류 위험이 줄어들고 코드의 가독성과 유지 관리성이 향상됩니다.

### 3.5. 데이터 중복 감소

워크플로에 많은 중복 데이터가 있습니다. 조인된 샘플의 각 항목은 `id` 및 `repeat` 필드를 반복합니다. 이 정보는 이미 그룹화 키에서 사용할 수 있으므로 이러한 중복을 피할 수 있습니다. 다시 말해서, 현재 데이터 구조는 다음과 같습니다:

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

`id` 및 `repeat` 필드는 그룹화 키에서 사용할 수 있으므로 중복을 방지하기 위해 각 채널 항목의 나머지 부분에서 제거해 보겠습니다. `subMap` 메서드를 사용하여 `type` 필드만 있는 새 맵을 생성하면 됩니다. 이 접근 방식을 통해 데이터 구조의 중복을 제거하면서 필요한 모든 정보를 유지할 수 있습니다.

=== "수정 후"

    ```groovy title="main.nf" linenums="8" hl_lines="1"
        getSampleIdAndReplicate = { meta, bam -> [ meta.subMap(['id', 'repeat']), meta.subMap(['type']), file(bam) ] }
    ```

=== "수정 전"

    ```groovy title="main.nf" linenums="8" hl_lines="1"
        getSampleIdAndReplicate = { meta, bam -> [ meta.subMap(['id', 'repeat']), meta, file(bam) ] }
    ```

이제 클로저는 첫 번째 요소에 `id` 및 `repeat` 필드가 포함되고 두 번째 요소에는 `type` 필드만 포함된 튜플을 반환합니다. 이는 `id` 및 `repeat` 정보를 그룹화 키에 한 번 저장하여 중복을 제거하면서 필요한 모든 정보를 유지합니다.

워크플로를 실행하여 어떻게 보이는지 확인하십시오:

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

그룹화 키에서 `id` 및 `repeat` 필드를 한 번만 명시하고 샘플 데이터에 `type` 필드가 있음을 알 수 있습니다. 정보를 잃지 않았지만 채널 내용을 더 간결하게 만들었습니다.

### 3.6. 중복 정보 제거

위에서 중복된 정보를 제거했지만 채널에는 여전히 다른 중복 정보가 있습니다.

처음에 우리는 `filter`를 사용하여 정상 및 종양 샘플을 분리한 다음 `id` 및 `repeat` 키를 기반으로 조인했습니다. `join` 연산자는 튜플이 병합되는 순서를 보존하므로 우리의 경우 왼쪽에 정상 샘플이 있고 오른쪽에 종양 샘플이 있는 경우 결과 채널은 이 구조를 유지합니다: `id, <정상 요소>, <종양 요소>`.

채널에서 각 요소의 위치를 알고 있으므로 `[type:normal]` 및 `[type:tumor]` 메타데이터를 삭제하여 구조를 더 단순화할 수 있습니다.

=== "수정 후"

    ```groovy title="main.nf" linenums="8" hl_lines="1"
        getSampleIdAndReplicate = { meta, file -> [ meta.subMap(['id', 'repeat']), file ] }
    ```

=== "수정 전"

    ```groovy title="main.nf" linenums="8" hl_lines="1"
        getSampleIdAndReplicate = { meta, file -> [ meta.subMap(['id', 'repeat']), meta.subMap(['type']), file ] }
    ```

다시 실행하여 결과를 확인하십시오:

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

### 핵심 요약

이 섹션에서는 다음을 배웠습니다:

- **튜플 조작**: `map`을 사용하여 튜플의 필드를 분리하는 방법
- **튜플 조인**: `join`을 사용하여 첫 번째 필드를 기반으로 튜플을 결합하는 방법
- **조인 키 생성**: `subMap`을 사용하여 새 조인 키를 생성하는 방법
- **명명된 클로저**: map에서 명명된 클로저를 사용하는 방법
- **여러 필드 조인**: 더 정확한 일치를 위해 여러 필드에서 조인하는 방법
- **데이터 구조 최적화**: 중복 정보를 제거하여 채널 구조를 간소화하는 방법

이제 샘플시트를 분할하고, 정상 및 종양 샘플을 필터링하고, 샘플 ID 및 복제 번호로 함께 조인한 다음 결과를 출력할 수 있는 워크플로가 있습니다.

이것은 독립적으로 처리한 후 샘플 또는 기타 유형의 데이터를 일치시켜야 하는 생물정보학 워크플로의 일반적인 패턴이므로 유용한 기술입니다. 다음으로 샘플을 여러 번 반복하는 방법을 살펴보겠습니다.

## 4. 샘플을 간격으로 분산

생물정보학 워크플로의 핵심 패턴은 게놈 영역에 분석을 분산하는 것입니다. 예를 들어, 변이 호출은 게놈을 간격(염색체 또는 더 작은 영역과 같은)으로 나누어 병렬화할 수 있습니다. 이 병렬화 전략은 계산 부하를 여러 코어 또는 노드에 분산시켜 전체 실행 시간을 줄임으로써 파이프라인 효율성을 크게 향상시킵니다.

다음 섹션에서는 샘플 데이터를 여러 게놈 간격에 분산하는 방법을 보여줍니다. 각 샘플을 모든 간격과 쌍으로 묶어 다른 게놈 영역의 병렬 처리를 허용합니다. 이렇게 하면 간격 수만큼 데이터 세트 크기가 배가되어 나중에 다시 가져올 수 있는 여러 독립적인 분석 단위가 생성됩니다.

### 4.1. `combine`을 사용하여 간격에 샘플 분산

간격 채널을 만드는 것부터 시작하겠습니다. 간단하게 유지하기 위해 수동으로 정의할 3개의 간격만 사용합니다. 실제
