# 메타데이터와 메타 맵

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } AI 지원 번역 - [자세히 알아보기 및 개선 제안](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

과학적 분석에서 우리는 원시 데이터 파일만으로 작업하는 경우가 거의 없습니다.
각 파일에는 고유한 추가 정보가 함께 제공됩니다: 그것이 무엇인지, 어디서 왔는지, 무엇이 특별한지 등입니다.
이러한 추가 정보를 우리는 메타데이터라고 부릅니다.

메타데이터는 다른 데이터를 설명하는 데이터입니다.
메타데이터는 파일과 실험 조건에 대한 중요한 세부 정보를 추적하고 각 데이터세트의 고유한 특성에 맞게 분석을 조정하는 데 도움이 됩니다.

도서관 목록과 같다고 생각하시면 됩니다: 책에는 실제 콘텐츠(원시 데이터)가 포함되어 있지만, 목록 카드는 각 책에 대한 필수 정보(언제 출판되었는지, 누가 썼는지, 어디서 찾을 수 있는지 등 메타데이터)를 제공합니다.
Nextflow 파이프라인에서 메타데이터는 다음과 같은 용도로 사용될 수 있습니다:

- 워크플로 전체에서 파일별 정보 추적
- 파일 특성에 따라 프로세스 구성
- 공동 분석을 위해 관련 파일 그룹화

### 학습 목표

이 사이드 퀘스트에서는 워크플로에서 메타데이터를 처리하는 방법을 탐구합니다.
기본 파일 정보가 포함된 간단한 데이터시트(생물정보학에서는 종종 샘플시트라고 함)부터 시작하여 다음을 학습합니다:

- CSV 파일에서 파일 메타데이터 읽기 및 파싱
- 메타데이터 맵 생성 및 조작
- 워크플로 실행 중 새 메타데이터 필드 추가
- 메타데이터를 사용하여 프로세스 동작 사용자 정의

이러한 기술은 복잡한 파일 관계와 처리 요구 사항을 처리할 수 있는 더 견고하고 유연한 파이프라인을 구축하는 데 도움이 됩니다.

### 사전 요구 사항

이 사이드 퀘스트를 시작하기 전에 다음을 충족해야 합니다:

- [Hello Nextflow](../hello_nextflow/README.md) 튜토리얼 또는 동등한 초급 과정을 완료했어야 합니다.
- 기본 Nextflow 개념과 메커니즘(프로세스, 채널, 연산자) 사용에 익숙해야 합니다.

---

## 0. 시작하기

#### 교육 코드스페이스 열기

아직 하지 않았다면 [환경 설정](../envsetup/index.md)에 설명된 대로 교육 환경을 열어야 합니다.

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

#### 프로젝트 디렉토리로 이동

이 튜토리얼의 파일이 있는 디렉토리로 이동하겠습니다.

```bash
cd side-quests/metadata
```

VSCode가 이 디렉토리에 집중하도록 설정할 수 있습니다:

```bash
code .
```

#### 자료 검토

메인 워크플로 파일과 데이터시트 및 여러 데이터 파일이 포함된 `data` 디렉토리를 찾을 수 있습니다.

??? abstract "디렉토리 내용"

    ```console
    .
    ├── data
    │   ├── bonjour.txt
    │   ├── ciao.txt
    │   ├── guten_tag.txt
    │   ├── hallo.txt
    │   ├── hello.txt
    │   ├── hola.txt
    │   ├── salut.txt
    │   └── datasheet.csv
    ├── main.nf
    └── nextflow.config
    ```

`main.nf` 파일의 워크플로는 단계적으로 완전히 작동하는 워크플로로 확장할 스텁입니다.

데이터시트는 데이터 파일의 경로와 관련 메타데이터를 3개의 열로 구성하여 나열합니다:

- `id`: 설명이 필요 없는 파일에 부여된 ID
- `character`: 캐릭터 이름으로, 나중에 다양한 생물을 그리는 데 사용됩니다
- `data`: 다양한 언어의 인사말이 포함된 `.txt` 파일 경로

```console title="datasheet.csv"
id,character,recording
sampleA,squirrel,/workspaces/training/side-quests/metadata/data/bonjour.txt
sampleB,tux,/workspaces/training/side-quests/metadata/data/guten_tag.txt
sampleC,sheep,/workspaces/training/side-quests/metadata/data/hallo.txt
sampleD,turkey,/workspaces/training/side-quests/metadata/data/hello.txt
sampleE,stegosaurus,/workspaces/training/side-quests/metadata/data/hola.txt
sampleF,moose,/workspaces/training/side-quests/metadata/data/salut.txt
sampleG,turtle,/workspaces/training/side-quests/metadata/data/ciao.txt
```

각 데이터 파일에는 5개 언어(fr: 프랑스어, de: 독일어, es: 스페인어, it: 이탈리아어, en: 영어) 중 하나로 된 인사말 텍스트가 포함되어 있습니다.

또한 `langid`라는 컨테이너화된 언어 분석 도구도 제공됩니다.

#### 과제 검토

귀하의 과제는 다음을 수행하는 Nextflow 워크플로를 작성하는 것입니다:

1. 각 파일의 언어를 자동으로 **식별**
2. 언어 계열(게르만어 대 로망스어)별로 파일 **그룹화**
3. 언어 및 메타데이터에 따라 각 파일의 처리 **사용자 정의**
4. 언어 그룹별로 출력 **구성**

이것은 파일별 메타데이터가 처리 결정을 주도하는 일반적인 워크플로 패턴을 나타내며, 메타데이터 맵이 우아하게 해결하는 문제 유형입니다.

#### 준비 체크리스트

시작할 준비가 되었다고 생각하십니까?

- [ ] 이 과정의 목표와 사전 요구 사항을 이해합니다
- [ ] 코드스페이스가 정상적으로 실행되고 있습니다
- [ ] 작업 디렉토리를 적절하게 설정했습니다
- [ ] 과제를 이해합니다

모든 항목을 체크할 수 있다면 시작할 준비가 된 것입니다.

---

## 1. 데이터시트에서 메타데이터 로드

`main.nf` 워크플로 파일을 열어 시작점으로 제공하는 워크플로 스텁을 확인하십시오.

```groovy title="main.nf" linenums="1"
#!/usr/bin/env nextflow

workflow  {

    ch_datasheet = channel.fromPath("./data/datasheet.csv")

}
```

예제 데이터시트를 파일로 로드하는 기본 채널 팩토리를 설정했지만, 아직 파일 내용을 읽지는 않습니다.
먼저 이것을 추가하겠습니다.

### 1.1. `splitCsv`로 내용 읽기

최소한의 노력으로 파일 내용을 적절하게 파싱하는 연산자를 선택해야 합니다.
데이터시트가 CSV 형식이므로 이것은 [`splitCsv`](https://www.nextflow.io/docs/latest/reference/operator.html#splitcsv) 연산자의 역할입니다. 이 연산자는 파일의 각 행을 채널의 요소로 로드합니다.

다음과 같이 변경하여 채널 구성 코드에 `splitCsv()` 작업을 추가하고, 파일 내용이 채널에 올바르게 로드되고 있는지 확인하기 위한 `view()` 작업을 추가하십시오.

=== "수정 후"

    ```groovy title="main.nf" linenums="3" hl_lines="4-5"
    workflow  {

        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .view()

    }
    ```

=== "수정 전"

    ```groovy title="main.nf" linenums="3"
    workflow  {

        ch_datasheet = channel.fromPath("./data/datasheet.csv")

    }
    ```

`header: true` 옵션을 사용하여 Nextflow에게 CSV 파일의 첫 번째 행을 헤더 행으로 읽도록 지시하고 있습니다.

결과를 확인해 볼까요?
워크플로를 실행합니다:

```bash
nextflow run main.nf
```

??? success "명령 출력"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [exotic_albattani] DSL2 - revision: c0d03cec83

    [id:sampleA, character:squirrel, recording:/workspaces/training/side-quests/metadata/data/bonjour.txt]
    [id:sampleB, character:tux, recording:/workspaces/training/side-quests/metadata/data/guten_tag.txt]
    [id:sampleC, character:sheep, recording:/workspaces/training/side-quests/metadata/data/hallo.txt]
    [id:sampleD, character:turkey, recording:/workspaces/training/side-quests/metadata/data/hello.txt]
    [id:sampleE, character:stegosaurus, recording:/workspaces/training/side-quests/metadata/data/hola.txt]
    [id:sampleF, character:moose, recording:/workspaces/training/side-quests/metadata/data/salut.txt]
    [id:sampleG, character:turtle, recording:/workspaces/training/side-quests/metadata/data/ciao.txt]
    ```

연산자가 CSV 파일의 각 행에 대한 키-값 쌍 맵을 구성했으며, 열 헤더를 해당 값의 키로 사용한 것을 볼 수 있습니다.

각 맵 항목은 데이터시트의 열에 해당합니다:

- `id`
- `character`
- `recording`

훌륭합니다! 이렇게 하면 각 파일의 특정 필드에 쉽게 액세스할 수 있습니다.
예를 들어 `id`로 파일 ID에 액세스하거나 `recording`으로 txt 파일 경로에 액세스할 수 있습니다.

??? info "(선택 사항) 맵에 대한 추가 정보"

    Groovy는 Nextflow가 구축된 프로그래밍 언어로, 맵은 Python의 딕셔너리, JavaScript의 객체 또는 Ruby의 해시와 유사한 키-값 데이터 구조입니다.

    다음은 실제로 맵을 정의하고 그 내용에 액세스하는 방법을 보여주는 실행 가능한 스크립트입니다:

    ```groovy title="examples/map_demo.nf"
    #!/usr/bin/env nextflow

    // 간단한 맵 생성
    def my_map = [id:'sampleA', character:'squirrel']

    // 전체 맵 출력
    println "map: ${my_map}"

    // 점 표기법을 사용하여 개별 값 액세스
    println "id: ${my_map.id}"
    println "character: ${my_map.character}"
    ```

    적절한 `workflow` 블록이 없더라도 Nextflow는 이것을 워크플로처럼 실행할 수 있습니다:

    ```bash
    nextflow run examples/map_demo.nf
    ```

    출력에서 다음을 볼 수 있습니다:

    ```console title="출력"
     N E X T F L O W   ~  version 25.10.2

    Launching `map_demo.nf` [cheesy_plateau] DSL2 - revision: fae5b8496e

    map: [id:sampleA, character:squirrel]
    id: sampleA
    character: squirrel
    ```

### 1.2. `map`으로 특정 필드 선택

데이터시트에서 `character` 열에 액세스하여 출력하고 싶다고 가정해 봅시다.
Nextflow `map` 연산자를 사용하여 채널의 각 항목을 반복하고 맵 객체에서 `character` 항목을 구체적으로 선택할 수 있습니다.

워크플로를 다음과 같이 편집하십시오:

=== "수정 후"

    ```groovy title="main.nf" linenums="3" hl_lines="5-7"
    workflow  {

        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .map{ row ->
                row.character
            }
            .view()

    }
    ```

=== "수정 전"

    ```groovy title="main.nf" linenums="3"
    workflow  {

        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .view()

    }
    ```

이제 워크플로를 다시 실행하십시오:

```bash
nextflow run main.nf
```

??? success "명령 출력"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [exotic_albattani] DSL2 - revision: c0d03cec83

    squirrel
    tux
    sheep
    turkey
    stegosaurus
    moose
    turtle
    ```

성공입니다! 데이터시트에서 파생된 맵 구조를 활용하여 각 행의 개별 열 값에 액세스했습니다.

이제 데이터시트를 성공적으로 읽었고 각 행의 데이터에 액세스할 수 있으므로 파이프라인 로직 구현을 시작할 수 있습니다.

### 1.3. 메타데이터를 '메타 맵'으로 구성

워크플로의 현재 상태에서 입력 파일(`recording` 키 아래)과 관련 메타데이터(`id`, `character`)는 모두 동일한 위치에 있습니다. 마치 모두 하나의 큰 가방에 있는 것과 같습니다.
실질적인 결과는 이 채널을 사용하는 모든 프로세스가 이 구조를 염두에 두고 구성되어야 한다는 것입니다:

```groovy
    input:
    tuple val(id), val(character), file(recording)
```

데이터시트의 열 수가 변경되지 않는 한 괜찮습니다.
그러나 데이터시트에 하나의 열만 추가해도 채널의 형태가 프로세스가 예상하는 것과 더 이상 일치하지 않아 워크플로에서 오류가 발생합니다.
또한 프로세스를 공유하기 어렵게 만들며, 스크립트 블록에서 필요하지 않은 변수를 프로세스에 하드코딩해야 할 수도 있습니다.

이 문제를 피하려면 데이터시트에 포함된 열 수에 관계없이 채널 구조를 일관되게 유지하는 방법을 찾아야 합니다.

튜플 내의 한 항목에 모든 메타데이터를 수집하여 이를 수행할 수 있습니다. 이를 메타데이터 맵 또는 더 간단히 '메타 맵'이라고 합니다.

`map` 작업을 다음과 같이 편집하십시오:

=== "수정 후"

    ```groovy title="main.nf" linenums="5" hl_lines="4"
        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .map { row ->
                    [ [id: row.id, character: row.character], row.recording ]
            }
            .view()
    ```

=== "수정 전"

    ```groovy title="main.nf" linenums="5" hl_lines="4"
        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .map{ row ->
                row.character
            }
            .view()
    ```

채널 요소를 메타 맵과 해당 파일 객체라는 두 요소로 구성된 튜플로 재구성했습니다.

워크플로를 실행해 봅시다:

```bash
nextflow run main.nf
```

??? success "명령 출력"

    ```console title="메타 맵 보기"
     N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [lethal_booth] DSL2 - revision: 0d8f844c07

    [[id:sampleA, character:squirrel], /workspaces/training/side-quests/metadata/data/bonjour.txt]
    [[id:sampleB, character:tux], /workspaces/training/side-quests/metadata/data/guten_tag.txt]
    [[id:sampleC, character:sheep], /workspaces/training/side-quests/metadata/data/hallo.txt]
    [[id:sampleD, character:turkey], /workspaces/training/side-quests/metadata/data/hello.txt]
    [[id:sampleE, character:stegosaurus], /workspaces/training/side-quests/metadata/data/hola.txt]
    [[id:sampleF, character:moose], /workspaces/training/side-quests/metadata/data/salut.txt]
    [[id:sampleG, character:turtle], /workspaces/training/side-quests/metadata/data/ciao.txt]
    ```

이제 채널의 각 요소에는 메타데이터 맵이 먼저 나오고 해당 파일 객체가 두 번째로 포함됩니다:

```console title="예제 출력 구조"
[
  [id:sampleA, character:squirrel],
  /workspaces/training/side-quests/metadata/data/bonjour.txt
]
```

결과적으로 데이터시트에 더 많은 열을 추가하면 `meta` 맵에서 더 많은 메타데이터를 사용할 수 있지만 채널 형태는 변경되지 않습니다.
이를 통해 메타데이터 항목을 입력 사양에 하드코딩하지 않고도 채널을 사용하는 프로세스를 작성할 수 있습니다:

```groovy title="구문 예제"
    input:
    tuple val(meta), file(recording)
```

이것은 Nextflow 워크플로에서 메타데이터를 구성하는 데 널리 사용되는 규칙입니다.

### 요점

이 섹션에서 다음을 학습했습니다:

- **메타데이터가 중요한 이유:** 데이터와 함께 메타데이터를 유지하면 워크플로 전체에서 중요한 파일 정보가 보존됩니다.
- **데이터시트 읽는 방법:** `splitCsv`를 사용하여 헤더 정보가 있는 CSV 파일을 읽고 행을 구조화된 데이터로 변환하는 방법
- **메타 맵 생성 방법:** 튜플 구조 `[ [id:value, ...], file ]`을 사용하여 메타데이터를 파일 데이터에서 분리하는 방법

---

## 2. 메타데이터 조작

이제 메타데이터를 로드했으니 활용해 봅시다!

각 생물의 녹음 파일에 포함된 언어를 식별하기 위해 [`langid`](https://github.com/saffsd/langid.py)라는 도구를 사용할 것입니다.
이 도구는 언어 세트에 대해 사전 학습되어 있으며, 텍스트 조각이 주어지면 언어 예측과 관련 확률 점수를 모두 `stdout`으로 출력합니다.

### 2.1. 프로세스 가져오기 및 코드 검토

`langid` 도구를 적용하는 `IDENTIFY_LANGUAGE`라는 사전 작성된 프로세스 모듈을 제공하므로 workflow 블록 앞에 include 문을 추가하기만 하면 됩니다.

워크플로를 다음과 같이 편집하십시오:

=== "수정 후"

    ```groovy title="main.nf" linenums="1" hl_lines="3"
    #!/usr/bin/env nextflow

    include { IDENTIFY_LANGUAGE } from './modules/langid.nf'

    workflow {
    ```

=== "수정 전"

    ```groovy title="main.nf" linenums="1"
    #!/usr/bin/env nextflow

    workflow {
    ```

모듈 파일을 열어 코드를 확인할 수 있습니다:

```groovy title="modules/langid.nf" linenums="1" hl_lines="9 12"
#!/usr/bin/env nextflow

// langid를 사용하여 각 입력 파일의 언어 예측
process IDENTIFY_LANGUAGE {

    container 'community.wave.seqera.io/library/pip_langid:b2269f456a5629ff'

    input:
    tuple val(meta), path(file)

    output:
    tuple val(meta), path(file), stdout

    script:
    """
    langid < ${file} -l en,de,fr,es,it | sed -E "s/.*\\('([a-z]+)'.*/\\1/" | tr -d '\\n'
    """
}
```

보시다시피 입력 정의는 입력 채널에 방금 적용한 것과 동일한 `tuple val(meta), path(file)` 구조를 사용합니다.

출력 정의는 입력과 유사한 구조의 튜플로 구성되지만 세 번째 요소로 `stdout`도 포함합니다.
이 `tuple val(meta), path(file), <output>` 패턴은 파이프라인을 통해 흐르는 메타데이터를 입력 데이터 및 출력과 연관시켜 유지합니다.

도구가 파일을 작성하지 않고 출력을 콘솔에 직접 출력하기 때문에 Nextflow의 [`stdout`](https://www.nextflow.io/docs/latest/process.html#outputs) 출력 한정자를 사용하고 있습니다. 그리고 명령줄에서 `sed`를 사용하여 확률 점수를 제거하고 개행 문자를 제거하여 문자열을 정리하고 언어 예측만 반환합니다.

### 2.2. `IDENTIFY_LANGUAGE` 호출 추가

이제 프로세스를 워크플로에서 사용할 수 있으므로 데이터 채널에서 실행하기 위해 `IDENTIFY_LANGUAGE` 프로세스 호출을 추가할 수 있습니다.

워크플로를 다음과 같이 편집하십시오:

=== "수정 후"

    ```groovy title="main.nf" linenums="7" hl_lines="7-9"
        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .map { row ->
                [[id: row.id, character: row.character], row.recording]
            }

        // langid를 실행하여 각 인사말의 언어 식별
        IDENTIFY_LANGUAGE(ch_datasheet)
        IDENTIFY_LANGUAGE.out.view()
    ```

=== "수정 전"

    ```groovy title="main.nf" linenums="7" hl_lines="6"
        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .map { row ->
                [[id: row.id, character: row.character], row.recording]
            }
            .view()
    ```

채널 구성의 원래 `.view()` 작업을 제거했습니다.

이제 워크플로를 실행할 수 있습니다:

```bash
nextflow run main.nf
```

??? success "명령 출력"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [voluminous_mcnulty] DSL2 - revision: f9bcfebabb

    executor >  local (7)
    [4e/f722fe] IDENTIFY_LANGUAGE (7) [100%] 7 of 7 ✔
    [[id:sampleA, character:squirrel], /workspaces/training/side-quests/metadata/work/eb/f7148ebdd898fbe1136bec6a714acb/bonjour.txt, fr]
    [[id:sampleB, character:tux], /workspaces/training/side-quests/metadata/work/16/71d72410952c22cd0086d9bca03680/guten_tag.txt, de]
    [[id:sampleD, character:turkey], /workspaces/training/side-quests/metadata/work/c4/b7562adddc1cc0b7d414ec45d436eb/hello.txt, en]
    [[id:sampleC, character:sheep], /workspaces/training/side-quests/metadata/work/ea/04f5d979429e4455e14b9242fb3b45/hallo.txt, de]
    [[id:sampleF, character:moose], /workspaces/training/side-quests/metadata/work/5a/6c2b84bf8fadb98e28e216426be079/salut.txt, fr]
    [[id:sampleE, character:stegosaurus], /workspaces/training/side-quests/metadata/work/af/ee7c69bcab891c40d0529305f6b9e7/hola.txt, es]
    [[id:sampleG, character:turtle], /workspaces/training/side-quests/metadata/work/4e/f722fe47271ba7ebcd69afa42964ca/ciao.txt, it]
    ```

훌륭합니다! 이제 각 캐릭터가 어떤 언어를 사용하는지 예측했습니다.

그리고 앞서 언급했듯이 출력에 입력 파일과 메타 맵도 포함했습니다. 즉, 방금 생성한 새 정보와 연결된 상태로 유지됩니다.
이것은 다음 단계에서 유용하게 사용될 것입니다.

!!! note

    더 일반적으로, 메타 맵을 결과와 연결된 상태로 유지하는 이 패턴은 동일한 식별자를 공유하는 관련 결과를 연결하기 쉽게 만듭니다.

    이미 배웠듯이 채널의 항목 순서에 의존하여 결과를 일치시킬 수 없습니다.
    대신 키를 사용하여 데이터를 올바르게 연결해야 하며, 메타 맵은 이 목적에 이상적인 구조를 제공합니다.

    이 사용 사례는 [분할 및 그룹화](./splitting_and_grouping.md) 사이드 퀘스트에서 자세히 탐구합니다.

### 2.3. 프로세스 출력으로 메타데이터 확장

방금 생성한 결과 자체가 파일 내용에 대한 메타데이터 형태이므로 메타 맵에 추가하면 유용합니다.

그러나 기존 메타 맵을 제자리에서 수정하고 싶지는 않습니다.
기술적 관점에서 그렇게 _할 수는_ 있지만 안전하지 않습니다.

따라서 대신 `+` 연산자(Groovy 기능)를 사용하여 기존 메타 맵의 내용과 새 정보를 보유하는 새 `lang: lang_id` 키-값 쌍을 포함하는 새 메타 맵을 만들 것입니다.
그리고 [`map`](https://www.nextflow.io/docs/latest/operator.html#map) 작업과 결합하여 기존 맵을 새 맵으로 교체합니다.

워크플로에 다음과 같은 편집을 해야 합니다:

=== "수정 후"

    ```groovy title="main.nf" linenums="13" hl_lines="3-7"
        // langid를 실행하여 각 인사말의 언어 식별
        IDENTIFY_LANGUAGE(ch_datasheet)
        IDENTIFY_LANGUAGE.out
            .map { meta, file, lang_id ->
                [meta + [lang: lang_id], file]
            }
            .view()
    ```

=== "수정 전"

    ```groovy title="main.nf" linenums="13" hl_lines="3"
        // langid를 실행하여 각 인사말의 언어 식별
        IDENTIFY_LANGUAGE(ch_datasheet)
        IDENTIFY_LANGUAGE.out.view()
    ```

`+` 연산자에 아직 익숙하지 않거나 혼란스럽다면 아래의 자세한 설명을 살펴보는 데 몇 분을 할애하십시오.

??? info "`+` 연산자를 사용한 새 메타 맵 생성"

    **먼저, Groovy 연산자 `+`를 사용하여 두 맵의 내용을 병합할 수 있다는 것을 알아야 합니다.**

    다음과 같은 맵이 있다고 가정해 봅시다:

    ```groovy
    map1 = [id: 'sampleA', character: 'squirrel']
    map2 = [lang: 'fr']
    ```

    이렇게 병합할 수 있습니다:

    ```groovy
    new_map = map1 + map2
    ```

    `new_map`의 내용은 다음과 같습니다:

    ```groovy
    [id: 'sampleA', character: 'squirrel', lang: 'fr']
    ```

    훌륭합니다!

    **그러나 맵의 일부가 아닌 필드를 추가해야 한다면 어떻게 해야 할까요?**

    `map1`에서 다시 시작하지만 언어 예측이 자체 맵에 없다고 가정해 봅시다(map2가 없습니다).
    대신 `lang_id`라는 변수에 보관되어 있으며 해당 값(`'fr'`)을 `lang` 키로 저장하려고 합니다.

    실제로 다음과 같이 할 수 있습니다:

    ```groovy
    new_map = [map1 + [lang: lang_id]]
    ```

    여기서 `[lang: new_info]`는 즉석에서 새로운 이름 없는 맵을 생성하고 `map1 +`는 `map1`을 새로운 이름 없는 맵과 병합하여 이전과 동일한 `new_map` 내용을 생성합니다.

    멋지지 않나요?

    **이제 이것을 Nextflow `channel.map()` 작업의 맥락으로 변환해 봅시다.**

    코드는 다음과 같이 됩니다:

    ```groovy
    .map { map1, lang_id ->
        [map1 + [lang: lang_id]]
    }
    ```

    이것은 다음을 수행합니다:

    - `map1, lang_id ->` 튜플의 두 항목을 가져옵니다
    - `[map1 + [lang: lang_id]]` 위에서 설명한 대로 새 맵을 생성합니다

    출력은 위 예제의 `new_map`과 동일한 내용을 가진 단일 이름 없는 맵입니다.
    따라서 효과적으로 다음을 변환했습니다:

    ```groovy
    [id: 'sampleA', character: 'squirrel'], 'it'
    ```

    다음으로:

    ```groovy
    [id: 'sampleA', character: 'squirrel', lang: 'fr']
    ```

    `map1`을 `meta`로 변경하면 워크플로의 메타 맵에 언어 예측을 추가하기 위해 필요한 거의 모든 것임을 알 수 있을 것입니다.

    한 가지를 제외하고!

    워크플로의 경우 **튜플에 `file` 객체의 존재도 고려해야 합니다**. 튜플은 `meta, file, lang_id`로 구성됩니다.

    따라서 여기의 코드는 다음과 같이 됩니다:

    ```groovy
    .map { meta, file, lang_id ->
        [meta + [lang: lang_id], file]
    }
    ```

    `file`이 `map` 작업에서 이동하는 것처럼 보이는 이유를 따라가기 어렵다면, `[meta + [lang: lang_id], file]` 대신 그 줄이 `[new_map, file]`로 읽힌다고 상상해 보십시오.
    이렇게 하면 튜플의 두 번째 위치에 원래 위치에 `file`을 그대로 두고 있다는 것이 더 명확해질 것입니다. `new_info` 값을 가져와 첫 번째 위치의 맵에 접었습니다.

    **그리고 이것이 우리를 `tuple val(meta), path(file)` 채널 구조로 다시 데려옵니다!**

이 코드가 무엇을 하는지 확신이 서면 워크플로를 실행하여 작동하는지 확인하십시오:

```bash
nextflow run main.nf -resume
```

??? success "명령 출력"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [cheeky_fermat] DSL2 - revision: d096281ee4

    [4e/f722fe] IDENTIFY_LANGUAGE (7) [100%] 7 of 7, cached: 7 ✔
    [[id:sampleA, character:squirrel, lang:fr], /workspaces/training/side-quests/metadata/work/eb/f7148ebdd898fbe1136bec6a714acb/bonjour.txt]
    [[id:sampleB, character:tux, lang:de], /workspaces/training/side-quests/metadata/work/16/71d72410952c22cd0086d9bca03680/guten_tag.txt]
    [[id:sampleC, character:sheep, lang:de], /workspaces/training/side-quests/metadata/work/ea/04f5d979429e4455e14b9242fb3b45/hallo.txt]
    [[id:sampleD, character:turkey, lang:en], /workspaces/training/side-quests/metadata/work/c4/b7562adddc1cc0b7d414ec45d436eb/hello.txt]
    [[id:sampleF, character:moose, lang:fr], /workspaces/training/side-quests/metadata/work/5a/6c2b84bf8fadb98e28e216426be079/salut.txt]
    [[id:sampleE, character:stegosaurus, lang:es], /workspaces/training/side-quests/metadata/work/af/ee7c69bcab891c40d0529305f6b9e7/hola.txt]
    [[id:sampleG, character:turtle, lang:it], /workspaces/training/side-quests/metadata/work/4e/f722fe47271ba7ebcd69afa42964ca/ciao.txt]
    ```

네, 확인되었습니다!
프로세스의 출력을 `meta, file, lang_id`에서 깔끔하게 재구성했으므로 `lang_id`가 이제 메타 맵의 키 중 하나이고 채널의 튜플이 `meta, file` 모델에 다시 맞습니다.

<!-- TODO (future) Should we also show how to remove a key using subMap?! Or note where to find that. -->

### 2.4. 조건문을 사용하여 언어 그룹 할당

이제 언어 예측이 있으므로 정보를 사용하여 새 그룹을 할당해 봅시다.

예제 데이터에서 캐릭터가 사용하는 언어는 게르만 언어(영어, 독일어)와 로망스 언어(프랑스어, 스페인어, 이탈리아어)로 그룹화할 수 있습니다.
파이프라인 후반에 해당 분류를 쉽게 사용할 수 있도록 메타 맵에 해당 정보를 추가해 봅시다.

그리고 좋은 소식은 이것이 `map` 연산자를 사용하기에 완벽한 또 다른 경우라는 것입니다!

구체적으로 `lang_group`이라는 변수를 정의하고 간단한 조건부 로직을 사용하여 각 데이터에 대해 `lang_group`에 할당할 값을 결정할 것입니다.

일반적인 구문은 다음과 같습니다:

```groovy
.map { meta, file ->

    // lang_group을 정의하는 조건부 로직이 여기에 들어갑니다

    [meta + [lang_group: lang_group], file]
}
```

이것이 이전 단계에서 사용한 즉석 맵 병합 작업과 매우 유사한 것을 볼 수 있습니다.
조건문만 작성하면 됩니다.

적용하려는 조건부 로직은 다음과 같습니다:

- 기본값이 `'unknown'`인 `lang_group`이라는 변수를 정의합니다.
- `lang`이 독일어(`'de'`) 또는 영어(`'en'`)이면 `lang_group`을 `germanic`으로 변경합니다.
- 그렇지 않고 `lang`이 프랑스어(`'fr'`), 스페인어(`'es'`) 및 이탈리아어(`'it'`)를 포함하는 목록에 포함되어 있으면 `lang_group`을 `romance`로 변경합니다.

Nextflow에서 조건문을 작성하는 방법을 이미 알고 있다면 직접 작성해 보십시오.

!!! tip

    맵 작업 내에서 `meta.lang`으로 `lang` 값에 액세스할 수 있습니다.

워크플로를 다음과 같이 변경해야 합니다:

=== "수정 후"

    ```groovy title="main.nf" linenums="13" hl_lines="7-19"
        // langid를 실행하여 각 인사말의 언어 식별
        IDENTIFY_LANGUAGE(ch_datasheet)
        IDENTIFY_LANGUAGE.out
            .map { meta, file, lang_id ->
                [meta + [lang: lang_id], file]
            }
            .map { meta, file ->

                def lang_group = "unknown"
                if (meta.lang.equals("de") || meta.lang.equals('en')) {
                    lang_group = "germanic"
                }
                else if (meta.lang in ["fr", "es", "it"]) {
                    lang_group = "romance"
                }

                [meta + [lang_group: lang_group], file]
            }
            .view()
    ```

=== "수정 전"

    ```groovy title="main.nf" linenums="13" hl_lines="7"
        // langid를 실행하여 각 인사말의 언어 식별
        IDENTIFY_LANGUAGE(ch_datasheet)
        IDENTIFY_LANGUAGE.out
            .map { meta, file, lang_id ->
                [meta + [lang: lang_id], file]
            }
            .view()
    ```

핵심 사항은 다음과 같습니다:

- `def lang_group = "unknown"`을 사용하여 기본값이 `unknown`으로 설정된 `lang_group` 변수를 생성합니다.
- 조건부 로직에 `if {} else if {}` 구조를 사용하며, 두 게르만 언어에 대한 대체 `.equals()` 테스트와 세 로망스 언어에 대한 목록 존재 테스트를 사용합니다.
- 이전과 같이 `meta + [lang_group:lang_group]` 병합 작업을 사용하여 업데이트된 메타 맵을 생성합니다.

<!-- TODO (future) Add note/links to relevant docs in additional resources section -->

모두 이해가 되면 워크플로를 다시 실행하여 결과를 확인하십시오:

```bash
nextflow run main.nf -resume
```

??? success "명령 출력"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [wise_almeida] DSL2 - revision: 46778c3cd0

    [da/652cc6] IDENTIFY_LANGUAGE (7) [100%] 7 of 7, cached: 7 ✔
    [[id:sampleA, character:squirrel, lang:fr, lang_group:romance], /workspaces/training/side-quests/metadata/data/bonjour.txt]
    [[id:sampleB, character:tux, lang:de, lang_group:germanic], /workspaces/training/side-quests/metadata/data/guten_tag.txt]
    [[id:sampleC, character:sheep, lang:de, lang_group:germanic], /workspaces/training/side-quests/metadata/data/hallo.txt]
    [[id:sampleD, character:turkey, lang:en, lang_group:germanic], /workspaces/training/side-quests/metadata/data/hello.txt]
    [[id:sampleE, character:stegosaurus, lang:es, lang_group:romance], /workspaces/training/side-quests/metadata/data/hola.txt]
    [[id:sampleF, character:moose, lang:fr, lang_group:romance], /workspaces/training/side-quests/metadata/data/salut.txt]
    [[id:sampleG, character:turtle, lang:it, lang_group:romance], /workspaces/training/side-quests/metadata/data/ciao.txt]
    ```

보시다시피 채널 요소는 `[meta, file]` 구조를 유지하지만 메타 맵에는 이제 이 새로운 분류가 포함됩니다.

### 요점

이 섹션에서 다음을 학습했습니다 :

- **출력 채널에 입력 메타데이터 적용**: 이러한 방식으로 메타데이터를 복사하면 나중에 메타데이터 내용을 기반으로 결과를 연결할 수 있습니다.
- **사용자 정의 키 생성**: 메타 맵에 두 개의 새 키를 생성하고 `meta + [new_key:value]`로 기존 메타 맵에 병합했습니다. 하나는 프로세스의 계산된 값을 기반으로 하고 다른 하나는 `map` 연산자에서 설정한 조건을 기반으로 합니다.

이를 통해 파이프라인을 진행하면서 파일과 새 메타데이터 및 기존 메타데이터를 연결할 수 있습니다.
프로세스의 일부로 메타데이터를 사용하지 않더라도 이와 같이 메타 맵을 데이터와 연결된 상태로 유지하면 모든 관련 정보를 함께 유지하기 쉽습니다.

---

## 3. 프로세스에서 메타 맵 정보 사용

이제 메타 맵을 만들고 업데이트하는 방법을 알았으므로 정말 재미있는 부분, 즉 프로세스에서 메타데이터를 실제로 사용하는 것으로 넘어갈 수 있습니다.

더 구체적으로, 각 동물을 ASCII 아트로 그리고 말풍선에 녹음된 텍스트를 말하게 하는 두 번째 단계를 워크플로에 추가할 것입니다.
[`cowpy`](https://github.com/jeffbuttars/cowpy)라는 도구를 사용하여 이 작업을 수행할 것입니다.

??? info "`cowpy`는 무엇을 하나요?"

    `cowpy`는 임의의 텍스트 입력을 재미있는 방식으로 표시하는 ASCII 아트를 생성하는 명령줄 도구입니다.
    Tony Monroe의 클래식 [cowsay](https://en.wikipedia.org/wiki/Cowsay) 도구의 Python 구현입니다.

    ```console
    cowpy "Hello Nextflow"
    ```

    ```console
    ______________________________________________________
    < Hello Nextflow >
    ------------------------------------------------------
        \   ^__^
          \  (oo)\_______
            (__)\       )\/\
              ||----w |
              ||     ||
    ```

    선택적으로 기본 소 대신 사용할 캐릭터(또는 'cowacter')를 선택할 수 있습니다.

    ```console
    cowpy "Hello Nextflow" -c tux
    ```

    ```console
    __________________
    < Hello Nextflow >
    ------------------
      \
        \
            .--.
          |o_o |
          |:_/ |
          //   \ \
        (|     | )
        /'\_   _/`\
        \___)=(___/
    ```

Hello Nextflow 과정을 진행했다면 이미 이 도구가 실행되는 것을 보았을 것입니다.
그렇지 않더라도 걱정하지 마십시오. 진행하면서 알아야 할 모든 것을 다룰 것입니다.

### 3.1. 프로세스 가져오기 및 코드 검토

`cowpy` 도구를 적용하는 `COWPY`라는 사전 작성된 프로세스 모듈을 제공하므로 workflow 블록 앞에 include 문을 추가하기만 하면 됩니다.

워크플로를 다음과 같이 편집하십시오:

=== "수정 후"

    ```groovy title="main.nf" linenums="1" hl_lines="4"
    #!/usr/bin/env nextflow

    include { IDENTIFY_LANGUAGE } from './modules/langid.nf'
    include { COWPY } from './modules/cowpy.nf'

    workflow {
    ```

=== "수정 전"

    ```groovy title="main.nf" linenums="1"
    #!/usr/bin/env nextflow

    include { IDENTIFY_LANGUAGE } from './modules/langid.nf'

    workflow {
    ```

모듈 파일을 열어 코드를 확인할 수 있습니다:

```groovy title="modules/cowpy.nf" linenums="1"
#!/usr/bin/env nextflow

// cowpy로 ASCII 아트 생성
process COWPY {

    publishDir "results/", mode: 'copy'

    container 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'

    input:
    path input_file
    val character

    output:
    path "cowpy-${input_file}"

    script:
    """
    cat ${input_file} | cowpy -c ${character} > cowpy-${input_file}
    """
}
```

보시다시피 이 프로세스는 현재 입력 파일(표시할 텍스트 포함)과 ASCII 아트로 그릴 캐릭터를 지정하는 값을 받도록 설계되어 있으며, 일반적으로 명령줄 매개변수에 의해 워크플로 수준에서 제공됩니다.

### 3.2. 메타 맵 필드를 입력으로 전달

Hello Nextflow 과정에서 `cowpy` 도구를 사용할 때 명령줄 매개변수를 사용하여 최종 이미지를 그리는 데 사용할 캐릭터를 결정했습니다.
파이프라인 실행당 하나의 이미지만 생성했기 때문에 그것이 타당했습니다.

그러나 이 튜토리얼에서는 처리하는 각 대상에 대해 적절한 이미지를 생성하려고 하므로 명령줄 매개변수를 사용하면 너무 제한적입니다.

좋은 소식: 데이터시트에 `character` 열이 있으므로 메타 맵에도 있습니다.
프로세스가 각 항목에 사용해야 하는 캐릭터를 설정하는 데 사용해 봅시다.

이를 위해 세 가지를 수행해야 합니다:

1. 이전 프로세스에서 나오는 출력 채널에 이름을 지정하여 더 편리하게 작업할 수 있도록 합니다.
2. 관심 있는 정보에 액세스하는 방법을 결정합니다.
3. 두 번째 프로세스 호출을 추가하고 정보를 적절하게 제공합니다.

시작하겠습니다.

#### 3.2.1. 이전 출력 채널에 이름 지정

첫 번째 프로세스 `IDENTIFY_LANGUAGE.out`의 출력 채널에 직접 이전 조작을 적용했습니다.
해당 채널의 내용을 다음 프로세스에 제공하기 위해(명확하고 읽기 쉬운 방식으로) 자체 이름 `ch_languages`를 지정하려고 합니다.

[`set`](https://www.nextflow.io/docs/latest/reference/operator.html#set) 연산자를 사용하여 그렇게 할 수 있습니다.

메인 워크플로에서 `.view()` 연산자를 `.set { ch_languages }`로 교체하고 이름으로 채널을 참조할 수 있는지 테스트하는 라인을 추가하십시오.

=== "수정 후"

    ```groovy title="main.nf" linenums="14" hl_lines="19 21 22"
        // langid를 실행하여 각 인사말의 언어 식별
        IDENTIFY_LANGUAGE(ch_datasheet)
        IDENTIFY_LANGUAGE.out
            .map { meta, file, lang_id ->
                [meta + [lang: lang_id], file]
            }
            .map { meta, file ->

                def lang_group = "unknown"
                if (meta.lang.equals("de") || meta.lang.equals('en')) {
                    lang_group = "germanic"
                }
                else if (meta.lang in ["fr", "es", "it"]) {
                    lang_group = "romance"
                }

                [meta + [lang_group: lang_group], file]
            }
            .set { ch_languages }

        // 임시: ch_languages 살펴보기
        ch_languages.view()
    ```

=== "수정 전"

    ```groovy title="main.nf" linenums="14" hl_lines="19"
        // langid를 실행하여 각 인사말의 언어 식별
        IDENTIFY_LANGUAGE(ch_datasheet)
        IDENTIFY_LANGUAGE.out
            .map { meta, file, lang_id ->
                [meta + [lang: lang_id], file]
            }
            .map { meta, file ->

                def lang_group = "unknown"
                if (meta.lang.equals("de") || meta.lang.equals('en')) {
                    lang_group = "germanic"
                }
                else if (meta.lang in ["fr", "es", "it"]) {
                    lang_group = "romance"
                }

                [meta + [lang_group: lang_group], file]
            }
            .view()
    ```

이것을 실행해 봅시다:

```bash
nextflow run main.nf
```

??? success "명령 출력"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `./main.nf` [friendly_austin] DSL2 - revision: 3dbe460fd6

    [36/cca6a7] IDENTIFY_LANGUAGE (7) | 7 of 7 ✔
    [[id:sampleB, character:tux, lang:de, lang_group:germanic], /workspaces/training/side-quests/metadata/work/e2/6db2402d83cf72081bcd2d11784714/guten_tag.txt]
    [[id:sampleA, character:squirrel, lang:fr, lang_group:romance], /workspaces/training/side-quests/metadata/work/6c/114c818317d169457d6e7336d5d55b/bonjour.txt]
    [[id:sampleC, character:sheep, lang:de, lang_group:germanic], /workspaces/training/side-quests/metadata/work/55/68c69c5efb527f3604ddb3daab8057/hallo.txt]
    [[id:sampleD, character:turkey, lang:en, lang_group:germanic], /workspaces/training/side-quests/metadata/work/2a/4752055ccb5d1370b0ef9da41d3993/hello.txt]
    [[id:sampleE, character:stegosaurus, lang:es, lang_group:romance], /workspaces/
