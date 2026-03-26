# 메타데이터와 메타 맵

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } AI 지원 번역 - [자세히 알아보기 및 개선 제안](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

과학적 분석에서 원시 데이터 파일만 다루는 경우는 거의 없습니다.
각 파일에는 고유한 추가 정보가 있습니다. 파일이 무엇인지, 어디서 왔는지, 무엇이 특별한지에 대한 정보입니다.
이러한 추가 정보를 메타데이터라고 합니다.

메타데이터는 다른 데이터를 설명하는 데이터입니다.
메타데이터는 파일과 실험 조건에 대한 중요한 세부 정보를 추적하며, 각 데이터셋의 고유한 특성에 맞게 분석을 조정하는 데 도움을 줍니다.

도서관 카탈로그를 생각해 보세요. 책에는 실제 내용(원시 데이터)이 담겨 있지만, 카탈로그 카드는 각 책에 대한 필수 정보(출판 시기, 저자, 위치 등)를 제공합니다(메타데이터).
Nextflow 파이프라인에서 메타데이터는 다음과 같은 용도로 사용됩니다.

- 워크플로우 전반에 걸쳐 파일별 정보 추적
- 파일 특성에 따른 프로세스 설정
- 공동 분석을 위한 관련 파일 그룹화

### 학습 목표

이 사이드 퀘스트에서는 워크플로우에서 메타데이터를 처리하는 방법을 살펴봅니다.
기본 파일 정보가 담긴 간단한 데이터시트(생물정보학에서는 흔히 samplesheet라고 함)를 시작으로 다음 내용을 학습합니다.

- CSV 파일에서 파일 메타데이터 읽기 및 파싱
- 메타데이터 맵 생성 및 조작
- 워크플로우 실행 중 새로운 메타데이터 필드 추가
- 메타데이터를 활용한 프로세스 동작 맞춤화

이러한 기술을 통해 복잡한 파일 관계와 처리 요구사항을 처리할 수 있는 더욱 견고하고 유연한 파이프라인을 구축할 수 있습니다.

### 사전 요구사항

이 사이드 퀘스트를 시작하기 전에 다음을 완료해야 합니다.

- [Hello Nextflow](../hello_nextflow/README.md) 튜토리얼 또는 동급의 입문 과정을 완료해야 합니다.
- 기본적인 Nextflow 개념과 메커니즘(프로세스, 채널, 연산자)에 익숙해야 합니다.

---

## 0. 시작하기

#### 교육 코드스페이스 열기

아직 열지 않았다면 [환경 설정](../envsetup/index.md)에 설명된 대로 교육 환경을 열어 주세요.

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

#### 프로젝트 디렉토리로 이동

이 튜토리얼의 파일이 있는 디렉토리로 이동합니다.

```bash
cd side-quests/metadata
```

VSCode에서 이 디렉토리에 포커스를 설정할 수 있습니다.

```bash
code .
```

#### 자료 검토

메인 워크플로우 파일과 데이터시트 및 몇 가지 데이터 파일이 포함된 `data` 디렉토리가 있습니다.

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

`main.nf` 파일의 워크플로우는 단계적으로 완전한 기능을 갖춘 워크플로우로 확장할 스텁입니다.

데이터시트에는 데이터 파일 경로와 관련 메타데이터가 3개의 열로 정리되어 있습니다.

- `id`: 파일에 부여된 ID
- `character`: 나중에 다양한 캐릭터를 그리는 데 사용할 캐릭터 이름
- `data`: 다양한 언어로 된 인사말이 담긴 `.txt` 파일 경로

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

여러분의 과제는 다음을 수행하는 Nextflow 워크플로우를 작성하는 것입니다.

1. 각 파일의 언어를 자동으로 **식별**합니다.
2. 언어 계열(게르만어 vs 로망스어)에 따라 파일을 **그룹화**합니다.
3. 언어와 메타데이터를 기반으로 각 파일의 처리를 **맞춤화**합니다.
4. 언어 그룹별로 출력을 **정리**합니다.

이는 파일별 메타데이터가 처리 결정을 이끄는 전형적인 워크플로우 패턴으로, 메타데이터 맵이 우아하게 해결하는 바로 그런 문제입니다.

#### 준비 체크리스트

시작할 준비가 되었나요?

- [ ] 이 과정의 목표와 사전 요구사항을 이해했습니다.
- [ ] 코드스페이스가 실행 중입니다.
- [ ] 작업 디렉토리를 적절히 설정했습니다.
- [ ] 과제를 이해했습니다.

모든 항목을 확인했다면 시작할 준비가 된 것입니다.

---

## 1. 데이터시트에서 메타데이터 불러오기

`main.nf` 워크플로우 파일을 열어 시작점으로 제공된 워크플로우 스텁을 확인합니다.

```groovy title="main.nf" linenums="1"
#!/usr/bin/env nextflow

workflow  {

    ch_datasheet = channel.fromPath("./data/datasheet.csv")

}
```

예제 데이터시트를 파일로 불러오는 기본 채널 팩토리를 설정했지만, 아직 파일의 내용을 읽어오지는 않습니다.
이 부분부터 시작합니다.

### 1.1. `splitCsv`로 내용 읽기

최소한의 노력으로 파일 내용을 적절히 파싱할 연산자를 선택해야 합니다.
데이터시트가 CSV 형식이므로, 파일의 각 행을 채널의 요소로 불러오는 [`splitCsv`](https://www.nextflow.io/docs/latest/reference/operator.html#splitcsv) 연산자를 사용합니다.

다음과 같이 채널 구성 코드에 `splitCsv()` 작업을 추가하고, 파일 내용이 채널에 올바르게 불러와지는지 확인하기 위한 `view()` 작업도 추가합니다.

=== "후"

    ```groovy title="main.nf" linenums="3" hl_lines="4-5"
    workflow  {

        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .view()

    }
    ```

=== "전"

    ```groovy title="main.nf" linenums="3"
    workflow  {

        ch_datasheet = channel.fromPath("./data/datasheet.csv")

    }
    ```

`header: true` 옵션을 사용하여 Nextflow가 CSV 파일의 첫 번째 행을 헤더 행으로 읽도록 지정합니다.

결과를 확인해 봅니다.
워크플로우를 실행합니다.

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

연산자가 CSV 파일의 각 행에 대해 키-값 쌍의 맵을 구성했으며, 열 헤더가 해당 값의 키로 사용된 것을 확인할 수 있습니다.

각 맵 항목은 데이터시트의 열에 해당합니다.

- `id`
- `character`
- `recording`

훌륭합니다! 각 파일의 특정 필드에 쉽게 접근할 수 있습니다.
예를 들어, `id`로 파일 ID에 접근하거나 `recording`으로 txt 파일 경로에 접근할 수 있습니다.

??? info "(선택 사항) 맵에 대한 추가 정보"

    Nextflow가 기반으로 하는 프로그래밍 언어인 Groovy에서 맵은 Python의 딕셔너리, JavaScript의 객체, Ruby의 해시와 유사한 키-값 데이터 구조입니다.

    다음은 맵을 정의하고 실제로 내용에 접근하는 방법을 보여주는 실행 가능한 스크립트입니다.

    ```groovy title="examples/map_demo.nf"
    #!/usr/bin/env nextflow

    // 간단한 맵 생성
    def my_map = [id:'sampleA', character:'squirrel']

    // 전체 맵 출력
    println "map: ${my_map}"

    // 점 표기법으로 개별 값 접근
    println "id: ${my_map.id}"
    println "character: ${my_map.character}"
    ```

    적절한 `workflow` 블록이 없더라도 Nextflow는 이를 워크플로우처럼 실행할 수 있습니다.

    ```bash
    nextflow run examples/map_demo.nf
    ```

    출력에서 다음과 같은 결과를 확인할 수 있습니다.

    ```console title="Output"
     N E X T F L O W   ~  version 25.10.2

    Launching `map_demo.nf` [cheesy_plateau] DSL2 - revision: fae5b8496e

    map: [id:sampleA, character:squirrel]
    id: sampleA
    character: squirrel
    ```

### 1.2. `map`으로 특정 필드 선택

데이터시트에서 `character` 열에 접근하여 출력하려고 합니다.
Nextflow `map` 연산자를 사용하여 채널의 각 항목을 순회하고 맵 객체에서 `character` 항목을 선택할 수 있습니다.

워크플로우를 다음과 같이 수정합니다.

=== "후"

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

=== "전"

    ```groovy title="main.nf" linenums="3"
    workflow  {

        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .view()

    }
    ```

워크플로우를 다시 실행합니다.

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

성공입니다! 데이터시트에서 파생된 맵 구조를 활용하여 각 행의 개별 열 값에 접근했습니다.

데이터시트를 성공적으로 읽고 각 행의 데이터에 접근할 수 있게 되었으므로, 이제 파이프라인 로직 구현을 시작할 수 있습니다.

### 1.3. 메타데이터를 '메타 맵'으로 구성

현재 워크플로우 상태에서는 입력 파일(`recording` 키 아래)과 관련 메타데이터(`id`, `character`)가 모두 동등한 위치에 있어, 마치 하나의 큰 가방에 담긴 것과 같습니다.
실질적인 결과로, 이 채널을 사용하는 모든 프로세스는 이 구조를 염두에 두고 설정해야 합니다.

```groovy
    input:
    tuple val(id), val(character), file(recording)
```

데이터시트의 열 수가 변경되지 않는 한 이 방식은 괜찮습니다.
그러나 데이터시트에 열을 하나만 추가해도 채널의 형태가 프로세스가 기대하는 것과 일치하지 않아 워크플로우에서 오류가 발생합니다.
또한 약간 다른 입력 데이터를 가진 다른 사람들과 프로세스를 공유하기 어렵게 만들고, 스크립트 블록에서 필요하지 않은 변수를 프로세스에 하드코딩해야 할 수도 있습니다.

이 문제를 피하려면 데이터시트에 열이 몇 개 있든 채널 구조를 일관되게 유지하는 방법이 필요합니다.

튜플 내의 항목에 모든 메타데이터를 수집하면 됩니다. 이를 메타데이터 맵, 또는 더 간단히 '메타 맵'이라고 합니다.

`map` 작업을 다음과 같이 수정합니다.

=== "후"

    ```groovy title="main.nf" linenums="5" hl_lines="4"
        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .map { row ->
                    [ [id: row.id, character: row.character], row.recording ]
            }
            .view()
    ```

=== "전"

    ```groovy title="main.nf" linenums="5" hl_lines="4"
        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .map{ row ->
                row.character
            }
            .view()
    ```

채널 요소를 메타 맵과 해당 파일 객체, 두 요소로 구성된 튜플로 재구성했습니다.

워크플로우를 실행합니다.

```bash
nextflow run main.nf
```

??? success "명령 출력"

    ```console title="View meta map"
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

이제 채널의 각 요소는 메타데이터 맵을 첫 번째로, 해당 파일 객체를 두 번째로 포함합니다.

```console title="Example output structure"
[
  [id:sampleA, character:squirrel],
  /workspaces/training/side-quests/metadata/data/bonjour.txt
]
```

결과적으로 데이터시트에 더 많은 열을 추가하면 `meta` 맵에서 더 많은 메타데이터를 사용할 수 있지만, 채널 형태는 변경되지 않습니다.
이를 통해 메타데이터 항목을 입력 사양에 하드코딩하지 않고도 채널을 사용하는 프로세스를 작성할 수 있습니다.

```groovy title="Syntax example"
    input:
    tuple val(meta), file(recording)
```

이는 Nextflow 워크플로우에서 메타데이터를 구성하는 데 널리 사용되는 관례입니다.

### 핵심 정리

이 섹션에서 다음 내용을 학습했습니다.

- **메타데이터가 중요한 이유:** 데이터와 함께 메타데이터를 유지하면 워크플로우 전반에 걸쳐 중요한 파일 정보가 보존됩니다.
- **데이터시트 읽는 방법:** `splitCsv`를 사용하여 헤더 정보가 있는 CSV 파일을 읽고 행을 구조화된 데이터로 변환합니다.
- **메타 맵 생성 방법:** `[ [id:value, ...], file ]` 튜플 구조를 사용하여 파일 데이터에서 메타데이터를 분리합니다.

---

## 2. 메타데이터 조작

메타데이터를 불러왔으니 이제 활용해 봅니다!

[`langid`](https://github.com/saffsd/langid.py)라는 도구를 사용하여 각 캐릭터의 녹음 파일에 포함된 언어를 식별합니다.
이 도구는 언어 세트로 사전 학습되어 있으며, 텍스트 조각이 주어지면 언어 예측과 관련 확률 점수를 `stdout`으로 출력합니다.

### 2.1. 프로세스 가져오기 및 코드 검토

`langid` 도구를 래핑하는 `IDENTIFY_LANGUAGE`라는 사전 작성된 프로세스 모듈을 제공합니다. 워크플로우 블록 앞에 include 문만 추가하면 됩니다.

워크플로우를 다음과 같이 수정합니다.

=== "후"

    ```groovy title="main.nf" linenums="1" hl_lines="3"
    #!/usr/bin/env nextflow

    include { IDENTIFY_LANGUAGE } from './modules/langid.nf'

    workflow {
    ```

=== "전"

    ```groovy title="main.nf" linenums="1"
    #!/usr/bin/env nextflow

    workflow {
    ```

모듈 파일을 열어 코드를 확인할 수 있습니다.

```groovy title="modules/langid.nf" linenums="1" hl_lines="9 12"
#!/usr/bin/env nextflow

// langid를 사용하여 각 입력 파일의 언어를 예측합니다
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

입력 정의가 입력 채널에 적용한 것과 동일한 `tuple val(meta), path(file)` 구조를 사용하는 것을 확인할 수 있습니다.

출력 정의는 입력과 유사한 구조의 튜플로 구성되어 있으며, 세 번째 요소로 `stdout`도 포함합니다.
이 `tuple val(meta), path(file), <output>` 패턴은 파이프라인을 통해 흐르면서 입력 데이터와 출력 모두에 메타데이터가 연결된 상태를 유지합니다.

도구가 파일을 작성하는 대신 콘솔에 직접 출력을 인쇄하기 때문에 Nextflow의 [`stdout`](https://www.nextflow.io/docs/latest/process.html#outputs) 출력 한정자를 사용합니다. 명령줄에서 `sed`를 사용하여 확률 점수를 제거하고, 개행 문자를 제거하여 문자열을 정리하고, 언어 예측만 반환합니다.

### 2.2. `IDENTIFY_LANGUAGE` 실행 추가

프로세스를 워크플로우에서 사용할 수 있게 되었으므로, 데이터 채널에서 `IDENTIFY_LANGUAGE` 프로세스를 실행하는 호출을 추가할 수 있습니다.

워크플로우를 다음과 같이 수정합니다.

=== "후"

    ```groovy title="main.nf" linenums="7" hl_lines="7-9"
        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .map { row ->
                [[id: row.id, character: row.character], row.recording]
            }

        // langid를 실행하여 각 인사말의 언어를 식별합니다
        IDENTIFY_LANGUAGE(ch_datasheet)
        IDENTIFY_LANGUAGE.out.view()
    ```

=== "전"

    ```groovy title="main.nf" linenums="7" hl_lines="6"
        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .map { row ->
                [[id: row.id, character: row.character], row.recording]
            }
            .view()
    ```

채널 구성에서 원래의 `.view()` 작업을 제거했습니다.

이제 워크플로우를 실행합니다.

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

훌륭합니다! 이제 각 캐릭터가 사용하는 언어에 대한 예측이 생성되었습니다.

앞서 언급했듯이, 출력에 입력 파일과 메타 맵도 포함되어 있어 방금 생성한 새로운 정보와 연결된 상태를 유지합니다.
이는 다음 단계에서 유용하게 활용됩니다.

!!! note "참고"

    더 일반적으로, 메타 맵을 결과와 연결된 상태로 유지하는 이 패턴은 동일한 식별자를 공유하는 관련 결과를 연결하기 쉽게 만듭니다.

    이미 학습했듯이, 채널의 항목 순서에 의존하여 결과를 매칭할 수 없습니다.
    대신 키를 사용하여 데이터를 올바르게 연결해야 하며, 메타 맵은 이 목적에 이상적인 구조를 제공합니다.

    이 사용 사례는 [Splitting & Grouping](../splitting_and_grouping/) 사이드 퀘스트에서 자세히 살펴봅니다.

### 2.3. 프로세스 출력으로 메타데이터 보강

방금 생성한 결과가 파일 내용에 대한 메타데이터의 한 형태이므로, 이를 메타 맵에 추가하면 유용합니다.

그러나 기존 메타 맵을 직접 수정하고 싶지는 않습니다.
기술적으로는 가능하지만, 안전하지 않습니다.

대신 `+` 연산자(Groovy 기능)를 사용하여 기존 메타 맵의 내용과 새로운 정보를 담는 `lang: lang_id` 키-값 쌍을 포함하는 새로운 메타 맵을 생성합니다.
그리고 [`map`](https://www.nextflow.io/docs/latest/operator.html#map) 작업과 결합하여 이전 맵을 새 맵으로 교체합니다.

워크플로우에 다음과 같은 수정을 합니다.

=== "후"

    ```groovy title="main.nf" linenums="13" hl_lines="3-7"
        // langid를 실행하여 각 인사말의 언어를 식별합니다
        IDENTIFY_LANGUAGE(ch_datasheet)
        IDENTIFY_LANGUAGE.out
            .map { meta, file, lang_id ->
                [meta + [lang: lang_id], file]
            }
            .view()
    ```

=== "전"

    ```groovy title="main.nf" linenums="13" hl_lines="3"
        // langid를 실행하여 각 인사말의 언어를 식별합니다
        IDENTIFY_LANGUAGE(ch_datasheet)
        IDENTIFY_LANGUAGE.out.view()
    ```

`+` 연산자에 익숙하지 않거나 혼란스럽다면, 아래의 자세한 설명을 잠시 살펴보세요.

??? info "`+` 연산자를 사용한 새 메타 맵 생성"

    **먼저, Groovy 연산자 `+`를 사용하여 두 맵의 내용을 병합할 수 있다는 것을 알아야 합니다.**

    다음과 같은 맵이 있다고 가정합니다.

    ```groovy
    map1 = [id: 'sampleA', character: 'squirrel']
    map2 = [lang: 'fr']
    ```

    다음과 같이 병합할 수 있습니다.

    ```groovy
    new_map = map1 + map2
    ```

    `new_map`의 내용은 다음과 같습니다.

    ```groovy
    [id: 'sampleA', character: 'squirrel', lang: 'fr']
    ```

    훌륭합니다!

    **그런데 맵에 아직 없는 필드를 추가해야 한다면 어떻게 할까요?**

    `map1`에서 다시 시작하지만, 언어 예측이 자체 맵에 없는 경우(`map2`가 없는 경우)를 가정합니다.
    대신 `lang_id`라는 변수에 저장되어 있고, 그 값(`'fr'`)을 `lang` 키로 저장하려고 합니다.

    다음과 같이 할 수 있습니다.

    ```groovy
    new_map = [map1 + [lang: lang_id]]
    ```

    여기서 `[lang: new_info]`는 즉석에서 새로운 이름 없는 맵을 생성하고, `map1 + `는 `map1`을 새로운 이름 없는 맵과 병합하여 이전과 동일한 `new_map` 내용을 생성합니다.

    멋지지 않나요?

    **이제 Nextflow `channel.map()` 작업의 맥락으로 적용해 봅니다.**

    코드는 다음과 같습니다.

    ```groovy
    .map { map1, lang_id ->
        [map1 + [lang: lang_id]]
    }
    ```

    이 코드는 다음을 수행합니다.

    - `map1, lang_id ->`는 튜플의 두 항목을 가져옵니다.
    - `[map1 + [lang: lang_id]]`는 위에서 설명한 대로 새 맵을 생성합니다.

    출력은 위 예시의 `new_map`과 동일한 내용을 가진 단일 이름 없는 맵입니다.
    즉, 다음을 변환한 것입니다.

    ```groovy
    [id: 'sampleA', character: 'squirrel'], 'fr'
    ```

    다음으로:

    ```groovy
    [id: 'sampleA', character: 'squirrel', lang: 'fr']
    ```

    `map1`을 `meta`로 변경하면, 워크플로우의 메타 맵에 언어 예측을 추가하는 데 필요한 것이 기본적으로 전부라는 것을 알 수 있습니다.

    단 한 가지만 제외하고요!

    워크플로우의 경우, **`meta, file, lang_id`로 구성된 튜플에서 `file` 객체의 존재도 고려해야 합니다.**

    따라서 코드는 다음과 같이 됩니다.

    ```groovy
    .map { meta, file, lang_id ->
        [meta + [lang: lang_id], file]
    }
    ```

    `map` 작업에서 `file`이 이동하는 것처럼 보이는 이유가 이해하기 어렵다면, `[meta + [lang: lang_id], file]` 대신 `[new_map, file]`로 읽어보세요.
    이렇게 하면 `file`을 튜플의 두 번째 위치에 그대로 두고, `new_info` 값을 첫 번째 위치의 맵에 통합했다는 것이 더 명확해집니다.

    **이것이 바로 `tuple val(meta), path(file)` 채널 구조로 돌아오는 것입니다!**

이 코드가 무엇을 하는지 이해했다면, 워크플로우를 실행하여 작동하는지 확인합니다.

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

확인되었습니다!
프로세스의 출력을 `meta, file, lang_id`에서 `lang_id`가 메타 맵의 키 중 하나가 되도록 깔끔하게 재구성했으며, 채널의 튜플이 다시 `meta, file` 모델에 맞게 되었습니다.

<!-- TODO (future) Should we also show how to remove a key using subMap?! Or note where to find that. -->

### 2.4. 조건문을 사용하여 언어 그룹 할당

언어 예측이 완료되었으니, 이 정보를 사용하여 새로운 그룹을 할당합니다.

예제 데이터에서 캐릭터들이 사용하는 언어는 게르만어(영어, 독일어)와 로망스어(프랑스어, 스페인어, 이탈리아어)로 그룹화할 수 있습니다.
나중에 파이프라인에서 이 분류를 쉽게 사용할 수 있으면 유용하므로, 메타 맵에 해당 정보를 추가합니다.

이 역시 `map` 연산자를 완벽하게 활용할 수 있는 경우입니다!

구체적으로, `lang_group`이라는 변수를 정의하고, 간단한 조건 로직을 사용하여 각 데이터에 `lang_group` 값을 결정합니다.

일반적인 문법은 다음과 같습니다.

```groovy
.map { meta, file ->

    // lang_group을 정의하는 조건 로직이 여기에 들어갑니다

    [meta + [lang_group: lang_group], file]
}
```

이전 단계에서 사용한 즉석 맵 병합 작업과 매우 유사합니다.
조건문만 작성하면 됩니다.

적용할 조건 로직은 다음과 같습니다.

- 기본값 `'unknown'`으로 `lang_group` 변수를 정의합니다.
- `lang`이 독일어(`'de'`) 또는 영어(`'en'`)이면 `lang_group`을 `germanic`으로 변경합니다.
- 그렇지 않고 `lang`이 프랑스어(`'fr'`), 스페인어(`'es'`), 이탈리아어(`'it'`)를 포함하는 목록에 있으면 `lang_group`을 `romance`로 변경합니다.

Nextflow에서 조건문 작성 방법을 이미 알고 있다면 직접 작성해 보세요.

!!! tip "팁"

    map 작업 내에서 `meta.lang`으로 `lang` 값에 접근할 수 있습니다.

워크플로우에 다음과 같은 변경을 해야 합니다.

=== "후"

    ```groovy title="main.nf" linenums="13" hl_lines="7-19"
        // langid를 실행하여 각 인사말의 언어를 식별합니다
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

=== "전"

    ```groovy title="main.nf" linenums="13" hl_lines="7"
        // langid를 실행하여 각 인사말의 언어를 식별합니다
        IDENTIFY_LANGUAGE(ch_datasheet)
        IDENTIFY_LANGUAGE.out
            .map { meta, file, lang_id ->
                [meta + [lang: lang_id], file]
            }
            .view()
    ```

주요 사항은 다음과 같습니다.

- `def lang_group = "unknown"`을 사용하여 기본값 `unknown`으로 `lang_group` 변수를 생성합니다.
- 두 게르만어에 대한 `.equals()` 테스트와 세 로망스어에 대한 목록 존재 테스트를 포함하는 `if {} else if {}` 구조를 조건 로직에 사용합니다.
- 이전과 같이 `meta + [lang_group:lang_group]` 병합 작업을 사용하여 업데이트된 메타 맵을 생성합니다.

<!-- TODO (future) Add note/links to relevant docs in additional resources section -->

모두 이해했다면, 워크플로우를 다시 실행하여 결과를 확인합니다.

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

채널 요소가 `[meta, file]` 구조를 유지하면서 메타 맵에 새로운 분류가 포함된 것을 확인할 수 있습니다.

### 핵심 정리

이 섹션에서 다음 내용을 학습했습니다.

- **입력 메타데이터를 출력 채널에 적용하기:** 이 방식으로 메타데이터를 복사하면 나중에 메타데이터 내용을 기반으로 결과를 연결할 수 있습니다.
- **사용자 정의 키 생성:** 메타 맵에 두 개의 새로운 키를 생성하고, `meta + [new_key:value]`로 기존 메타 맵에 병합했습니다. 하나는 프로세스의 계산된 값을 기반으로, 다른 하나는 `map` 연산자에서 설정한 조건을 기반으로 합니다.

이를 통해 파이프라인을 진행하면서 새로운 메타데이터와 기존 메타데이터를 파일과 연결할 수 있습니다.
프로세스의 일부로 메타데이터를 사용하지 않더라도, 이처럼 메타 맵을 데이터와 연결된 상태로 유지하면 모든 관련 정보를 함께 쉽게 관리할 수 있습니다.

---

## 3. 프로세스에서 메타 맵 정보 활용

메타 맵을 생성하고 업데이트하는 방법을 알았으니, 이제 정말 재미있는 부분으로 넘어갑니다. 실제로 프로세스에서 메타데이터를 활용하는 것입니다.

구체적으로, 워크플로우에 두 번째 단계를 추가하여 각 동물을 ASCII 아트로 그리고 말풍선에 녹음된 텍스트를 표시합니다.
[`cowpy`](https://github.com/jeffbuttars/cowpy)라는 도구를 사용합니다.

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

Hello Nextflow 과정을 진행했다면 이 도구를 이미 사용해 보셨을 것입니다.
그렇지 않더라도 걱정하지 마세요. 진행하면서 필요한 모든 내용을 다룹니다.

### 3.1. 프로세스 가져오기 및 코드 검토

`cowpy` 도구를 래핑하는 `COWPY`라는 사전 작성된 프로세스 모듈을 제공합니다. 워크플로우 블록 앞에 include 문만 추가하면 됩니다.

워크플로우를 다음과 같이 수정합니다.

=== "후"

    ```groovy title="main.nf" linenums="1" hl_lines="4"
    #!/usr/bin/env nextflow

    include { IDENTIFY_LANGUAGE } from './modules/langid.nf'
    include { COWPY } from './modules/cowpy.nf'

    workflow {
    ```

=== "전"

    ```groovy title="main.nf" linenums="1"
    #!/usr/bin/env nextflow

    include { IDENTIFY_LANGUAGE } from './modules/langid.nf'

    workflow {
    ```

모듈 파일을 열어 코드를 확인할 수 있습니다.

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

이 프로세스는 현재 입력 파일(표시할 텍스트 포함)과 ASCII 아트로 그릴 캐릭터를 지정하는 값을 받도록 설계되어 있으며, 일반적으로 명령줄 매개변수를 통해 워크플로우 수준에서 제공됩니다.

### 3.2. 메타 맵 필드를 입력으로 전달

Hello Nextflow 과정에서 `cowpy` 도구를 사용할 때는 명령줄 매개변수를 사용하여 최종 이미지를 그릴 캐릭터를 결정했습니다.
파이프라인 실행당 하나의 이미지만 생성했기 때문에 그 방식이 적합했습니다.

그러나 이 튜토리얼에서는 처리하는 각 대상에 적합한 이미지를 생성하려고 하므로, 명령줄 매개변수를 사용하는 것은 너무 제한적입니다.

다행히 데이터시트에 `character` 열이 있으므로 메타 맵에도 있습니다.
이를 사용하여 각 항목에 대해 프로세스가 사용할 캐릭터를 설정합니다.

이를 위해 세 가지를 해야 합니다.

1. 이전 프로세스에서 나오는 출력 채널에 이름을 지정하여 더 편리하게 작업할 수 있도록 합니다.
2. 관심 있는 정보에 접근하는 방법을 결정합니다.
3. 두 번째 프로세스 실행을 추가하고 정보를 적절히 전달합니다.

시작합니다.

#### 3.2.1. 이전 출력 채널에 이름 지정

이전 조작을 첫 번째 프로세스 `IDENTIFY_LANGUAGE.out`의 출력 채널에 직접 적용했습니다.
해당 채널의 내용을 다음 프로세스에 전달하기 위해(명확하고 읽기 쉬운 방식으로) 자체 이름 `ch_languages`를 부여합니다.

[`set`](https://www.nextflow.io/docs/latest/reference/operator.html#set) 연산자를 사용하면 됩니다.

메인 워크플로우에서 `.view()` 연산자를 `.set { ch_languages }`로 교체하고, 채널을 이름으로 참조할 수 있는지 테스트하는 줄을 추가합니다.

=== "후"

    ```groovy title="main.nf" linenums="14" hl_lines="19 21 22"
        // langid를 실행하여 각 인사말의 언어를 식별합니다
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

        // 임시: ch_languages 내용 확인
        ch_languages.view()
    ```

=== "전"

    ```groovy title="main.nf" linenums="14" hl_lines="19"
        // langid를 실행하여 각 인사말의 언어를 식별합니다
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

실행합니다.

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
    [[id:sampleE, character:stegosaurus, lang:es, lang_group:romance], /workspaces/training/side-quests/metadata/work/f4/fcd3186dc666d5d239ffa6c37d125d/hola.txt]
    [[id:sampleF, character:moose, lang:fr, lang_group:romance], /workspaces/training/side-quests/metadata/work/c3/3b2627f733f278a7088332a5806108/salut.txt]
    [[id:sampleG, character:turtle, lang:it, lang_group:romance], /workspaces/training/side-quests/metadata/work/36/cca6a7dbfa26ac24f9329787a32e9d/ciao.txt]
    ```

이제 채널을 이름으로 참조할 수 있음이 확인되었습니다.

#### 3.2.2. 파일 및 캐릭터 메타데이터 접근

모듈 코드를 보면 `COWPY` 프로세스가 텍스트 파일과 `character` 값을 받아야 한다는 것을 알 수 있습니다.
`COWPY` 프로세스 실행을 작성하려면 채널의 각 요소에서 해당 파일 객체와 메타데이터를 추출하는 방법만 알면 됩니다.

가장 간단한 방법은 `map` 작업을 사용하는 것입니다.

채널에는 `[meta, file]` 구조의 튜플이 포함되어 있으므로 `file` 객체에 직접 접근할 수 있고, 메타 맵 내부에 저장된 `character` 값은 `meta.character`로 참조할 수 있습니다.

메인 워크플로우에서 다음과 같이 코드를 변경합니다.

=== "후"

    ```groovy title="main.nf" linenums="34"
        // 임시: 파일과 캐릭터에 접근
        ch_languages.map { meta, file -> file }.view { file -> "File: " + file }
        ch_languages.map { meta, file -> meta.character }.view { character -> "Character: " + character }
    ```

=== "전"

    ```groovy title="main.nf" linenums="34"
        // 임시: ch_languages 내용 확인
        ch_languages.view()
    ```

`.view` 작업의 출력을 더 읽기 쉽게 만들기 위해 closure(예: `{ file -> "File: " + file }`)를 사용합니다.

실행합니다.

```bash
nextflow run main.nf -resume
```

??? success "명령 출력"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `./main.nf` [cheesy_cantor] DSL2 - revision: 15af9c1ec7

    [43/05df08] IDENTIFY_LANGUAGE (7) [100%] 7 of 7, cached: 7 ✔
    Character: squirrel
    File: /workspaces/training/side-quests/metadata/work/8d/4b9498bbccb7a74f04e41877cdc3e5/bonjour.txt
    File: /workspaces/training/side-quests/metadata/work/d3/604274985406e40d79021dea658e60/guten_tag.txt
    Character: tux
    Character: turkey
    File: /workspaces/training/side-quests/metadata/work/d4/fafcc9415b61d2b0fea872e6a05e8a/hello.txt
    File: /workspaces/training/side-quests/metadata/work/02/468ac9efb27f636715e8144b37e9a7/hallo.txt
    Character: sheep
    Character: moose
    Character: stegosaurus
    File: /workspaces/training/side-quests/metadata/work/d4/61a7e1188b4f2742bc72004e226eca/salut.txt
    File: /workspaces/training/side-quests/metadata/work/ae/68364be238c11149c588bf6fc858b1/hola.txt
    File: /workspaces/training/side-quests/metadata/work/43/05df081af5d879ab52e5828fa0357e/ciao.txt
    Character: turtle
    ```

_파일 경로와 캐릭터 값은 출력에서 다른 순서로 나올 수 있습니다._

채널의 각 요소에서 파일과 캐릭터에 접근할 수 있음이 확인되었습니다.

#### 3.2.3. `COWPY` 프로세스 실행

이제 모든 것을 합쳐 `ch_languages` 채널에서 실제로 `COWPY` 프로세스를 실행합니다.

메인 워크플로우에서 다음과 같이 코드를 변경합니다.

=== "후"

    ```groovy title="main.nf" linenums="34"
        // cowpy를 실행하여 ASCII 아트 생성
        COWPY(
            ch_languages.map { meta, file -> file },
            ch_languages.map { meta, file -> meta.character }
        )
    ```

=== "전"

    ```groovy title="main.nf" linenums="34"
        // 임시: 파일과 캐릭터에 접근
        ch_languages.map { meta, file -> [file, meta.character] }
            .view()
    ```

두 map 작업(`.view()` 문 제외)을 프로세스 실행의 입력으로 그대로 복사합니다.
두 입력 사이의 쉼표를 잊지 마세요!

다소 번거롭지만, 다음 섹션에서 더 나은 방법을 살펴봅니다.

실행합니다.

```bash
nextflow run main.nf -resume
```

??? success "명령 출력"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [suspicious_crick] DSL2 - revision: 25541014c5

    executor >  local (7)
    [43/05df08] IDENTIFY_LANGUAGE (7) [100%] 7 of 7, cached: 7 ✔
    [e7/317c18] COWPY (6)             [100%] 7 of 7 ✔
    ```

results 디렉토리를 확인하면 해당 캐릭터가 말하는 각 인사말의 ASCII 아트가 담긴 개별 파일을 볼 수 있습니다.

??? abstract "디렉토리 및 예제 파일 내용"

    ```console
    results/
    ├── cowpy-bonjour.txt
    ├── cowpy-ciao.txt
    ├── cowpy-guten_tag.txt
    ├── cowpy-hallo.txt
    ├── cowpy-hello.txt
    ├── cowpy-hola.txt
    └── cowpy-salut.txt
    ```

    ```text title="results/cowpy-bonjour.txt"
     _________________
    / Bonjour         \
    \ Salut, à demain /
    -----------------
      \
        \
                      _ _
          | \__/|  .~    ~.
          /oo `./      .'
          {o__,   \    {
            / .  . )    \
            `-` '-' \    }
          .(   _(   )_.'
          '---.~_ _ _|
    ```

메타 맵의 정보를 사용하여 파이프라인의 두 번째 단계에서 명령을 매개변수화할 수 있었습니다.

그러나 앞서 언급했듯이, 워크플로우 본문의 맥락에서 메타데이터를 언패킹해야 했기 때문에 일부 코드가 다소 번거로웠습니다.
이 방식은 메타 맵에서 소수의 필드를 사용할 때는 잘 작동하지만, 더 많은 필드를 사용하려면 확장성이 떨어집니다.

`multiMap()`이라는 또 다른 연산자를 사용하면 이를 약간 간소화할 수 있지만, 그래도 이상적이지는 않습니다.

??? info "(선택 사항) `multiMap()`을 사용한 대안 버전"

    참고로, `file`과 `character`를 모두 출력하는 단일 `map()` 작업을 작성할 수 없었습니다. 그렇게 하면 튜플로 반환되기 때문입니다.
    `file`과 `character` 요소를 프로세스에 별도로 전달하기 위해 두 개의 별도 `map()` 작업을 작성해야 했습니다.

    기술적으로는 여러 채널을 내보낼 수 있는 `multiMap()` 연산자를 사용하여 단일 매핑 작업으로 이를 수행하는 다른 방법이 있습니다.
    예를 들어, 위의 `COWPY` 실행을 다음 코드로 교체할 수 있습니다.

    === "후"

        ```groovy title="main.nf" linenums="34"
            // cowpy를 실행하여 ASCII 아트 생성
            COWPY(
                ch_languages.multiMap { meta, file ->
                    file: file
                    character: meta.character
                }
            )
        ```

    === "전"

        ```groovy title="main.nf" linenums="34"
            // cowpy를 실행하여 ASCII 아트 생성
            COWPY(
                ch_languages.map { meta, file -> file },
                ch_languages.map { meta, file -> meta.character }
            )
        ```

    이는 정확히 동일한 결과를 생성합니다.

어느 경우든, 워크플로우 수준에서 언패킹을 해야 한다는 것이 번거롭습니다.

전체 메타 맵을 프로세스에 전달하고 필요한 것을 그 안에서 선택할 수 있다면 더 좋을 것입니다.

### 3.3. 전체 메타 맵 전달 및 활용

메타 맵의 목적은 결국 모든 메타데이터를 하나의 묶음으로 함께 전달하는 것입니다.
위에서 그렇게 할 수 없었던 유일한 이유는 프로세스가 메타 맵을 받도록 설정되어 있지 않았기 때문입니다.
하지만 프로세스 코드를 제어하므로 변경할 수 있습니다.

`COWPY` 프로세스를 수정하여 첫 번째 프로세스에서 사용한 `[meta, file]` 튜플 구조를 받도록 하여 워크플로우를 간소화합니다.

이를 위해 세 가지를 해야 합니다.

1. `COWPY` 프로세스 모듈의 입력 정의 수정
2. 메타 맵을 사용하도록 프로세스 명령 업데이트
3. 워크플로우 본문에서 프로세스 실행 업데이트

준비되셨나요? 시작합니다!

#### 3.3.1. `COWPY` 모듈 입력 수정

`cowpy.nf` 모듈 파일을 다음과 같이 수정합니다.

=== "후"

    ```groovy title="cowpy.nf" linenums="10" hl_lines="2"
    input:
    tuple val(meta), path(input_file)
    ```

=== "전"

    ```groovy title="cowpy.nf" linenums="10" hl_lines="2-3"
    input:
    path(input_file)
    val character
    ```

이를 통해 튜토리얼 앞부분에서 다룬 `[meta, file]` 튜플 구조를 사용할 수 있습니다.

튜토리얼을 간소화하기 위해 프로세스 출력 정의를 메타 맵을 출력하도록 업데이트하지 않았지만, `IDENTIFY_LANGUAGE` 프로세스를 모델로 직접 연습해 보세요.

#### 3.3.2. 메타 맵 필드를 사용하도록 명령 업데이트

이제 전체 메타 맵을 프로세스 내부에서 사용할 수 있으므로, 명령 블록 내부에서 직접 포함된 정보를 참조할 수 있습니다.

`cowpy.nf` 모듈 파일을 다음과 같이 수정합니다.

=== "후"

    ```groovy title="cowpy.nf" linenums="16" hl_lines="3"
    script:
    """
    cat ${input_file} | cowpy -c ${meta.character} > cowpy-${input_file}
    """
    ```

=== "전"

    ```groovy title="cowpy.nf" linenums="16" hl_lines="3"
    script:
    """
    cat ${input_file} | cowpy -c ${character} > cowpy-${input_file}
    """
    ```

이전에 독립 입력으로 전달된 `character` 값에 대한 참조를 `meta.character`를 사용하여 메타 맵에 저장된 값으로 교체했습니다.

이제 프로세스 실행을 업데이트합니다.

#### 3.3.3. 프로세스 실행 업데이트 및 실행

프로세스는 이제 `[meta, file]` 튜플 구조를 입력으로 받도록 설정되어 있으며, 이는 이전 프로세스가 출력하는 것과 동일하므로 전체 `ch_languages` 채널을 `COWPY` 프로세스에 전달할 수 있습니다.

메인 워크플로우를 다음과 같이 수정합니다.

=== "후"

    ```groovy title="main.nf" linenums="34" hl_lines="2"
    // cowpy를 실행하여 ASCII 아트 생성
    COWPY(ch_languages)
    ```

=== "전"

    ```groovy title="main.nf" linenums="34" hl_lines="3-4"
    // cowpy를 실행하여 ASCII 아트 생성
    COWPY(
        ch_languages.map { meta, file -> file },
        ch_languages.map { meta, file -> meta.character }
    )
    ```

실행이 훨씬 간단해졌습니다!

이전 실행 결과를 삭제하고 실행합니다.

```bash
rm -r results
nextflow run main.nf
```

??? success "명령 출력"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [wise_sammet] DSL2 - revision: 99797b1e92

    executor >  local (14)
    [5d/dffd4e] process > IDENTIFY_LANGUAGE (7) [100%] 7 of 7 ✔
    [25/9243df] process > COWPY (7)             [100%] 7 of 7 ✔
    ```

results 디렉토리를 확인하면 이전과 동일한 출력, 즉 해당 캐릭터가 말하는 각 인사말의 ASCII 아트가 담긴 개별 파일을 볼 수 있습니다.

??? abstract "디렉토리 내용"

    ```console
    ./results/
    ├── cowpy-bonjour.txt
    ├── cowpy-ciao.txt
    ├── cowpy-guten_tag.txt
    ├── cowpy-hallo.txt
    ├── cowpy-hello.txt
    ├── cowpy-hola.txt
    └── cowpy-salut.txt
    ```

더 간단한 코드로 이전과 동일한 결과를 생성합니다.

물론, 이는 프로세스 코드를 수정할 수 있다고 가정합니다.
경우에 따라 수정할 수 없는 기존 프로세스에 의존해야 할 수도 있으며, 이는 선택지를 제한합니다.
[nf-core](https://nf-co.re/) 프로젝트의 모듈을 사용할 계획이라면 좋은 소식이 있습니다. nf-core 모듈은 모두 `[meta, file]` 튜플 구조를 표준으로 사용하도록 설정되어 있습니다.

### 3.4. 누락된 필수 입력 문제 해결

`character` 값은 `COWPY` 프로세스가 성공적으로 실행되기 위해 필요합니다.
설정 파일에서 기본값을 설정하지 않으면 데이터시트에서 반드시 값을 제공해야 합니다.

**제공하지 않으면 어떻게 될까요?**
입력 데이터시트의 내용과 실행 중인 워크플로우 버전에 따라 다릅니다.

#### 3.4.1. character 열이 존재하지만 비어 있는 경우

데이터 수집 오류를 시뮬레이션하기 위해 데이터시트의 한 항목에서 character 값을 삭제한다고 가정합니다.

```csv title="datasheet.csv" linenums="1" hl_lines="2"
id,character,recording
sampleA,,/workspaces/training/side-quests/metadata/data/bonjour.txt
sampleB,tux,/workspaces/training/side-quests/metadata/data/guten_tag.txt
sampleC,sheep,/workspaces/training/side-quests/metadata/data/hallo.txt
sampleD,turkey,/workspaces/training/side-quests/metadata/data/hello.txt
sampleE,stegosaurus,/workspaces/training/side-quests/metadata/data/hola.txt
sampleF,moose,/workspaces/training/side-quests/metadata/data/salut.txt
sampleG,turtle,/workspaces/training/side-quests/metadata/data/ciao.txt
```

위에서 사용한 워크플로우 버전 모두에서, 데이터시트를 읽을 때 모든 항목에 대해 `character` 키가 생성되지만 `sampleA`의 경우 값이 빈 문자열이 됩니다.

이로 인해 오류가 발생합니다.

??? failure "명령 출력"

    ```console hl_lines="8 11 16 28"
     N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [marvelous_hirsch] DSL2 - revision: 0dfeee3cc1

    executor >  local (9)
    [c1/c5dd4f] process > IDENTIFY_LANGUAGE (7) [ 85%] 6 of 7
    [d3/b7c415] process > COWPY (2)             [  0%] 0 of 6
    ERROR ~ Error executing process > 'COWPY (1)'

    Caused by:
      Process `COWPY (1)` terminated with an error exit status (2)


    Command executed:

      cat bonjour.txt | cowpy -c  > cowpy-bonjour.txt

    Command exit status:
      2

    Command output:
      (empty)

    Command error:
      usage: cowpy [-h] [-l] [-L] [-t] [-u] [-e EYES] [-c COWACTER] [-E] [-r] [-x]
                  [-C]
                  [msg ...]
      cowpy: error: argument -c/--cowacter: expected one argument

    Work dir:
      /workspaces/training/side-quests/metadata/work/ca/9d49796612a54dec5ed466063c809b

    Container:
      community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273

    Tip: you can try to figure out what's wrong by changing to the process work dir and showing the script file named `.command.sh`

    -- Check '.nextflow.log' file for details
    ```

Nextflow가 해당 샘플에 대해 `cowpy` 명령줄을 실행할 때 `${meta.character}`가 빈 문자열로 채워지므로, `cowpy` 도구는 `-c` 인자에 값이 제공되지 않았다는 오류를 발생시킵니다.

#### 3.4.2. 데이터시트에 character 열이 없는 경우

이번에는 데이터시트에서 `character` 열을 완전히 삭제한다고 가정합니다.

```csv title="datasheet.csv" linenums="1"
id,recording
sampleA,/workspaces/training/side-quests/metadata/data/bonjour.txt
sampleB,/workspaces/training/side-quests/metadata/data/guten_tag.txt
sampleC,/workspaces/training/side-quests/metadata/data/hallo.txt
sampleD,/workspaces/training/side-quests/metadata/data/hello.txt
sampleE,/workspaces/training/side-quests/metadata/data/hola.txt
sampleF,/workspaces/training/side-quests/metadata/data/salut.txt
sampleG,/workspaces/training/side-quests/metadata/data/ciao.txt
```

이 경우 데이터시트를 읽을 때 `character` 키가 전혀 생성되지 않습니다.

##### 3.4.2.1. 워크플로우 수준에서 값에 접근하는 경우

섹션 3.2에서 작성한 코드 버전을 사용하는 경우, Nextflow는 `COWPY` 프로세스를 실행하기 전에 메타 맵에서 `character` 키에 접근하려고 시도합니다.

해당 명령과 일치하는 요소를 찾지 못하므로 `COWPY`를 전혀 실행하지 않습니다.

??? success "명령 출력"

    ```console hl_lines="7"
     N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [desperate_montalcini] DSL2 - revision: 0dfeee3cc1

    executor >  local (7)
    [1a/df2544] process > IDENTIFY_LANGUAGE (7) [100%] 7 of 7 ✔
    [-        ] process > COWPY                 -
    ```

Nextflow 입장에서는 이 워크플로우가 성공적으로 실행된 것입니다!
그러나 원하는 출력이 전혀 생성되지 않습니다.

##### 3.4.2.2. 프로세스 수준에서 값에 접근하는 경우

섹션 3.3의 버전을 사용하는 경우, Nextflow는 전체 메타 맵을 `COWPY` 프로세스에 전달하고 명령 실행을 시도합니다.

이로 인해 오류가 발생하지만, 첫 번째 경우와는 다른 오류입니다.

??? failure "명령 출력"

    ```console hl_lines="8 11 16"
     N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [jovial_bohr] DSL2 - revision: eaaf375827

    executor >  local (9)
    [0d/ada9db] process > IDENTIFY_LANGUAGE (5) [ 85%] 6 of 7
    [06/28065f] process > COWPY (2)             [  0%] 0 of 6
    ERROR ~ Error executing process > 'COWPY (2)'

    Caused by:
      Process `COWPY (2)` terminated with an error exit status (1)


    Command executed:

      cat guten_tag.txt | cowpy -c null > cowpy-guten_tag.txt

    Command exit status:
      1

    Command output:
      (empty)

    Command error:
      Traceback (most recent call last):
        File "/opt/conda/bin/cowpy", line 10, in <module>
          sys.exit(main())
                  ~~~~^^
        File "/opt/conda/lib/python3.13/site-packages/cowpy/cow.py", line 1215, in main
          print(cow(eyes=args.eyes,
                ~~~^^^^^^^^^^^^^^^^
                tongue=args.tongue,
                ^^^^^^^^^^^^^^^^^^^
                thoughts=args.thoughts
                ^^^^^^^^^^^^^^^^^^^^^^
                    ).milk(msg)
                    ^
      TypeError: 'str' object is not callable

    Work dir:
      /workspaces/training/side-quests/metadata/work/06/28065f7d9fd7d22bba084aa941b6d6

    Container:
      community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273

    Tip: you can replicate the issue by changing to the process work dir and entering the command `bash .command.run`

    -- Check '.nextflow.log' file for details
    ```

`meta.character`가 존재하지 않으므로 접근 시도가 `null`을 반환합니다. 결과적으로 Nextflow는 명령줄에 `null`을 그대로 삽입하며, 이는 당연히 `cowpy` 도구에서 인식되지 않습니다.

#### 3.4.3. 해결책

워크플로우 설정의 일부로 기본값을 제공하는 것 외에도, 더 견고하게 처리하기 위해 두 가지를 할 수 있습니다.

1. 워크플로우에 입력 유효성 검사를 구현하여 데이터시트에 필요한 모든 정보가 포함되어 있는지 확인합니다. Hello nf-core 교육 과정의 [입력 유효성 검사 소개](../hello_nf-core/05_input_validation.md)를 참조하세요. <!-- TODO (future) pending a proper Validation side quest -->

2. 프로세스 모듈을 사용하는 누구나 필수 입력을 즉시 식별할 수 있도록 하려면, 필수 메타데이터 속성을 명시적 입력으로 만들 수도 있습니다.

다음은 그 방법의 예시입니다.

먼저, 프로세스 수준에서 입력 정의를 다음과 같이 업데이트합니다.

=== "후"

    ```groovy title="cowpy.nf" linenums="12" hl_lines="2"
        input:
        tuple val(meta), val(character), path(input_file)
    ```

=== "전"

    ```groovy title="cowpy.nf" linenums="12" hl_lines="2"
        input:
        tuple val(meta), path(input_file)
    ```

그런 다음, 워크플로우 수준에서 매핑 작업을 사용하여 메타데이터에서 `character` 속성을 추출하고 입력 튜플의 명시적 구성 요소로 만듭니다.

=== "후"

    ```groovy title="main.nf" linenums="37" hl_lines="1"
        COWPY(ch_languages.map{meta, file -> [meta, meta.character, file]})
    ```

=== "전"

    ```groovy title="main.nf" linenums="37" hl_lines="1"
        COWPY(ch_languages)
    ```

이 방식은 `character`가 필수임을 명시적으로 보여주고, 다른 맥락에서 프로세스를 재배포하기 쉽게 만드는 장점이 있습니다.

이는 중요한 설계 원칙을 강조합니다.

**메타 맵은 선택적이고 설명적인 정보에 사용하되, 필수 값은 명시적 입력으로 추출합니다.**

메타 맵은 채널 구조를 깔끔하게 유지하고 임의적인 채널 구조를 방지하는 데 탁월하지만, 프로세스에서 직접 참조되는 필수 요소의 경우 명시적 입력으로 추출하면 더 견고하고 유지 관리하기 쉬운 코드가 됩니다.

### 핵심 정리

이 섹션에서는 메타데이터를 활용하여 프로세스 실행을 맞춤화하는 방법을 학습했습니다. 워크플로우 수준 또는 프로세스 수준에서 메타데이터에 접근할 수 있습니다.

---

## 보충 연습

프로세스 내부에서 메타 맵 정보를 사용하는 연습을 하고 싶다면, `lang` 및 `lang_group`과 같은 메타 맵의 다른 정보를 사용하여 출력의 이름 지정 및/또는 구성을 맞춤화해 보세요.

예를 들어, 다음 결과를 생성하도록 코드를 수정해 보세요.

```console title="Results directory contents"
results/
├── germanic
│   ├── de-guten_tag.txt
│   ├── de-hallo.txt
│   └── en-hello.txt
└── romance
    ├── es-hola.txt
    ├── fr-bonjour.txt
    ├── fr-salut.txt
    └── it-ciao.txt
```

<!-- TODO (future) Provide worked out solution -->
<!-- the renaming should use the meta inside the process -->
<!-- the output org should use the meta in the workflow outputs -->

---

## 요약

이 사이드 퀘스트에서는 Nextflow 워크플로우에서 메타데이터를 효과적으로 다루는 방법을 살펴보았습니다.

메타데이터를 명시적으로 유지하고 데이터에 연결하는 이 패턴은 Nextflow의 핵심 모범 사례로, 파일 정보를 하드코딩하는 것보다 여러 가지 장점을 제공합니다.

- 파일 메타데이터가 워크플로우 전반에 걸쳐 파일과 연결된 상태를 유지합니다.
- 파일별로 프로세스 동작을 맞춤화할 수 있습니다.
- 출력 구성이 파일 메타데이터를 반영할 수 있습니다.
- 파이프라인 실행 중에 파일 정보를 확장할 수 있습니다.

이 패턴을 자신의 작업에 적용하면 견고하고 유지 관리하기 쉬운 생물정보학 워크플로우를 구축할 수 있습니다.

### 핵심 패턴

1.  **메타데이터 읽기 및 구조화:** CSV 파일을 읽고 데이터 파일과 연결된 상태를 유지하는 구조화된 메타데이터 맵을 생성합니다.

    ```groovy
    channel.fromPath('datasheet.csv')
      .splitCsv(header: true)
      .map { row ->
          [ [id:row.id, character:row.character], row.recording ]
      }
    ```

2.  **워크플로우 중 메타데이터 확장:** 프로세스 출력을 추가하고 조건 로직을 통해 값을 도출하여 파이프라인이 진행되면서 메타데이터에 새로운 정보를 추가합니다.

    - 프로세스 출력을 기반으로 새 키 추가

    ```groovy
    .map { meta, file, lang ->
      [ meta + [lang:lang], file ]
    }
    ```

    - 조건절을 사용하여 새 키 추가

    ```groovy
    .map{ meta, file ->
        if ( meta.lang.equals("de") || meta.lang.equals('en') ){
            lang_group = "germanic"
        } else if ( meta.lang in ["fr", "es", "it"] ) {
            lang_group = "romance"
        } else {
            lang_group = "unknown"
        }
    }
    ```

3.  **프로세스 동작 맞춤화:** 프로세스 내부에서 메타데이터를 활용합니다.

    ```groovy
    cat $input_file | cowpy -c ${meta.character} > cowpy-${input_file}
    ```

### 추가 자료

- [map](https://www.nextflow.io/docs/latest/operator.html#map)
- [stdout](https://www.nextflow.io/docs/latest/process.html#outputs)

---

## 다음 단계

[사이드 퀘스트 메뉴](../)로 돌아가거나 페이지 오른쪽 하단의 버튼을 클릭하여 목록의 다음 주제로 이동하세요.
