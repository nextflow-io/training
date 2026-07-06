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
- "메타 맵 + 데이터 파일" 인터페이스가 널리 사용되는 관례인 이유 이해
- 워크플로우 실행 중 새로운 메타데이터 필드 추가
- 메타데이터를 활용한 프로세스 동작 맞춤화 및 출력 구성

이러한 기술을 통해 복잡한 파일 관계와 처리 요구사항을 처리할 수 있는 더욱 견고하고 유연한 파이프라인을 구축할 수 있습니다.

### 사전 요구사항

이 사이드 퀘스트를 시작하기 전에 다음을 완료해야 합니다.

- [Hello Nextflow](../../hello_nextflow/index.md) 튜토리얼 또는 동급의 입문 과정을 완료해야 합니다.
- 기본적인 Nextflow 개념과 메커니즘(프로세스, 채널, 연산자)에 익숙해야 합니다.

---

## 0. 시작하기

#### 교육 코드스페이스 열기

아직 열지 않았다면 [환경 설정](../../envsetup/index.md)에 설명된 대로 교육 환경을 열어 주세요.

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

편집기가 해당 프로젝트 디렉토리에 포커스된 상태로 열립니다.

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

[`COWPY`](https://github.com/jeffbuttars/cowpy)라는 도구를 사용하여 각 캐릭터가 녹음된 인사말을 말하는 ASCII 아트를 생성합니다.

??? info "`COWPY`는 무엇을 하나요?"

    `COWPY`는 임의의 텍스트 입력을 재미있는 방식으로 표시하는 ASCII 아트를 생성하는 명령줄 도구입니다.
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

또한 `langid`라는 언어 분석 도구를 사용하여 각 캐릭터가 사용하는 언어를 식별하고, 그에 따라 파이프라인의 출력을 구성합니다.

#### 과제 검토

여러분의 과제는 다음을 수행하는 Nextflow 워크플로우를 작성하는 것입니다.

1. 각 캐릭터의 **ASCII 아트를 생성**합니다.
2. 언어 계열(게르만어 vs 로망스어)에 따라 출력을 **구성**합니다.

이는 파일별 메타데이터가 처리 결정을 이끄는 전형적인 워크플로우 패턴으로, 메타데이터 맵이 우아하게 해결하는 바로 그런 문제입니다.

#### 준비 체크리스트

시작할 준비가 되었나요?

- [ ] 이 과정의 목표와 사전 요구사항을 이해했습니다.
- [ ] 코드스페이스가 실행 중입니다.
- [ ] 작업 디렉토리를 적절히 설정했습니다.
- [ ] 과제를 이해했습니다.

모든 항목을 확인했다면 시작할 준비가 된 것입니다.

---

## 1. 메타데이터를 불러오고 사용하는 기본 방법

`main.nf` 워크플로우 파일을 열어 시작점으로 제공된 워크플로우 스텁을 확인합니다.

```groovy title="main.nf" linenums="1"
#!/usr/bin/env nextflow

workflow  {
    main:
    ch_datasheet = channel.fromPath("./data/datasheet.csv")
        .splitCsv(header: true)
        .view()

    publish:
    cowpy_art = channel.empty()
}

output {
    cowpy_art {
    }
}
```

[`splitCsv`](https://www.nextflow.io/docs/latest/reference/operator.html#splitcsv) 연산자는 파일의 각 행을 채널 요소로 읽어옵니다.
이는 입문 과정인 Hello Nextflow에서 CSV 데이터를 불러올 때 사용하는 방식과 동일합니다.
작동 방식이 기억나지 않는다면 [이 섹션](../../hello_nextflow/02_hello_channels.md#4-read-input-values-from-a-csv-file)을 참조하세요.

`header: true` 옵션을 사용하면 첫 번째 행이 열 헤더로 처리되어, 각 요소가 열 이름을 키로 하는 키-값 쌍의 맵이 됩니다.

아직 데이터에 대해 어떤 프로세스도 실행하지 않으므로, `publish`와 `output` 블록은 스텁 상태입니다.

### 1.1. 워크플로우 실행

워크플로우를 실행하여 모든 데이터가 불러와진 후 채널 내용이 어떻게 구성되는지 확인합니다.

```bash
nextflow run main.nf
```

??? success "명령 출력"

    ```console
     N E X T F L O W   ~  version 25.10.4

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

이를 통해 각 행의 특정 필드에 쉽게 접근할 수 있습니다.
예를 들어, `id`로 파일 ID에 접근하거나 `recording`으로 txt 파일 경로에 접근할 수 있습니다.

??? info "(선택 사항) Groovy 맵에 대한 추가 정보"

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
     N E X T F L O W   ~  version 25.10.4

    Launching `map_demo.nf` [cheesy_plateau] DSL2 - revision: fae5b8496e

    map: [id:sampleA, character:squirrel]
    id: sampleA
    character: squirrel
    ```

### 1.2. `map`으로 특정 필드 선택

`map` 연산자를 사용하여 채널의 각 요소를 순회하고, 점 표기법으로 이름에 접근하여 `character` 필드만 선택합니다.

#### 1.2.1. map 작업 추가

`character` 열에 접근하기 위해 `.view()` 작업 앞에 다음과 같이 `map` 작업을 추가합니다.

=== "후"

    ```groovy title="main.nf" linenums="5" hl_lines="3-5"
        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .map { row ->
                row.character
            }
            .view()
    ```

=== "전"

    ```groovy title="main.nf" linenums="5"
        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .view()
    ```

특정 필드에 접근하는 이 방식은 Hello Nextflow의 [이 섹션](../../hello_nextflow/02_hello_channels.md#43-use-the-map-operator-to-extract-the-greetings)에서 더 자세히 설명합니다.

#### 1.2.2. 워크플로우 실행

워크플로우를 실행하여 추출된 캐릭터 이름을 확인합니다.

```bash
nextflow run main.nf
```

??? success "명령 출력"

    ```console
     N E X T F L O W   ~  version 25.10.4

    Launching `main.nf` [exotic_albattani] DSL2 - revision: c0d03cec83

    squirrel
    tux
    sheep
    turkey
    stegosaurus
    moose
    turtle
    ```

각 행의 `character` 열 값에 접근할 수 있음이 확인되었습니다.

이제 이 데이터를 활용하여 `character`와 `recording` 필드를 함께 사용해 `COWPY`로 ASCII 아트를 생성합니다.

### 1.3. `multiMap`으로 서브 채널 내보내기

사전 작성된 `COWPY` 프로세스 모듈을 제공하므로, 먼저 프로세스의 입력 요구사항을 확인합니다.

파일을 열어 프로세스 코드를 확인합니다.

```groovy title="modules/cowpy.nf" linenums="1" hl_lines="7 8"
// cowpy로 ASCII 아트 생성
process COWPY {

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

프로세스가 두 개의 별도 입력(녹음 파일과 캐릭터 이름)을 받는 것을 확인할 수 있습니다.
두 값 모두 존재하지만, 현재 채널의 각 요소 안에 묶여 있습니다.

여러 필드를 별도의 채널로 추출하는 방법 중 하나는 [`multiMap`](https://www.nextflow.io/docs/latest/reference/operator.html#multimap) 연산자입니다. 이 연산자는 하나의 채널을 단일 작업으로 여러 개의 이름 있는 서브 채널로 분리합니다.

#### 1.3.1. multiMap 작업 추가

`map` 작업을 `multiMap`으로 교체합니다.

=== "후"

    ```groovy title="main.nf" linenums="5" hl_lines="3-6"
        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .multiMap { row ->
                file: row.recording
                character: row.character
            }
    ```

=== "전"

    ```groovy title="main.nf" linenums="5" hl_lines="3-6"
        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .map { row ->
                row.character
            }
            .view()
    ```

`multiMap` 블록은 각 행에서 두 개의 이름 있는 서브 채널(`file`과 `character`)을 정의하며, `ch_datasheet.file`과 `ch_datasheet.character`로 접근할 수 있습니다.

#### 1.3.2. 서브 채널에서 COWPY 실행

`COWPY` 프로세스를 포함하고 각 서브 채널을 별도의 인자로 전달합니다.

=== "후"

    ```groovy title="main.nf" linenums="1" hl_lines="3 14"
    #!/usr/bin/env nextflow

    include { COWPY } from './modules/cowpy.nf'

    workflow {
        main:
        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .multiMap { row ->
                file: row.recording
                character: row.character
            }

        COWPY(ch_datasheet.file, ch_datasheet.character)

        publish:
        cowpy_art = channel.empty()
    }

    output {
        cowpy_art {
        }
    }
    ```

=== "전"

    ```groovy title="main.nf" linenums="1"
    #!/usr/bin/env nextflow

    workflow {
        main:
        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .multiMap { row ->
                file: row.recording
                character: row.character
            }

        publish:
        cowpy_art = channel.empty()
    }

    output {
        cowpy_art {
        }
    }
    ```

이를 통해 `COWPY`가 요구하는 대로 두 필드를 별도로 전달할 수 있습니다.

#### 1.3.3. 출력 게시 설정

마지막으로 `COWPY`의 출력을 `publish:` 블록에 추가합니다.

=== "후"

    ```groovy title="main.nf" linenums="14" hl_lines="4"
        COWPY(ch_datasheet.file, ch_datasheet.character)

        publish:
        cowpy_art = COWPY.out
    ```

=== "전"

    ```groovy title="main.nf" linenums="14" hl_lines="4"
        COWPY(ch_datasheet.file, ch_datasheet.character)

        publish:
        cowpy_art = channel.empty()
    ```

이를 통해 워크플로우가 생성한 출력을 쉽게 확인할 수 있습니다.

#### 1.3.4. 워크플로우 실행

워크플로우를 실행하여 `COWPY`가 제공한 입력에서 올바르게 실행되는지 확인합니다.

```bash
nextflow run main.nf
```

??? success "명령 출력"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [clever_dijkstra] DSL2 - revision: a1b2c3d4e5

    executor >  local (7)
    [3a/f1c290] COWPY (7) [100%] 7 of 7 ✔
    ```

`COWPY`가 각 파일에 대해 올바른 캐릭터를 사용하여 실행된 것을 확인할 수 있습니다.

??? abstract "결과 디렉토리 내용"

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

??? example "results/cowpy-guten_tag.txt 내용"

    ```console
    $ cat results/cowpy-guten_tag.txt
     _____________________________
    / Guten Tag, wie geht es dir? \
    \ Auf Wiedersehen, bis morgen /
     -----------------------------
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

이 방식은 작동하지만 한 가지 제한이 있습니다. 채널을 두 개의 별도 서브 채널로 분리해야 했습니다.
프로세스에 더 많은 필드를 전달하려면 더 많은 서브 채널로 분리해야 하며, 이는 번거롭고 복잡해질 수 있습니다.

더 간단한 방법이 있습니다.

### 1.4. 모든 것을 단일 입력으로 프로세스에 전달

필드를 별도의 채널로 분리하는 대신, 프로세스가 모든 입력을 단일 튜플로 받도록 업데이트하면 프로세스 실행이 간소화됩니다.

#### 1.4.1. COWPY 프로세스 업데이트

각 행의 세 요소에 해당하는 튜플을 받도록 `COWPY`를 업데이트합니다.

=== "후"

    ```groovy title="modules/cowpy.nf" linenums="1" hl_lines="7 10 14"
    // cowpy로 ASCII 아트 생성
    process COWPY {

        container 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'

        input:
        tuple val(id), val(character), path(recording)

        output:
        path "cowpy-${recording}"

        script:
        """
        cat ${recording} | cowpy -c ${character} > cowpy-${recording}
        """
    }
    ```

=== "전"

    ```groovy title="modules/cowpy.nf" linenums="1" hl_lines="7-8 11 15"
    // cowpy로 ASCII 아트 생성
    process COWPY {

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

이제 프로세스는 필요한 모든 것을 포함하는 단일 입력을 받습니다.

#### 1.4.2. `map()`으로 입력 튜플 생성

프로세스에 전달할 튜플의 요소를 열거하기 위해 매핑 작업을 사용합니다.

=== "후"

    ```groovy title="main.nf" linenums="5" hl_lines="3-5"
        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .map { row ->
                [row.id, row.character, row.recording]
            }
    ```

=== "전"

    ```groovy title="main.nf" linenums="5" hl_lines="3-6"
        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .multiMap { row ->
                file: row.recording
                character: row.character
            }
    ```

`splitCsv`에서 나오는 Groovy 맵 전체를 그대로 전달할 수 없는 이유가 궁금할 수 있습니다.
Nextflow에게 녹음 파일이 경로로 처리되어야 한다는 것(즉, 올바르게 스테이징되어야 한다는 것)을 명시적으로 알려야 하기 때문입니다.
이는 `COWPY`의 입력 인터페이스 수준에서 `recording` 요소가 `path`로 명시적으로 지정될 때 처리됩니다.

#### 1.4.3. 프로세스 실행 업데이트

프로세스 실행에서 두 개의 별도 입력을 방금 생성한 단일 튜플로 교체합니다.

=== "후"

    ```groovy title="main.nf" linenums="5" hl_lines="7"
        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .map { row ->
                [row.id, row.character, row.recording]
            }

        COWPY(ch_datasheet)
    ```

=== "전"

    ```groovy title="main.nf" linenums="5" hl_lines="7"
        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .map { row ->
                [row.id, row.character, row.recording]
            }

        COWPY(ch_datasheet.file, ch_datasheet.character)
    ```

프로세스 실행이 약간 간소화되었습니다.

#### 1.4.4. 워크플로우 실행

워크플로우를 실행하여 `COWPY`가 데이터를 올바르게 처리하는지 확인합니다.

```bash
nextflow run main.nf
```

??? success "명령 출력"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [pedantic_lovelace] DSL2 - revision: b2c3d4e5f6

    executor >  local (7)
    [5e/2a1b34] COWPY (7) [100%] 7 of 7 ✔
    ```

이전과 동일한 7개의 `cowpy-*.txt` 파일이 출력되며, 이제 더 간단한 `COWPY` 실행으로 생성됩니다.

??? abstract "결과 디렉토리 내용"

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

??? example "results/cowpy-guten_tag.txt 내용"

    ```console
    $ cat results/cowpy-guten_tag.txt
     _____________________________
    / Guten Tag, wie geht es dir? \
    \ Auf Wiedersehen, bis morgen /
     -----------------------------
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

`multiMap` 방식보다 약간 개선되었습니다.
그러나 입력 튜플을 생성하기 위해 원래의 Groovy 맵을 언패킹해야 했으며, 프로세스와 데이터시트 간의 결합이 여전히 강합니다. `COWPY` 입력 정의가 열 이름 `id`, `character`, `recording`을 직접 참조합니다.

```groovy
input:
tuple val(id), val(character), path(recording)
```

협업자가 추가 열이 있거나 열 순서가 다른 데이터시트를 사용한다면, 이 프로세스는 수정 없이는 작동하지 않습니다.
프로세스의 입력 구조가 데이터시트의 정확한 구성에 종속되어 있어 취약합니다.

이를 해결하려면 프로세스 인터페이스에 정확한 구조를 하드코딩하지 않고 모든 메타데이터를 묶음으로 전달하는 방법이 필요합니다.

### 1.5. 메타 맵 + 파일 인터페이스 사용

해결책은 채널에서 두 가지 관심사를 분리하는 것입니다. **샘플에 대한 메타데이터**와 **데이터 파일** 자체입니다.
모든 메타데이터를 단일 맵인 "메타 맵"으로 묶으면, 데이터시트에 메타데이터 열이 몇 개 있든 일관된 두 요소 튜플을 얻을 수 있습니다.

```groovy title="Syntax example"
input:
tuple val(meta), path(file)
```

데이터시트에서 열을 추가하거나 제거하면 `meta` 내부의 내용이 변경되지만, 튜플 형태 `[meta, file]`은 그대로 유지됩니다.
이 구조를 받는 프로세스는 메타데이터 필드가 몇 개 있는지 알 필요가 없습니다.

#### 1.5.1. 튜플 내용을 메타 맵으로 재구성

`map` 작업을 `[meta, file]` 튜플을 생성하도록 재구성합니다.

=== "후"

    ```groovy title="main.nf" linenums="5" hl_lines="4 6 8 11"
        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .map { row ->
                [[id: row.id, character: row.character], row.recording]
            }
            .view()

        // COWPY(ch_datasheet)  // 다음 단계에서 업데이트 예정

        publish:
        cowpy_art = channel.empty() // COWPY.out
    ```

=== "전"

    ```groovy title="main.nf" linenums="5" hl_lines="4 7"
        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .map { row ->
                [row.id, row.character, row.recording]
            }

        COWPY(ch_datasheet)

        publish:
        cowpy_art = COWPY.out
    ```

`view()` 문을 추가하고, `COWPY` 실행을 주석 처리하고, `COWPY.out`을 `channel.empty()`로 교체했습니다. 프로세스 입력 정의가 아직 새로운 구조와 일치하지 않기 때문입니다.

#### 1.5.2. 워크플로우를 실행하여 재구성된 내용 확인

워크플로우를 실행하여 새로운 채널 형태를 확인합니다.

```bash
nextflow run main.nf
```

??? success "명령 출력"

    ```console title="View meta map"
     N E X T F L O W   ~  version 25.10.4

    Launching `main.nf` [lethal_booth] DSL2 - revision: 0d8f844c07

    [[id:sampleA, character:squirrel], /workspaces/training/side-quests/metadata/data/bonjour.txt]
    [[id:sampleB, character:tux], /workspaces/training/side-quests/metadata/data/guten_tag.txt]
    [[id:sampleC, character:sheep], /workspaces/training/side-quests/metadata/data/hallo.txt]
    [[id:sampleD, character:turkey], /workspaces/training/side-quests/metadata/data/hello.txt]
    [[id:sampleE, character:stegosaurus], /workspaces/training/side-quests/metadata/data/hola.txt]
    [[id:sampleF, character:moose], /workspaces/training/side-quests/metadata/data/salut.txt]
    [[id:sampleG, character:turtle], /workspaces/training/side-quests/metadata/data/ciao.txt]
    ```

채널의 각 요소는 이제 두 요소 튜플입니다. 첫 번째는 메타 맵, 두 번째는 파일입니다.

```console title="Example element structure"
[
  [id:sampleA, character:squirrel],
  /workspaces/training/side-quests/metadata/data/bonjour.txt
]
```

나중에 데이터시트에 `language` 열을 추가하면 프로세스 입력 정의를 변경하지 않고도 `meta.language`로 접근할 수 있습니다.

#### 1.5.3. 메타 맵을 사용하도록 `COWPY` 프로세스 업데이트

`[meta, file]` 튜플 구조를 받도록 `COWPY`를 업데이트합니다.

=== "후"

    ```groovy title="modules/cowpy.nf" linenums="1" hl_lines="7 10 14"
    // cowpy로 ASCII 아트 생성
    process COWPY {

        container 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'

        input:
        tuple val(meta), path(input_file)

        output:
        path "cowpy-${input_file}"

        script:
        """
        cat ${input_file} | cowpy -c ${meta.character} > cowpy-${input_file}
        """
    }
    ```

=== "전"

    ```groovy title="modules/cowpy.nf" linenums="1" hl_lines="7 10 14"
    // cowpy로 ASCII 아트 생성
    process COWPY {

        container 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'

        input:
        tuple val(id), val(character), path(recording)

        output:
        path "cowpy-${recording}"

        script:
        """
        cat ${recording} | cowpy -c ${character} > cowpy-${recording}
        """
    }
    ```

스크립트 블록 내부에서 `meta.character`는 메타 맵의 `character` 필드에 접근합니다.
메타 맵의 모든 필드는 동일한 방식으로 접근할 수 있습니다.

#### 1.5.4. 프로세스 실행 업데이트

`COWPY` 실행을 복원하고 출력을 게시에 연결합니다.

=== "후"

    ```groovy title="main.nf" linenums="5" hl_lines="7 10"
        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .map { row ->
                [[id: row.id, character: row.character], row.recording]
            }

        COWPY(ch_datasheet)

        publish:
        cowpy_art = COWPY.out
    ```

=== "전"

    ```groovy title="main.nf" linenums="5" hl_lines="6 8 11"
        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .map { row ->
                [[id: row.id, character: row.character], row.recording]
            }
            .view()

        // COWPY(ch_datasheet)  // 다음 단계에서 업데이트 예정

        publish:
        cowpy_art = channel.empty() // COWPY.out
    ```

출력 게시도 복원했습니다.

#### 1.5.5. 워크플로우 실행

워크플로우를 실행하여 모든 것이 올바르게 작동하는지 확인합니다.

```bash
nextflow run main.nf
```

??? success "명령 출력"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [wise_sammet] DSL2 - revision: 99797b1e92

    executor >  local (7)
    [5d/dffd4e] COWPY (7) [100%] 7 of 7 ✔
    ```

results 디렉토리에 ASCII 아트 파일이 생성되었습니다.

??? abstract "디렉토리 내용"

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

??? example "results/cowpy-guten_tag.txt 내용"

    ```console
    $ cat results/cowpy-guten_tag.txt
     _____________________________
    / Guten Tag, wie geht es dir? \
    \ Auf Wiedersehen, bis morgen /
     -----------------------------
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

프로세스는 이제 `meta`를 통해 모든 메타데이터를 묶음으로 받고, 필요한 것(`meta.character`)을 사용하며, 나머지는 무시합니다.

이것이 모든 [nf-core](https://nf-co.re/) 모듈에서 사용하는 표준 인터페이스입니다.
`tuple val(meta), path(file)` 패턴은 nf-core 모듈 라이브러리 전반에 걸쳐 일관되게 사용되므로, 이 관례를 채택한 워크플로우는 최소한의 수정으로 nf-core 모듈을 교체하여 사용할 수 있습니다.

### 핵심 정리

이 섹션에서 다음 내용을 학습했습니다.

- **데이터시트 읽는 방법:** `splitCsv`를 사용하여 헤더 정보가 있는 CSV 파일을 파싱합니다.
- **메타 맵 관례가 존재하는 이유:** 메타데이터와 데이터 파일을 `[meta, file]` 튜플로 분리하면 데이터시트가 변경되어도 채널 구조가 안정적으로 유지됩니다.
- **프로세스 내부에서 메타 맵 필드를 사용하는 방법:** 메타 맵의 모든 필드는 스크립트 블록에서 점 표기법으로 접근할 수 있습니다.

---

## 2. 추가적인 메타데이터 조작

메타 맵 인터페이스가 갖춰졌으니, 데이터가 파이프라인을 흐르면서 메타 맵을 더욱 풍부하게 만들 수 있습니다.

[`langid`](https://github.com/saffsd/langid.py)라는 도구를 사용하여 각 녹음 파일의 언어를 식별합니다.
텍스트 조각이 주어지면 언어 예측과 확률 점수를 `stdout`으로 출력합니다.

### 2.1. 언어 식별 단계 추가

`langid` 도구를 적용하는 `IDENTIFY_LANGUAGE`라는 사전 작성된 프로세스 모듈을 제공합니다.

모듈 파일을 열어 코드를 확인합니다.

```groovy title="modules/langid.nf" linenums="1" hl_lines="7 10"
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

입력 정의가 섹션 1에서 구축한 것과 동일한 `tuple val(meta), path(file)` 구조를 사용하므로, `ch_datasheet`를 별도의 수정 없이 이 프로세스에 바로 전달할 수 있습니다.

출력은 세 번째 요소로 `stdout`을 추가합니다. 이는 `langid`가 콘솔에 출력하는 언어 예측을 캡처합니다.
`sed` 명령은 확률 점수와 후행 개행 문자를 제거하여 두 글자 언어 코드만 남깁니다.

#### 2.1.1. `IDENTIFY_LANGUAGE` 실행 추가

`IDENTIFY_LANGUAGE` 프로세스 모듈을 포함하고 데이터시트 채널에서 실행합니다.

=== "후"

    ```groovy title="main.nf" linenums="1" hl_lines="4 14-16"
    #!/usr/bin/env nextflow

    include { COWPY } from './modules/cowpy.nf'
    include { IDENTIFY_LANGUAGE } from './modules/langid.nf'

    workflow {
        main:
        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .map { row ->
                [[id: row.id, character: row.character], row.recording]
            }

        // langid를 실행하여 각 인사말의 언어를 식별합니다
        IDENTIFY_LANGUAGE(ch_datasheet)
        IDENTIFY_LANGUAGE.out.view()

        COWPY(ch_datasheet)
    ```

=== "전"

    ```groovy title="main.nf" linenums="1"
    #!/usr/bin/env nextflow

    include { COWPY } from './modules/cowpy.nf'

    workflow {
        main:
        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .map { row ->
                [[id: row.id, character: row.character], row.recording]
            }

        COWPY(ch_datasheet)
    ```

이 프로세스의 주요 출력은 문자열이므로 게시할 출력 파일이 없습니다.
대신 `IDENTIFY_LANGUAGE.out.view()`를 사용하여 작업 결과를 확인합니다.

#### 2.1.2. 워크플로우 실행

`-resume`을 사용하여 `COWPY` 작업을 다시 실행하지 않고 언어 식별을 실행합니다.

```bash
nextflow run main.nf -resume
```

??? success "명령 출력"

    ```console
     N E X T F L O W   ~  version 25.10.4

    Launching `main.nf` [voluminous_mcnulty] DSL2 - revision: f9bcfebabb

    executor >  local (14)
    [5d/dffd4e] COWPY (7)             [100%] 7 of 7, cached: 7 ✔
    [4e/f722fe] IDENTIFY_LANGUAGE (7) [100%] 7 of 7 ✔
    [[id:sampleA, character:squirrel], /workspaces/training/side-quests/metadata/work/eb/f7148ebdd898fbe1136bec6a714acb/bonjour.txt, fr]
    [[id:sampleB, character:tux], /workspaces/training/side-quests/metadata/work/16/71d72410952c22cd0086d9bca03680/guten_tag.txt, de]
    [[id:sampleD, character:turkey], /workspaces/training/side-quests/metadata/work/c4/b7562adddc1cc0b7d414ec45d436eb/hello.txt, en]
    [[id:sampleC, character:sheep], /workspaces/training/side-quests/metadata/work/ea/04f5d979429e4455e14b9242fb3b45/hallo.txt, de]
    [[id:sampleF, character:moose], /workspaces/training/side-quests/metadata/work/5a/6c2b84bf8fadb98e28e216426be079/salut.txt, fr]
    [[id:sampleE, character:stegosaurus], /workspaces/training/side-quests/metadata/work/af/ee7c69bcab891c40d0529305f6b9e7/hola.txt, es]
    [[id:sampleG, character:turtle], /workspaces/training/side-quests/metadata/work/4e/f722fe47271ba7ebcd69afa42964ca/ciao.txt, it]
    ```

데이터셋의 각 파일에 대한 언어 예측이 생성되었습니다.

출력 튜플이 `[meta, file, lang_id]`로 구성되어 있어, 메타 맵과 파일이 새로운 결과와 함께 전달되는 것을 확인할 수 있습니다.

!!! note "참고"

    메타 맵을 결과와 연결된 상태로 유지하는 이 패턴은 나중에 채널 간 결과를 연결하기 쉽게 만듭니다.
    채널의 항목 순서에 의존하여 데이터를 올바르게 연결할 수 없습니다.
    대신 키를 사용해야 합니다.
    메타 맵은 이 목적에 이상적인 구조를 제공합니다.

    이 사용 사례는 [Splitting & Grouping](../splitting_and_grouping/index.md) 사이드 퀘스트에서 자세히 살펴봅니다.

### 2.2. 프로세스 출력으로 메타데이터 보강

언어 예측은 파일 내용에 대한 메타데이터의 한 형태입니다.
별도의 요소로 유지하는 대신, 메타 맵에 다시 통합합니다.

#### 2.2.1. 새롭고 확장된 메타 맵 생성

Groovy `+` 연산자를 사용하여 원래 메타 맵을 대체하는 새로운 메타 맵을 생성합니다.

=== "후"

    ```groovy title="main.nf" linenums="14" hl_lines="3-7"
        // langid를 실행하여 각 인사말의 언어를 식별합니다
        IDENTIFY_LANGUAGE(ch_datasheet)
        IDENTIFY_LANGUAGE.out
            .map { meta, file, lang_id ->
                [meta + [lang: lang_id], file]
            }
            .view()
    ```

=== "전"

    ```groovy title="main.nf" linenums="14" hl_lines="3"
        // langid를 실행하여 각 인사말의 언어를 식별합니다
        IDENTIFY_LANGUAGE(ch_datasheet)
        IDENTIFY_LANGUAGE.out.view()
    ```

이 작업의 핵심은 `#!groovy meta + [lang: lang_id]`입니다.

이 코드는 언어 코드를 담은 단일 키-값 쌍의 임시 맵(`[lang: lang_id]`)을 생성한 다음, Groovy `+` 연산자를 사용하여 기존 메타데이터가 담긴 원래 `meta` 맵과 결합하여 새롭고 확장된 메타 맵을 생성합니다.

더 자세한 설명은 아래 박스를 참조하세요.

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
    new_map = map1 + [lang: lang_id]
    ```

    여기서 `[lang: lang_id]`는 즉석에서 새로운 이름 없는 맵을 생성하고, `map1 + `는 `map1`을 새로운 이름 없는 맵과 병합하여 이전과 동일한 `new_map` 내용을 생성합니다.

    멋지지 않나요?

    **이제 Nextflow `channel.map()` 작업의 맥락으로 적용해 봅니다.**

    코드는 다음과 같습니다.

    ```groovy
    .map { map1, lang_id ->
        map1 + [lang: lang_id]
    }
    ```

    이 코드는 다음을 수행합니다.

    - `#!groovy map1, lang_id ->`는 튜플의 두 항목을 가져옵니다.
    - `#!groovy map1 + [lang: lang_id]`는 위에서 설명한 대로 새 맵을 생성합니다.

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

    `map` 작업에서 `file`이 이동하는 것처럼 보이는 이유가 이해하기 어렵다면, `#!groovy [meta + [lang: lang_id], file]` 대신 `[new_map, file]`로 읽어보세요.
    이렇게 하면 `file`을 튜플의 두 번째 위치에 그대로 두고, `new_info` 값을 첫 번째 위치의 맵에 통합했다는 것이 더 명확해집니다.

    **이것이 바로 `tuple val(meta), path(file)` 채널 구조로 돌아오는 것입니다!**

#### 2.2.2. 워크플로우 실행

코드가 무엇을 하는지 이해했다면, 워크플로우를 실행하여 작동하는지 확인합니다.

```bash
nextflow run main.nf -resume
```

??? success "명령 출력"

    ```console
     N E X T F L O W   ~  version 25.10.4

    Launching `main.nf` [cheeky_fermat] DSL2 - revision: d096281ee4

    [5d/dffd4e] COWPY (7)             [100%] 7 of 7, cached: 7 ✔
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

!!! tip "메타 맵에서 키 제거"

    Groovy [`subMap`](https://docs.groovy-lang.org/latest/html/groovy-jdk/java/util/Map.html#subMap(java.util.Collection)) 메서드를 사용하여 메타 맵에서 키를 제거할 수 있습니다. 이 메서드는 지정한 키만 포함하는 새 맵을 반환합니다.

    ```groovy
    meta.subMap(['id', 'character'])  // 'id'와 'character'만 포함하는 맵 반환
    ```

    이는 다운스트림 프로세스나 모듈이 메타 맵에 누적된 모든 필드를 필요로 하지 않을 때 유용합니다.

### 2.3. 조건문을 사용하여 언어 그룹 할당

메타 맵에 언어 예측이 포함되었으니, 이를 바탕으로 추가 메타데이터를 도출할 수 있습니다.
데이터셋의 언어는 게르만어(영어, 독일어)와 로망스어(프랑스어, 스페인어, 이탈리아어) 두 계열로 나뉩니다.
`lang_group` 필드를 추가하면 이 분류를 다운스트림에서 쉽게 활용할 수 있습니다.

#### 2.3.1. 조건 로직이 포함된 `map` 작업 추가

두 번째 `map` 작업과 조건 로직을 사용하여 언어 계열을 할당합니다.

```groovy
.map { meta, file ->

    // lang_group을 정의하는 조건 로직이 여기에 들어갑니다

    [meta + [lang_group: lang_group], file]
}
```

적용할 로직은 다음과 같습니다.

- 기본값으로 `lang_group = 'unknown'`을 설정합니다.
- `meta.lang`이 `'de'` 또는 `'en'`이면 `lang_group`을 `'germanic'`으로 설정합니다.
- 그렇지 않고 `meta.lang`이 `['fr', 'es', 'it']`에 포함되면 `lang_group`을 `'romance'`로 설정합니다.

!!! tip "팁"

    map 작업 내에서 `meta.lang`으로 `lang` 값에 접근할 수 있습니다.

워크플로우에 다음과 같은 변경을 합니다.

=== "후"

    ```groovy title="main.nf" linenums="14" hl_lines="7-19 21"
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

        ch_languages.view()
    ```

=== "전"

    ```groovy title="main.nf" linenums="14" hl_lines="7"
        // langid를 실행하여 각 인사말의 언어를 식별합니다
        IDENTIFY_LANGUAGE(ch_datasheet)
        IDENTIFY_LANGUAGE.out
            .map { meta, file, lang_id ->
                [meta + [lang: lang_id], file]
            }
            .view()
    ```

주요 사항은 다음과 같습니다.

- `def lang_group = "unknown"`으로 안전한 기본값을 가진 변수를 초기화합니다.
- `if / else if` 구조로 두 언어 계열을 처리하며, 그 외의 경우는 `'unknown'`으로 유지됩니다.
- `#!groovy .set { ch_languages }`는 결과 채널에 이름을 부여하여 다음 단계에서 사용할 수 있게 합니다.

<!-- TODO (future) Add note/links to relevant docs in additional resources section -->

#### 2.3.2. 워크플로우 실행

워크플로우를 실행하여 올바르게 작동하는지 확인합니다.

```bash
nextflow run main.nf -resume
```

??? success "명령 출력"

    ```console
     N E X T F L O W   ~  version 25.10.4

    Launching `main.nf` [wise_almeida] DSL2 - revision: 46778c3cd0

    [5d/dffd4e] COWPY (7)             [100%] 7 of 7, cached: 7 ✔
    [da/652cc6] IDENTIFY_LANGUAGE (7) [100%] 7 of 7, cached: 7 ✔
    [[id:sampleA, character:squirrel, lang:fr, lang_group:romance], /workspaces/training/side-quests/metadata/data/bonjour.txt]
    [[id:sampleB, character:tux, lang:de, lang_group:germanic], /workspaces/training/side-quests/metadata/data/guten_tag.txt]
    [[id:sampleC, character:sheep, lang:de, lang_group:germanic], /workspaces/training/side-quests/metadata/data/hallo.txt]
    [[id:sampleD, character:turkey, lang:en, lang_group:germanic], /workspaces/training/side-quests/metadata/data/hello.txt]
    [[id:sampleE, character:stegosaurus, lang:es, lang_group:romance], /workspaces/training/side-quests/metadata/data/hola.txt]
    [[id:sampleF, character:moose, lang:fr, lang_group:romance], /workspaces/training/side-quests/metadata/data/salut.txt]
    [[id:sampleG, character:turtle, lang:it, lang_group:romance], /workspaces/training/side-quests/metadata/data/ciao.txt]
    ```

메타 맵에 이제 `id`, `character`, `lang`, `lang_group` 네 개의 필드가 포함되었습니다.
채널 구조는 여전히 `[meta, file]`입니다.

### 2.4. 메타데이터를 사용하여 출력 이름 지정 및 구성

메타 맵에 `lang`과 `lang_group`이 포함되었으니, 이를 사용하여 출력 파일 이름에 언어 코드를 추가하고 언어 계열별로 서브디렉토리에 구성합니다.

세 가지 변경이 필요합니다. `COWPY` 프로세스를 업데이트하여 출력 이름을 변경하고 `meta`를 내보내도록 하고, `COWPY` 실행을 `ch_languages`에서 실행하도록 업데이트하고, 서브디렉토리 경로를 지정하도록 output 블록을 업데이트합니다.

#### 2.4.1. `COWPY` 프로세스 업데이트

메타 맵의 언어 코드를 사용하여 출력 파일 이름을 변경하고, output 블록이 서브디렉토리 라우팅을 위해 `lang_group`에 접근할 수 있도록 출력에 `meta`를 추가합니다.

=== "후"

    ```groovy title="modules/cowpy.nf" linenums="9" hl_lines="2 6"
        output:
        tuple val(meta), path("${meta.lang}-${input_file}")

        script:
        """
        cat ${input_file} | cowpy -c ${meta.character} > ${meta.lang}-${input_file}
        """
    ```

=== "전"

    ```groovy title="modules/cowpy.nf" linenums="9" hl_lines="2 6"
        output:
        path "cowpy-${input_file}"

        script:
        """
        cat ${input_file} | cowpy -c ${meta.character} > cowpy-${input_file}
        """
    ```

이는 입력 정의를 전혀 수정하지 않고도 다른 메타데이터 필드를 활용하여 프로세스 동작을 맞춤화하는 방법을 보여줍니다.

#### 2.4.2. `ch_languages`에서 `COWPY` 실행하도록 업데이트

`COWPY(ch_datasheet)`를 `COWPY(ch_languages)`로 교체합니다.

=== "후"

    ```groovy title="main.nf" linenums="32" hl_lines="3"
        .set { ch_languages }

        COWPY(ch_languages)

        publish:
        cowpy_art = COWPY.out
    }
    ```

=== "전"

    ```groovy title="main.nf" linenums="32" hl_lines="3 5"
        .set { ch_languages }

        ch_languages.view()

        COWPY(ch_datasheet)

        publish:
        cowpy_art = COWPY.out
    }
    ```

채널 내용을 더 이상 확인할 필요가 없으므로 `ch_languages.view()` 줄도 제거합니다.

#### 2.4.3. output 블록 업데이트

`output {}` 블록에 `path` closure를 추가하여 각 파일을 언어 그룹 서브디렉토리로 라우팅합니다.

=== "후"

    ```groovy title="main.nf" linenums="40" hl_lines="3"
    output {
        cowpy_art {
            path { meta, file -> meta.lang_group }
        }
    }
    ```

=== "전"

    ```groovy title="main.nf" linenums="40" hl_lines="2 3"
    output {
        cowpy_art {
        }
    }
    ```

이는 메타데이터를 사용하여 출력을 유연하게 구성하는 방법을 보여줍니다.

#### 2.4.4. 전체 파이프라인 실행

이전 결과를 삭제하고 전체 파이프라인을 실행합니다.

```bash
rm -r results
nextflow run main.nf
```

??? success "명령 출력"

    ```console
     N E X T F L O W   ~  version 25.10.4

    Launching `main.nf` [suspicious_crick] DSL2 - revision: 25541014c5

    executor >  local (14)
    [5d/dffd4e] IDENTIFY_LANGUAGE (7) [100%] 7 of 7 ✔
    [e7/317c18] COWPY (7)             [100%] 7 of 7 ✔
    ```

results 디렉토리가 이제 언어 계열별로 구성되었으며, 각 파일은 감지된 언어에 따라 이름이 지정되었습니다.

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

`output {}` 블록의 `path` closure는 각 `[meta, file]` 튜플을 받아 `meta.lang_group`을 서브디렉토리 이름으로 반환합니다.
파일 이름 자체는 프로세스가 출력하는 것(`#!groovy "${meta.lang}-${input_file}"`)에서 옵니다.
두 메타데이터(언어 코드와 언어 그룹) 모두 이 섹션에서 구축한 풍부해진 메타 맵에서 옵니다.

### 핵심 정리

이 섹션에서 다음 내용을 학습했습니다.

- **프로세스 출력으로 메타 맵을 보강하는 방법:** `#!groovy meta + [key: value]`로 새 키를 추가하면 `[meta, file]` 채널 구조를 유지하면서 메타데이터를 풍부하게 만들 수 있습니다.
- **메타데이터에서 메타데이터를 도출하는 방법:** `map` 작업 내의 조건 로직으로 기존 필드에서 새 필드를 계산할 수 있습니다.
- **출력 구성에 메타데이터를 사용하는 방법:** `output {}` 블록의 `path` closure가 메타 맵을 읽어 파일을 서브디렉토리로 라우팅할 수 있습니다.

---

## 3. 견고성 고려사항

메타데이터 값이 프로세스 동작을 이끌 때, 누락되거나 불완전한 데이터는 진단하기 어려운 문제를 일으킬 수 있습니다.
예상되는 상황과 처리 방법을 살펴봅니다.

### 3.1. 필수 메타데이터 필드가 누락된 경우

`character` 값은 `COWPY` 프로세스가 유효한 결과를 생성하기 위해 필요합니다.
오류 발생 방식은 열이 데이터시트에 존재하지만 비어 있는 경우와 완전히 없는 경우에 따라 다릅니다.

#### 3.1.1. 열은 존재하지만 값이 비어 있는 경우

데이터시트의 한 항목에서 `character` 필드가 비어 있다고 가정합니다.

```csv title="datasheet.csv" linenums="1" hl_lines="2"
id,character,recording
sampleA,,/workspaces/training/side-quests/metadata/data/bonjour.txt
sampleB,tux,/workspaces/training/side-quests/metadata/data/guten_tag.txt
...
```

데이터시트를 파싱할 때 모든 항목에 대해 `character` 키가 생성되지만, `sampleA`의 `meta.character`는 빈 문자열이 됩니다.
Nextflow가 명령에 `#!groovy ${meta.character}`를 대입할 때 `COWPY` 도구는 `-c`에 빈 인자를 받아 오류를 발생시킵니다.

??? failure "명령 출력"

    ```console hl_lines="8 11 16 28"
     N E X T F L O W   ~  version 25.10.4

    Launching `main.nf` [marvelous_hirsch] DSL2 - revision: 0dfeee3cc1

    executor >  local (9)
    [c1/c5dd4f] process > IDENTIFY_LANGUAGE (7) [ 85%] 6 of 7
    [d3/b7c415] process > COWPY (2)             [  0%] 0 of 6
    ERROR ~ Error executing process > 'COWPY (1)'

    Caused by:
      Process `COWPY (1)` terminated with an error exit status (2)


    Command executed:

      cat bonjour.txt | cowpy -c  > fr-bonjour.txt

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

오류 메시지(`expected one argument`)는 빈 `-c` 플래그를 가리킵니다.
작업 디렉토리의 `.command.sh` 파일을 확인하면 명령이 빈 값으로 실행되었음을 확인할 수 있습니다.

#### 3.1.2. 데이터시트에 열이 없는 경우

`character` 열이 완전히 없는 경우를 가정합니다.

```csv title="datasheet.csv" linenums="1"
id,recording
sampleA,/workspaces/training/side-quests/metadata/data/bonjour.txt
sampleB,/workspaces/training/side-quests/metadata/data/guten_tag.txt
...
```

메타 맵에 `character` 키가 전혀 생성되지 않습니다.
프로세스 스크립트가 `#!groovy ${meta.character}`를 평가할 때 누락된 키는 `null`을 반환하며, Nextflow는 명령에 문자열 `null`을 그대로 대입합니다.

??? failure "명령 출력"

    ```console hl_lines="8 11 16"
     N E X T F L O W   ~  version 25.10.4

    Launching `main.nf` [jovial_bohr] DSL2 - revision: eaaf375827

    executor >  local (9)
    [0d/ada9db] process > IDENTIFY_LANGUAGE (5) [ 85%] 6 of 7
    [06/28065f] process > COWPY (2)             [  0%] 0 of 6
    ERROR ~ Error executing process > 'COWPY (2)'

    Caused by:
      Process `COWPY (2)` terminated with an error exit status (1)


    Command executed:

      cat guten_tag.txt | cowpy -c null > de-guten_tag.txt

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

실행된 명령의 `cowpy -c null`이 진단 단서입니다.

### 3.2. 누락된 메타데이터 처리 전략

워크플로우를 누락된 메타데이터에 더 견고하게 만드는 두 가지 보완적인 접근 방식이 있습니다.

**1. 입력 유효성 검사**

가장 신뢰할 수 있는 해결책은 처리가 시작되기 전에 데이터시트를 검증하는 것입니다. 이를 통해 실행 중간에 알 수 없는 프로세스 오류로 나타나는 대신, 명확한 오류 메시지와 함께 문제를 조기에 발견할 수 있습니다.
[Hello nf-core](../../hello_nf-core/05_input_validation.md) 교육 과정에서 nf-schema 플러그인을 사용하여 입력 유효성 검사를 추가하는 방법을 다룹니다. <!-- TODO (future) pending a proper Validation side quest -->

**2. 필수 값에 대한 명시적 프로세스 입력**

특정 값이 필수임을 프로세스 인터페이스 자체에서 전달하려면, 메타 맵에서 해당 값을 명시적 입력으로 추출하는 것을 고려하세요.

=== "프로세스 정의"

    ```groovy title="modules/cowpy.nf" linenums="6"
    input:
    tuple val(meta), val(character), path(input_file)
    ```

=== "워크플로우 실행"

    ```groovy title="main.nf"
    COWPY(ch_languages.map { meta, file -> [meta, meta.character, file] })
    ```

이 방식은 `character`를 프로세스 계약의 가시적이고 필수적인 부분으로 만듭니다.
모듈을 읽는 누구나 캐릭터 값이 반드시 제공되어야 한다는 것을 즉시 알 수 있습니다.
필드가 없으면 프로세스가 실행되기 전에 채널 수준에서 명확하게 오류가 발생합니다.

이는 유용한 설계 원칙을 강조합니다.

**메타 맵은 선택적이고 설명적인 정보에 사용하되, 필수 값은 명시적 입력으로 추출합니다.**

메타 맵은 채널 구조를 깔끔하고 안정적으로 유지하는 데 탁월하지만, 프로세스에서 실제로 필요한 값의 경우 명명된 입력으로 표면화하면 명확성이 향상되고 다른 맥락에서 모듈을 올바르게 사용하기 쉬워집니다.

### 핵심 정리

이 섹션에서 다음 내용을 확인했습니다.

- **누락된 메타데이터가 나타나는 방식:** 빈 필드는 빈 인자를 생성하고, 없는 필드는 명령에 `null`이 그대로 대입됩니다.
- **두 가지 보완적 전략:** 문제를 조기에 발견하기 위한 입력 유효성 검사와, 요구사항을 명확하게 전달하기 위한 명시적 프로세스 입력.

---

## 요약

이 사이드 퀘스트에서는 Nextflow 워크플로우에서 메타데이터를 효과적으로 다루는 방법을 살펴보았습니다.

"메타 맵 + 데이터 파일" 튜플 패턴은 Nextflow의 핵심 관례로, 메타데이터를 개별 값으로 전달하는 것보다 여러 가지 장점을 제공합니다.

- 데이터시트가 변경되어도 채널 구조가 안정적으로 유지됩니다.
- 필드 이름을 하드코딩하지 않고 샘플별로 프로세스 동작을 맞춤화할 수 있습니다.
- 출력의 이름 지정, 그룹화, 구성을 위해 파이프라인 전반에 걸쳐 메타데이터를 활용할 수 있습니다.
- 이 인터페이스로 작성된 모듈은 nf-core 모듈을 포함하여 상호 교환 가능합니다.

### 핵심 패턴

1.  **메타데이터 읽기 및 구조화:** CSV 데이터시트를 파싱하고 메타 맵을 생성합니다.

    ```groovy
    channel.fromPath('datasheet.csv')
        .splitCsv(header: true)
        .map { row ->
            [ [id: row.id, character: row.character], row.recording ]
        }
    ```

2.  **워크플로우 중 메타데이터 확장:** 프로세스 출력이나 도출된 로직으로 새 키를 추가합니다.

    ```groovy
    // 프로세스 출력에서
    .map { meta, file, lang ->
        [ meta + [lang: lang], file ]
    }

    // 조건 로직에서
    .map { meta, file ->
        def lang_group = "unknown"
        if (meta.lang in ["de", "en"]) { lang_group = "germanic" }
        else if (meta.lang in ["fr", "es", "it"]) { lang_group = "romance" }
        [ meta + [lang_group: lang_group], file ]
    }
    ```

3.  **프로세스 내부에서 메타데이터 사용:** 스크립트 블록에서 점 표기법으로 모든 필드에 접근합니다.

    ```groovy
    cat ${input_file} | cowpy -c ${meta.character} > ${meta.lang}-${input_file}
    ```

4.  **메타데이터 값으로 출력 구성:** `output {}` 블록의 `path` closure를 사용합니다.

    ```groovy
    output {
        cowpy_art {
            path { meta, file -> meta.lang_group }
        }
    }
    ```

### 추가 자료

- [map 연산자](https://www.nextflow.io/docs/latest/operator.html#map)
- [multiMap 연산자](https://www.nextflow.io/docs/latest/reference/operator.html#multimap)
- [stdout 출력 한정자](https://www.nextflow.io/docs/latest/process.html#outputs)

---

## 다음 단계

[사이드 퀘스트 메뉴](../index.md)로 돌아가거나 페이지 오른쪽 하단의 버튼을 클릭하여 목록의 다음 주제로 이동하세요.
