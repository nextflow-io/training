# 파트 1: Nextflow 실행

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } AI 지원 번역 - [자세히 알아보기 및 개선 제안](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

이 파트에서는 Nextflow 파이프라인 실행의 핵심 개념을 학습합니다.
간단한 Hello World 워크플로우부터 시작하여, 컨테이너를 사용해 여러 입력을 병렬로 처리하는 완전한 다단계 파이프라인까지 단계적으로 진행합니다.

---

## 1. Hello World

워크플로우 `1-hello.nf`는 명령줄 인자를 통해 인사말을 받아 파일에 작성합니다.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello_world.svg"
</figure>

### 1.1. 워크플로우 실행

터미널에서 다음 명령을 실행합니다.

```bash
nextflow run 1-hello.nf --input 'Hello World!'
```

??? success "명령 출력"

    ```console hl_lines="6"
    N E X T F L O W   ~  version 26.04.4

    Launching `1-hello.nf` [infallible_volhard] revision: 82d40dbf94

    executor >  local (1)
    [6d/740edd] sayHello | 1 of 1 ✔

    Outputs:

      /workspaces/training/nextflow-run/results

      first_output: 1-hello/Hello World!-output.txt
    ```

출력에서 핵심이 되는 줄은 프로세스 상태 줄입니다:

```console
[6d/740edd] sayHello | 1 of 1 ✔
```

이 줄은 `sayHello` 프로세스가 성공적으로 한 번 실행되었음을 나타냅니다.
`[6d/740edd]` 접두사는 작업의 work 디렉토리 경로를 축약한 것으로, 자세한 내용은 아래에서 다룹니다.
이어지는 `Outputs:` 블록에는 파이프라인이 게시한 모든 파일이 나열되며, 아래 [1.4](#14-optional-code-walkthrough)에서 다루는 `output` 블록에 따라 레이블이 지정됩니다.

### 1.2. 출력 파일 확인

이 워크플로우는 출력을 `results` 디렉토리에 게시하도록 설정되어 있습니다.
실행 후 해당 디렉토리에서 출력 파일을 확인할 수 있습니다:

```console title="results/"
results
└── 1-hello
    └── Hello World!-output.txt
```

파일을 열어 `Hello World!`가 포함되어 있는지 확인합니다.

### 1.3. `work/` 디렉토리 살펴보기

Nextflow는 내부적으로 `work/`라는 디렉토리 안에 모든 process 실행마다 고유한 작업 디렉토리를 생성합니다.
콘솔 출력에 표시된 해시값(`[6d/740edd]`)이 해당 디렉토리의 경로입니다.

```bash
ls work/6d/740edd*
```

디렉토리 안에는 출력 파일과 함께 여러 숨김 로그 파일이 있습니다:

- **`.command.sh`**: Nextflow가 실행한 정확한 명령
- **`.command.out`** / **`.command.err`**: 프로세스의 표준 출력 (stdout) 및 표준 오류 (stderr)
- **`.command.log`**: 통합 로그 출력
- **`.exitcode`**: 프로세스 종료 코드

`.command.sh` 파일은 디버깅 시 특히 유용합니다. 실제로 실행된 내용을 정확하게 보여줍니다.

### 1.4. 선택 사항: 코드 살펴보기

파이프라인 실행만이 목적이라면 코드를 이해하는 것이 필수는 아니지만, 궁금하다면 살펴볼 가치가 있습니다.

??? optional "클릭하여 이 연습과 관련된 코드 살펴보기"

    `1-hello.nf`를 열어 주요 구성 요소를 확인합니다.

    ```groovy title="1-hello.nf" linenums="1"
    #!/usr/bin/env nextflow

    include { sayHello } from './modules/sayHello.nf'

    /*
     * 파이프라인 매개변수
     */
    params {
        input: String
    }

    workflow {

        main:
        // 인사말을 내보냅니다
        sayHello(params.input)

        publish:
        first_output = sayHello.out
    }

    output {
        first_output {
            path '1-hello'
            mode 'copy'
        }
    }
    ```

    다음 구성 요소를 확인할 수 있습니다:

    - `process` 모듈을 가리키는 `include` 구문
    - 파이프라인 매개변수를 정의하는 `params` 블록
    - 수행할 작업을 기술하는 `workflow` 블록
    - 출력 처리 방법을 기술하는 `output` 블록

    각 구성 요소를 순서대로 살펴봅니다.

    ### `process` 모듈

    `include` 구문은 Nextflow에게 별도의 코드 파일에서 `sayHello`라는 항목을 불러오도록 지시합니다.

    ```groovy title="1-hello.nf" linenums="3"
    include { sayHello } from './modules/sayHello.nf'
    ```

    해당 파일에는 `sayHello`라는 프로세스의 정의가 있습니다:

    ```groovy title="modules/sayHello.nf" linenums="4"
    process sayHello {

        input:
        val greeting

        output:
        path "${greeting}-output.txt"

        script:
        """
        echo '${greeting}' > '${greeting}-output.txt'
        """
    }
    ```

    **process**는 파이프라인의 단일 단계를 정의합니다.
    입력, 출력, 실행할 스크립트를 선언합니다.
    `val` 한정자는 입력이 일반 값(문자열, 숫자 등)임을 의미합니다.
    `path` 한정자는 출력이 파일 경로임을 의미합니다.

    프로세스 정의를 메인 워크플로우 파일에 작성할 수도 있지만, 별도의 모듈 파일로 분리하면 재사용성이 높아집니다. 동일한 모듈을 여러 워크플로우 스크립트에서 가져올 수 있습니다.

    ### `params` 블록

    `params` 블록은 워크플로우가 허용하는 명령줄 매개변수를 선언합니다:

    ```groovy title="1-hello.nf" linenums="8"
    params {
        input: String
    }
    ```

    여기에 선언된 매개변수는 이중 대시(`--input`)를 사용하여 명령줄에서 사용할 수 있습니다.
    지원되는 타입으로는 `String`, `Integer`, `Float`, `Boolean`, `Path`가 있습니다.

    !!! tip "팁"

        워크플로우 매개변수는 항상 이중 대시(`--input`)를 사용합니다. 이는 단일 대시를 사용하는 Nextflow 자체 CLI 플래그(예: `-resume`)와 구분하기 위함입니다.

    ### `workflow` 블록

    **workflow** 블록은 데이터 흐름 로직, 즉 어떤 프로세스를 어떤 순서로 실행할지를 정의합니다.

    ```groovy title="1-hello.nf" linenums="12"
    workflow {

        main:
        // 인사말을 내보냅니다
        sayHello(params.input)

        publish:
        first_output = sayHello.out
    }
    ```

    여기서는 하나의 프로세스만 실행되므로 매우 단순합니다. 보다 현실적인 예제는 이후에 다룹니다.

    `main:` 섹션은 `--input` 값을 사용하여 `sayHello` 프로세스를 호출합니다.
    `publish:` 섹션은 results 디렉토리에 복사할 출력을 나열합니다.

    ### `output` 블록

    파일 하단의 `output` 블록은 대상 경로와 복사 모드를 지정합니다.

    ```groovy title="1-hello.nf" linenums="22"
    output {
        first_output {
            path '1-hello'
            mode 'copy'
        }
    }
    ```

    각 이름이 지정된 항목은 워크플로우의 `publish:` 레이블에 대응하며, `results/` 아래의 하위 디렉토리에 매핑됩니다.

### 핵심 정리

Nextflow 파이프라인을 실행하고 출력을 찾는 방법을 알게 되었으며, 작업이 `work/` 아래의 작업 디렉토리에서 실행된다는 것을 이해했습니다.

### 다음 단계

Nextflow가 여러 입력을 효율적으로 처리하는 방법을 확인합니다.

---

## 2. 여러 입력 처리

실제 파이프라인은 일반적으로 하나가 아닌 많은 데이터를 처리합니다.
워크플로우 `2-inputs.nf`는 CSV 파일을 읽어 각 행마다 `sayHello`를 병렬로 실행합니다.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-pipeline-multi-inputs-csv.svg"
</figure>

먼저 워크플로우를 실행한 후, Nextflow가 여러 입력을 처리하는 데 사용하는 메커니즘을 살펴봅니다.

### 2.1. 워크플로우 실행

터미널에서 다음 명령을 실행합니다.

```bash
nextflow run 2-inputs.nf --input data/greetings.csv
```

??? success "명령 출력"

    ```console hl_lines="6"
    N E X T F L O W   ~  version 26.04.4

    Launching `2-inputs.nf` [nauseous_babbage] revision: b90778224d

    executor >  local (3)
    [66/de7844] sayHello (3) | 3 of 3 ✔

    Outputs:

      /workspaces/training/nextflow-run/results

      first_output:
        - 2-inputs/Bonjour-output.txt
        - 2-inputs/Hello-output.txt
        - 2-inputs/Hola-output.txt
    ```

`3 of 3`은 CSV의 각 행마다 한 번씩, `sayHello` 프로세스가 세 번 호출되었음을 나타냅니다.

`results` 디렉토리에는 인사말마다 하나씩, 총 세 개의 출력 파일이 생성됩니다:

```console title="results/2-inputs/"
2-inputs
├── Bonjour-output.txt
├── Hello-output.txt
└── Hola-output.txt
```

출력 파일 중 하나를 열어 각 파일에 인사말이 포함되어 있는지 확인합니다.

위의 요약 출력은 `sayHello`에 대한 단일 요약 줄을 보여주지만, Nextflow는 실제로 CSV의 각 행마다 하나씩 세 개의 별도 작업 실행을 시작하고, 시스템 리소스가 허용하는 즉시 병렬로 실행합니다.

[1.3](#13-explore-the-work-directory)에서 살펴본 단일 작업과 마찬가지로, 이 세 번의 실행 각각은 `work/` 아래에 완전히 격리된 자체 작업 디렉토리를 갖습니다:

```console title="work/"
work
├── 2d/276c63.../
│   ├── .command.sh
│   └── Hola-output.txt
├── ab/007682.../
│   ├── .command.sh
│   └── Bonjour-output.txt
└── d9/2476082.../
    ├── .command.sh
    └── Hello-output.txt
```

각 `.command.sh`에는 해당 인사말에 대한 명령만 포함됩니다:

```console title="work/2d/276c63.../.command.sh"
#!/bin/bash -ue
echo 'Hola' > 'Hola-output.txt'
```

이러한 격리 덕분에 병렬 실행이 안전합니다. 동시에 실행되는 세 작업은 작업 디렉토리를 공유하지 않으므로, 같은 이름의 파일을 생성하더라도 한 작업이 다른 작업의 파일을 덮어쓰거나 충돌하는 일이 발생하지 않습니다.
또한 이것이 `-resume`(다음에 다룸)이 개별 작업을 독립적으로 캐시하고 재사용할 수 있는 이유이기도 합니다. 각 작업의 입력, 출력, 로그는 모두 자체 디렉토리 안에 있으며, 작업 간에 공유되는 것이 없어 동기화 문제가 발생하지 않습니다.

### 2.2. `-ansi-log false` 옵션으로 워크플로우 재실행

기본적으로 Nextflow는 출력을 프로세스당 단일 요약 줄로 압축합니다.
각 process 실행을 개별적으로 나열하려면 `-ansi-log false`를 추가합니다:

```bash
nextflow run 2-inputs.nf --input data/greetings.csv -ansi-log false
```

??? success "명령 출력"

    ```console
    N E X T F L O W  ~  version 26.04.4
    Launching `2-inputs.nf` [extravagant_bardeen] - revision: b90778224d
    [43/0bac1c] Submitted process > sayHello (1)
    [2d/99f604] Submitted process > sayHello (2)
    [6d/7578d7] Submitted process > sayHello (3)

    Outputs:

      /workspaces/training/nextflow-run/results

      first_output:
        - 2-inputs/Hello-output.txt
        - 2-inputs/Bonjour-output.txt
        - 2-inputs/Hola-output.txt
    ```

세 번의 process 실행과 각각 생성된 고유한 work 하위 디렉토리가 표시됩니다.

### 2.3. `-resume`으로 완료된 작업 건너뛰기

이제 두 개의 인사말이 추가된 확장 입력 파일로 전환하고, 명령줄에 `-resume`을 추가합니다:

```bash
nextflow run 2-inputs.nf --input data/greetings-extended.csv -resume
```

??? success "명령 출력"

    ```console hl_lines="6"
    N E X T F L O W   ~  version 26.04.4

    Launching `2-inputs.nf` [adoring_mayer] revision: b90778224d

    executor >  local (2)
    [84/2f3067] sayHello (5) | 5 of 5, cached: 3 ✔

    Outputs:

      /workspaces/training/nextflow-run/results

      first_output:
        - 2-inputs/Bonjour-output.txt
        - 2-inputs/Ciao-output.txt
        - 2-inputs/Hello-output.txt
        - 2-inputs/Hola-output.txt
        - 2-inputs/Ola-output.txt
    ```

Nextflow는 새로운 두 입력만 실행했습니다.
이전 실행에서 처리된 세 개의 인사말은 캐시되어 자동으로 재사용되었습니다.

이 기능은 다단계 파이프라인에서 이미 성공적으로 완료된 프로세스의 실행을 건너뛰는 데도 활용됩니다.
예를 들어, 시스템 오류로 파이프라인 실행이 중단되었거나, 개발 중인 파이프라인에 새로운 단계를 추가한 경우에 유용합니다.

`-resume` 기능은 실패 복구 시 중요한 시간과 리소스를 절약할 수 있어, 긴 파이프라인에서 특히 유용합니다.

### 2.4. 선택 사항: 코드 살펴보기

파이프라인 실행만이 목적이라면 코드를 이해하는 것이 필수는 아니지만, 궁금하다면 살펴볼 가치가 있습니다.

??? optional "클릭하여 이 연습과 관련된 코드 살펴보기"

    `2-inputs.nf`의 핵심 변경 사항은 워크플로우의 `main:` 섹션에 있습니다:

    ```groovy title="2-inputs.nf" linenums="14" hl_lines="3 4 5"
        main:
        // CSV 파일에서 입력을 위한 채널 생성
        greeting_ch = channel.fromPath(params.input)
                            .splitCsv()
                            .map { line -> line[0] }
        // 인사말을 내보냅니다
        sayHello(greeting_ch)
    ```

    여기서 볼 수 있는 것은 **channel**이라고 불리는 구조입니다. 이는 입력 데이터를 처리하여 작업을 쉽게 병렬화할 수 있게 해주는 큐 구조입니다.

    - `channel.fromPath(params.input)`은 `--input`으로 지정된 파일 경로에서 채널을 생성합니다
    - `.splitCsv()`는 CSV를 행으로 파싱합니다
    - `#!groovy .map { line -> line[0] }`은 각 행에서 첫 번째 열을 추출합니다

    결과는 `Hello`, `Bonjour`, `Hola`를 포함하는 채널입니다.
    `sayHello(greeting_ch)`에 전달되면, Nextflow는 항목마다 자동으로 프로세스를 한 번씩 호출하고, 리소스가 허용하는 경우 병렬로 실행합니다.

### 핵심 정리

CSV 파일에서 여러 입력을 병렬로 처리하는 방법과, `-resume`을 사용하여 완료된 작업을 반복하지 않는 방법을 알게 되었습니다.

### 다음 단계

완전한 다단계 파이프라인이 채널을 사용하여 프로세스를 연결하는 방법과, 컨테이너를 사용하여 분석 도구와 의존성을 관리하는 방법을 학습합니다.

---

## 3. 다단계 파이프라인 실행

지금까지 단일 프로세스를 실행하고, 여러 입력에 대해 병렬로 여러 번 실행해 보았습니다.
실제 파이프라인은 보통 더 나아가, 여러 프로세스를 연결하여 한 프로세스의 출력을 다음 프로세스의 입력으로 전달하며, 종종 여러 소프트웨어를 함께 사용합니다.
워크플로우 `main.nf`는 이 두 가지를 결합하여 완전한 파이프라인을 구성합니다.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello_pipeline_complete.svg"
</figure>

각 입력 인사말은 네 단계를 모두 거칩니다. `sayHello`가 파일에 인사말을 작성하고, `convertToUpper`가 텍스트를 대문자로 변환하고, `collectGreetings`가 모든 결과를 하나의 파일로 병합하고, `cowpy`가 컨테이너화된 도구를 사용하여 병합된 출력에서 ASCII 아트를 생성합니다.
Nextflow는 채널을 통해 이 단계들을 연결합니다. 한 프로세스의 출력이 다음 프로세스의 입력이 되므로, 데이터가 준비되는 즉시 전체 체인이 자동으로 실행됩니다. 각 단계를 직접 조율할 필요가 없습니다.

이 워크플로우는 모듈을 사용합니다. 각 프로세스는 `modules/` 아래의 자체 파일에 정의되어 있으며, `main.nf`는 인라인으로 정의하는 대신 `include` 구문으로 가져옵니다.
이를 통해 코드를 중복하지 않고 각 프로세스를 여러 워크플로우에서 재사용할 수 있습니다. 자세한 내용은 아래의 코드 살펴보기 섹션을 참조하세요.

### 3.1. 워크플로우 실행

터미널에서 다음 명령을 실행합니다.

```bash
nextflow run main.nf --input data/greetings.csv
```

`character` 매개변수는 `nextflow.config`에서 기본값이 `turkey`로 설정되어 있으므로, 재정의하지 않으면 ASCII 아트에 칠면조가 사용됩니다(`--character tux`를 추가해 보세요).

??? success "명령 출력"

    ```console hl_lines="6 7 8 9"
    N E X T F L O W   ~  version 26.04.4

    Launching `main.nf` [nostalgic_brahmagupta] revision: ce74f81996

    executor >  local (8)
    [56/8499f6] sayHello (3)       | 3 of 3 ✔
    [cc/0ee42a] convertToUpper (3) | 3 of 3 ✔
    [eb/0f2e24] collectGreetings   | 1 of 1 ✔
    [b5/34e07f] cowpy              | 1 of 1 ✔

    Outputs:

      /workspaces/training/nextflow-run/results

      first_output:
        - batch/intermediates/Bonjour-output.txt
        - batch/intermediates/Hello-output.txt
        - batch/intermediates/Hola-output.txt

      uppercased:
        - batch/intermediates/UPPER-Bonjour-output.txt
        - batch/intermediates/UPPER-Hello-output.txt
        - batch/intermediates/UPPER-Hola-output.txt

      collected: batch/intermediates/COLLECTED-batch-output.txt

      batch_report: batch/batch-report.txt

      cowpy_art: batch/cowpy-COLLECTED-batch-output.txt
    ```

네 개의 프로세스가 실행되었지만, 실행 횟수는 각각 다릅니다.
`sayHello`와 `convertToUpper`는 각 입력마다 한 번씩 실행되었습니다(3 of 3). 각 인사말을 개별적으로 파일에 작성하고 대문자로 변환해야 하기 때문입니다.
`collectGreetings`와 `cowpy`는 각각 한 번만 실행되었습니다(1 of 1). 인사말 병합과 ASCII 아트 생성은 모든 개별 결과가 준비된 후에만 의미가 있기 때문입니다.
여러 병렬 작업이 더 적은 수의 다운스트림 작업으로 수렴하는 이러한 팬아웃-팬인 구조는 실제 파이프라인에서 흔히 볼 수 있습니다.

Nextflow는 전체 단계가 완료될 때까지 기다리지 않고 다음 단계를 시작합니다.
하나의 `sayHello` 출력이 준비되는 즉시 해당 `convertToUpper` 작업이 시작될 수 있으므로, 서로 다른 프로세스의 작업이 엄격한 배치 순서가 아닌 동시에 실행됩니다.
`collectGreetings`와 `cowpy`는 각각 모든 업스트림 결과가 준비될 때까지 기다려야 합니다.

`results` 디렉토리는 이러한 팬인 구조와 파이프라인 작성자가 선택한 게시 위치를 반영합니다. 이 구조를 정의하는 것은 1.4의 코드 살펴보기에서 다룬 `output` 블록입니다.

```console title="results/"
results
└── batch
    ├── batch-report.txt
    ├── cowpy-COLLECTED-batch-output.txt
    └── intermediates
        ├── Bonjour-output.txt
        ├── COLLECTED-batch-output.txt
        ├── Hello-output.txt
        ├── Hola-output.txt
        ├── UPPER-Bonjour-output.txt
        ├── UPPER-Hello-output.txt
        └── UPPER-Hola-output.txt
```

최상위 디렉토리는 `batch` 매개변수의 이름을 따르며, 기본값은 `batch`입니다. 이후 연습에서 변경되는 것을 확인할 수 있습니다.

ASCII 아트 파일은 `cowpy-COLLECTED-batch-output.txt`에서 확인합니다.

??? abstract "파일 내용"

    ```console title="results/batch/cowpy-COLLECTED-batch-output.txt"
     _________
    / HELLO   \
    | BONJOUR |
    \ HOLA    /
     ---------
      \                                  ,+*^^*+___+++_
       \                           ,*^^^^              )
        \                       _+*                     ^**+_
         \                    +^       _ _++*+_+++_,         )
                  _+^^*+_    (     ,+*^ ^          \+_        )
                 {       )  (    ,(    ,_+--+--,      ^)      ^\
                { (\@)    } f   ,(  ,+-^ __*_*_  ^^\_   ^\       )
               {:;-/    (_+*-+^^^^^+*+*<_ _++_)_    )    )      /
              ( /  (    (        ,___    ^*+_+* )   <    <      \
               U _/     )    *--<  ) ^\-----++__)   )    )       )
                (      )  _(^)^^))  )  )\^^^^^))^*+/    /       /
              (      /  (_))_^)) )  )  ))^^^^^))^^^)__/     +^^
             (     ,/    (^))^))  )  ) ))^^^^^^^))^^)       _)
              *+__+*       (_))^)  ) ) ))^^^^^^))^^^^^)____*^
              \             \_)^)_)) ))^^^^^^^^^^))^^^^)
               (_             ^\__^^^^^^^^^^^^))^^^^^^^)
                 ^\___            ^\__^^^^^^))^^^^^^^^)\\
                      ^^^^^\uuu/^^\uuu/^^^^\^\^\^\^\^\^\^\
                         ___) >____) >___   ^\_\_\_\_\_\_\)
                        ^^^//\\_^^//\\_^       ^(\_\_\_\)
                          ^^^ ^^ ^^^ ^
    ```

[2.1](#21-run-the-workflow)에서와 마찬가지로, 네 개의 프로세스에 걸친 8번의 작업 실행 각각은 `work/` 아래에 완전히 격리된 자체 디렉토리를 갖습니다.
`collectGreetings`는 이것이 왜 중요한지를 잘 보여줍니다. 이 프로세스는 세 개의 서로 다른 작업 디렉토리에 있는 세 `convertToUpper` 작업의 출력에 의존합니다. 따라서 Nextflow는 업스트림 작업 디렉토리에서 직접 읽는 대신, `collectGreetings` 자체 디렉토리 안에 해당 파일들의 심볼릭 링크를 생성합니다:

```console title="work/eb/0f2e24.../"
COLLECTED-batch-output.txt
UPPER-Bonjour-output.txt -> ../../69/e6c057.../UPPER-Bonjour-output.txt
UPPER-Hello-output.txt   -> ../../cc/ba7d19.../UPPER-Hello-output.txt
UPPER-Hola-output.txt    -> ../../cc/0ee42a.../UPPER-Hola-output.txt
batch-report.txt
.command.sh
```

각 작업은 출처에 관계없이 필요한 특정 파일만 볼 수 있으며, 다른 작업 디렉토리의 내부 내용은 볼 수 없습니다.
전체 파이프라인에 걸쳐, [2.1](#21-run-the-workflow)에서 단일 프로세스로 확인한 것과 동일한 격리 덕분에 Nextflow는 모든 프로세스의 모든 작업을 안전하게 동시에 실행할 수 있습니다.

!!! note "참고"

    `cowpy` 단계는 로컬에 설치된 소프트웨어에 의존하지 않고 Docker 컨테이너 안에서 실행됩니다.
    컨테이너는 애플리케이션과 실행에 필요한 모든 것을 함께 패키징하므로, 의존성을 직접 설치하고 관리할 필요가 없으며, 컨테이너를 실행할 수 있는 모든 시스템에서 파이프라인이 동일하게 동작합니다.
    Nextflow는 컨테이너의 대안으로 Conda도 지원합니다. 전환 방법은 [파트 2](./02_configure_pipeline.md)를 참조하세요.

### 3.2. 선택 사항: 코드 살펴보기

파이프라인 실행만이 목적이라면 코드를 이해하는 것이 필수는 아니지만, 궁금하다면 살펴볼 가치가 있습니다.

??? optional "클릭하여 이 연습과 관련된 코드 살펴보기"

    ### 한 단계에서 다음 단계로 데이터가 흐르는 방식

    각 프로세스는 출력 채널을 다음 프로세스에 전달합니다:

    ```groovy title="main.nf" linenums="19" hl_lines="7 8 9"
        main:
        // CSV 파일에서 입력을 위한 채널 생성
        greeting_ch = channel.fromPath(params.input)
                            .splitCsv()
                            .map { line -> line[0] }
        sayHello(greeting_ch)
        convertToUpper(sayHello.out)
        collectGreetings(convertToUpper.out.collect(), params.batch)
        cowpy(collectGreetings.out.outfile, params.character)
    ```

    `processName.out` 패턴은 프로세스의 출력 채널을 참조합니다.

    `.collect()` 연산자는 `convertToUpper`의 모든 개별 출력을 단일 채널 항목으로 수집한 후 `collectGreetings`에 전달합니다.

    ### 프로세스 모듈 사용

    `main.nf`는 프로세스 코드를 직접 정의하지 않습니다.
    대신, `modules/` 아래의 각 파일에서 프로세스를 가져옵니다:

    ```groovy title="main.nf" linenums="3"
    include { sayHello } from './modules/sayHello.nf'
    include { convertToUpper } from './modules/convertToUpper.nf'
    include { collectGreetings } from './modules/collectGreetings.nf'
    include { cowpy } from './modules/cowpy.nf'
    ```

    각 모듈 파일에는 [1.4](#14-optional-code-walkthrough)의 `sayHello` 모듈과 동일한 구조로 단일 프로세스 정의가 포함되어 있습니다.
    프로세스를 별도의 파일로 분리하면 코드를 중복하지 않고 여러 워크플로우에서 재사용할 수 있습니다.

    <figure class="excalidraw">
    --8<-- "docs/en/docs/hello_nextflow/img/modules.svg"
    </figure>

    ### 컨테이너화된 소프트웨어 사용

    `cowpy` 프로세스는 모듈 파일에 지정된 Docker 컨테이너 안에서 실행됩니다:

    ```groovy title="modules/cowpy.nf" linenums="2" hl_lines="3"
    process cowpy {

        container 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'

        input:
        path input_file
        val character

        output:
        path "cowpy-${input_file}"

        script:
        """
        cat ${input_file} | cowpy -c "${character}" > cowpy-${input_file}
        """
    }
    ```

    Nextflow는 자동으로 이미지를 가져오고, 컨테이너 안에서 스크립트를 실행한 후 정리합니다.
    Docker는 `nextflow.config`에서 이 프로젝트에 대해 활성화되어 있습니다:

    ```groovy title="nextflow.config"
    docker.enabled = true
    ```

    이 한 줄로 파이프라인에서 컨테이너가 지정된 모든 프로세스에 Docker가 활성화됩니다.

### 핵심 정리

컨테이너화된 도구를 사용하여 여러 입력을 병렬로 처리하는 완전한 다단계 파이프라인을 실행했습니다.

### 다음 단계

[파트 2](./02_configure_pipeline.md)로 이동하여 `nextflow.config`를 사용하여 파이프라인 동작을 설정하는 방법을 학습합니다.

---

## 요약

이 파트에서 학습한 내용:

- Nextflow 워크플로우를 실행하고 출력을 찾는 방법
- `work/` 디렉토리와 로그 파일 살펴보기
- CSV 파일에서 여러 입력을 병렬로 처리하는 방법
- 새로운 입력을 추가할 때 `-resume`으로 완료된 작업을 건너뛰는 방법
- 컨테이너화된 도구를 사용하는 다단계 파이프라인 실행
