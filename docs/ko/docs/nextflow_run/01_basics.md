# 파트 1: 기본 작업 실행하기

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } AI 지원 번역 - [자세히 알아보기 및 개선 사항 제안하기](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Nextflow Run 교육 과정의 첫 번째 파트에서는 매우 기본적인 도메인 독립적 Hello World 예제로 시작합니다. 이를 통해 필수 작업을 시연하고 해당하는 Nextflow 코드 구성 요소를 설명합니다.

??? info "Hello World 예제란 무엇인가요?"

    "Hello World!"는 프로그래밍 언어나 소프트웨어 프레임워크의 기본 구문과 구조를 보여주기 위한 최소한의 예제입니다.
    일반적으로 "Hello, World!" 문구를 콘솔이나 터미널과 같은 출력 장치에 출력하거나 파일에 작성하는 것으로 구성됩니다.

---

## 1. Hello World를 직접 실행하기

Nextflow로 적용하기 전에 어떤 작업을 하는지 보여주기 위해 터미널에서 직접 실행하는 간단한 명령으로 이 개념을 시연하겠습니다.

!!! tip

    [시작하기](00_orientation.md) 페이지에 설명된 대로 `nextflow-run/` 디렉토리 안에 있어야 합니다.

### 1.1. 터미널에서 인사말 출력하기

터미널에서 다음 명령을 실행하세요.

```bash
echo 'Hello World!'
```

??? success "명령 출력"

    ```console
    Hello World!
    ```

이 명령은 터미널에 바로 'Hello World' 텍스트를 출력합니다.

### 1.2. 출력을 파일에 작성하기

파이프라인 실행은 주로 파일에서 데이터를 읽고 결과를 다른 파일에 작성하는 것과 관련이 있으므로, 예제를 좀 더 실용적으로 만들기 위해 텍스트 출력을 파일에 작성하도록 명령을 수정하겠습니다.

```bash
echo 'Hello World!' > output.txt
```

??? success "명령 출력"

    ```console

    ```

이 명령은 터미널에 아무것도 출력하지 않습니다.

### 1.3. 출력 찾기

'Hello World' 텍스트는 이제 우리가 지정한 출력 파일인 `output.txt`에 있어야 합니다.
파일 탐색기에서 열거나 명령줄에서 `cat` 유틸리티를 사용하여 열 수 있습니다.

??? abstract "파일 내용"

    ```console title="output.txt" linenums="1"
    Hello World!
    ```

이것이 우리가 첫 번째 Nextflow 워크플로우로 재현하려는 것입니다.

### 핵심 정리

터미널에서 텍스트를 출력하는 간단한 명령을 실행하는 방법과 선택적으로 출력을 파일에 작성하는 방법을 배웠습니다.

### 다음 단계

동일한 결과를 달성하는 Nextflow 워크플로우를 실행하는 데 필요한 것이 무엇인지 알아보세요.

---

## 2. 워크플로우 실행하기

`--input`이라는 명령줄 인수를 통해 입력 인사말을 받아 해당 인사말이 포함된 텍스트 파일을 생성하는 `1-hello.nf`라는 워크플로우 스크립트를 제공합니다.

아직 코드를 살펴보지 않을 것입니다. 먼저 실행하는 것이 어떤 모습인지 살펴보겠습니다.

### 2.1. 워크플로우 실행 및 모니터링

터미널에서 다음 명령을 실행하세요.

```bash
nextflow run 1-hello.nf --input 'Hello World!'
```

??? success "명령 출력"

    ```console hl_lines="6"
    N E X T F L O W   ~  version 25.10.2

    Launching `1-hello.nf` [goofy_torvalds] DSL2 - revision: c33d41f479

    executor >  local (1)
    [a3/7be2fa] sayHello | 1 of 1 ✔
    ```

콘솔 출력이 이와 비슷하다면, 축하합니다. 첫 번째 Nextflow 워크플로우를 실행했습니다!

여기서 가장 중요한 출력은 위에서 강조 표시된 마지막 줄입니다:

```console
[a3/7be2fa] sayHello | 1 of 1 ✔
```

이는 `sayHello` 프로세스가 한 번 성공적으로 실행되었음을 알려줍니다(`1 of 1 ✔`).

좋습니다만, 출력은 어디에 있을까요?

### 2.2. `results` 디렉토리에서 출력 파일 찾기

이 워크플로우는 출력을 results 디렉토리에 게시하도록 구성되어 있습니다.
현재 디렉토리를 보면 워크플로우를 실행했을 때 Nextflow가 `results`라는 새 디렉토리와 그 아래에 `1-hello`라는 하위 디렉토리를 생성했으며, 그 안에 `output.txt`라는 파일이 있는 것을 볼 수 있습니다.

```console title="results/"
results
└── 1-hello
    └── output.txt
```

파일을 열어보세요. 내용이 명령줄에서 지정한 문자열과 일치해야 합니다.

```console title="results/1-hello/output.txt" linenums="1"
Hello World!
```

좋습니다. 워크플로우가 의도한 대로 작동했습니다!

### 2.3. 결과를 다른 디렉토리에 저장하기

기본적으로 Nextflow는 파이프라인 출력을 현재 경로의 `results`라는 디렉토리에 저장합니다.
파일이 게시되는 위치를 변경하려면 `-output-dir` CLI 플래그(또는 줄여서 `-o`)를 사용하세요.

!!! danger

    `--input`은 하이픈이 두 개이고 `-output-dir`은 하나입니다!
    이는 `--input`이 파이프라인 _매개변수_이고 `-output-dir`은 Nextflow 코어 CLI 플래그이기 때문입니다.
    이에 대해서는 나중에 자세히 설명하겠습니다.

```bash
nextflow run 1-hello.nf --input 'Hello World!' -output-dir hello_results
```

??? success "명령 출력"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `1-hello.nf` [hungry_celsius] DSL2 - revision: f048d6ea78

    executor >  local (1)
    [a3/1e1535] sayHello [100%] 1 of 1 ✔
    ```

이제 출력이 `results` 대신 `hello_results`라는 디렉토리에 게시된 것을 볼 수 있습니다:

```console title="hello_results/"
hello_results
└── 1-hello
    └── output.txt
```

이 디렉토리 내의 파일은 이전과 동일하며, 최상위 디렉토리만 다릅니다.
그러나 두 경우 모두 '게시된' 결과는 Nextflow가 워크플로우를 실행할 때 생성한 실제 출력의 복사본(또는 경우에 따라 심볼릭 링크)임을 유의하세요.

이제 Nextflow가 실제로 작업을 실행한 위치를 확인하기 위해 내부를 살펴보겠습니다.

!!! Warning

    모든 워크플로우가 results 디렉토리에 출력을 게시하도록 설정되어 있는 것은 아니며, 디렉토리 이름과 구조가 다를 수 있습니다.
    이 섹션의 뒷부분에서 이 동작이 어디에 지정되어 있는지 확인하는 방법을 보여드리겠습니다.

### 2.4. `work/` 디렉토리에서 원본 출력 및 로그 찾기

워크플로우를 실행하면 Nextflow는 워크플로우의 각 프로세스 호출마다(=파이프라인의 각 단계마다) 고유한 '작업 디렉토리'를 생성합니다.
각 작업 디렉토리에서 필요한 입력을 준비하고, 관련 명령을 실행하며, 출력과 로그 파일을 작성합니다. 이 디렉토리는 고유성을 위해 해시를 사용하여 자동으로 이름이 지정됩니다.

이러한 모든 작업 디렉토리는 현재 디렉토리(명령을 실행하는 위치) 아래의 `work`라는 디렉토리에 있습니다.

혼란스러울 수 있으니 실제로 어떻게 보이는지 살펴보겠습니다.

앞서 실행한 워크플로우의 콘솔 출력으로 돌아가면 다음과 같은 줄이 있었습니다:

```console
[a3/1e1535] sayHello [100%] 1 of 1 ✔
```

줄이 `[a3/1e1535]`로 시작하는 것을 보셨나요?
이것은 해당 프로세스 호출의 작업 디렉토리 경로를 축약한 형태이며, `work/` 디렉토리 경로 내에서 `sayHello` 프로세스 호출의 출력을 찾을 수 있는 위치를 알려줍니다.

다음 명령을 입력하고(터미널에 표시된 것으로 `a3/1e1535`를 대체) 탭 키를 눌러 경로를 자동 완성하거나 별표를 추가하여 전체 경로를 찾을 수 있습니다:

```bash
ls work/a3/1e1535*
```

이렇게 하면 전체 디렉토리 경로가 나타납니다: `work/a3/1e153543b0a7f9d2c4735ddb4ab231`

그 안에 무엇이 있는지 살펴보겠습니다.

??? abstract "디렉토리 내용"

    ```console
    work
    ├── a3
    │   └── 1e153543b0a7f9d2c4735ddb4ab231
    │       ├── .command.begin
    │       ├── .command.err
    │       ├── .command.log
    │       ├── .command.out
    │       ├── .command.run
    │       ├── .command.sh
    │       ├── .exitcode
    │       └── output.txt
    └── a4
        └── aa3694b8808bdcc1135ef4a1187a4d
            ├── .command.begin
            ├── .command.err
            ├── .command.log
            ├── .command.out
            ├── .command.run
            ├── .command.sh
            ├── .exitcode
            └── output.txt
    ```

??? question "같은 내용이 보이지 않나요?"

    정확한 하위 디렉토리 이름은 시스템마다 다릅니다.

    VSCode 파일 탐색기에서 작업 하위 디렉토리의 내용을 탐색하면 모든 파일이 바로 표시됩니다.
    그러나 로그 파일은 터미널에서 보이지 않도록 설정되어 있으므로 `ls` 또는 `tree`를 사용하여 보려면 숨김 파일을 표시하는 관련 옵션을 설정해야 합니다.

    ```bash
    tree -a work
    ```

`work/`에는 우리가 수행한 두 번의 파이프라인 실행에서 생성된 두 세트의 디렉토리가 있습니다.
각 작업 실행은 작업할 고유하고 격리된 디렉토리를 갖습니다.
이 경우 파이프라인이 두 번 모두 동일한 작업을 수행했으므로 각 작업 디렉토리의 내용은 동일합니다.

`output.txt` 파일을 바로 알아볼 수 있을 것입니다. 이것은 실제로 `results` 디렉토리에 게시된 `sayHello` 프로세스의 원본 출력입니다.
파일을 열면 다시 `Hello World!` 인사말을 찾을 수 있습니다.

```console title="work/a3/1e153543b0a7f9d2c4735ddb4ab231/output.txt"
Hello World!
```

그렇다면 다른 모든 파일은 무엇일까요?

이것들은 Nextflow가 작업 실행의 일부로 작성한 도우미 및 로그 파일입니다:

- **`.command.begin`**: 작업이 시작되자마자 생성되는 센티널 파일
- **`.command.err`**: 프로세스 호출에서 발생한 오류 메시지(`stderr`)
- **`.command.log`**: 프로세스 호출에서 발생한 전체 로그 출력
- **`.command.out`**: 프로세스 호출의 일반 출력(`stdout`)
- **`.command.run`**: Nextflow가 프로세스 호출을 실행하기 위해 실행한 전체 스크립트
- **`.command.sh`**: 프로세스 호출에서 실제로 실행된 명령
- **`.exitcode`**: 명령 실행 결과로 나온 종료 코드

`.command.sh` 파일은 모든 부기 및 작업/환경 설정을 제외하고 Nextflow가 실행한 주요 명령을 보여주기 때문에 특히 유용합니다.

```console title="work/a3/1e153543b0a7f9d2c4735ddb4ab231/.command.sh"
#!/bin/bash -ue
echo 'Hello World!' > output.txt

```

이것은 워크플로우가 앞서 명령줄에서 직접 실행한 것과 동일한 명령을 구성했음을 확인시켜 줍니다.

문제가 발생하여 무슨 일이 일어났는지 해결해야 할 때, `command.sh` 스크립트를 보고 워크플로우 지침, 변수 보간 등을 기반으로 Nextflow가 정확히 어떤 명령을 구성했는지 확인하는 것이 유용할 수 있습니다.

### 2.5. 다른 인사말로 워크플로우 다시 실행하기

`--input` 인수에 다른 값을 사용하여 워크플로우를 몇 번 다시 실행한 다음 작업 디렉토리를 살펴보세요.

??? abstract "디렉토리 내용"

    ```console
    work/
    ├── 09
    │   └── 5ea8665939daf6f04724286c9b3c8a
    │       ├── .command.begin
    │       ├── .command.err
    │       ├── .command.log
    │       ├── .command.out
    │       ├── .command.run
    │       ├── .command.sh
    │       ├── .exitcode
    │       └── output.txt
    ├── 92
    │   └── ceb95e05d87621c92a399da9bd2067
    │       ├── .command.begin
    │       ├── .command.err
    │       ├── .command.log
    │       ├── .command.out
    │       ├── .command.run
    │       ├── .command.sh
    │       ├── .exitcode
    │       └── output.txt
    ├── 93
    │   └── 6708dbc20c7efdc6769cbe477061ec
    │       ├── .command.begin
    │       ├── .command.err
    │       ├── .command.log
    │       ├── .command.out
    │       ├── .command.run
    │       ├── .command.sh
    │       ├── .exitcode
    │       └── output.txt
    ├── a3
    │   └── 1e153543b0a7f9d2c4735ddb4ab231
    │       ├── .command.begin
    │       ├── .command.err
    │       ├── .command.log
    │       ├── .command.out
    │       ├── .command.run
    │       ├── .command.sh
    │       ├── .exitcode
    │       └── output.txt
    └── a4
        └── aa3694b8808bdcc1135ef4a1187a4d
            ├── .command.begin
            ├── .command.err
            ├── .command.log
            ├── .command.out
            ├── .command.run
            ├── .command.sh
            ├── .exitcode
            └── output.txt
    ```

각 실행마다 완전한 출력 및 로그 파일 세트가 포함된 새 하위 디렉토리가 생성된 것을 볼 수 있습니다.

반면 `results` 디렉토리를 보면 여전히 하나의 결과 세트만 있으며, 출력 파일의 내용은 마지막으로 실행한 것에 해당합니다.

??? abstract "디렉토리 내용"

    ```console title="results/"
    results
    └── 1-hello
        └── output.txt
    ```

이는 게시된 결과가 후속 실행에 의해 덮어쓰여지는 반면, `work/` 아래의 작업 디렉토리는 보존됨을 보여줍니다.

### 핵심 정리

간단한 Nextflow 스크립트를 실행하고, 실행을 모니터링하며, 출력을 찾는 방법을 배웠습니다.

### 다음 단계

기본 Nextflow 스크립트를 읽고 구성 요소가 기능과 어떻게 관련되는지 파악하는 방법을 배워보세요.

---

## 3. Hello World 워크플로우 시작 스크립트 살펴보기

지금까지 워크플로우 스크립트를 블랙박스처럼 다뤘습니다.
이제 무엇을 하는지 보았으니 상자를 열고 내부를 살펴보겠습니다.

여기서 우리의 목표는 Nextflow 코드의 구문을 암기하는 것이 아니라, 주요 구성 요소가 무엇이고 어떻게 구성되어 있는지에 대한 기본적인 직관을 형성하는 것입니다.

### 3.1. 전체 코드 구조 살펴보기

현재 디렉토리인 `nextflow-run`에서 `1-hello.nf` 스크립트를 찾을 수 있습니다. 편집기 창에서 열어보세요.

??? full-code "전체 코드 파일"

    ```groovy title="1-hello.nf" linenums="1"
    #!/usr/bin/env nextflow

    /*
    * echo를 사용하여 'Hello World!'를 파일에 출력
    */
    process sayHello {

        input:
        val greeting

        output:
        path 'output.txt'

        script:
        """
        echo '${greeting}' > output.txt
        """
    }

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

Nextflow 워크플로우 스크립트는 일반적으로 하나 이상의 **process** 정의, **workflow** 자체, 그리고 **params** 및 **output**과 같은 몇 가지 선택적 블록을 포함합니다.

각 **process**는 파이프라인의 해당 단계가 수행해야 하는 작업을 설명하고, **workflow**는 다양한 단계를 연결하는 데이터 흐름 로직을 설명합니다.

먼저 **process** 블록을 자세히 살펴본 다음 **workflow** 블록을 살펴보겠습니다.

### 3.2. `process` 정의

첫 번째 코드 블록은 [**process**](https://nextflow.io/docs/latest/process.html)를 설명합니다.
프로세스 정의는 `process` 키워드로 시작하고, 프로세스 이름이 이어지며, 마지막으로 중괄호로 구분된 프로세스 본문이 옵니다.
프로세스 본문에는 실행할 명령을 지정하는 script 블록이 포함되어야 하며, 이는 명령줄 터미널에서 실행할 수 있는 모든 것이 될 수 있습니다.

```groovy title="1-hello.nf" linenums="3"
/*
* echo를 사용하여 인사말을 파일에 출력
*/
process sayHello {

    input:
    val greeting

    output:
    path 'output.txt'

    script:
    """
    echo '${greeting}' > output.txt
    """
}
```

여기에는 `greeting`이라는 **input** 변수를 받아 `output.txt`라는 파일에 **output**을 작성하는 `sayHello`라는 **process**가 있습니다.

<figure class="excalidraw">
--8<-- "docs/en/docs/nextflow_run/img/sayhello_with_input.svg"
</figure>

이것은 `input` 정의, `output` 정의 및 실행할 `script`만 포함하는 매우 최소한의 프로세스 정의입니다.

`input` 정의에는 `val` 한정자가 포함되어 있으며, 이는 Nextflow에게 어떤 종류의 값(문자열, 숫자 등)을 기대하도록 지시합니다.

`output` 정의에는 `path` 한정자가 포함되어 있으며, 이는 Nextflow에게 이것을 경로(디렉토리 경로와 파일 모두 포함)로 처리하도록 지시합니다.

### 3.3. `workflow` 정의

두 번째 코드 블록은 [**workflow**](https://nextflow.io/docs/latest/workflow.html) 자체를 설명합니다.
워크플로우 정의는 `workflow` 키워드로 시작하고, 선택적 이름이 이어지며, 중괄호로 구분된 워크플로우 본문이 옵니다.

여기에는 `main:` 블록과 `publish:` 블록으로 구성된 **workflow**가 있습니다.
`main:` 블록은 워크플로우의 주요 본문이고 `publish:` 블록은 `results` 디렉토리에 게시되어야 하는 출력을 나열합니다.

```groovy title="1-hello.nf" linenums="27"
workflow {

    main:
    // 인사말을 내보냅니다
    sayHello(params.input)

    publish:
    first_output = sayHello.out
}
```

이 경우 `main:` 블록에는 `sayHello` 프로세스 호출이 포함되어 있으며 인사말로 사용할 `params.input`이라는 입력을 제공합니다.

곧 더 자세히 논의하겠지만, `params.input`은 명령줄에서 `--input` 매개변수에 제공한 값을 보유합니다.

`publish:` 블록은 `sayHello()` 프로세스 호출의 출력을 나열하며, 이를 `sayHello.out`으로 참조하고 `first_output`이라는 이름을 부여합니다(워크플로우 작성자가 원하는 대로 지정할 수 있습니다).

이것은 매우 최소한의 **workflow** 정의입니다.
실제 파이프라인에서 워크플로우는 일반적으로 **channel**로 연결된 여러 **process** 호출을 포함하며, 변수 입력에 대한 기본값이 설정되어 있을 수 있습니다.

이에 대해서는 과정의 파트 2에서 다루겠습니다.
지금은 워크플로우가 입력과 출력을 어떻게 처리하는지 자세히 살펴보겠습니다.

### 3.4. 명령줄 매개변수의 `params` 시스템

`sayHello()` 프로세스 호출에 제공하는 `params.input`은 Nextflow 코드의 멋진 부분이며 추가로 시간을 할애할 가치가 있습니다.

위에서 언급했듯이, 이것이 `--input` 명령줄 매개변수의 값을 `sayHello()` 프로세스 호출에 전달하는 방법입니다.
실제로 `params.someParameterName`을 선언하는 것만으로도 워크플로우에 명령줄에서 `--someParameterName`이라는 매개변수를 제공할 수 있습니다.

여기서는 워크플로우가 기대하는 입력 유형을 지정하는 `params` 블록을 설정하여 해당 매개변수 선언을 공식화했습니다(Nextflow 25.10.2 이상).

```groovy title="1-hello.nf" linenums="20"
/*
 * 파이프라인 매개변수
 */
params {
    input: String
}
```

지원되는 유형에는 `String`, `Integer`, `Float`, `Boolean` 및 `Path`가 포함됩니다.
자세한 내용은 Nextflow 참조 문서의 [워크플로우 매개변수](https://nextflow.io/docs/latest/config.html#workflow-parameters)를 참조하세요.

!!! tip

    `params` 시스템을 사용하여 선언된 _워크플로우_ 매개변수는 명령줄에서 항상 하이픈 두 개(`--`)를 사용합니다.
    이는 하이픈 하나(`-`)만 사용하는 _Nextflow 수준_ CLI 플래그와 구별됩니다.

### 3.5. `publish` 지시문

워크플로우의 다른 쪽 끝에서 `publish:` 블록을 이미 살펴봤습니다.
이것은 출력 처리 시스템의 절반이며, 나머지 절반은 아래에 있는 `output` 블록입니다.

```groovy title="1-hello.nf" linenums="37"
output {
    first_output {
        path '1-hello'
        mode 'copy'
    }
}
```

이것은 `publish:` 블록에 나열된 `first_output` 출력이 기본 `results` 출력 디렉토리 아래의 `1-hello`라는 하위 디렉토리에 복사되어야 함을 지정합니다.

`mode 'copy'` 줄은 시스템의 기본 동작을 재정의하는데, 기본 동작은 적절한 복사본 대신 `work/` 디렉토리의 원본 파일에 대한 심볼릭 링크(또는 symlink)를 만드는 것입니다.

게시 동작을 제어하기 위한 옵션이 여기에 표시된 것보다 더 많이 있습니다. 나중에 몇 가지를 다루겠습니다.
또한 워크플로우가 여러 출력을 생성할 때 각 출력이 `output` 블록에 이런 식으로 나열되는 것을 볼 수 있습니다.

자세한 내용은 Nextflow 참조 문서의 [출력 게시](https://nextflow.io/docs/latest/workflow.html#publishing-outputs)를 참조하세요.

??? info "`publishDir`을 사용한 출력 게시의 이전 구문"

    최근까지 출력을 게시하는 확립된 방법은 `publishDir` 지시문을 사용하여 각 개별 프로세스 수준에서 수행하는 것이었습니다.

    이전 Nextflow 파이프라인과 프로세스 모듈에서 이 코드 패턴을 여전히 많이 볼 수 있으므로 이를 인식하는 것이 중요합니다.

    워크플로우에 `publish:` 블록과 최상위 수준의 `output` 블록을 갖는 대신, `sayHello` 프로세스 정의에 `publishDir` 줄이 있는 것을 볼 수 있습니다:

    ```groovy title="구문 예제" linenums="1" hl_lines="3"
    process sayHello {

        publishDir 'results/1-hello', mode: 'copy'

        output:
        path 'output.txt'

        script:
        """
        echo 'Hello World!' > output.txt
        """
    }
    ```

    그러나 향후 Nextflow 언어 버전에서 결국 허용되지 않을 것이므로 새로운 작업에서는 이것을 사용하지 않는 것이 좋습니다.

### 핵심 정리

간단한 Nextflow 워크플로우가 어떻게 구성되어 있는지, 그리고 기본 구성 요소가 기능과 어떻게 관련되는지 배웠습니다.

### 다음 단계

워크플로우 실행을 편리하게 관리하는 방법을 배워보세요.

---

## 4. 워크플로우 실행 관리하기

워크플로우를 실행하고 출력을 검색하는 방법을 아는 것은 좋지만, 삶을 더 쉽게 만들어줄 워크플로우 관리의 몇 가지 다른 측면이 있다는 것을 곧 알게 될 것입니다.

여기서는 동일한 워크플로우를 다시 실행해야 할 때 `resume` 기능을 활용하는 방법, `nextflow log`로 실행 로그를 검사하는 방법, `nextflow clean`으로 이전 작업 디렉토리를 삭제하는 방법을 보여드립니다.

### 4.1. `-resume`으로 워크플로우 다시 실행하기

때때로 이전에 실행한 파이프라인을 이미 성공적으로 완료된 작업을 다시 수행하지 않고 다시 실행하고 싶을 것입니다.

Nextflow에는 이를 수행할 수 있는 `-resume`이라는 옵션이 있습니다.
구체적으로, 이 모드에서는 정확히 동일한 코드, 설정 및 입력으로 이미 실행된 프로세스는 건너뜁니다.
즉, Nextflow는 마지막 실행 이후 추가하거나 수정한 프로세스 또는 새로운 설정이나 입력을 제공하는 프로세스만 실행합니다.

이렇게 하면 두 가지 주요 이점이 있습니다:

- 파이프라인을 개발하는 중이라면 변경 사항을 테스트하기 위해 적극적으로 작업 중인 프로세스만 실행하면 되므로 더 빠르게 반복할 수 있습니다.
- 프로덕션에서 파이프라인을 실행하다가 문제가 발생하면 많은 경우 문제를 해결하고 파이프라인을 다시 실행할 수 있으며, 실패 지점부터 실행을 재개하여 많은 시간과 컴퓨팅 리소스를 절약할 수 있습니다.

사용하려면 명령에 `-resume`을 추가하고 실행하세요:

```bash
nextflow run 1-hello.nf --input 'Hello World!' -resume
```

??? success "명령 출력"

    ```console linenums="1"
    N E X T F L O W   ~  version 25.10.2

    Launching `1-hello.nf` [tiny_noyce] DSL2 - revision: c33d41f479

    [a3/7be2fa] sayHello | 1 of 1, cached: 1 ✔
    ```

콘솔 출력은 익숙해 보이지만 이전과 약간 다른 점이 하나 있습니다.

프로세스 상태 줄(5번째 줄)에 추가된 `cached:` 부분을 찾아보세요. 이는 Nextflow가 이미 이 작업을 수행했음을 인식하고 이전의 성공적인 실행 결과를 단순히 재사용했음을 의미합니다.

또한 작업 하위 디렉토리 해시가 이전 실행과 동일한 것을 볼 수 있습니다.
Nextflow는 말 그대로 이전 실행을 가리키며 "저기서 이미 했어요"라고 말하고 있습니다.

!!! tip

    `resume`으로 파이프라인을 다시 실행할 때 Nextflow는 이전에 성공적으로 실행된 실행에 의해 work 디렉토리 외부에 게시된 파일을 덮어쓰지 않습니다.

    자세한 내용은 Nextflow 참조 문서의 [캐시 및 재개](https://nextflow.io/docs/latest/cache-and-resume.html)를 참조하세요.

### 4.2. 과거 실행 로그 검사하기

nextflow 워크플로우를 실행할 때마다 현재 작업 디렉토리의 `.nextflow`라는 숨김 디렉토리 아래에 있는 `history`라는 로그 파일에 한 줄이 기록됩니다.

??? abstract "파일 내용"

    ```txt title=".nextflow/history" linenums="1"
    2025-07-04 19:27:09	1.8s	wise_watson	OK	3539118582ccde68dde471cc2c66295c	a02c9c46-c3c7-4085-9139-d1b9b5b194c8	nextflow run 1-hello.nf --input 'Hello World'
    2025-07-04 19:27:20	2.9s	spontaneous_blackwell	OK	3539118582ccde68dde471cc2c66295c	59a5db23-d83c-4c02-a54e-37ddb73a337e	nextflow run 1-hello.nf --input Bonjour
    2025-07-04 19:27:31	1.8s	gigantic_yonath	OK	3539118582ccde68dde471cc2c66295c	5acaa83a-6ad6-4509-bebc-cb25d5d7ddd0	nextflow run 1-hello.nf --input 'Dobry den'
    2025-07-04 19:27:45	2.4s	backstabbing_swartz	OK	3539118582ccde68dde471cc2c66295c	5f4b3269-5b53-404a-956c-cac915fbb74e	nextflow run 1-hello.nf --input Konnichiwa
    2025-07-04 19:27:57	2.1s	goofy_wilson	OK	3539118582ccde68dde471cc2c66295c	5f4b3269-5b53-404a-956c-cac915fbb74e	nextflow run 1-hello.nf --input Konnichiwa -resume
    ```

이 파일은 현재 작업 디렉토리 내에서 실행된 모든 Nextflow 실행에 대한 타임스탬프, 실행 이름, 상태, 리비전 ID, 세션 ID 및 전체 명령줄을 제공합니다.

이 정보에 액세스하는 더 편리한 방법은 [`nextflow log`](https://nextflow.io/docs/latest/reference/cli.html#log) 명령을 사용하는 것입니다.

```bash
nextflow log
```

??? success "명령 출력"

    ```console linenums="1"
    TIMESTAMP               DURATION        RUN NAME                STATUS  REVISION ID     SESSION ID                              COMMAND
    2025-07-04 19:27:09     1.8s            wise_watson             OK       3539118582     a02c9c46-c3c7-4085-9139-d1b9b5b194c8    nextflow run 1-hello.nf --input 'Hello World'
    2025-07-04 19:27:20     2.9s            spontaneous_blackwell   OK       3539118582     59a5db23-d83c-4c02-a54e-37ddb73a337e    nextflow run 1-hello.nf --input Bonjour
    2025-07-04 19:27:31     1.8s            gigantic_yonath         OK       3539118582     5acaa83a-6ad6-4509-bebc-cb25d5d7ddd0    nextflow run 1-hello.nf --input 'Dobry den'
    2025-07-04 19:27:45     2.4s            backstabbing_swartz     OK       3539118582     5f4b3269-5b53-404a-956c-cac915fbb74e    nextflow run 1-hello.nf --input Konnichiwa
    2025-07-04 19:27:57     2.1s            goofy_wilson            OK       3539118582     5f4b3269-5b53-404a-956c-cac915fbb74e    nextflow run 1-hello.nf --input Konnichiwa -resume
    ```

이렇게 하면 로그 파일의 내용이 헤더 줄과 함께 터미널에 출력됩니다.

새로운 `nextflow run` 명령을 실행할 때마다 세션 ID가 변경되는 것을 알 수 있습니다. 단, `-resume` 옵션을 사용하는 경우는 예외입니다.
이 경우 세션 ID는 동일하게 유지됩니다.

Nextflow는 세션 ID를 사용하여 `.nextflow` 아래에 있는 `cache` 디렉토리에 실행 캐싱 정보를 그룹화합니다.

### 4.3. 이전 작업 디렉토리 삭제하기

많은 파이프라인을 실행하면 많은 하위 디렉토리에 걸쳐 매우 많은 파일이 축적될 수 있습니다.
하위 디렉토리 이름이 무작위로 지정되므로 이름만으로는 이전 실행과 최근 실행을 구별하기 어렵습니다.

다행히 Nextflow에는 더 이상 관심이 없는 과거 실행의 작업 하위 디렉토리를 자동으로 삭제할 수 있는 [`nextflow clean`](https://www.nextflow.io/docs/latest/reference/cli.html#clean)이라는 유용한 명령이 포함되어 있습니다.

#### 4.3.1. 삭제 기준 결정하기

삭제할 항목을 결정하는 여러 옵션이 있으며, 위에 링크된 문서에서 탐색할 수 있습니다.
여기서는 실행 이름을 사용하여 지정된 특정 실행 이전의 모든 실행에서 모든 하위 디렉토리를 삭제하는 예제를 보여드립니다.

`-resume`을 사용하지 않은 가장 최근의 성공적인 실행을 찾아보세요. 우리의 경우 실행 이름은 `backstabbing_swartz`였습니다.

실행 이름은 `Launching (...)` 콘솔 출력 줄의 대괄호에 표시된 기계 생성 두 부분 문자열입니다.
Nextflow 로그를 사용하여 타임스탬프 및/또는 명령줄을 기반으로 실행을 조회할 수도 있습니다.

#### 4.3.2. 드라이 런 수행하기

먼저 드라이 런 플래그 `-n`을 사용하여 명령이 주어졌을 때 무엇이 삭제될지 확인합니다:

```bash
nextflow clean -before backstabbing_swartz -n
```

??? success "명령 출력"

    ```console
    Would remove /workspaces/training/hello-nextflow/work/eb/1a5de36637b475afd88fca7f79e024
    Would remove /workspaces/training/hello-nextflow/work/6b/19b0e002ea13486d3a0344c336c1d0
    Would remove /workspaces/training/hello-nextflow/work/45/9a6dd7ab771f93003d040956282883
    ```

출력에는 다른 작업 디렉토리 이름이 있고 줄 수가 다를 수 있지만 예제와 비슷해야 합니다.

출력된 줄이 없으면 유효한 실행 이름을 제공하지 않았거나 삭제할 과거 실행이 없는 것입니다. 예제 명령의 `backstabbing_swartz`를 로그에서 해당하는 최신 실행 이름으로 변경해야 합니다.

#### 4.3.3. 삭제 진행하기

출력이 예상대로 보이고 삭제를 진행하려면 `-n` 대신 `-f` 플래그를 사용하여 명령을 다시 실행하세요:

```bash
nextflow clean -before backstabbing_swartz -f
```

??? success "명령 출력"

    ```console
    Removed /workspaces/training/hello-nextflow/work/eb/1a5de36637b475afd88fca7f79e024
    Removed /workspaces/training/hello-nextflow/work/6b/19b0e002ea13486d3a0344c336c1d0
    Removed /workspaces/training/hello-nextflow/work/45/9a6dd7ab771f93003d040956282883
    ```

출력은 이전과 비슷해야 하지만 이제 'Would remove' 대신 'Removed'라고 표시됩니다.
이것은 두 문자 하위 디렉토리(위의 `eb/`와 같은)를 제거하지는 않지만 내용을 비웁니다.

!!! Warning

    과거 실행의 작업 하위 디렉토리를 삭제하면 Nextflow의 캐시에서 제거되고 해당 디렉토리에 저장된 모든 출력이 삭제됩니다.
    즉, 해당 프로세스를 다시 실행하지 않고 실행을 재개하는 Nextflow의 기능이 중단됩니다.

    관심 있는 출력을 저장하는 것은 사용자의 책임입니다! 이것이 `publish` 지시문에 `symlink` 모드보다 `copy` 모드를 사용하는 것을 선호하는 주된 이유입니다.

### 핵심 정리

이미 동일한 방식으로 실행된 단계를 반복하지 않고 파이프라인을 다시 실행하는 방법, 실행 로그를 검사하는 방법, `nextflow clean` 명령을 사용하여 이전 작업 디렉토리를 정리하는 방법을 배웠습니다.

### 다음 단계

잠시 휴식을 취하세요! Nextflow 구문과 기본 사용 지침의 구성 요소를 방금 습득했습니다.

이 교육의 다음 섹션에서는 Nextflow가 여러 입력을 효율적으로 처리하고, 여러 단계가 연결된 워크플로우를 실행하며, 모듈식 코드 구성 요소를 활용하고, 더 나은 재현성과 이식성을 위해 컨테이너를 활용할 수 있는 방법을 보여주는 Hello World 파이프라인의 점진적으로 더 현실적인 네 가지 버전을 살펴보겠습니다.

---

## 퀴즈

<quiz>
콘솔 출력 줄 `[a3/7be2fa] SAYHELLO | 1 of 1 ✔`에서 `[a3/7be2fa]`는 무엇을 나타냅니까?
- [ ] 프로세스 버전 번호
- [ ] 고유한 실행 식별자
- [x] 작업의 work 디렉토리에 대한 축약된 경로
- [ ] 출력 파일의 체크섬

자세히 알아보기: [2.3. `work/` 디렉토리에서 원본 출력 및 로그 찾기](#23-work-디렉토리에서-원본-출력-및-로그-찾기)
</quiz>

<quiz>
작업 디렉토리의 `.command.sh` 파일의 목적은 무엇입니까?
- [ ] 작업의 구성 설정을 저장합니다
- [x] 프로세스에서 실행된 실제 명령을 보여줍니다
- [ ] 실패한 작업의 오류 메시지를 포함합니다
- [ ] 작업을 위해 준비된 입력 파일을 나열합니다

자세히 알아보기: [2.3. `work/` 디렉토리에서 원본 출력 및 로그 찾기](#23-work-디렉토리에서-원본-출력-및-로그-찾기)
</quiz>

<quiz>
`-resume` 없이 워크플로우를 다시 실행하면 게시된 결과는 어떻게 됩니까?
- [ ] 별도의 타임스탬프가 지정된 디렉토리에 보존됩니다
- [x] 새 실행으로 덮어쓰여집니다
- [ ] Nextflow가 덮어쓰기를 방지하고 실패합니다
- [ ] 자동으로 백업됩니다

자세히 알아보기: [2.4. 다른 인사말로 워크플로우 다시 실행하기](#24-다른-인사말로-워크플로우-다시-실행하기)
</quiz>

<quiz>
이 콘솔 출력은 무엇을 나타냅니까?

```console
[skipped  ] process > sayHello (1) [100%] 1 of 1, cached: 1 ✔
```

- [ ] 작업이 실패하여 건너뛰었습니다
- [ ] 작업이 대기열에서 대기 중입니다
- [x] Nextflow가 이전의 동일한 실행 결과를 재사용했습니다
- [ ] 작업이 수동으로 취소되었습니다

자세히 알아보기: [4.1. `-resume`으로 워크플로우 다시 실행하기](#41--resume으로-워크플로우-다시-실행하기)
</quiz>

<quiz>
`nextflow log` 명령이 표시하는 실행 기록은 어디에 저장됩니까?
- [ ] results 디렉토리에
- [ ] work 디렉토리에
- [x] `.nextflow/history` 파일에
- [ ] `nextflow.config`에

자세히 알아보기: [4.2. 과거 실행 로그 검사하기](#42-과거-실행-로그-검사하기)
</quiz>

<quiz>
워크플로우 파일의 `params` 블록의 목적은 무엇입니까?
- [ ] 프로세스 리소스 요구 사항을 정의합니다
- [ ] executor를 구성합니다
- [x] 워크플로우 입력 매개변수를 선언하고 유형을 지정합니다
- [ ] 출력 게시 옵션을 지정합니다

자세히 알아보기: [3.4. 명령줄 매개변수의 `params` 시스템](#34-명령줄-매개변수의-params-시스템)
</quiz>

<quiz>
워크플로우의 `output` 블록에서 `mode 'copy'`는 무엇을 합니까?
- [ ] work 디렉토리의 백업을 생성합니다
- [x] 심볼릭 링크 대신 파일의 전체 복사본을 만듭니다
- [ ] 워크플로우 스크립트를 results에 복사합니다
- [ ] 증분 파일 복사를 활성화합니다

자세히 알아보기: [3.5. `publish` 지시문](#35-publish-지시문)
</quiz>

<quiz>
실제로 파일을 삭제하기 전에 `nextflow clean` 명령과 함께 사용하는 권장 플래그는 무엇입니까?
- [x] `-n` (드라이 런)으로 삭제될 항목을 미리 봅니다
- [ ] `-v` (상세)로 자세한 출력을 봅니다
- [ ] `-a` (모두)로 모든 디렉토리를 선택합니다
- [ ] `-q` (조용히)로 경고를 억제합니다

자세히 알아보기: [4.3. 이전 작업 디렉토리 삭제하기](#43-이전-작업-디렉토리-삭제하기)
</quiz>
