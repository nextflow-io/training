# 파트 1: 기본 작업 실행

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } AI 지원 번역 - [자세히 알아보기 및 개선 제안](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Nextflow Run 교육 과정의 첫 번째 파트에서는 매우 기본적인 도메인에 구애받지 않는 Hello World 예제로 주제를 시작하며, 이를 사용하여 필수 작업을 시연하고 해당 Nextflow 코드 구성 요소를 설명합니다.

??? info "Hello World 예제란 무엇인가요?"

    "Hello World!"는 프로그래밍 언어나 소프트웨어 프레임워크의 기본 구문과 구조를 보여주기 위한 최소한의 예제입니다.
    일반적으로 "Hello, World!"라는 문구를 콘솔이나 터미널과 같은 출력 장치에 출력하거나 파일에 쓰는 것으로 구성됩니다.

---

## 1. Hello World 직접 실행

Nextflow로 적용하기 전에 터미널에서 직접 실행하는 간단한 명령으로 이 개념을 시연하여 무엇을 하는지 보여드리겠습니다.

!!! tip

    [시작하기](00_orientation.md) 페이지에 설명된 대로 지금 `nextflow-run/` 디렉토리 안에 있어야 합니다.

### 1.1. 터미널에 hello 출력하기

터미널에서 다음 명령을 실행하세요.

```bash
echo 'Hello World!'
```

??? success "명령 출력"

    ```console
    Hello World!
    ```

이렇게 하면 터미널에 'Hello World' 텍스트가 바로 출력됩니다.

### 1.2. 출력을 파일에 쓰기

pipeline 실행은 대부분 파일에서 데이터를 읽고 다른 파일에 결과를 쓰는 것을 포함하므로, 예제를 좀 더 관련성 있게 만들기 위해 텍스트 출력을 파일에 쓰도록 명령을 수정해 보겠습니다.

```bash
echo 'Hello World!' > output.txt
```

??? success "명령 출력"

    ```console

    ```

이것은 터미널에 아무것도 출력하지 않습니다.

### 1.3. 출력 찾기

'Hello World' 텍스트는 이제 지정한 출력 파일인 `output.txt`에 있어야 합니다.
파일 탐색기에서 열거나 예를 들어 `cat` 유틸리티를 사용하여 명령줄에서 열 수 있습니다.

??? abstract "파일 내용"

    ```console title="output.txt" linenums="1"
    Hello World!
    ```

이것이 우리의 첫 번째 Nextflow workflow로 복제하려는 것입니다.

### 요약

이제 텍스트를 출력하는 간단한 명령을 터미널에서 실행하는 방법과 선택적으로 출력을 파일에 쓰는 방법을 알게 되었습니다.

### 다음 단계

동일한 결과를 달성하는 Nextflow workflow를 실행하는 데 무엇이 필요한지 알아봅니다.

---

## 2. Workflow 실행

`--input` 명령줄 인수를 통해 입력 인사말을 받아 해당 인사말이 포함된 텍스트 파일을 생성하는 `1-hello.nf`라는 workflow 스크립트를 제공합니다.

코드를 아직 살펴보지 않겠습니다. 먼저 실행하면 어떻게 되는지 살펴보겠습니다.

### 2.1. Workflow 시작 및 실행 모니터링

터미널에서 다음 명령을 실행하세요:

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

콘솔 출력이 이와 비슷하다면 축하합니다. 첫 번째 Nextflow workflow를 실행했습니다!

여기서 가장 중요한 출력은 위 출력에서 강조된 마지막 줄입니다:

```console
[a3/7be2fa] sayHello | 1 of 1 ✔
```

이것은 `sayHello` process가 한 번 성공적으로 실행되었음을 알려줍니다(`1 of 1 ✔`).

훌륭합니다. 하지만 출력은 어디에 있을까요?

### 2.2. `results` 디렉토리에서 출력 파일 찾기

이 workflow는 출력을 results 디렉토리에 게시하도록 구성되어 있습니다.
현재 디렉토리를 보면 workflow를 실행할 때 Nextflow가 `results`라는 새 디렉토리와 그 아래에 `1-hello`라는 하위 디렉토리를 생성하고 `output.txt`라는 파일을 포함하는 것을 볼 수 있습니다.

```console title="results/"
results
└── 1-hello
    └── output.txt
```

파일을 열어보세요. 내용은 명령줄에서 지정한 문자열과 일치해야 합니다.

```console title="results/1-hello/output.txt" linenums="1"
Hello World!
```

훌륭합니다. workflow가 해야 할 일을 했습니다!

그러나 '게시된' 결과는 Nextflow가 workflow를 실행할 때 생성한 실제 출력의 사본(또는 경우에 따라 심볼릭 링크)임을 알아두세요.

이제 Nextflow가 실제로 작업을 실행한 위치를 확인하기 위해 내부를 살펴보겠습니다.

!!! Warning

    모든 workflow가 results 디렉토리에 출력을 게시하도록 설정되어 있는 것은 아니며, 디렉토리 이름과 구조가 다를 수 있습니다.
    이 섹션의 조금 더 뒤에서 이 동작이 어디에 지정되어 있는지 확인하는 방법을 보여드리겠습니다.

### 2.3. `work/` 디렉토리에서 원본 출력 및 로그 찾기

workflow를 실행할 때 Nextflow는 workflow의 각 process 호출(=pipeline의 각 단계)에 대해 별도의 '작업 디렉토리'를 생성합니다.
각각에 대해 필요한 입력을 스테이징하고, 관련 명령을 실행하고, 출력과 로그 파일을 해당 디렉토리 내에 작성합니다. 이 디렉토리는 고유하게 만들기 위해 해시를 사용하여 자동으로 이름이 지정됩니다.

이러한 모든 작업 디렉토리는 현재 디렉토리(명령을 실행하는 위치) 내의 `work`라는 디렉토리 아래에 위치합니다.

혼란스럽게 들릴 수 있으므로 실제로 어떻게 보이는지 살펴보겠습니다.

앞서 실행한 workflow의 콘솔 출력으로 돌아가면 다음 줄이 있었습니다:

```console
[a3/7be2fa] sayHello | 1 of 1 ✔
```

줄이 `[a3/7be2fa]`로 시작하는 것을 보셨나요?
이것은 해당 process 호출의 작업 디렉토리 경로의 축약된 형태이며, `work/` 디렉토리 경로 내에서 `sayHello` process 호출의 출력을 찾을 위치를 알려줍니다.

다음 명령을 입력하고(자신의 터미널에 표시된 것으로 `a3/7be2fa`를 대체) Tab 키를 눌러 경로를 자동 완성하거나 별표를 추가하여 전체 경로를 찾을 수 있습니다:

```bash
ls work/a3/7be2fa*
```

이렇게 하면 전체 경로 디렉토리 경로가 표시됩니다: `work/a3/7be2fa7be2fad5e71e5f49998f795677fd68`

그 안에 무엇이 있는지 살펴보겠습니다.

??? abstract "디렉토리 내용"

    ```console
    work
    └── a3
        └── 7be2fad5e71e5f49998f795677fd68
            ├── .command.begin
            ├── .command.err
            ├── .command.log
            ├── .command.out
            ├── .command.run
            ├── .command.sh
            ├── .exitcode
            └── output.txt
    ```

??? question "같은 것이 보이지 않나요?"

    정확한 하위 디렉토리 이름은 시스템에서 다를 것입니다.

    VSCode 파일 탐색기에서 작업 하위 디렉토리의 내용을 찾아보면 모든 파일이 바로 표시됩니다.
    그러나 로그 파일은 터미널에서 보이지 않도록 설정되어 있으므로 `ls` 또는 `tree`를 사용하여 보려면 보이지 않는 파일을 표시하는 관련 옵션을 설정해야 합니다.

    ```bash
    tree -a work
    ```

`output.txt` 파일을 바로 알아볼 수 있을 것입니다. 이 파일은 실제로 `results` 디렉토리에 게시된 `sayHello` process의 원본 출력입니다.
열어보면 `Hello World!` 인사말을 다시 찾을 수 있습니다.

```console title="work/a3/7be2fa7be2fad5e71e5f49998f795677fd68/output.txt"
Hello World!
```

그렇다면 다른 모든 파일들은 무엇일까요?

이것들은 Nextflow가 작업 실행의 일부로 작성한 도우미 및 로그 파일입니다:

- **`.command.begin`**: 작업이 시작되자마자 생성되는 센티넬 파일.
- **`.command.err`**: process 호출에서 발생한 오류 메시지(`stderr`)
- **`.command.log`**: process 호출에서 발생한 전체 로그 출력
- **`.command.out`**: process 호출의 일반 출력(`stdout`)
- **`.command.run`**: Nextflow가 process 호출을 실행하기 위해 실행한 전체 스크립트
- **`.command.sh`**: process 호출에서 실제로 실행된 명령
- **`.exitcode`**: 명령에서 발생한 종료 코드

`.command.sh` 파일은 모든 부기 및 작업/환경 설정을 포함하지 않고 Nextflow가 실행한 기본 명령을 보여주기 때문에 특히 유용합니다.

```console title="work/a3/7be2fa7be2fad5e71e5f49998f795677fd68/command.sh"
#!/bin/bash -ue
echo 'Hello World!' > output.txt

```

이것은 workflow가 이전에 명령줄에서 직접 실행한 것과 동일한 명령을 구성했음을 확인합니다.

문제가 발생하여 무슨 일이 있었는지 해결해야 할 때, `command.sh` 스크립트를 보고 Nextflow가 workflow 지침, 변수 보간 등을 기반으로 어떤 명령을 구성했는지 정확히 확인하는 것이 유용할 수 있습니다.

### 2.4. 다른 인사말로 workflow 다시 실행

`--input` 인수에 다른 값을 사용하여 workflow를 몇 번 더 실행한 다음 작업 디렉토리를 살펴보세요.

??? abstract "디렉토리 내용"

    ```console
    work
    ├── 0f
    │   └── 52b7e07b0e274a80843fca48ed21b8
    │       ├── .command.begin
    │       ├── .command.err
    │       ├── .command.log
    │       ├── .command.out
    │       ├── .command.run
    │       ├── .command.sh
    │       ├── .exitcode
    │       └── output.txt
    ├── 67
    │   ├── 134e6317f90726c6c17ad53234a32b
    │   │   ├── .command.begin
    │   │   ├── .command.err
    │   │   ├── .command.log
    │   │   ├── .command.out
    │   │   ├── .command.run
    │   │   ├── .command.sh
    │   │   ├── .exitcode
    │   │   └── output.txt
    │   └── e029f2e75305874a9ab263d21ebc2c
    │       ├── .command.begin
    │       ├── .command.err
    │       ├── .command.log
    │       ├── .command.out
    │       ├── .command.run
    │       ├── .command.sh
    │       ├── .exitcode
    │       └── output.txt
    ├── 6c
    │   └── d4fd787e0b01b3c82e85696c297500
    │       ├── .command.begin
    │       ├── .command.err
    │       ├── .command.log
    │       ├── .command.out
    │       ├── .command.run
    │       ├── .command.sh
    │       ├── .exitcode
    │       └── output.txt
    └── e8
        └── ab99fad46ade52905ec973ff39bb80
            ├── .command.begin
            ├── .command.err
            ├── .command.log
            ├── .command.out
            ├── .command.run
            ├── .command.sh
            ├── .exitcode
            └── output.txt
    ```

각 실행에 대해 완전한 출력 및 로그 파일 세트가 포함된 새 하위 디렉토리가 생성된 것을 볼 수 있습니다.

반면에 `results` 디렉토리를 보면 여전히 결과 세트가 하나뿐이며 출력 파일의 내용은 마지막으로 실행한 것에 해당합니다.

??? abstract "디렉토리 내용"

    ```console title="results/"
    results
    └── 1-hello
        └── output.txt
    ```

이것은 게시된 결과가 후속 실행에 의해 덮어쓰여지는 반면, `work/` 아래의 작업 디렉토리는 보존된다는 것을 보여줍니다.

### 요약

간단한 Nextflow 스크립트를 실행하고 실행을 모니터링하고 출력을 찾는 방법을 알게 되었습니다.

### 다음 단계

기본 Nextflow 스크립트를 읽고 구성 요소가 기능과 어떻게 관련되는지 식별하는 방법을 배웁니다.

---

## 3. Hello World workflow 시작 스크립트 검토

방금 한 것은 기본적으로 workflow 스크립트를 블랙 박스로 취급한 것입니다.
이제 무엇을 하는지 보았으니 상자를 열고 내부를 살펴보겠습니다.

여기서 목표는 Nextflow 코드의 구문을 암기하는 것이 아니라 주요 구성 요소가 무엇이고 어떻게 구성되어 있는지에 대한 기본적인 직관을 형성하는 것입니다.

### 3.1. 전체 코드 구조 검토

`1-hello.nf` 스크립트는 현재 디렉토리인 `nextflow-run`에서 찾을 수 있습니다. 편집기 창에서 열어보세요.

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
        // 인사말 출력
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

Nextflow workflow 스크립트는 일반적으로 하나 이상의 **process** 정의, **workflow** 자체, 그리고 **params** 및 **output**과 같은 몇 가지 선택적 블록을 포함합니다.

각 **process**는 pipeline의 해당 단계가 수행해야 할 작업을 설명하고, **workflow**는 다양한 단계를 연결하는 데이터 흐름 로직을 설명합니다.

먼저 **process** 블록을 자세히 살펴본 다음 **workflow** 블록을 살펴보겠습니다.

### 3.2. `process` 정의

첫 번째 코드 블록은 **process**를 설명합니다.
process 정의는 `process` 키워드로 시작하고 그 뒤에 process 이름이 오고 마지막으로 중괄호로 구분된 process 본문이 옵니다.
process 본문에는 실행할 명령을 지정하는 script 블록이 포함되어야 하며, 이 명령은 명령줄 터미널에서 실행할 수 있는 모든 것이 될 수 있습니다.

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

여기서는 `greeting`이라는 **input** 변수를 받아 `output.txt`라는 파일에 **output**을 쓰는 `sayHello`라는 **process**가 있습니다.

<figure class="excalidraw">
--8<-- "docs/en/docs/nextflow_run/img/sayhello_with_input.svg"
</figure>

이것은 `input` 정의, `output` 정의 및 실행할 `script`만 포함하는 매우 최소한의 process 정의입니다.

`input` 정의에는 Nextflow에게 어떤 종류의 값(문자열, 숫자 등)을 기대하도록 알려주는 `val` 한정자가 포함되어 있습니다.

`output` 정의에는 이것이 경로로 처리되어야 함을 Nextflow에게 알려주는 `path` 한정자가 포함되어 있습니다(디렉토리 경로와 파일 모두 포함).

### 3.3. `workflow` 정의

두 번째 코드 블록은 **workflow** 자체를 설명합니다.
workflow 정의는 `workflow` 키워드로 시작하고 그 뒤에 선택적 이름이 오고, 그 다음 중괄호로 구분된 workflow 본문이 옵니다.

여기서는 `main:` 블록과 `publish:` 블록으로 구성된 **workflow**가 있습니다.
`main:` 블록은 workflow의 기본 본문이고 `publish:` 블록은 `results` 디렉토리에 게시되어야 하는 출력을 나열합니다.

```groovy title="1-hello.nf" linenums="27"
workflow {

    main:
    // 인사말 출력
    sayHello(params.input)

    publish:
    first_output = sayHello.out
}
```

이 경우 `main:` 블록에는 `sayHello` process 호출이 포함되어 있으며 인사말로 사용할 `params.input`이라는 입력을 제공합니다.

잠시 후 더 자세히 논의하겠지만, `params.input`은 명령줄에서 `--input` 매개변수에 제공한 값을 보유합니다.

`publish:` 블록은 `sayHello()` process 호출의 출력을 나열하며, 이를 `sayHello.out`으로 참조하고 `first_output`이라는 이름을 부여합니다(이것은 workflow 작성자가 원하는 대로 지정할 수 있습니다).

이것은 매우 최소한의 **workflow** 정의입니다.
실제 pipeline에서 workflow는 일반적으로 **channels**로 연결된 **processes**에 대한 여러 호출을 포함하며, 변수 입력에 대한 기본값이 설정될 수 있습니다.

이것은 과정의 Part 2에서 다룰 것입니다.
지금은 workflow가 입력과 출력을 어떻게 처리하는지 자세히 살펴보겠습니다.

### 3.4. 명령줄 매개변수의 `params` 시스템

`sayHello()` process 호출에 제공하는 `params.input`은 Nextflow 코드의 깔끔한 부분이며 잠시 더 살펴볼 가치가 있습니다.

위에서 언급했듯이, 이것이 `--input` 명령줄 매개변수의 값을 `sayHello()` process 호출에 전달하는 방법입니다.
실제로 `params.someParameterName`을 선언하는 것만으로 workflow에 명령줄에서 `--someParameterName`이라는 매개변수를 부여할 수 있습니다.

여기서는 workflow가 기대하는 입력 유형을 지정하는 `params` 블록을 설정하여 해당 매개변수 선언을 공식화했습니다(Nextflow 25.10.2 이상).

```groovy title="1-hello.nf" linenums="20"
/*
 * Pipeline parameters
 */
params {
    input: String
}
```

지원되는 유형에는 `String`, `Integer`, `Float`, `Boolean` 및 `Path`가 포함됩니다.

!!! tip

    `params` 시스템을 사용하여 선언된 workflow 매개변수는 항상 명령줄에서 두 개의 대시(`--`)를 사용합니다.
    이것은 하나의 대시(`-`)만 사용하는 Nextflow 수준 매개변수와 구분됩니다.

### 3.5. `publish` 지시문

workflow의 다른 쪽에서 `publish:` 블록을 이미 살펴보았습니다.
그것은 출력 처리 시스템의 절반이며, 나머지 절반은 아래에 있는 `output` 블록입니다.

```groovy title="1-hello.nf" linenums="37"
output {
    first_output {
        path '1-hello'
        mode 'copy'
    }
}
```

이것은 `publish:` 블록에 나열된 `first_output` 출력이 기본 `results` 출력 디렉토리 아래의 `1-hello`라는 하위 디렉토리에 복사되어야 함을 지정합니다.

`mode 'copy'` 줄은 `work/` 디렉토리의 원본 파일에 대한 심볼릭 링크(또는 symlink)를 만드는 대신 적절한 복사본을 만드는 시스템의 기본 동작을 재정의합니다.

게시 동작을 제어하기 위한 옵션이 여기에 표시된 것보다 더 많이 있습니다. 나중에 몇 가지를 다룰 것입니다.
또한 workflow가 여러 출력을 생성할 때 각각이 이 방식으로 `output` 블록에 나열되는 것을 볼 수 있습니다.

??? info "`publishDir`을 사용한 출력 게시의 이전 구문"

    매우 최근까지 출력을 게시하는 확립된 방법은 `publishDir` 지시문을 사용하여 각 개별 process 수준에서 수행하는 것이었습니다.

    이 코드 패턴은 여전히 이전 Nextflow pipeline 및 process 모듈 전체에서 찾을 수 있으므로 이를 알고 있는 것이 중요합니다.

    workflow의 `publish:` 블록과 최상위 수준의 `output` 블록 대신 `sayHello` process 정의에 `publishDir` 줄이 표시됩니다:

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

    그러나 향후 Nextflow 언어 버전에서 결국 허용되지 않을 것이므로 새 작업에서 이것을 사용하는 것은 권장하지 않습니다.

### 요약

이제 간단한 Nextflow workflow가 어떻게 구조화되어 있고 기본 구성 요소가 기능과 어떻게 관련되는지 알게 되었습니다.

### 다음 단계

workflow 실행을 편리하게 관리하는 방법을 배웁니다.

---

## 4. Workflow 실행 관리

workflow를 시작하고 출력을 검색하는 방법을 아는 것은 훌륭하지만, 작업을 더 쉽게 만들어 주는 workflow 관리의 몇 가지 다른 측면이 있습니다.

여기서는 동일한 workflow를 다시 시작해야 할 때 `resume` 기능을 활용하는 방법, `nextflow log`로 실행 로그를 검사하는 방법, `nextflow clean`으로 이전 작업 디렉토리를 삭제하는 방법을 보여드립니다.

### 4.1. `-resume`으로 workflow 다시 시작

때때로 이전에 이미 시작한 pipeline을 이미 성공적으로 완료된 작업을 다시 수행하지 않고 다시 실행하고 싶을 것입니다.

Nextflow에는 이를 수행할 수 있는 `-resume`이라는 옵션이 있습니다.
구체적으로 이 모드에서는 정확히 동일한 코드, 설정 및 입력으로 이미 실행된 모든 process가 건너뛰어집니다.
이것은 Nextflow가 마지막 실행 이후 추가하거나 수정한 process 또는 새 설정이나 입력을 제공하는 process만 실행한다는 것을 의미합니다.

이렇게 하면 두 가지 주요 이점이 있습니다:

- pipeline을 개발하는 중이라면 변경 사항을 테스트하기 위해 현재 작업 중인 process만 실행하면 되므로 더 빠르게 반복할 수 있습니다.
- 프로덕션에서 pipeline을 실행 중이고 문제가 발생하면 많은 경우 문제를 수정하고 pipeline을 다시 시작할 수 있으며, 실패 지점부터 실행을 재개하여 많은 시간과 컴퓨팅을 절약할 수 있습니다.

사용하려면 명령에 `-resume`을 추가하고 실행하면 됩니다:

```bash
nextflow run 1-hello.nf --input 'Hello World!' -resume
```

??? success "명령 출력"

    ```console linenums="1"
    N E X T F L O W   ~  version 25.10.2

    Launching `1-hello.nf` [tiny_noyce] DSL2 - revision: c33d41f479

    [a3/7be2fa] sayHello | 1 of 1, cached: 1 ✔
    ```

콘솔 출력은 익숙해 보이지만 이전과 약간 다른 점이 있습니다.

process 상태 줄(5번째 줄)에 추가된 `cached:` 부분을 찾아보세요. 이것은 Nextflow가 이 작업을 이미 수행했음을 인식하고 이전 성공적인 실행의 결과를 단순히 재사용했음을 의미합니다.

작업 하위 디렉토리 해시가 이전 실행과 동일한 것도 볼 수 있습니다.
Nextflow는 문자 그대로 이전 실행을 가리키며 "저기서 이미 했어요"라고 말하는 것입니다.

!!! tip

    `resume`으로 pipeline을 다시 실행할 때 Nextflow는 이전에 성공적으로 실행된 실행에 의해 작업 디렉토리 외부에 게시된 파일을 덮어쓰지 않습니다.

### 4.2. 과거 실행 로그 검사

nextflow workflow를 시작할 때마다 현재 작업 디렉토리의 `.nextflow`라는 숨겨진 디렉토리 아래에 `history`라는 로그 파일에 줄이 기록됩니다.

??? abstract "파일 내용"

    ```txt title=".nextflow/history" linenums="1"
    2025-07-04 19:27:09	1.8s	wise_watson	OK	3539118582ccde68dde471cc2c66295c	a02c9c46-c3c7-4085-9139-d1b9b5b194c8	nextflow run 1-hello.nf --input 'Hello World'
    2025-07-04 19:27:20	2.9s	spontaneous_blackwell	OK	3539118582ccde68dde471cc2c66295c	59a5db23-d83c-4c02-a54e-37ddb73a337e	nextflow run 1-hello.nf --input Bonjour
    2025-07-04 19:27:31	1.8s	gigantic_yonath	OK	3539118582ccde68dde471cc2c66295c	5acaa83a-6ad6-4509-bebc-cb25d5d7ddd0	nextflow run 1-hello.nf --input 'Dobry den'
    2025-07-04 19:27:45	2.4s	backstabbing_swartz	OK	3539118582ccde68dde471cc2c66295c	5f4b3269-5b53-404a-956c-cac915fbb74e	nextflow run 1-hello.nf --input Konnichiwa
    2025-07-04 19:27:57	2.1s	goofy_wilson	OK	3539118582ccde68dde471cc2c66295c	5f4b3269-5b53-404a-956c-cac915fbb74e	nextflow run 1-hello.nf --input Konnichiwa -resume
    ```

이 파일은 현재 작업 디렉토리 내에서 시작된 모든 Nextflow 실행에 대한 타임스탬프, 실행 이름, 상태, 리비전 ID, 세션 ID 및 전체 명령줄을 제공합니다.

이 정보에 액세스하는 더 편리한 방법은 `nextflow log` 명령을 사용하는 것입니다.

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

이렇게 하면 헤더 줄이 추가된 로그 파일의 내용이 터미널에 출력됩니다.

새 `nextflow run` 명령을 실행할 때마다 세션 ID가 변경되지만, `-resume` 옵션을 사용하는 경우에는 예외입니다.
그 경우 세션 ID가 동일하게 유지됩니다.

Nextflow는 세션 ID를 사용하여 `.nextflow` 아래에 있는 `cache` 디렉토리 아래에 실행 캐싱 정보를 그룹화합니다.

### 4.3. 이전 작업 디렉토리 삭제

많은 pipeline을 실행하면 많은 하위 디렉토리에 매우 많은 파일이 축적될 수 있습니다.
하위 디렉토리는 무작위로 이름이 지정되므로 이름만으로는 이전 실행과 최근 실행을 구분하기 어렵습니다.

다행히 Nextflow에는 더 이상 관심이 없는 과거 실행에 대한 작업 하위 디렉토리를 자동으로 삭제할 수 있는 유용한 `clean` 하위 명령이 포함되어 있습니다.

#### 4.3.1. 삭제 기준 결정

삭제할 항목을 결정하는 여러 [옵션](https://www.nextflow.io/docs/latest/reference/cli.html#clean)이 있습니다.

여기서는 실행 이름을 사용하여 지정된 실행 이전의 모든 하위 디렉토리를 삭제하는 예를 보여드립니다.

`-resume`을 사용하지 않은 가장 최근의 성공적인 실행을 찾아보세요. 우리의 경우 실행 이름은 `backstabbing_swartz`였습니다.

실행 이름은 `Launching (...)` 콘솔 출력 줄에서 대괄호 안에 표시되는 기계 생성 두 부분 문자열입니다.
Nextflow 로그를 사용하여 타임스탬프 및/또는 명령줄을 기반으로 실행을 찾을 수도 있습니다.

#### 4.3.2. 드라이 런 수행

먼저 드라이 런 플래그 `-n`을 사용하여 명령에 따라 삭제될 항목을 확인합니다:

```bash
nextflow clean -before backstabbing_swartz -n
```

??? success "명령 출력"

    ```console
    Would remove /workspaces/training/hello-nextflow/work/eb/1a5de36637b475afd88fca7f79e024
    Would remove /workspaces/training/hello-nextflow/work/6b/19b0e002ea13486d3a0344c336c1d0
    Would remove /workspaces/training/hello-nextflow/work/45/9a6dd7ab771f93003d040956282883
    ```

출력에는 다른 작업 디렉토리 이름이 있고 줄 수가 다를 수 있지만 예제와 비슷해 보여야 합니다.

줄이 출력되지 않으면 유효한 실행 이름을 제공하지 않았거나 삭제할 과거 실행이 없는 것입니다. 예제 명령의 `backstabbing_swartz`를 로그의 해당 최신 실행 이름으로 변경해야 합니다.

#### 4.3.3. 삭제 진행

출력이 예상대로 보이고 삭제를 진행하려면 `-n` 대신 `-f` 플래그로 명령을 다시 실행하세요:

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
이것은 두 문자 하위 디렉토리(위의 `eb/`와 같은)를 제거하지 않지만 그 내용을 비웁니다.

!!! Warning

    과거 실행에서 작업 하위 디렉토리를 삭제하면 Nextflow의 캐시에서 제거되고 해당 디렉토리에 저장된 모든 출력이 삭제됩니다.
    이것은 해당 process를 다시 실행하지 않고 실행을 재개하는 Nextflow의 기능을 손상시킵니다.

    관심 있는 출력을 저장하는 것은 사용자의 책임입니다! 이것이 `publish` 지시문에 `symlink` 모드보다 `copy` 모드를 선호하는 주된 이유입니다.

### 요약

이미 동일한 방식으로 실행된 단계를 반복하지 않고 pipeline을 다시 시작하는 방법, 실행 로그를 검사하는 방법, `nextflow clean` 명령을 사용하여 이전 작업 디렉토리를 정리하는 방법을 알게 되었습니다.

### 다음 단계

잠시 휴식을 취하세요! Nextflow 구문과 기본 사용 지침의 구성 요소를 방금 흡수했습니다.

교육의 다음 섹션에서는 Hello World pipeline의 점점 더 현실적인 네 가지 버전을 살펴보며 Nextflow가 여러 입력을 효율적으로 처리하고, 함께 연결된 여러 단계로 구성된 workflow를 실행하고, 모듈식 코드 구성 요소를 활용하고, 더 큰 재현성과 이식성을 위해 컨테이너를 활용하는 방법을 시연합니다.

---

## 퀴즈

<quiz>
콘솔 출력 줄 `[a3/7be2fa] SAYHELLO | 1 of 1 ✔`에서 `[a3/7be2fa]`는 무엇을 나타내나요?
- [ ] process 버전 번호
- [ ] 고유 실행 식별자
- [x] 작업의 작업 디렉토리에 대한 축약된 경로
- [ ] 출력 파일의 체크섬

자세히 알아보기: [2.3. `work/` 디렉토리에서 원본 출력 및 로그 찾기](#23-work-디렉토리에서-원본-출력-및-로그-찾기)
</quiz>

<quiz>
작업 디렉토리의 `.command.sh` 파일의 목적은 무엇인가요?
- [ ] 작업의 구성 설정을 저장합니다
- [x] process에서 실행된 실제 명령을 보여줍니다
- [ ] 실패한 작업의 오류 메시지를 포함합니다
- [ ] 작업을 위해 스테이징된 입력 파일을 나열합니다

자세히 알아보기: [2.3. `work/` 디렉토리에서 원본 출력 및 로그 찾기](#23-work-디렉토리에서-원본-출력-및-로그-찾기)
</quiz>

<quiz>
`-resume` 없이 workflow를 다시 실행하면 게시된 결과는 어떻게 되나요?
- [ ] 별도의 타임스탬프가 있는 디렉토리에 보존됩니다
- [x] 새 실행에 의해 덮어쓰여집니다
- [ ] Nextflow가 덮어쓰기를 방지하고 실패합니다
- [ ] 자동으로 백업됩니다

자세히 알아보기: [2.4. 다른 인사말로 workflow 다시 실행](#24-다른-인사말로-workflow-다시-실행)
</quiz>

<quiz>
이 콘솔 출력은 무엇을 나타내나요?

```console
[skipped  ] process > sayHello (1) [100%] 1 of 1, cached: 1 ✔
```

- [ ] 작업이 실패하여 건너뛰었습니다
- [ ] 작업이 대기열에서 대기 중입니다
- [x] Nextflow가 이전 동일한 실행의 결과를 재사용했습니다
- [ ] 작업이 수동으로 취소되었습니다

자세히 알아보기: [4.1. `-resume`으로 workflow 다시 시작](#41--resume으로-workflow-다시-시작)
</quiz>

<quiz>
`nextflow log` 명령이 표시하는 실행 히스토리는 어디에 저장되나요?
- [ ] results 디렉토리에
- [ ] work 디렉토리에
- [x] `.nextflow/history` 파일에
- [ ] `nextflow.config`에

자세히 알아보기: [4.2. 과거 실행 로그 검사](#42-과거-실행-로그-검사)
</quiz>

<quiz>
workflow 파일에서 `params` 블록의 목적은 무엇인가요?
- [ ] process 리소스 요구 사항을 정의합니다
- [ ] executor를 구성합니다
- [x] workflow 입력 매개변수를 선언하고 유형을 지정합니다
- [ ] 출력 게시 옵션을 지정합니다

자세히 알아보기: [3.4. 명령줄 매개변수의 params 시스템](#34-명령줄-매개변수의-params-시스템)
</quiz>

<quiz>
workflow의 `output` 블록에서 `mode 'copy'`는 무엇을 하나요?
- [ ] 작업 디렉토리의 백업을 생성합니다
- [x] 심볼릭 링크 대신 파일의 전체 복사본을 만듭니다
- [ ] workflow 스크립트를 results에 복사합니다
- [ ] 증분 파일 복사를 활성화합니다

자세히 알아보기: [3.5. publish 지시문](#35-publish-지시문)
</quiz>

<quiz>
파일을 실제로 삭제하기 전에 `nextflow clean` 명령과 함께 사용하는 것이 권장되는 플래그는 무엇인가요?
- [x] `-n` (드라이 런)으로 삭제될 항목을 미리 봅니다
- [ ] `-v` (상세)로 자세한 출력을 봅니다
- [ ] `-a` (전체)로 모든 디렉토리를 선택합니다
- [ ] `-q` (조용)으로 경고를 숨깁니다

자세히 알아보기: [4.3. 이전 작업 디렉토리 삭제](#43-이전-작업-디렉토리-삭제)
</quiz>
