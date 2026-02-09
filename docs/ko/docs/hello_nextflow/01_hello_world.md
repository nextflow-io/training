# 파트 1: Hello World

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/tOukLxWCHiA?si=y8lAedhEHWaTV4zd&amp;list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&amp;cc_load_policy=1" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

/// caption
:fontawesome-brands-youtube:{ .youtube } [전체 재생목록](https://youtube.com/playlist?list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&si=eF7cLR62goy-lc6n)은 Nextflow YouTube 채널에서 확인하세요.

:green_book: 비디오 스크립트는 [여기](./transcripts/01_hello_world.md)에서 확인할 수 있습니다.
///

Hello Nextflow 교육 과정의 첫 번째 파트에서는 매우 기본적이고 도메인에 구애받지 않는 Hello World 예제로 시작합니다. 이 예제를 단계적으로 확장하면서 Nextflow의 기본 로직과 구성 요소 사용법을 학습합니다.

??? info "Hello World 예제란 무엇인가요?"

    "Hello World!"는 프로그래밍 언어나 소프트웨어 프레임워크의 기본 구문과 구조를 보여주기 위한 최소한의 예제입니다.
    일반적으로 콘솔이나 터미널과 같은 출력 장치에 "Hello, World!" 문구를 출력하거나 파일에 작성하는 것으로 구성됩니다.

---

## 0. 준비 운동: Hello World 예제를 직접 실행하기

Nextflow로 감싸기 전에 터미널에서 직접 간단한 명령을 실행하여 어떻게 작동하는지 보여드리겠습니다.

!!! tip

    [시작하기](00_orientation.md) 페이지에 설명된 대로 `hello-nextflow/` 디렉토리 안에 있어야 합니다.

### 0.1. 터미널에서 인사말 출력하기

터미널에서 다음 명령을 실행하세요.

```bash
echo 'Hello World!'
```

??? success "명령 출력"

    ```console
    Hello World!
    ```

이 명령은 터미널에 바로 'Hello World' 텍스트를 출력합니다.

### 0.2. 출력을 파일에 작성하기

파이프라인 실행은 주로 파일에서 데이터를 읽고 결과를 다른 파일에 작성하는 것과 관련이 있으므로, 예제를 좀 더 실용적으로 만들기 위해 텍스트 출력을 파일에 작성하도록 명령을 수정하겠습니다.

```bash
echo 'Hello World!' > output.txt
```

??? success "명령 출력"

    ```console

    ```

이 명령은 터미널에 아무것도 출력하지 않습니다.

### 0.3. 출력 찾기

'Hello World' 텍스트는 이제 우리가 지정한 `output.txt`라는 출력 파일에 있어야 합니다.
파일 탐색기에서 열거나 명령줄에서 `cat` 유틸리티 등을 사용하여 열 수 있습니다.

??? abstract "파일 내용"

    ```console title="output.txt" linenums="1"
    Hello World!
    ```

이것이 우리가 첫 번째 Nextflow 워크플로우로 재현하려는 것입니다.

### 핵심 정리

터미널에서 텍스트를 출력하는 간단한 명령을 실행하는 방법과 선택적으로 출력을 파일에 작성하는 방법을 배웠습니다.

### 다음 단계

Nextflow 워크플로우로 작성하면 어떻게 보이는지 알아보세요.

---

## 1. 스크립트 살펴보고 실행하기

이전과 동일한 작업(즉, 'Hello World!' 작성)을 수행하지만 Nextflow를 사용하는 완전히 작동하는 최소한의 워크플로우 스크립트 `hello-world.nf`를 제공합니다.

시작하기 위해 워크플로우 스크립트를 열어 구조를 파악해 보겠습니다.
그런 다음 실행하고 출력을 찾아보겠습니다.

### 1.1. 코드 살펴보기

현재 디렉토리인 `hello-nextflow`에서 `hello-world.nf` 스크립트를 찾을 수 있습니다. 편집기 창에서 열어보세요.

??? full-code "전체 코드 파일"

    ```groovy title="hello-world.nf" linenums="1"
    #!/usr/bin/env nextflow

    /*
    * echo를 사용하여 'Hello World!'를 파일에 출력합니다
    */
    process sayHello {

        output:
        path 'output.txt'

        script:
        """
        echo 'Hello World!' > output.txt
        """
    }

    workflow {

        main:
        // emit a greeting
        sayHello()
    }
    ```

Nextflow 워크플로우 스크립트는 일반적으로 하나 이상의 [**process**](https://nextflow.io/docs/latest/process.html) 정의와 [**workflow**](https://nextflow.io/docs/latest/workflow.html) 자체, 그리고 몇 가지 선택적 블록(여기에는 없음)을 포함하며, 이는 나중에 소개하겠습니다.

각 **process**는 파이프라인의 해당 단계가 수행해야 하는 작업을 설명하고, **workflow**는 다양한 단계를 연결하는 데이터 흐름 로직을 설명합니다.

먼저 **process** 블록을 자세히 살펴본 다음 **workflow** 블록을 살펴보겠습니다.

#### 1.1.1. `process` 정의

첫 번째 코드 블록은 **process**를 설명합니다.

프로세스 정의는 `process` 키워드로 시작하고, 프로세스 이름이 이어지며, 마지막으로 중괄호로 구분된 프로세스 본문이 옵니다.
프로세스 본문에는 실행할 명령을 지정하는 script 블록이 포함되어야 하며, 이는 명령줄 터미널에서 실행할 수 있는 모든 것이 될 수 있습니다.

```groovy title="hello-world.nf" linenums="3"
/*
* echo를 사용하여 'Hello World!'를 파일에 출력합니다
*/
process sayHello {

    output:
    path 'output.txt'

    script:
    """
    echo 'Hello World!' > output.txt
    """
}
```

여기에는 `output.txt`라는 파일에 **output**을 작성하는 `sayHello`라는 **process**가 있습니다.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello_world.svg"
</figure>

이것은 `output` 정의와 실행할 `script`만 포함하는 매우 최소한의 프로세스 정의입니다.

`output` 정의에는 `path` 한정자가 포함되어 있으며, 이는 Nextflow에게 이것을 경로로 처리해야 한다고 알려줍니다(디렉토리 경로와 파일 모두 포함).
또 다른 일반적인 한정자는 `val`입니다.

중요한 점은 output 정의가 어떤 출력이 생성될지를 _결정하지_ 않는다는 것입니다.
단순히 예상되는 출력이 무엇인지를 _선언_하여 Nextflow가 실행이 완료되면 이를 찾을 수 있도록 합니다.
이는 명령이 성공적으로 실행되었는지 확인하고 필요한 경우 출력을 다운스트림 프로세스로 전달하는 데 필요합니다. output 블록에 선언된 것과 일치하지 않는 생성된 출력은 다운스트림 프로세스로 전달되지 않습니다.

!!! warning

    이 예제는 출력 파일명을 두 곳(script와 output 블록)에 하드코딩했기 때문에 취약합니다.
    하나를 변경하고 다른 하나를 변경하지 않으면 스크립트가 중단됩니다.
    나중에 이 문제를 완화하기 위해 변수를 사용하는 방법을 배우게 됩니다.

실제 파이프라인에서 프로세스는 일반적으로 지시문 및 입력과 같은 추가 블록을 포함하며, 이는 곧 소개하겠습니다.

#### 1.1.2. `workflow` 정의

두 번째 코드 블록은 **workflow** 자체를 설명합니다.
워크플로우 정의는 `workflow` 키워드로 시작하고, 선택적 이름이 이어지며, 중괄호로 구분된 워크플로우 본문이 옵니다.

여기에는 `sayHello` 프로세스 호출을 포함하는 `main:` 블록('이것이 워크플로우의 메인 본문입니다'라는 의미)으로 구성된 **workflow**가 있습니다.

```groovy title="hello-world.nf" linenums="17"
workflow {

    main:
    // emit a greeting
    sayHello()
}
```

이것은 매우 최소한의 **workflow** 정의입니다.
실제 파이프라인에서 워크플로우는 일반적으로 **채널**로 연결된 여러 **프로세스** 호출을 포함하며, 프로세스는 하나 이상의 변수 **입력**을 예상합니다.

이 교육 모듈의 뒷부분에서 변수 입력을 추가하는 방법을 배우게 되며, 이 과정의 파트 3에서 더 많은 프로세스를 추가하고 채널로 연결하는 방법을 배우게 됩니다.

!!! tip

    기술적으로 이와 같은 간단한 워크플로우에는 `main:` 줄이 필요하지 않으므로 이것이 없는 워크플로우를 접할 수 있습니다.
    하지만 워크플로우 수준 출력을 활용하려면 필요하므로 처음부터 포함하는 것이 좋습니다.

### 1.2. 워크플로우 실행하기

코드를 보는 것만큼 실행하는 것이 재미있지 않으므로 실제로 시도해 보겠습니다.

#### 1.2.1. 워크플로우 시작 및 실행 모니터링

터미널에서 다음 명령을 실행하세요:

```bash
nextflow run hello-world.nf
```

??? success "명령 출력"

    ```console hl_lines="7"
    N E X T F L O W   ~  version 25.10.2

    Launching `hello-world.nf` [goofy_torvalds] DSL2 - revision: c33d41f479

    executor >  local (1)
    [65/7be2fa] sayHello | 1 of 1 ✔
    ```

콘솔 출력이 이와 비슷하다면 축하합니다. 첫 번째 Nextflow 워크플로우를 실행했습니다!

여기서 가장 중요한 출력은 위에서 강조 표시된 마지막 줄입니다:

```console
[65/7be2fa] sayHello | 1 of 1 ✔
```

이것은 `sayHello` 프로세스가 한 번 성공적으로 실행되었음을 알려줍니다(`1 of 1 ✔`).

중요한 점은 이 줄이 `sayHello` 프로세스 호출의 출력을 찾을 수 있는 위치도 알려준다는 것입니다.
이제 그것을 살펴보겠습니다.

#### 1.2.2. `work` 디렉토리에서 출력 및 로그 찾기

특정 디렉토리에서 Nextflow를 처음 실행하면 `work`라는 디렉토리가 생성되며, 여기에 실행 과정에서 생성된 모든 파일(및 모든 심볼릭 링크)이 작성됩니다.

`work` 디렉토리 내에서 Nextflow는 프로세스 호출별로 출력과 로그를 구성합니다.
각 프로세스 호출에 대해 Nextflow는 고유하게 만들기 위해 해시로 이름이 지정된 중첩된 하위 디렉토리를 생성하며, 여기에 필요한 모든 입력을 준비하고(기본적으로 심볼릭 링크 사용), 헬퍼 파일을 작성하고, 로그와 프로세스의 모든 출력을 작성합니다.

해당 하위 디렉토리의 경로는 콘솔 출력에 축약된 형태로 대괄호 안에 표시됩니다.
위에 표시된 실행에서 얻은 것을 보면 sayHello 프로세스의 콘솔 로그 줄이 `[65/7be2fa]`로 시작합니다. 이는 다음 디렉토리 경로에 해당합니다: `work/65/7be2fad5e71e5f49998f795677fd68`

그 안에 무엇이 있는지 살펴보겠습니다.

??? abstract "디렉토리 내용"

    ```console
    work
    └── 65
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

    정확한 하위 디렉토리 이름은 시스템마다 다릅니다.

    VSCode 파일 탐색기에서 작업 하위 디렉토리의 내용을 탐색하면 모든 파일을 바로 볼 수 있습니다.
    그러나 로그 파일은 터미널에서 보이지 않도록 설정되어 있으므로 `ls` 또는 `tree`를 사용하여 보려면 보이지 않는 파일을 표시하는 관련 옵션을 설정해야 합니다.

    ```bash
    tree -a work
    ```

가장 먼저 살펴볼 것은 워크플로우의 실제 출력, 즉 `sayHello` 프로세스에서 생성된 `output.txt` 파일입니다.
열어보면 최소한의 워크플로우의 목적이었던 `Hello World!` 인사말을 찾을 수 있습니다.

??? abstract "파일 내용"

    ```console title="output.txt"
    Hello World!
    ```

작동했습니다!

물론 이렇게 작은 결과에 대해 많은 래퍼 코드처럼 보일 수 있지만, 입력 파일을 읽고 여러 단계를 연결하기 시작하면 모든 래퍼 코드의 가치가 더 명확해질 것입니다.

그렇긴 하지만 해당 디렉토리의 다른 파일도 살펴보겠습니다. 이것들은 Nextflow가 작업 실행의 일부로 생성한 헬퍼 및 로그 파일입니다.

- **`.command.begin`**: 프로세스 호출 실행 시작과 관련된 메타데이터
- **`.command.err`**: 프로세스 호출에서 발생한 오류 메시지(`stderr`)
- **`.command.log`**: 프로세스 호출에서 발생한 전체 로그 출력
- **`.command.out`**: 프로세스 호출의 일반 출력(`stdout`)
- **`.command.run`**: 프로세스 호출을 실행하기 위해 Nextflow가 실행한 전체 스크립트
- **`.command.sh`**: 프로세스 호출에서 실제로 실행된 명령
- **`.exitcode`**: 명령에서 발생한 종료 코드

`.command.sh` 파일은 Nextflow가 실행한 주요 명령을 알려주기 때문에 특히 유용하며, 모든 부기 및 작업/환경 설정은 포함하지 않습니다.

??? abstract "파일 내용"

    ```console title=".command.sh"
    #!/bin/bash -ue
    echo 'Hello World!' > output.txt
    ```

이것은 우리가 이전에 수동으로 실행한 것과 일치합니다.

이 경우 프로세스 명령이 하드코딩되어 있어서 매우 간단하지만, 과정 후반부에서 일부 변수 보간을 포함하는 프로세스 명령을 보게 될 것입니다.
실패한 실행 문제를 해결할 때 Nextflow가 코드를 어떻게 해석했는지, 어떤 명령이 생성되었는지 정확히 볼 수 있다는 것은 특히 유용합니다.

### 1.3. 워크플로우 다시 실행하기

워크플로우를 몇 번 다시 실행한 다음 `work/` 아래의 작업 디렉토리를 살펴보세요.

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
    ├── 65
        └── 7be2fad5e71e5f49998f795677fd68
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

각 실행에 대해 전체 출력 및 로그 파일 세트가 포함된 새 하위 디렉토리가 생성된 것을 볼 수 있습니다.
이것은 동일한 워크플로우를 여러 번 실행해도 이전 실행의 결과를 덮어쓰지 않는다는 것을 보여줍니다.

### 핵심 정리

간단한 Nextflow 스크립트를 해독하고 실행하며 work 디렉토리에서 출력 및 관련 로그 파일을 찾는 방법을 배웠습니다.

### 다음 단계

워크플로우 출력을 더 편리한 위치에 게시하는 방법을 배우세요.

---

## 2. 출력 게시하기

방금 배운 것처럼 파이프라인에서 생성된 출력은 여러 계층 깊이의 작업 디렉토리에 묻혀 있습니다.
이것은 의도적으로 수행됩니다. Nextflow가 이 디렉토리를 제어하며 우리는 이와 상호 작용하지 않아야 합니다.
그러나 이로 인해 우리가 관심 있는 출력을 검색하기가 불편합니다.

다행히 Nextflow는 [워크플로우 출력 정의](https://nextflow.io/docs/latest/workflow.html#workflow-outputs)를 사용하여 지정된 디렉토리에 출력을 게시하는 방법을 제공합니다.

### 2.1. 기본 사용법

이것은 두 가지 새로운 코드 조각을 포함합니다:

1. `workflow` 본문 내부의 `publish:` 블록으로 프로세스 출력을 선언합니다.
2. 모드 및 위치와 같은 출력 옵션을 지정하는 스크립트에 `output` 블록을 추가합니다.

#### 2.1.1. `sayHello` 프로세스의 출력 선언하기

워크플로우 본문에 `publish:` 블록(`main:` 블록과 동일한 종류의 코드 요소)을 추가하고 `sayHello()` 프로세스의 출력을 나열해야 합니다.

워크플로우 스크립트 파일 `hello-world.nf`에 다음 코드 줄을 추가하세요:

=== "후"

    ```groovy title="hello-world.nf" linenums="17" hl_lines="7-8"
    workflow {

        main:
        // emit a greeting
        sayHello()

        publish:
        first_output = sayHello.out
    }
    ```

=== "전"

    ```groovy title="hello-world.nf" linenums="17"
    workflow {

        main:
        // emit a greeting
        sayHello()
    }
    ```

`sayHello().out`을 수행하여 프로세스의 출력을 참조하고 임의의 이름인 `first_output`을 할당할 수 있습니다.

#### 2.1.2. 스크립트에 `output:` 블록 추가하기

이제 출력 디렉토리 경로가 지정될 `output:` 블록을 추가하기만 하면 됩니다. 이 새 블록은 스크립트 내에서 `workflow` 블록 **외부** 및 **아래**에 위치합니다.

워크플로우 스크립트 파일 `hello-world.nf`에 다음 코드 줄을 추가하세요:

=== "후"

    ```groovy title="hello-world.nf" linenums="17" hl_lines="11-15"
    workflow {

        main:
        // emit a greeting
        sayHello()

        publish:
        first_output = sayHello.out
    }

    output {
        first_output {
            path '.'
        }
    }
    ```

=== "전"

    ```groovy title="hello-world.nf" linenums="17"
    workflow {

        main:
        // emit a greeting
        sayHello()

        publish:
        first_output = sayHello.out
    }
    ```

이것을 사용하여 `workflow` 블록에 선언된 모든 프로세스 출력에 특정 경로를 할당할 수 있습니다.
나중에 정교한 출력 디렉토리 구조를 생성하는 방법을 배우게 되지만, 지금은 단순성을 위해 최소한의 경로만 하드코딩하고 있습니다.

#### 2.1.3. 워크플로우 실행하기

이제 수정된 워크플로우 스크립트를 실행하세요:

```bash
nextflow run hello-world.nf
```

??? success "명령 출력"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `hello-world.nf` [jovial_mayer] DSL2 - revision: 35bd3425e5

    executor >  local (1)
    [9f/48ef97] sayHello | 1 of 1 ✔
    ```

터미널 출력은 익숙해 보일 것입니다. 외부적으로는 아무것도 변경되지 않았습니다.

그러나 파일 탐색기를 확인하세요: 이번에는 Nextflow가 `results/`라는 새 디렉토리를 생성했습니다.

??? abstract "디렉토리 내용"

    ```console hl_lines="10-11 22"
    .
    ├── greetings.csv
    ├── hello-channels.nf
    ├── hello-config.nf
    ├── hello-containers.nf
    ├── hello-modules.nf
    ├── hello-workflow.nf
    ├── hello-world.nf
    ├── nextflow.config
    ├── results
    │   └── output.txt -> /workspaces/training/hello-nextflow/work/9f/48ef97f110b0dbd83635d7cbe288d2/output.txt
    ├── solutions
    │   ├── 1-hello-world
    │   ├── 2-hello-channels
    │   ├── 3-hello-workflow
    │   ├── 4-hello-modules
    │   ├── 5-hello-containers
    │   └── 6-hello-config
    ├── test-params.json
    └── work
        ├── 65
        └── 9f
    ```

`results` 디렉토리 내부에서 방금 실행한 명령으로 work 디렉토리에 생성된 `output.txt`에 대한 심볼릭 링크를 찾을 수 있습니다.

이를 통해 work 하위 디렉토리를 파헤칠 필요 없이 출력 파일을 쉽게 검색할 수 있습니다.

### 2.2. 사용자 지정 위치 설정하기

기본 위치가 있는 것은 좋지만 결과가 저장되는 위치와 구성 방법을 사용자 지정할 수 있습니다.

예를 들어 출력을 하위 디렉토리로 구성할 수 있습니다.
이를 수행하는 가장 간단한 방법은 출력별로 특정 출력 경로를 할당하는 것입니다.

#### 2.2.1. 출력 경로 수정하기

다시 한 번, 특정 출력에 대한 게시 동작을 수정하는 것은 정말 간단합니다.
사용자 지정 위치를 설정하려면 `path`를 적절히 편집하기만 하면 됩니다:

=== "후"

    ```groovy title="hello-world.nf" linenums="27" hl_lines="3"
    output {
        first_output {
            path 'hello_world'
        }
    }
    ```

=== "전"

    ```groovy title="hello-world.nf" linenums="27" hl_lines="3"
    output {
        first_output {
            path '.'
        }
    }
    ```

이것은 개별 출력 수준에서 설정되므로 필요에 맞게 다른 위치와 하위 디렉토리를 지정할 수 있습니다.

#### 2.2.2. 워크플로우 다시 실행하기

시도해 보겠습니다.

```bash
nextflow run hello-world.nf
```

??? success "명령 출력"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `hello-world.nf` [tiny_shaw] DSL2 - revision: 757723adc1

    executor >  local (1)
    [8c/79499c] process > sayHello [100%] 1 of 1 ✔
    ```

이번에는 결과가 지정된 하위 디렉토리 아래에 작성됩니다.

??? abstract "디렉토리 내용"

    ```console hl_lines="2-3"
    results/
    ├── hello_world
    │   └── output.txt -> /workspaces/training/hello-nextflow/work/8c/79499c2e506b79e2e01acb808d9d12/output.txt
    └── output.txt -> /workspaces/training/hello-nextflow/work/65/f56f2cd75df1352e106fcdd084b97b/output.txt
    ```

이전 실행의 결과가 여전히 있는 것을 볼 수 있습니다.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello_world_output.svg"
</figure>

원하는 만큼 많은 수준의 중첩을 사용할 수 있습니다.
프로세스 이름이나 다른 변수를 사용하여 결과를 구성하는 데 사용되는 디렉토리의 이름을 지정할 수도 있으며, 최상위 출력 디렉토리의 기본 이름을 변경할 수도 있습니다(`-o` CLI 플래그 또는 config 변수 `outputDir`로 제어됨).
이러한 옵션은 교육 후반부에서 다룰 것입니다.

### 2.3. 게시 모드를 copy로 설정하기

기본적으로 출력은 `work` 디렉토리에서 심볼릭 링크로 게시됩니다.
즉, 파일 시스템에 단일 파일만 있습니다.

이것은 여러 복사본을 저장하고 싶지 않은 매우 큰 파일을 다룰 때 좋습니다.
그러나 어느 시점에서 work 디렉토리를 삭제하면(정리 작업은 곧 다룰 것입니다) 파일에 대한 액세스 권한을 잃게 됩니다.
따라서 중요한 파일의 복사본을 안전한 장소에 저장하기 위한 계획이 필요합니다.

한 가지 쉬운 옵션은 관심 있는 출력에 대해 게시 모드를 copy로 전환하는 것입니다.

#### 2.3.1. mode 지시문 추가하기

이 부분은 정말 간단합니다.
관련 워크플로우 수준 출력 정의에 `mode 'copy'`를 추가하기만 하면 됩니다:

=== "후"

    ```groovy title="hello-world.nf" linenums="27" hl_lines="4"
    output {
        first_output {
            path 'hello_world'
            mode 'copy'
        }
    }
    ```

=== "전"

    ```groovy title="hello-world.nf" linenums="27"
    output {
        first_output {
            path 'hello_world'
        }
    }
    ```

이것은 해당 특정 출력에 대한 게시 모드를 설정합니다.

#### 2.3.2. 워크플로우 다시 실행하기

시도해 보겠습니다.

```bash
nextflow run hello-world.nf
```

??? success "명령 출력"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `hello-world.nf` [tiny_shaw] DSL2 - revision: 757723adc1

    executor >  local (1)
    [df/521638] process > sayHello [100%] 1 of 1 ✔
    ```

이번에는 결과를 보면 파일이 단순한 심볼릭 링크가 아닌 적절한 복사본입니다.

??? abstract "디렉토리 내용"

    ```console hl_lines="3"
    results/
    ├── hello_world
    │   └── output.txt
    └── output.txt -> /workspaces/training/hello-nextflow/work/65/f56f2cd75df1352e106fcdd084b97b/output.txt
    ```

이것도 개별 출력 수준에서 설정되므로 세밀한 방식으로 게시 모드를 설정할 수 있습니다.
이것은 나중에 다단계 파이프라인으로 이동할 때 특히 유용할 것입니다. 예를 들어 최종 출력만 복사하고 중간 출력은 심볼릭 링크로 남겨둘 수 있습니다.

앞서 언급했듯이 출력 게시 방법을 제어하기 위한 다른 더 정교한 옵션이 있습니다.
Nextflow 여정에서 적절한 시기에 이를 사용하는 방법을 보여드리겠습니다.

### 2.4. 프로세스 수준 `publishDir` 지시문에 대한 참고 사항

최근까지 출력을 게시하는 확립된 방법은 `publishDir` 지시문을 사용하여 각 개별 프로세스 수준에서 수행하는 것이었습니다.

`sayHello` 프로세스의 출력에 대해 방금 수행한 작업을 달성하려면 대신 프로세스 정의에 다음 줄을 추가했을 것입니다:

```groovy title="hello-world.nf" linenums="6" hl_lines="3"
process sayHello {

    publishDir 'results/hello_world', mode: 'copy'

    output:
    path 'output.txt'

    script:
    """
    echo 'Hello World!' > output.txt
    """
}
```

이전 Nextflow 파이프라인과 프로세스 모듈에서 이 코드 패턴을 여전히 많이 찾을 수 있으므로 이를 인식하는 것이 중요합니다.
그러나 향후 Nextflow 언어 버전에서는 결국 허용되지 않을 것이므로 새로운 작업에서는 이를 사용하지 않는 것이 좋습니다.

### 핵심 정리

워크플로우 출력을 더 편리한 위치에 게시하는 방법을 배웠습니다.

### 다음 단계

명령줄 매개변수를 통해 제공되는 변수 입력을 사용하고 기본값을 효과적으로 활용하는 방법을 배우세요.

---

## 3. 명령줄에서 전달된 변수 입력 사용하기

현재 상태에서 워크플로우는 프로세스 명령에 하드코딩된 인사말을 사용합니다.
런타임에 인사말을 더 쉽게 변경할 수 있도록 입력 변수를 사용하여 유연성을 추가하려고 합니다.

이를 위해 스크립트를 세 가지 변경해야 합니다:

1. 프로세스가 변수 입력을 예상하도록 변경
2. 사용자 입력을 캡처하기 위한 명령줄 매개변수 설정
3. 워크플로우 본문에서 프로세스에 입력 전달

이러한 변경 사항을 한 번에 하나씩 수행하겠습니다.

### 3.1. `sayHello` 프로세스가 변수 입력을 예상하도록 변경하기

프로세스 정의를 편집하여 (1) 입력 변수를 받아들이고 (2) 명령줄에서 해당 변수를 사용해야 합니다.

#### 3.1.1. 프로세스 정의에 input 블록 추가하기

먼저 `greeting`이라는 입력을 받아들이도록 프로세스 정의를 조정하겠습니다.

프로세스 블록에서 다음 코드 변경을 수행하세요:

=== "후"

    ```groovy title="hello-world.nf" linenums="6" hl_lines="3-4"
    process sayHello {

        input:
        val greeting

        output:
        path 'output.txt'
    ```

=== "전"

    ```groovy title="hello-world.nf" linenums="6"
    process sayHello {

        output:
        path 'output.txt'
    ```

`greeting` 변수 앞에 `val`이 붙어 Nextflow에게 이것이 값(경로가 아님)임을 알려줍니다.

#### 3.1.2. 입력 변수를 사용하도록 프로세스 명령 편집하기

이제 원래 하드코딩된 값을 받을 것으로 예상되는 입력 변수의 값으로 교체합니다.

프로세스 블록에서 다음 코드 변경을 수행하세요:

=== "후"

    ```groovy title="hello-world.nf" linenums="14" hl_lines="3"
    script:
    """
    echo '${greeting}' > output.txt
    """
    ```

=== "전"

    ```groovy title="hello-world.nf" linenums="14" hl_lines="3"
    script:
    """
    echo 'Hello World!' > output.txt
    """
    ```

`$` 기호와 중괄호(`{ }`)는 Nextflow에게 이것이 실제 입력 값으로 대체되어야 하는 변수 이름임을 알려줍니다(=보간됨).

!!! tip

    중괄호(`{ }`)는 이전 버전의 Nextflow에서는 기술적으로 선택 사항이었으므로 `echo '$greeting' > output.txt`로 작성된 이전 워크플로우를 볼 수 있습니다.

이제 `sayHello()` 프로세스가 변수 입력을 받을 준비가 되었으므로 워크플로우 수준에서 프로세스 호출에 입력 값을 제공하는 방법이 필요합니다.

### 3.2. 사용자 입력을 캡처하기 위한 명령줄 매개변수 설정하기

프로세스 호출을 `sayHello('Hello World!')`로 만들어 입력을 직접 하드코딩할 수 있습니다.
그러나 워크플로우로 실제 작업을 수행할 때는 명령줄에서 입력을 제어할 수 있기를 원할 것이므로 다음과 같이 할 수 있습니다:

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello_world_input.svg"
</figure>

다행히 Nextflow에는 CLI 매개변수를 쉽게 선언하고 사용할 수 있는 [`params`](https://nextflow.io/docs/latest/config.html#params)라는 내장 워크플로우 매개변수 시스템이 있습니다.

일반적인 구문은 `params.<parameter_name>`을 선언하여 Nextflow에게 명령줄에서 `--<parameter_name>` 매개변수를 예상하도록 알려주는 것입니다.

여기서는 `--input`이라는 매개변수를 만들고 싶으므로 워크플로우 어딘가에 `params.input`을 선언해야 합니다.
원칙적으로 어디에나 작성할 수 있지만 `sayHello()` 프로세스 호출에 제공하려고 하므로 `sayHello(params.input)`을 작성하여 직접 연결할 수 있습니다.

워크플로우 블록에서 다음 코드 변경을 수행하세요:

=== "후"

    ```groovy title="hello-world.nf" linenums="23" hl_lines="2"
    // emit a greeting
    sayHello(params.input)
    ```

=== "전"

    ```groovy title="hello-world.nf" linenums="23" hl_lines="2"
    // emit a greeting
    sayHello()
    ```

이것은 Nextflow에게 `--input` 매개변수를 통해 제공된 값에 대해 `sayHello` 프로세스를 실행하도록 알려줍니다.

실제로 섹션 시작 부분에 설명된 단계 (2)와 (3)을 한 번에 수행했습니다.

### 3.3. 워크플로우 명령 실행하기

실행해 보겠습니다!

```bash
nextflow run hello-world.nf --input 'Bonjour le monde!'
```

??? success "명령 출력"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `hello-world.nf` [elated_lavoisier] DSL2 - revision: 7c031b42ea

    executor >  local (1)
    [4b/654319] sayHello | 1 of 1 ✔
    ```

이러한 편집을 모두 올바르게 수행했다면 또 다른 성공적인 실행을 얻어야 합니다.

이제 새 버전의 인사말이 있는지 확인하기 위해 출력 파일을 열어보세요.

??? abstract "파일 내용"

    ```console title="results/hello_world/output.txt"
    Bonjour le monde!
    ```

Et voilà!

새 실행이 `results` 디렉토리에 게시된 출력 파일을 덮어썼습니다.
그러나 이전 실행의 결과는 여전히 `work` 아래의 작업 디렉토리에 보존됩니다.

!!! tip

    Nextflow 수준 매개변수와 파이프라인 수준 매개변수를 쉽게 구별할 수 있습니다.

    - 파이프라인에 적용되는 매개변수는 항상 이중 하이픈(`--`)을 사용합니다.
    - 이전에 사용한 `-resume` 기능과 같이 Nextflow 설정을 수정하는 매개변수는 단일 하이픈(`-`)을 사용합니다.

### 3.4. 명령줄 매개변수에 기본값 사용하기

좋습니다. 편리했지만 많은 경우 모든 실행에 대해 지정할 필요가 없도록 주어진 매개변수에 기본값을 제공하는 것이 합리적입니다.

#### 3.4.1. CLI 매개변수에 기본값 설정하기

워크플로우 정의 전에 선언하여 `input` 매개변수에 기본값을 지정하겠습니다.

```groovy title="hello-world.nf" linenums="20"
/*
 * 파이프라인 매개변수
 */
params {
    input: String = 'Holà mundo!'
}
```

보시다시피 워크플로우가 예상하는 입력 유형을 지정할 수 있습니다(Nextflow 25.10.2 이상).
구문은 `name: Type = default_value`입니다.
지원되는 유형에는 `String`, `Integer`, `Float`, `Boolean` 및 `Path`가 포함됩니다.

!!! info

    이전 워크플로우에서는 전체 `params` 블록이 단순히 `input = 'Holà mundo!'`로 작성된 것을 볼 수 있습니다.

파이프라인에 더 많은 매개변수를 추가할 때 기본값을 제공해야 하는지 여부에 관계없이 이 블록에 모두 추가해야 합니다.
이렇게 하면 구성 가능한 모든 매개변수를 한눈에 쉽게 찾을 수 있습니다.

#### 3.4.2. 매개변수를 지정하지 않고 워크플로우 다시 실행하기

이제 기본값이 설정되었으므로 명령줄에서 값을 지정하지 않고도 워크플로우를 다시 실행할 수 있습니다.

```bash
nextflow run hello-world.nf
```

??? success "명령 출력"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `hello-world.nf` [determined_edison] DSL2 - revision: 3539118582

    executor >  local (1)
    [72/394147] sayHello | 1 of 1 ✔
    ```

출력은 이전과 같은 위치에 있지만 내용은 새 텍스트로 업데이트되어야 합니다.

??? abstract "파일 내용"

    ```console title="results/hello_world/output.txt"
    Holà mundo!
    ```

Nextflow는 인사말 매개변수의 기본값을 사용하여 출력을 생성했습니다.

#### 3.4.3. 기본값 재정의하기

명령줄에서 매개변수를 제공하면 CLI 값이 기본값을 재정의합니다.

시도해 보세요:

```bash
nextflow run hello-world.nf --input 'Konnichiwa!'
```

??? success "명령 출력"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `hello-world.nf` [elegant_faraday] DSL2 - revision: 3539118582

    executor >  local (1)
    [6f/a12a91] sayHello | 1 of 1 ✔
    ```

다시 한 번 결과 디렉토리에서 해당하는 업데이트된 출력을 찾아야 합니다.

??? abstract "파일 내용"

    ```console title="results/hello_world/output.txt"
    Konnichiwa!
    ```

!!! note

    Nextflow에서는 매개변수 값을 지정할 수 있는 여러 위치가 있습니다.
    동일한 매개변수가 여러 위치에서 다른 값으로 설정된 경우 Nextflow는 [여기](https://www.nextflow.io/docs/latest/config.html)에 설명된 우선 순위에 따라 사용할 값을 결정합니다.

    이것은 파트 6(구성)에서 더 자세히 다룰 것입니다.

### 핵심 정리

명령줄 매개변수를 통해 런타임에 제공되는 간단한 변수 입력을 사용하는 방법과 기본값을 설정, 사용 및 재정의하는 방법을 배웠습니다.

### 다음 단계

워크플로우 실행을 더 편리하게 관리하는 방법을 배우세요.

---

## 4. 워크플로우 실행 관리하기

워크플로우를 시작하고 출력을 검색하는 방법을 아는 것은 좋지만, 특히 자신의 워크플로우를 개발하는 경우 삶을 더 쉽게 만들어줄 워크플로우 관리의 몇 가지 다른 측면이 있다는 것을 곧 알게 될 것입니다.

여기서는 동일한 워크플로우를 다시 시작해야 할 때 [`-resume`](https://nextflow.io/docs/latest/cache-and-resume.html) 기능을 사용하는 방법, [`nextflow log`](https://nextflow.io/docs/latest/reference/cli.html#log)로 과거 실행 로그를 검사하는 방법, [`nextflow clean`](https://nextflow.io/docs/latest/reference/cli.html#clean)으로 이전 work 디렉토리를 삭제하는 방법을 보여드립니다.

### 4.1. `-resume`으로 워크플로우 다시 시작하기

때때로 이미 성공적으로 완료된 단계를 다시 수행하지 않고 이전에 시작한 파이프라인을 다시 실행하고 싶을 것입니다.

Nextflow에는 이를 수행할 수 있는 [`-resume`](https://nextflow.io/docs/latest/cache-and-resume.html)이라는 옵션이 있습니다.
구체적으로 이 모드에서는 정확히 동일한 코드, 설정 및 입력으로 이미 실행된 프로세스는 건너뜁니다.
즉, Nextflow는 마지막 실행 이후 추가하거나 수정한 프로세스 또는 새 설정이나 입력을 제공하는 프로세스만 실행합니다.

이렇게 하면 두 가지 주요 이점이 있습니다:

- 파이프라인을 개발하는 중이라면 변경 사항을 테스트하기 위해 적극적으로 작업 중인 프로세스만 실행하면 되므로 더 빠르게 반복할 수 있습니다.
- 프로덕션에서 파이프라인을 실행하고 문제가 발생하면 많은 경우 문제를 해결하고 파이프라인을 다시 시작할 수 있으며, 실패 지점부터 실행을 재개하여 많은 시간과 컴퓨팅을 절약할 수 있습니다.

사용하려면 명령에 `-resume`을 추가하고 실행하기만 하면 됩니다:

```bash
nextflow run hello-world.nf -resume
```

??? success "명령 출력"

    ```console hl_lines="5"
    N E X T F L O W   ~  version 25.10.2

    Launching `hello-world.nf` [golden_cantor] DSL2 - revision: 35bd3425e5

    [62/49a1f8] sayHello | 1 of 1, cached: 1 ✔
    ```

콘솔 출력은 익숙해 보이지만 이전과 약간 다른 점이 하나 있습니다.

프로세스 상태 줄(5번째 줄)에 추가된 `cached:` 부분을 찾아보세요. 이는 Nextflow가 이미 이 작업을 수행했음을 인식하고 이전 성공적인 실행의 결과를 단순히 재사용했음을 의미합니다.

work 하위 디렉토리 해시가 이전 실행과 동일한 것도 볼 수 있습니다.
Nextflow는 말 그대로 이전 실행을 가리키며 "저기서 이미 했어요"라고 말하고 있습니다.

!!! tip

    `resume`으로 파이프라인을 다시 실행할 때 Nextflow는 이전에 성공적으로 실행된 실행에 의해 work 디렉토리 외부에 게시된 파일을 덮어쓰지 않습니다.

### 4.2. 과거 실행 로그 검사하기

새 파이프라인을 개발하든 프로덕션에서 파이프라인을 실행하든 어느 시점에서 과거 실행에 대한 정보를 조회해야 할 것입니다.
방법은 다음과 같습니다.

nextflow 워크플로우를 시작할 때마다 현재 작업 디렉토리의 `.nextflow`라는 숨겨진 디렉토리 아래에 `history`라는 로그 파일에 줄이 작성됩니다.

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

이것은 로그 파일의 내용을 헤더 줄과 함께 터미널에 출력합니다.

새 `nextflow run` 명령을 실행할 때마다 세션 ID가 변경되는 것을 알 수 있습니다. 단, `-resume` 옵션을 사용하는 경우는 예외입니다.
이 경우 세션 ID는 동일하게 유지됩니다.

Nextflow는 세션 ID를 사용하여 `.nextflow` 아래에 있는 `cache` 디렉토리에 실행 캐싱 정보를 그룹화합니다.

### 4.3. 이전 work 디렉토리 삭제하기

개발 프로세스 중에 일반적으로 초안 파이프라인을 여러 번 실행하게 되며, 이로 인해 많은 하위 디렉토리에 많은 파일이 축적될 수 있습니다.

다행히 Nextflow에는 더 이상 관심이 없는 과거 실행의 work 하위 디렉토리를 자동으로 삭제할 수 있는 유용한 `clean` 하위 명령이 포함되어 있습니다.

#### 4.3.1. 삭제 기준 결정하기

삭제할 항목을 결정하는 여러 [옵션](https://www.nextflow.io/docs/latest/reference/cli.html#clean)이 있습니다.

여기서는 실행 이름을 사용하여 지정된 특정 실행 이전의 실행에서 모든 하위 디렉토리를 삭제하는 예를 보여드립니다.

`-resume`을 사용하지 않은 가장 최근의 성공적인 실행을 찾아보세요. 우리의 경우 실행 이름은 `golden_cantor`였습니다.

실행 이름은 `Launching (...)` 콘솔 출력 줄의 대괄호 안에 표시된 기계 생성 두 부분 문자열입니다.
Nextflow 로그를 사용하여 타임스탬프 및/또는 명령줄을 기반으로 실행을 조회할 수도 있습니다.

#### 4.3.2. 드라이 런 수행하기

먼저 드라이 런 플래그 `-n`을 사용하여 명령이 주어졌을 때 무엇이 삭제될지 확인합니다:

```bash
nextflow clean -before golden_cantor -n
```

??? success "명령 출력"

    ```console
    Would remove /workspaces/training/hello-nextflow/work/a3/7be2fad5e71e5f49998f795677fd68
    ```

출력에는 다른 작업 디렉토리 이름이 있고 줄 수가 다를 수 있지만 예제와 비슷해야 합니다.

출력된 줄이 없으면 유효한 실행 이름을 제공하지 않았거나 삭제할 과거 실행이 없는 것입니다. 예제 명령의 `golden_cantor`를 로그에서 해당하는 최신 실행 이름으로 변경해야 합니다.

#### 4.3.3. 삭제 진행하기

출력이 예상대로 보이고 삭제를 진행하려면 `-n` 대신 `-f` 플래그를 사용하여 명령을 다시 실행하세요:

```bash
nextflow clean -before golden_cantor -f
```

??? success "명령 출력"

    ```console
    Removed /workspaces/training/hello-nextflow/work/a3/7be2fad5e71e5f49998f795677fd68
    ```

출력은 이전과 비슷해야 하지만 이제 'Would remove' 대신 'Removed'라고 표시됩니다.
이것은 두 문자 하위 디렉토리(위의 `a3/`와 같은)를 제거하지 않지만 내용을 비웁니다.

!!! Warning

    과거 실행의 work 하위 디렉토리를 삭제하면 Nextflow의 캐시에서 제거되고 해당 디렉토리에 저장된 모든 출력이 삭제됩니다.
    즉, 해당 프로세스를 다시 실행하지 않고 실행을 재개하는 Nextflow의 기능이 중단됩니다.

    관심 있거나 의존할 계획인 출력을 저장하는 것은 귀하의 책임입니다! 이것이 `publish` 지시문에 대해 `symlink` 모드보다 `copy` 모드를 사용하는 것을 선호하는 주된 이유입니다.

### 핵심 정리

특정 디렉토리에 출력을 게시하고, 이미 동일한 방식으로 실행된 단계를 반복하지 않고 파이프라인을 다시 시작하고, `nextflow clean` 명령을 사용하여 이전 work 디렉토리를 정리하는 방법을 배웠습니다.

더 일반적으로 간단한 Nextflow 워크플로우를 해석하고 실행을 관리하며 출력을 검색하는 방법을 배웠습니다.

### 다음 단계

잠시 휴식을 취하세요. 충분히 받을 자격이 있습니다!

준비가 되면 [**파트 2: Hello Channels**](./02_hello_channels.md)로 이동하여 채널을 사용하여 워크플로우에 입력을 공급하는 방법을 배우세요. 이를 통해 Nextflow의 내장 데이터 흐름 병렬 처리 및 기타 강력한 기능을 활용할 수 있습니다.

---

## 퀴즈

<quiz>
Nextflow 프로세스의 최소 필수 구성 요소는 무엇입니까?
- [ ] 입력 및 출력 블록만
- [x] 출력 및 script 블록
- [ ] 입력, 출력 및 script 블록
- [ ] script 블록만

자세히 알아보기: [1.1.1. 프로세스 정의](#111-the-process-definition)
</quiz>

<quiz>
프로세스에서 output 블록의 목적은 무엇입니까?
- [ ] 결과를 콘솔에 출력하기 위해
- [ ] work 디렉토리에 파일을 저장하기 위해
- [x] 프로세스에서 예상되는 출력을 선언하기 위해
- [ ] 환경 변수를 정의하기 위해

자세히 알아보기: [1.1.1. 프로세스 정의](#111-the-process-definition)
</quiz>

<quiz>
Nextflow 워크플로우를 실행하는 데 사용되는 명령은 무엇입니까?
- [ ] `nextflow start`
- [ ] `nextflow execute`
- [x] `nextflow run`
- [ ] `nextflow launch`
</quiz>

<quiz>
작업의 work 디렉토리를 보면 실제로 실행된 명령이 포함된 파일은 무엇입니까?

```
work/a3/7be2fa.../
├── .command.begin
├── .command.err
├── .command.log
├── .command.out
├── .command.run
├── .command.sh
├── .exitcode
└── output.txt
```

- [ ] `.command.run`
- [x] `.command.sh`
- [ ] `.command.log`
- [ ] `.command.out`

자세히 알아보기: [1.2.2. `work` 디렉토리에서 출력 및 로그 찾기](#122-find-the-output-and-logs-in-the-work-directory)
</quiz>

<quiz>
`-resume` 플래그는 무엇을 합니까?
- [ ] 워크플로우를 처음부터 다시 시작합니다
- [ ] 워크플로우를 일시 중지합니다
- [x] 이미 성공적으로 완료된 프로세스를 건너뜁니다
- [ ] 워크플로우의 백업을 생성합니다

자세히 알아보기: [4.1. `-resume`으로 워크플로우 다시 시작하기](#41-re-launch-a-workflow-with--resume)
</quiz>

<quiz>
워크플로우 출력을 게시하는 기본 모드는 무엇입니까?
- [ ] 출력 디렉토리에 파일 복사
- [x] 출력 디렉토리에 심볼릭 링크 생성
- [ ] 출력 디렉토리로 파일 이동
- [ ] 출력 디렉토리에서 파일 압축

자세히 알아보기: [2.3. 게시 모드를 copy로 설정하기](#23-set-the-publish-mode-to-copy)
</quiz>

<quiz>
명령줄에서 Nextflow 워크플로우에 매개변수 값을 전달하는 방법은 무엇입니까?
- [ ] `-parameter value`
- [ ] `--parameter:value`
- [x] `--parameter value`
- [ ] `-p parameter=value`

자세히 알아보기: [3.2. 사용자 입력을 캡처하기 위한 명령줄 매개변수 설정하기](#32-set-up-a-command-line-parameter-to-capture-user-input)
</quiz>

<quiz>
Nextflow script 블록 내에서 변수를 참조하는 방법은 무엇입니까?
- [ ] `%variable%` 구문 사용
- [x] `#!groovy ${variable}` 구문 사용
- [ ] `{{variable}}` 구문 사용
- [ ] `[variable]` 구문 사용
</quiz>
