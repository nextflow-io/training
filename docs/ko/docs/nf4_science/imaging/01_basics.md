# 파트 1: 기본 작업 실행하기

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } AI 지원 번역 - [자세히 알아보기 및 개선 제안](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Nextflow for Bioimaging 교육 과정의 첫 번째 파트에서는 매우 기본적인 도메인 독립적 Hello World 예제를 사용하여 필수 작업을 시연하고 해당하는 Nextflow 코드 구성 요소를 설명하겠습니다.

## 1. 워크플로우 실행하기

`--greeting`이라는 명령줄 인수를 통해 입력을 받아 해당 인사말이 포함된 텍스트 파일을 생성하는 `hello-world.nf`라는 워크플로우 스크립트를 제공합니다.
아직 코드를 살펴보지는 않을 것이며, 먼저 실행하는 것이 어떤 모습인지 확인해 보겠습니다.

### 1.1. 워크플로우 시작 및 실행 모니터링하기

터미널에서 다음 명령을 실행하십시오:

```bash
nextflow run hello-world.nf --greeting 'Hello World!'
```

콘솔 출력은 다음과 같이 표시되어야 합니다:

```console title="출력" linenums="1"
 N E X T F L O W   ~  version 25.04.3

Launching `hello-world.nf` [goofy_torvalds] DSL2 - revision: c33d41f479

executor >  local (1)
[a3/7be2fa] sayHello | 1 of 1 ✔
```

축하합니다. 첫 번째 Nextflow 워크플로우를 실행하셨습니다!

여기서 가장 중요한 출력은 마지막 줄(6번 줄)입니다:

```console title="출력" linenums="6"
[a3/7be2fa] sayHello | 1 of 1 ✔
```

이는 `sayHello` process가 한 번 성공적으로 실행되었음을 알려줍니다(`1 of 1 ✔`).

좋습니다만, 출력이 어디에 있는지 궁금하실 것입니다.

### 1.2. `results` 디렉토리에서 출력 파일 찾기

이 워크플로우는 출력을 `results`라는 디렉토리에 게시하도록 구성되어 있습니다.
현재 디렉토리를 살펴보면 워크플로우를 실행할 때 Nextflow가 `results`라는 새 디렉토리를 생성했으며, 이 디렉토리에는 `output.txt`라는 파일이 포함되어 있습니다.

```console title="results/" linenums="1"
results
└── output.txt
```

파일을 열어보십시오. 내용은 명령줄에서 지정한 인사말과 일치해야 합니다.

<details>
  <summary>파일 내용</summary>

```console title="results/output.txt" linenums="1"
Hello World!
```

</details>

훌륭합니다. 워크플로우가 의도한 대로 작동했습니다!

그러나 '게시된' 결과는 Nextflow가 워크플로우를 실행할 때 생성한 실제 출력의 사본(또는 경우에 따라 심볼릭 링크)이라는 점을 알아두십시오.

이제 Nextflow가 실제로 작업을 실행한 위치를 확인하기 위해 내부를 살펴보겠습니다.

!!! warning "경고"

    모든 워크플로우가 results 디렉토리에 출력을 게시하도록 설정되어 있는 것은 아니며, 디렉토리 이름이 다를 수 있습니다.
    이 섹션 후반부에서 이 동작이 어디에 지정되어 있는지 확인하는 방법을 보여드리겠습니다.

### 1.3. `work/` 디렉토리에서 원본 출력 및 로그 찾기

워크플로우를 실행하면 Nextflow는 워크플로우의 각 process 호출마다 별도의 'task 디렉토리'를 생성합니다(=파이프라인의 각 단계).
각각에 대해 필요한 입력을 준비하고, 관련 명령을 실행하며, 해당 디렉토리 내에 출력 및 로그 파일을 작성합니다. 이 디렉토리는 고유성을 위해 해시를 사용하여 자동으로 이름이 지정됩니다.

이러한 모든 task 디렉토리는 현재 디렉토리(명령을 실행하는 위치) 내의 `work`라는 디렉토리 아래에 위치합니다.

혼란스러울 수 있으니 실제로 어떤 모습인지 살펴보겠습니다.

앞서 실행한 워크플로우의 콘솔 출력으로 돌아가면 다음과 같은 줄이 있었습니다:

```console title="명령 출력 발췌" linenums="6"
[a3/7be2fa] sayHello | 1 of 1 ✔
```

줄이 `[a3/7be2fa]`로 시작하는 것이 보이시나요?
이것은 해당 process 호출에 대한 task 디렉토리 경로의 축약된 형태이며, `work/` 디렉토리 경로 내에서 `sayHello` process 호출의 출력을 어디에서 찾을 수 있는지 알려줍니다.

다음 명령을 입력하고(`a3/7be2fa`를 터미널에 표시된 것으로 교체) tab 키를 눌러 경로를 자동 완성하거나 별표를 추가하여 전체 경로를 찾을 수 있습니다:

```bash
tree work/a3/7be2fa*
```

전체 디렉토리 경로가 표시됩니다: `work/a3/7be2fa7be2fad5e71e5f49998f795677fd68`

그 안에 무엇이 있는지 살펴보겠습니다.

!!! Tip "팁"

    VSCode 파일 탐색기에서 task 하위 디렉토리의 내용을 탐색하면 모든 파일을 바로 볼 수 있습니다.
    그러나 로그 파일은 터미널에서 보이지 않도록 설정되어 있으므로 `ls` 또는 `tree`를 사용하여 보려면 보이지 않는 파일을 표시하는 관련 옵션을 설정해야 합니다.

    ```bash
    tree -a work
    ```

시스템에서는 정확한 하위 디렉토리 이름이 다를 것입니다.

<details>
  <summary>디렉토리 내용</summary>

```console title="work/"
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

</details>

`output.txt` 파일을 바로 확인할 수 있을 것입니다. 이것은 실제로 `results` 디렉토리에 게시된 `sayHello` process의 원본 출력입니다.
파일을 열면 다시 `Hello World!` 인사말을 찾을 수 있습니다.

<details>
  <summary>output.txt의 파일 내용</summary>

```console title="work/a3/7be2fa7be2fad5e71e5f49998f795677fd68/output.txt" linenums="1"
Hello World!
```

</details>

그렇다면 다른 모든 파일은 무엇일까요?

이것들은 Nextflow가 task 실행의 일부로 작성한 도우미 및 로그 파일입니다:

- **`.command.begin`**: task가 시작되자마자 생성되는 센티넬 파일
- **`.command.err`**: process 호출에서 발생한 오류 메시지(`stderr`)
- **`.command.log`**: process 호출에서 발생한 전체 로그 출력
- **`.command.out`**: process 호출의 일반 출력(`stdout`)
- **`.command.run`**: Nextflow가 process 호출을 실행하기 위해 실행한 전체 스크립트
- **`.command.sh`**: process 호출에서 실제로 실행된 명령
- **`.exitcode`**: 명령 실행 결과로 생성된 종료 코드

`.command.sh` 파일은 Nextflow가 모든 장부 처리와 task/환경 설정을 제외하고 실제로 실행한 주요 명령을 보여주기 때문에 특히 유용합니다.

<details>
  <summary>파일 내용</summary>

```console title="work/a3/7be2fa7be2fad5e71e5f49998f795677fd68/.command.sh" linenums="1"
#!/bin/bash -ue
echo 'Hello World!' > output.txt

```

</details>

!!! Tip "팁"

    문제가 발생하여 무엇이 잘못되었는지 문제를 해결해야 할 때, `command.sh` 스크립트를 살펴보면 워크플로우 지시사항, 변수 보간 등을 기반으로 Nextflow가 구성한 정확한 명령을 확인하는 데 유용할 수 있습니다.

### 1.4. 선택 연습: 다른 인사말로 다시 실행하기

`--greeting` 인수에 다른 값을 사용하여 워크플로우를 여러 번 다시 실행한 다음 `results/` 디렉토리와 task 디렉토리의 내용을 모두 살펴보십시오.

격리된 task 디렉토리의 출력 및 로그가 어떻게 보존되는지 관찰하고, `results` 디렉토리의 내용은 후속 실행의 출력으로 덮어쓰이는 것을 확인하십시오.

### 핵심 사항

간단한 Nextflow 스크립트를 실행하고, 실행을 모니터링하며, 출력을 찾는 방법을 알게 되었습니다.

### 다음 단계

기본 Nextflow 스크립트를 읽고 그 구성 요소가 기능과 어떻게 연관되는지 식별하는 방법을 배웁니다.

---

## 2. Hello World 워크플로우 스타터 스크립트 살펴보기

앞서 우리가 한 것은 기본적으로 워크플로우 스크립트를 블랙박스처럼 취급한 것입니다.
이제 무엇을 하는지 확인했으니 상자를 열어 내부를 살펴보겠습니다.

_여기서의 목표는 Nextflow 코드의 구문을 암기하는 것이 아니라 주요 구성 요소가 무엇이며 어떻게 구성되는지에 대한 기본적인 직관을 형성하는 것입니다._

### 2.1. 전체 코드 구조 살펴보기

편집기 창에서 `hello-world.nf` 스크립트를 열어보겠습니다.

<details>
  <summary>코드</summary>

```groovy title="hello-world.nf" linenums="1"
#!/usr/bin/env nextflow

/*
 * Use echo to print a greeting to a file
 */
process sayHello {

    publishDir 'results', mode: 'copy'

    input:
    val greeting

    output:
    path 'output.txt'

    script:
    """
    echo '$greeting' > output.txt
    """
}

workflow {

    // 인사말 출력
    sayHello(params.greeting)
}
```

</details>

Nextflow 스크립트는 두 가지 주요 핵심 구성 요소를 포함합니다: 하나 이상의 **process**와 **workflow** 자체입니다.
각 **process**는 파이프라인의 해당 단계가 수행해야 할 작업을 설명하며, **workflow**는 다양한 단계를 연결하는 데이터 흐름 논리를 설명합니다.

먼저 **process** 블록을 자세히 살펴본 다음 **workflow** 블록을 살펴보겠습니다.

### 2.2. `process` 정의

첫 번째 코드 블록은 **process**를 설명합니다.
process 정의는 `process` 키워드로 시작하고 process 이름이 뒤따르며 마지막으로 중괄호로 구분된 process 본문이 옵니다.
process 본문에는 실행할 명령을 지정하는 script 블록이 포함되어야 하며, 명령줄 터미널에서 실행할 수 있는 모든 것이 가능합니다.

여기에는 `greeting`이라는 **input** 변수를 받아 `output.txt`라는 파일에 **output**을 작성하는 `sayHello`라는 **process**가 있습니다.

<details>
  <summary>코드</summary>

```groovy title="hello-world.nf" linenums="3"
/*
 * Use echo to print a greeting to a file
 */
process sayHello {

    publishDir 'results', mode: 'copy'

    input:
    val greeting

    output:
    path 'output.txt'

    script:
    """
    echo '$greeting' > output.txt
    """
}
```

</details>

이것은 `input` 정의, `output` 정의 및 실행할 `script`만 포함하는 매우 최소한의 process 정의입니다.

`input` 정의에는 `val` 한정자가 포함되어 있으며, 이는 Nextflow에게 어떤 종류의 값(문자열, 숫자 등)을 기대하도록 지시합니다.

`output` 정의에는 `path` 한정자가 포함되어 있으며, 이는 Nextflow에게 이것을 경로로 처리하도록 지시합니다(디렉토리 경로와 파일 모두 포함).

!!! Tip "팁"

    output 정의는 어떤 출력이 생성될지 _결정_하지 않습니다.
    단순히 예상되는 출력 파일을 어디에서 찾을 수 있는지 _선언_하므로 Nextflow가 실행이 완료되면 이를 찾을 수 있습니다.

    이는 명령이 성공적으로 실행되었는지 확인하고 필요한 경우 출력을 다운스트림 process에 전달하는 데 필요합니다.
    output 블록에 선언된 것과 일치하지 않는 출력은 다운스트림 process로 전달되지 않습니다.

실제 파이프라인에서 process는 일반적으로 추가 정보(예: process 지시문)를 포함하며, 이는 곧 소개하겠습니다.

### 2.3. `workflow` 정의

두 번째 코드 블록은 **workflow** 자체를 설명합니다.
workflow 정의는 `workflow` 키워드로 시작하고 선택적 이름이 뒤따르며 중괄호로 구분된 workflow 본문이 옵니다.

여기에는 `--greeting` 매개변수에 제공한 값을 담고 있는 `params.greeting`이라는 입력을 받는 `sayHello` process를 한 번 호출하는 **workflow**가 있습니다.

```groovy title="hello-world.nf" linenums="22"
workflow {

    // 인사말 출력
    sayHello(params.greeting)
}
```

이것은 매우 최소한의 **workflow** 정의입니다.
실제 파이프라인에서 workflow는 일반적으로 **channel**로 연결된 여러 **process** 호출을 포함하며, 변수 입력에 대한 기본값이 설정될 수 있습니다.

파트 2에서 nf-core/molkart를 실행할 때 이를 실제로 볼 수 있습니다.

### 2.4. 명령줄 매개변수의 `params` 시스템

`sayHello()` process 호출에 제공하는 `params.greeting`은 Nextflow 코드의 깔끔한 부분이며 추가로 1분 정도 시간을 투자할 가치가 있습니다.

위에서 언급한 것처럼 이것이 `--greeting` 명령줄 매개변수의 값을 `sayHello()` process 호출에 전달하는 방법입니다.
실제로 `params.someParameterName`을 선언하기만 하면 명령줄에서 `--someParameterName`이라는 매개변수를 워크플로우에 제공할 수 있게 됩니다.

!!! Tip "팁"

    `params` 시스템을 사용하여 선언된 이러한 워크플로우 매개변수는 항상 두 개의 대시(`--`)를 사용합니다.
    이는 한 개의 대시(`-`)만 사용하는 Nextflow 수준 매개변수와 구별됩니다.

### 핵심 사항

이제 간단한 Nextflow 워크플로우가 어떻게 구조화되어 있으며 기본 구성 요소가 기능과 어떻게 연관되는지 알게 되었습니다.

### 다음 단계

워크플로우 실행을 편리하게 관리하는 방법을 배웁니다.

---

## 3. 워크플로우 실행 관리하기

워크플로우를 시작하고 출력을 검색하는 방법을 아는 것은 훌륭하지만, 삶을 더 쉽게 만들어줄 워크플로우 관리의 몇 가지 다른 측면이 있다는 것을 곧 알게 될 것입니다.

여기서는 동일한 워크플로우를 다시 시작해야 할 때 `resume` 기능을 활용하는 방법, `nextflow log`로 실행 로그를 검사하는 방법, `nextflow clean`으로 이전 작업 디렉토리를 삭제하는 방법을 보여드립니다.

### 3.1. `-resume`으로 워크플로우 다시 시작하기

때때로 이전에 이미 시작한 파이프라인을 성공적으로 완료된 작업을 다시 수행하지 않고 다시 실행하고 싶을 것입니다.

Nextflow에는 이를 수행할 수 있는 `-resume`이라는 옵션이 있습니다.
구체적으로 이 모드에서는 정확히 동일한 코드, 설정 및 입력으로 이미 실행된 process는 건너뜁니다.
즉, Nextflow는 마지막 실행 이후 추가하거나 수정한 process 또는 새로운 설정이나 입력을 제공하는 process만 실행합니다.

이렇게 하면 두 가지 주요 이점이 있습니다:

- 파이프라인을 개발 중인 경우 변경 사항을 테스트하기 위해 작업 중인 process만 실행하면 되므로 더 빠르게 반복할 수 있습니다.
- 프로덕션에서 파이프라인을 실행하다가 문제가 발생한 경우 대부분의 경우 문제를 해결하고 파이프라인을 다시 시작할 수 있으며, 실패 지점부터 실행이 재개되므로 많은 시간과 컴퓨팅 리소스를 절약할 수 있습니다.

사용하려면 명령에 `-resume`을 추가하고 실행하십시오:

```bash
nextflow run hello-world.nf --greeting 'Hello World!' -resume
```

??? success "명령 출력"

    ```console
    N E X T F L O W   ~  version 25.04.3

    Launching `hello-world.nf` [tiny_noyce] DSL2 - revision: c33d41f479

    [a3/7be2fa] process > sayHello [100%] 1 of 1, cached: 1 ✔
    ```

process 상태 줄(5번 줄)에 추가된 `cached:` 부분을 찾아보십시오. 이는 Nextflow가 이미 이 작업을 수행했음을 인식하고 이전에 성공한 실행의 결과를 단순히 재사용했음을 의미합니다.

작업 하위 디렉토리 해시가 이전 실행과 동일하다는 것도 확인할 수 있습니다.
Nextflow는 말 그대로 이전 실행을 가리키며 "저기에서 이미 했습니다"라고 말하고 있습니다.

!!! Tip "팁"

    `resume`으로 파이프라인을 다시 실행할 때 Nextflow는 이전에 성공적으로 실행된 process 호출에 의해 `publishDir` 디렉토리에 작성된 파일을 덮어쓰지 않습니다.

### 3.2. 과거 실행 로그 검사하기

nextflow 워크플로우를 시작할 때마다 현재 작업 디렉토리의 `.nextflow`라는 숨겨진 디렉토리 아래에 있는 `history`라는 로그 파일에 한 줄이 기록됩니다.

이 정보에 액세스하는 더 편리한 방법은 `nextflow log` 명령을 사용하는 것입니다.

```bash
nextflow log
```

이렇게 하면 로그 파일의 내용이 터미널에 출력되어 현재 작업 디렉토리 내에서 시작된 모든 Nextflow 실행의 타임스탬프, 실행 이름, 상태 및 전체 명령줄을 보여줍니다.

### 3.3. 이전 작업 디렉토리 삭제하기

개발 프로세스 중에 일반적으로 초안 파이프라인을 여러 번 실행하게 되며, 이로 인해 많은 하위 디렉토리에 걸쳐 매우 많은 파일이 축적될 수 있습니다.
하위 디렉토리의 이름이 무작위로 지정되므로 이름만으로는 이전 실행과 최근 실행을 구별하기 어렵습니다.

Nextflow에는 더 이상 필요하지 않은 과거 실행의 작업 하위 디렉토리를 자동으로 삭제할 수 있는 편리한 `clean` 하위 명령이 포함되어 있으며, 삭제할 항목을 제어하는 여러 [옵션](https://www.nextflow.io/docs/latest/reference/cli.html#clean)이 있습니다.

Nextflow 로그를 사용하여 타임스탬프 및/또는 명령줄을 기반으로 실행을 조회한 다음 `nextflow clean -before <run_name> -f`를 사용하여 이전 실행의 작업 디렉토리를 삭제할 수 있습니다.

!!! Warning "경고"

    과거 실행의 작업 하위 디렉토리를 삭제하면 Nextflow의 캐시에서 제거되고 해당 디렉토리에 저장된 모든 출력이 삭제됩니다.
    즉, 해당 process를 다시 실행하지 않고 Nextflow가 실행을 재개할 수 있는 기능이 손상됩니다.

    보관하거나 의존할 계획인 모든 출력을 저장하는 것은 귀하의 책임입니다! 이 목적으로 `publishDir` 지시문을 사용하는 경우 `symlink` 모드가 아닌 `copy` 모드를 사용해야 합니다.

### 핵심 사항

이미 동일한 방식으로 실행된 단계를 반복하지 않고 파이프라인을 다시 시작하고, 실행 로그를 검사하며, `nextflow clean` 명령을 사용하여 이전 작업 디렉토리를 정리하는 방법을 알게 되었습니다.

### 다음 단계

이제 기본 Nextflow 작업을 이해했으므로 nf-core/molkart로 실제 생물 이미징 파이프라인을 실행할 준비가 되었습니다.
