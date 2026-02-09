# 워크플로우의 워크플로우

파이프라인을 개발하다 보면 서로 다른 데이터 유형이나 분석 단계에 대해 유사한 프로세스 시퀀스를 만드는 경우가 많습니다. 이러한 프로세스 시퀀스를 복사하여 붙여넣게 되면 유지 관리가 어려운 중복 코드가 생기거나, 이해하고 수정하기 어려운 하나의 거대한 워크플로우를 만들게 될 수 있습니다.

Nextflow의 가장 강력한 기능 중 하나는 더 작고 재사용 가능한 워크플로우 모듈로부터 복잡한 파이프라인을 구성할 수 있다는 점입니다. 이러한 모듈식 접근 방식은 파이프라인을 더 쉽게 개발하고, 테스트하고, 유지 관리할 수 있게 만듭니다.

### 학습 목표

이 사이드 퀘스트에서는 독립적으로 테스트하고 사용할 수 있는 워크플로우 모듈을 개발하는 방법, 이러한 모듈을 더 큰 파이프라인으로 구성하는 방법, 그리고 모듈 간 데이터 흐름을 관리하는 방법을 학습합니다.

이 사이드 퀘스트를 마치면 다음을 수행할 수 있습니다:

- 복잡한 파이프라인을 논리적이고 재사용 가능한 단위로 분해하기
- 각 워크플로우 모듈을 독립적으로 테스트하기
- 워크플로우를 조합하여 새로운 파이프라인 만들기
- 여러 파이프라인에서 공통 워크플로우 모듈 공유하기
- 코드를 더 유지 관리하기 쉽고 이해하기 쉽게 만들기

이러한 기술은 깔끔하고 유지 관리 가능한 코드 구조를 유지하면서 복잡한 파이프라인을 구축하는 데 도움이 됩니다.

### 사전 요구 사항

이 사이드 퀘스트를 시작하기 전에 다음을 완료해야 합니다:

- [Hello Nextflow](../hello_nextflow/README.md) 튜토리얼 또는 이에 상응하는 초급 과정을 완료했어야 합니다.
- 기본적인 Nextflow 개념과 메커니즘(프로세스, 채널, 연산자, 모듈)을 편안하게 사용할 수 있어야 합니다.

---

## 0. 시작하기

#### 교육 코드스페이스 열기

아직 열지 않았다면 [환경 설정](../envsetup/index.md)에 설명된 대로 교육 환경을 열어주세요.

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

#### 프로젝트 디렉토리로 이동

이 튜토리얼의 파일이 있는 디렉토리로 이동하겠습니다.

```bash
cd side-quests/workflows_of_workflows
```

VSCode가 이 디렉토리에 집중하도록 설정할 수 있습니다:

```bash
code .
```

#### 자료 검토

'Hello Nextflow'에서 학습한 내용을 기반으로 하는 여러 프로세스 정의가 포함된 `modules` 디렉토리를 찾을 수 있습니다:

```console title="Directory contents"
modules/
├── say_hello.nf             # Creates a greeting (from Hello Nextflow)
├── say_hello_upper.nf       # Converts to uppercase (from Hello Nextflow)
├── timestamp_greeting.nf    # Adds timestamps to greetings
├── validate_name.nf         # Validates input names
└── reverse_text.nf          # Reverses text content
```

#### 과제 검토

여러분의 과제는 이러한 모듈을 두 개의 별도 워크플로우로 조합한 다음, 이를 메인 워크플로우로 구성하는 것입니다:

- 이름을 검증하고, 인사말을 생성하고, 타임스탬프를 추가하는 `GREETING_WORKFLOW`
- 텍스트를 대문자로 변환하고 반전시키는 `TRANSFORM_WORKFLOW`

#### 준비 체크리스트

시작할 준비가 되었나요?

- [ ] 이 과정의 목표와 사전 요구 사항을 이해했습니다
- [ ] 코드스페이스가 실행 중입니다
- [ ] 작업 디렉토리를 적절하게 설정했습니다
- [ ] 과제를 이해했습니다

모든 항목을 체크할 수 있다면 시작할 준비가 된 것입니다.

---

## 1. 인사말 워크플로우 생성

이름을 검증하고 타임스탬프가 있는 인사말을 생성하는 워크플로우를 만들어 보겠습니다.

### 1.1. 워크플로우 구조 생성

```bash title="Create workflow directory and file"
mkdir -p workflows
touch workflows/greeting.nf
```

### 1.2. 첫 번째 (서브)워크플로우 코드 추가

`workflows/greeting.nf`에 다음 코드를 추가하세요:

```groovy title="workflows/greeting.nf" linenums="1"
include { VALIDATE_NAME } from '../modules/validate_name'
include { SAY_HELLO } from '../modules/say_hello'
include { TIMESTAMP_GREETING } from '../modules/timestamp_greeting'

workflow {

    names_ch = channel.of('Alice', 'Bob', 'Charlie')

    // Chain processes: validate -> create greeting -> add timestamp
    validated_ch = VALIDATE_NAME(names_ch)
    greetings_ch = SAY_HELLO(validated_ch)
    timestamped_ch = TIMESTAMP_GREETING(greetings_ch)
}
```

이것은 'Hello Nextflow' 튜토리얼에서 본 것과 유사한 구조를 가진 완전한 워크플로우이며, 독립적으로 테스트할 수 있습니다. 지금 시도해 보겠습니다:

```bash
nextflow run workflows/greeting.nf
```

??? success "명령 출력"

    ```console
    N E X T F L O W  ~  version 24.10.0
    Launching `workflows/greeting.nf` [peaceful_montalcini] DSL2 - revision: 90f61b7093
    executor >  local (9)
    [51/4f980f] process > VALIDATE_NAME (validating Bob)                    [100%] 3 of 3 ✔
    [2b/dd8dc2] process > SAY_HELLO (greeting Bob)                          [100%] 3 of 3 ✔
    [8e/882565] process > TIMESTAMP_GREETING (adding timestamp to greeting) [100%] 3 of 3 ✔
    ```

예상대로 작동하지만, 구성 가능하게 만들려면 몇 가지를 변경해야 합니다.

### 1.3. 워크플로우를 구성 가능하게 만들기

구성 가능한 워크플로우는 'Hello Nextflow' 튜토리얼에서 본 것과 몇 가지 차이점이 있습니다:

- workflow 블록에 이름이 필요합니다
- 입력은 `take:` 키워드를 사용하여 선언됩니다
- 워크플로우 내용은 `main:` 블록 안에 배치됩니다
- 출력은 `emit:` 키워드를 사용하여 선언됩니다

인사말 워크플로우를 이 구조에 맞게 업데이트해 보겠습니다. 코드를 다음과 같이 변경하세요:

```groovy title="workflows/greeting.nf" linenums="1" hl_lines="6 7 9 15 16 17"
include { VALIDATE_NAME } from '../modules/validate_name'
include { SAY_HELLO } from '../modules/say_hello'
include { TIMESTAMP_GREETING } from '../modules/timestamp_greeting'

workflow GREETING_WORKFLOW {
    take:
        names_ch        // Input channel with names

    main:
        // Chain processes: validate -> create greeting -> add timestamp
        validated_ch = VALIDATE_NAME(names_ch)
        greetings_ch = SAY_HELLO(validated_ch)
        timestamped_ch = TIMESTAMP_GREETING(greetings_ch)

    emit:
        greetings = greetings_ch      // Original greetings
        timestamped = timestamped_ch  // Timestamped greetings
}
```

워크플로우에 이제 이름이 있고 `take:` 및 `emit:` 블록이 있으며, 이것들이 상위 레벨 워크플로우를 구성하는 데 사용할 연결점입니다.
워크플로우 내용도 `main:` 블록 안에 배치되었습니다. 또한 `names_ch` 입력 채널 선언을 제거했는데, 이제 워크플로우에 인수로 전달되기 때문입니다.

워크플로우가 예상대로 작동하는지 다시 테스트해 보겠습니다:

```bash
nextflow run workflows/greeting.nf
```

??? failure "명령 출력"

    ```console
    N E X T F L O W  ~  version 24.10.0
    Launching `workflows/greeting.nf` [high_brahmagupta] DSL2 - revision: 8f5857af25
    No entry workflow specified
    ```

이것은 '엔트리 워크플로우'라는 또 다른 새로운 개념에 대해 알려줍니다. 엔트리 워크플로우는 Nextflow 스크립트를 실행할 때 호출되는 워크플로우입니다. 기본적으로 Nextflow는 이름이 없는 워크플로우가 있을 때 이를 엔트리 워크플로우로 사용하며, 지금까지 다음과 같이 시작하는 workflow 블록으로 이를 수행해 왔습니다:

```groovy title="hello.nf" linenums="1"
workflow {
```

하지만 우리의 인사말 워크플로우에는 이름이 없는 워크플로우가 없고, 대신 이름이 있는 워크플로우가 있습니다:

```groovy title="workflows/greeting.nf" linenums="1"
workflow GREETING_WORKFLOW {
```

그래서 Nextflow가 오류를 발생시키고 우리가 원하는 작업을 수행하지 않은 것입니다.

우리가 `take:`/`emit:` 구문을 추가한 것은 워크플로우를 직접 호출하기 위해서가 아니라 다른 워크플로우와 구성할 수 있도록 하기 위해서입니다. 해결책은 이름이 있는 워크플로우를 가져와서 호출하는 이름이 없는 엔트리 워크플로우가 있는 메인 스크립트를 만드는 것입니다.

### 1.4. 메인 워크플로우 생성 및 테스트

이제 `greeting` 워크플로우를 가져와서 사용하는 메인 워크플로우를 만들겠습니다.

`main.nf`를 생성하세요:

```groovy title="main.nf" linenums="1"
include { GREETING_WORKFLOW } from './workflows/greeting'

workflow {
    names = channel.of('Alice', 'Bob', 'Charlie')
    GREETING_WORKFLOW(names)

    GREETING_WORKFLOW.out.greetings.view { "Original: $it" }
    GREETING_WORKFLOW.out.timestamped.view { "Timestamped: $it" }
}

```

이 파일의 워크플로우 엔트리는 이름이 없으며, 이는 엔트리 워크플로우로 사용할 것이기 때문입니다.

이것을 실행하고 출력을 확인하세요:

```bash
nextflow run main.nf
```

??? success "명령 출력"

    ```console
    N E X T F L O W  ~  version 24.10.0
    Launching `main.nf` [goofy_mayer] DSL2 - revision: 543f8742fe
    executor >  local (9)
    [05/3cc752] process > GREETING_WORKFLOW:VALIDATE_NAME (validating Char... [100%] 3 of 3 ✔
    [b1/b56ecf] process > GREETING_WORKFLOW:SAY_HELLO (greeting Charlie)      [100%] 3 of 3 ✔
    [ea/342168] process > GREETING_WORKFLOW:TIMESTAMP_GREETING (adding tim... [100%] 3 of 3 ✔
    Original: /workspaces/training/side_quests/workflows_of_workflows/work/bb/c8aff3df0ebc15a4d7d35f736db44c/Alice-output.txt
    Original: /workspaces/training/side_quests/workflows_of_workflows/work/fb/fa877776e8a5d90b537b1bcd3b6f5b/Bob-output.txt
    Original: /workspaces/training/side_quests/workflows_of_workflows/work/b1/b56ecf938fda8bcbec211847c8f0be/Charlie-output.txt
    Timestamped: /workspaces/training/side_quests/workflows_of_workflows/work/06/877bc909f140bbf8223343450cea36/timestamped_Alice-output.txt
    Timestamped: /workspaces/training/side_quests/workflows_of_workflows/work/aa/bd31b71cdb745b7c155ca7f8837b8a/timestamped_Bob-output.txt
    Timestamped: /workspaces/training/side_quests/workflows_of_workflows/work/ea/342168d4ba04cc899a89c56cbfd9b0/timestamped_Charlie-output.txt
    ```

작동합니다! 이름이 있는 인사말 워크플로우를 이름이 없는 엔트리 `workflow` 블록이 있는 메인 워크플로우로 감쌌습니다. 메인 워크플로우는 `GREETING_WORKFLOW` 워크플로우를 프로세스처럼(완전히 같지는 않지만) 사용하고 있으며, `names` 채널을 인수로 전달하고 있습니다.

### 핵심 정리

이 섹션에서는 몇 가지 중요한 개념을 학습했습니다:

- **이름이 있는 워크플로우**: 가져와서 재사용할 수 있는 이름이 있는 워크플로우(`GREETING_WORKFLOW`) 생성하기
- **워크플로우 인터페이스**: 구성 가능한 워크플로우를 만들기 위해 `take:`로 명확한 입력을 정의하고 `emit:`로 출력을 정의하기
- **엔트리 포인트**: Nextflow가 스크립트를 실행하려면 이름이 없는 엔트리 워크플로우가 필요하다는 것 이해하기
- **워크플로우 구성**: 다른 워크플로우 내에서 이름이 있는 워크플로우를 가져와서 사용하기
- **워크플로우 네임스페이스**: `.out` 네임스페이스를 사용하여 워크플로우 출력에 접근하기(`GREETING_WORKFLOW.out.greetings`)

이제 다음과 같은 작동하는 인사말 워크플로우가 있습니다:

- 이름 채널을 입력으로 받습니다
- 각 이름을 검증합니다
- 각 유효한 이름에 대한 인사말을 생성합니다
- 인사말에 타임스탬프를 추가합니다
- 원본 및 타임스탬프가 있는 인사말을 모두 출력으로 노출합니다

이러한 모듈식 접근 방식을 사용하면 인사말 워크플로우를 독립적으로 테스트하거나 더 큰 파이프라인의 구성 요소로 사용할 수 있습니다.

---

## 2. 변환 워크플로우 추가

이제 인사말에 텍스트 변환을 적용하는 워크플로우를 만들어 보겠습니다.

### 2.1. 워크플로우 파일 생성

```bash
touch workflows/transform.nf
```

### 2.2. 워크플로우 코드 추가

`workflows/transform.nf`에 다음 코드를 추가하세요:

```groovy title="workflows/transform.nf" linenums="1"
include { SAY_HELLO_UPPER } from '../modules/say_hello_upper'
include { REVERSE_TEXT } from '../modules/reverse_text'

workflow TRANSFORM_WORKFLOW {
    take:
        input_ch         // Input channel with messages

    main:
        // Apply transformations in sequence
        upper_ch = SAY_HELLO_UPPER(input_ch)
        reversed_ch = REVERSE_TEXT(upper_ch)

    emit:
        upper = upper_ch        // Uppercase greetings
        reversed = reversed_ch  // Reversed uppercase greetings
}
```

여기서는 구성 가능한 구문에 대한 설명을 반복하지 않겠지만, 이름이 있는 워크플로우가 다시 `take:` 및 `emit:` 블록으로 선언되었고, 워크플로우 내용이 `main:` 블록 안에 배치되었다는 점에 주목하세요.

### 2.3. 메인 워크플로우 업데이트

두 워크플로우를 모두 사용하도록 `main.nf`를 업데이트하세요:

```groovy title="main.nf" linenums="1"
include { GREETING_WORKFLOW } from './workflows/greeting'
include { TRANSFORM_WORKFLOW } from './workflows/transform'

workflow {
    names = channel.of('Alice', 'Bob', 'Charlie')

    // Run the greeting workflow
    GREETING_WORKFLOW(names)

    // Run the transform workflow
    TRANSFORM_WORKFLOW(GREETING_WORKFLOW.out.timestamped)

    // View results
    TRANSFORM_WORKFLOW.out.upper.view { "Uppercase: $it" }
    TRANSFORM_WORKFLOW.out.reversed.view { "Reversed: $it" }
}
```

전체 파이프라인을 실행하세요:

```bash
nextflow run main.nf
```

??? success "명령 출력"

    ```console
    N E X T F L O W  ~  version 24.10.0
    Launching `main.nf` [sick_kimura] DSL2 - revision: 8dc45fc6a8
    executor >  local (13)
    executor >  local (15)
    [83/1b51f4] process > GREETING_WORKFLOW:VALIDATE_NAME (validating Alice)  [100%] 3 of 3 ✔
    [68/556150] process > GREETING_WORKFLOW:SAY_HELLO (greeting Alice)        [100%] 3 of 3 ✔
    [de/511abd] process > GREETING_WORKFLOW:TIMESTAMP_GREETING (adding tim... [100%] 3 of 3 ✔
    [cd/e6a7e0] process > TRANSFORM_WORKFLOW:SAY_HELLO_UPPER (converting t... [100%] 3 of 3 ✔
    [f0/74ba4a] process > TRANSFORM_WORKFLOW:REVERSE_TEXT (reversing UPPER... [100%] 3 of 3 ✔
    Uppercase: /workspaces/training/side_quests/workflows_of_workflows/work/a0/d4f5df4d6344604498fa47a6084a11/UPPER-timestamped_Bob-output.txt
    Uppercase: /workspaces/training/side_quests/workflows_of_workflows/work/69/b5e37f6c79c2fd38adb75d0eca8f87/UPPER-timestamped_Charlie-output.txt
    Uppercase: /workspaces/training/side_quests/workflows_of_workflows/work/cd/e6a7e0b17e7d5a2f71bb8123cd53a7/UPPER-timestamped_Alice-output.txt
    Reversed: /workspaces/training/side_quests/workflows_of_workflows/work/7a/7a222f7957b35d1d121338566a24ac/REVERSED-UPPER-timestamped_Bob-output.txt
    Reversed: /workspaces/training/side_quests/workflows_of_workflows/work/46/8d19af62e33a5a6417c773496e0f90/REVERSED-UPPER-timestamped_Charlie-output.txt
    Reversed: /workspaces/training/side_quests/workflows_of_workflows/work/f0/74ba4a10d9ef5c82f829d1c154d0f6/REVERSED-UPPER-timestamped_Alice-output.txt
    ```

반전된 파일 중 하나를 살펴보면, 대문자 버전의 인사말이 반전된 것을 볼 수 있습니다:

```bash
cat /workspaces/training/side_quests/workflows_of_workflows/work/f0/74ba4a10d9ef5c82f829d1c154d0f6/REVERSED-UPPER-timestamped_Alice-output.txt
```

```console title="Reversed file content"
!ECILA ,OLLEH ]04:50:71 60-30-5202[
```

### 핵심 정리

이제 다음과 같은 완전한 파이프라인이 있어야 합니다:

- 인사말 워크플로우를 통해 이름을 처리합니다
- 타임스탬프가 있는 인사말을 변환 워크플로우에 전달합니다
- 인사말의 대문자 및 반전 버전을 모두 생성합니다

---

## 요약

이 사이드 퀘스트에서는 더 작고 재사용 가능한 구성 요소로부터 복잡한 파이프라인을 구축할 수 있게 해주는 Nextflow의 강력한 워크플로우 구성 개념을 학습했습니다.

이러한 모듈식 접근 방식은 모놀리식 파이프라인에 비해 여러 가지 장점을 제공합니다:

- 각 워크플로우를 독립적으로 개발하고, 테스트하고, 디버그할 수 있습니다
- 워크플로우를 여러 파이프라인에서 재사용할 수 있습니다
- 전체 파이프라인 구조가 더 읽기 쉽고 유지 관리하기 쉬워집니다
- 인터페이스가 일관성을 유지하는 한 한 워크플로우의 변경이 다른 워크플로우에 반드시 영향을 미치지 않습니다
- 필요에 따라 파이프라인의 다른 부분을 실행하도록 엔트리 포인트를 구성할 수 있습니다

_하지만 워크플로우를 호출하는 것이 프로세스를 호출하는 것과 약간 비슷하지만, 실제로는 같은 것이 아니라는 점을 주목하는 것이 중요합니다. 예를 들어, 크기가 N인 채널로 워크플로우를 호출하여 워크플로우를 N번 실행할 수는 없습니다 - 크기가 N인 채널을 워크플로우에 전달하고 내부적으로 반복해야 합니다._

이러한 기술을 여러분의 작업에 적용하면 유지 관리 가능하고 확장 가능한 상태를 유지하면서 복잡한 생물정보학 작업을 처리할 수 있는 더 정교한 Nextflow 파이프라인을 구축할 수 있습니다.

### 주요 패턴

1.  **워크플로우 구조**: `take:` 및 `emit:` 구문을 사용하여 각 워크플로우에 대한 명확한 입력과 출력을 정의하여 구성 요소 간에 잘 정의된 인터페이스를 만들고, `main:` 블록 내에 워크플로우 로직을 감쌌습니다.

    ```groovy
    workflow EXAMPLE_WORKFLOW {
        take:
            // Input channels are declared here
            input_ch

        main:
            // Workflow logic goes here
            // This is where processes are called and channels are manipulated
            result_ch = SOME_PROCESS(input_ch)

        emit:
            // Output channels are declared here
            output_ch = result_ch
    }
    ```

2.  **워크플로우 가져오기:** 두 개의 독립적인 워크플로우 모듈을 구축하고 include 문으로 메인 파이프라인에 가져왔습니다.

    - 단일 워크플로우 포함

    ```groovy
    include { WORKFLOW_NAME } from './path/to/workflow'
    ```

    - 여러 워크플로우 포함

    ```groovy
    include { WORKFLOW_A; WORKFLOW_B } from './path/to/workflows'
    ```

    - 이름 충돌을 피하기 위해 별칭으로 포함

    ```groovy
    include { WORKFLOW_A as WORKFLOW_A_ALIAS } from './path/to/workflow'
    ```

3.  **엔트리 포인트**: Nextflow는 실행을 시작할 위치를 알기 위해 이름이 없는 엔트리 워크플로우가 필요합니다. 이 엔트리 워크플로우가 이름이 있는 워크플로우를 호출합니다.

    - 이름이 없는 워크플로우 (엔트리 포인트)

    ```groovy
    workflow {
        // This is the entry point when the script is run
        NAMED_WORKFLOW(input_ch)
    }
    ```

    - 이름이 있는 워크플로우 (엔트리 워크플로우에서 호출됨)

    ```groovy
    workflow NAMED_WORKFLOW {
        // Must be called from the entry workflow
    }
    ```

4.  **데이터 흐름 관리:** 네임스페이스 표기법(`WORKFLOW_NAME.out.channel_name`)을 사용하여 워크플로우 출력에 접근하고 이를 다른 워크플로우에 전달하는 방법을 학습했습니다.

    ```nextflow
    WORKFLOW_A(input_ch)
    WORKFLOW_B(WORKFLOW_A.out.some_channel)
    ```

### 추가 자료

- [Nextflow Workflow 문서](https://www.nextflow.io/docs/latest/workflow.html)
- [Channel 연산자 참조](https://www.nextflow.io/docs/latest/operator.html)
- [Error Strategy 문서](https://www.nextflow.io/docs/latest/process.html#errorstrategy)

---

## 다음 단계

[사이드 퀘스트 메뉴](./index.md)로 돌아가거나 페이지 오른쪽 하단의 버튼을 클릭하여 목록의 다음 주제로 이동하세요.
