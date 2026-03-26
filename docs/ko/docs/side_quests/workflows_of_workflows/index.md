# 워크플로우의 워크플로우

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } AI 지원 번역 - [자세히 알아보기 및 개선 제안](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

파이프라인을 개발하다 보면, 서로 다른 데이터 유형이나 분석 단계에 대해 유사한 프로세스 시퀀스를 반복적으로 만들게 되는 경우가 많습니다. 이러한 프로세스 시퀀스를 복사하여 붙여넣다 보면 유지 관리하기 어려운 중복 코드가 생기거나, 이해하고 수정하기 어려운 거대한 워크플로우 하나가 만들어질 수 있습니다.

Nextflow의 가장 강력한 기능 중 하나는 더 작고 재사용 가능한 워크플로우 모듈로 복잡한 파이프라인을 구성할 수 있다는 점입니다. 이러한 모듈식 접근 방식은 파이프라인의 개발, 테스트, 유지 관리를 더 쉽게 만들어 줍니다.

### 학습 목표

이 사이드 퀘스트에서는 독립적으로 테스트하고 사용할 수 있는 워크플로우 모듈을 개발하는 방법, 해당 모듈들을 더 큰 파이프라인으로 구성하는 방법, 그리고 모듈 간 데이터 흐름을 관리하는 방법을 살펴봅니다.

이 사이드 퀘스트를 마치면 다음을 수행할 수 있습니다:

- 복잡한 파이프라인을 논리적이고 재사용 가능한 단위로 분리하기
- 각 워크플로우 모듈을 독립적으로 테스트하기
- 워크플로우를 조합하여 새로운 파이프라인 만들기
- 서로 다른 파이프라인에서 공통 워크플로우 모듈 공유하기
- 코드를 더 유지 관리하기 쉽고 이해하기 쉽게 만들기

이러한 기술은 깔끔하고 유지 관리하기 쉬운 코드 구조를 유지하면서 복잡한 파이프라인을 구축하는 데 도움이 됩니다.

### 사전 요구 사항

이 사이드 퀘스트를 시작하기 전에 다음을 완료해야 합니다:

- [Hello Nextflow](../hello_nextflow/README.md) 튜토리얼 또는 동급의 입문 과정을 완료해야 합니다.
- 기본적인 Nextflow 개념과 메커니즘(프로세스, 채널, 연산자, 모듈)을 편안하게 사용할 수 있어야 합니다.

---

## 0. 시작하기

#### 교육 코드스페이스 열기

아직 열지 않았다면, [환경 설정](../envsetup/index.md)에 설명된 대로 교육 환경을 열어야 합니다.

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

#### 프로젝트 디렉토리로 이동

이 튜토리얼의 파일이 있는 디렉토리로 이동합니다.

```bash
cd side-quests/workflows_of_workflows
```

VSCode가 이 디렉토리에 집중하도록 설정할 수 있습니다:

```bash
code .
```

#### 자료 검토

'Hello Nextflow'에서 학습한 내용을 기반으로 하는 여러 프로세스 정의가 포함된 `modules` 디렉토리를 확인할 수 있습니다:

```console title="Directory contents"
modules/
├── say_hello.nf             # Creates a greeting (from Hello Nextflow)
├── say_hello_upper.nf       # Converts to uppercase (from Hello Nextflow)
├── timestamp_greeting.nf    # Adds timestamps to greetings
├── validate_name.nf         # Validates input names
└── reverse_text.nf          # Reverses text content
```

#### 과제 검토

여러분의 과제는 이 모듈들을 두 개의 별도 워크플로우로 조합하고, 이를 하나의 메인 워크플로우로 구성하는 것입니다:

- 이름을 유효성 검사하고, 인사말을 생성하며, 타임스탬프를 추가하는 `GREETING_WORKFLOW`
- 텍스트를 대문자로 변환하고 뒤집는 `TRANSFORM_WORKFLOW`

<!-- TODO: give a bit more details, similar to how it's done in the Metadata side quest -->

#### 준비 체크리스트

시작할 준비가 되었나요?

- [ ] 이 과정의 목표와 사전 요구 사항을 이해했습니다
- [ ] 코드스페이스가 실행 중입니다
- [ ] 작업 디렉토리를 적절하게 설정했습니다
- [ ] 과제를 이해했습니다

모든 항목을 확인했다면 시작할 준비가 된 것입니다.

---

## 1. Greeting 워크플로우 만들기

이름을 유효성 검사하고 타임스탬프가 포함된 인사말을 생성하는 워크플로우를 만들어 봅니다.

### 1.1. 워크플로우 구조 만들기

```bash title="Create workflow directory and file"
mkdir -p workflows
touch workflows/greeting.nf
```

### 1.2. 첫 번째 (서브)워크플로우 코드 추가

`workflows/greeting.nf`에 다음 코드를 추가합니다:

```groovy title="workflows/greeting.nf" linenums="1"
include { VALIDATE_NAME } from '../modules/validate_name'
include { SAY_HELLO } from '../modules/say_hello'
include { TIMESTAMP_GREETING } from '../modules/timestamp_greeting'

workflow {

    names_ch = channel.of('Alice', 'Bob', 'Charlie')

    // 프로세스 연결: 유효성 검사 -> 인사말 생성 -> 타임스탬프 추가
    validated_ch = VALIDATE_NAME(names_ch)
    greetings_ch = SAY_HELLO(validated_ch)
    timestamped_ch = TIMESTAMP_GREETING(greetings_ch)
}
```

이것은 'Hello Nextflow' 튜토리얼에서 본 것과 유사한 구조를 가진 완전한 워크플로우로, 독립적으로 테스트할 수 있습니다. 지금 바로 테스트해 봅니다:

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

예상대로 작동하지만, 구성 가능하게 만들기 위해 몇 가지를 변경해야 합니다.

### 1.3. 워크플로우를 구성 가능하게 만들기

구성 가능한 워크플로우는 'Hello Nextflow' 튜토리얼에서 본 것과 몇 가지 차이점이 있습니다:

- 워크플로우 블록에 이름이 있어야 합니다
- 입력은 `take:` 키워드를 사용하여 선언합니다
- 워크플로우 내용은 `main:` 블록 안에 배치합니다
- 출력은 `emit:` 키워드를 사용하여 선언합니다

이 구조에 맞게 greeting 워크플로우를 업데이트합니다. 코드를 다음과 같이 변경합니다:

<!-- TODO: switch to before/after tabs -->

```groovy title="workflows/greeting.nf" linenums="1" hl_lines="6 7 9 15 16 17"
include { VALIDATE_NAME } from '../modules/validate_name'
include { SAY_HELLO } from '../modules/say_hello'
include { TIMESTAMP_GREETING } from '../modules/timestamp_greeting'

workflow GREETING_WORKFLOW {
    take:
        names_ch        // 이름이 담긴 입력 채널

    main:
        // 프로세스 연결: 유효성 검사 -> 인사말 생성 -> 타임스탬프 추가
        validated_ch = VALIDATE_NAME(names_ch)
        greetings_ch = SAY_HELLO(validated_ch)
        timestamped_ch = TIMESTAMP_GREETING(greetings_ch)

    emit:
        greetings = greetings_ch      // 원본 인사말
        timestamped = timestamped_ch  // 타임스탬프가 추가된 인사말
}
```

워크플로우에 이름이 생겼고 `take:` 및 `emit:` 블록이 추가된 것을 확인할 수 있습니다. 이것이 상위 레벨 워크플로우를 구성하는 데 사용할 연결 지점입니다.
워크플로우 내용도 `main:` 블록 안에 배치되었습니다. 또한 `names_ch` 입력 채널 선언이 제거된 것을 확인할 수 있는데, 이제 워크플로우의 인자로 전달되기 때문입니다.

워크플로우가 예상대로 작동하는지 다시 테스트해 봅니다:

```bash
nextflow run workflows/greeting.nf
```

??? failure "명령 출력"

    ```console
    N E X T F L O W  ~  version 24.10.0
    Launching `workflows/greeting.nf` [high_brahmagupta] DSL2 - revision: 8f5857af25
    No entry workflow specified
    ```

이것은 '진입 워크플로우(entry workflow)'라는 새로운 개념을 알려줍니다. 진입 워크플로우는 Nextflow 스크립트를 실행할 때 호출되는 워크플로우입니다. 기본적으로 Nextflow는 이름 없는 워크플로우가 있을 경우 이를 진입 워크플로우로 사용하며, 지금까지는 다음과 같이 시작하는 워크플로우 블록으로 이 방식을 사용해 왔습니다:

```groovy title="hello.nf" linenums="1"
workflow {
```

하지만 greeting 워크플로우에는 이름 없는 워크플로우가 없고, 대신 이름이 있는 워크플로우가 있습니다:

```groovy title="workflows/greeting.nf" linenums="1"
workflow GREETING_WORKFLOW {
```

그래서 Nextflow가 오류를 발생시키고 원하는 대로 실행되지 않은 것입니다.

`take:`/`emit:` 구문을 추가한 것은 워크플로우를 직접 호출하기 위해서가 아니라, 다른 워크플로우와 구성하기 위해서입니다. 해결책은 이름이 있는 워크플로우를 가져와서 호출하는 이름 없는 진입 워크플로우가 있는 메인 스크립트를 만드는 것입니다.

### 1.4. 메인 워크플로우 만들기 및 테스트

이제 `greeting` 워크플로우를 가져와서 사용하는 메인 워크플로우를 만들겠습니다.

`main.nf`를 만듭니다:

```groovy title="main.nf" linenums="1"
include { GREETING_WORKFLOW } from './workflows/greeting'

workflow {
    names = channel.of('Alice', 'Bob', 'Charlie')
    GREETING_WORKFLOW(names)

    GREETING_WORKFLOW.out.greetings.view { "Original: $it" }
    GREETING_WORKFLOW.out.timestamped.view { "Timestamped: $it" }
}

```

이 파일의 워크플로우 진입점은 이름이 없으며, 진입 워크플로우로 사용할 것이기 때문입니다.

실행하여 출력을 확인합니다:

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

작동합니다! 이름이 있는 greeting 워크플로우를 이름 없는 진입 `workflow` 블록이 있는 메인 워크플로우로 감쌌습니다. 메인 워크플로우는 `GREETING_WORKFLOW` 워크플로우를 프로세스와 거의 유사하게(완전히 같지는 않음) 사용하며, `names` 채널을 인자로 전달합니다.

### 핵심 정리

이 섹션에서는 다음과 같은 중요한 개념들을 학습했습니다:

- **이름 있는 워크플로우**: 가져와서 재사용할 수 있는 이름 있는 워크플로우(`GREETING_WORKFLOW`) 만들기
- **워크플로우 인터페이스**: `take:`로 입력을, `emit:`으로 출력을 명확하게 정의하여 구성 가능한 워크플로우 만들기
- **진입점**: Nextflow가 스크립트를 실행하려면 이름 없는 진입 워크플로우가 필요하다는 것 이해하기
- **워크플로우 구성**: 다른 워크플로우 안에서 이름 있는 워크플로우를 가져와 사용하기
- **워크플로우 네임스페이스**: `.out` 네임스페이스를 사용하여 워크플로우 출력에 접근하기(`GREETING_WORKFLOW.out.greetings`)

이제 다음을 수행하는 작동하는 greeting 워크플로우가 생겼습니다:

- 이름 채널을 입력으로 받기
- 각 이름의 유효성 검사
- 유효한 각 이름에 대한 인사말 생성
- 인사말에 타임스탬프 추가
- 원본 인사말과 타임스탬프가 추가된 인사말을 모두 출력으로 노출

이 모듈식 접근 방식을 통해 greeting 워크플로우를 독립적으로 테스트하거나 더 큰 파이프라인의 구성 요소로 사용할 수 있습니다.

---

## 2. Transform 워크플로우 추가하기

이제 인사말에 텍스트 변환을 적용하는 워크플로우를 만들어 봅니다.

### 2.1. 워크플로우 파일 만들기

```bash
touch workflows/transform.nf
```

### 2.2. 워크플로우 코드 추가

`workflows/transform.nf`에 다음 코드를 추가합니다:

```groovy title="workflows/transform.nf" linenums="1"
include { SAY_HELLO_UPPER } from '../modules/say_hello_upper'
include { REVERSE_TEXT } from '../modules/reverse_text'

workflow TRANSFORM_WORKFLOW {
    take:
        input_ch         // 메시지가 담긴 입력 채널

    main:
        // 순서대로 변환 적용
        upper_ch = SAY_HELLO_UPPER(input_ch)
        reversed_ch = REVERSE_TEXT(upper_ch)

    emit:
        upper = upper_ch        // 대문자 인사말
        reversed = reversed_ch  // 뒤집힌 대문자 인사말
}
```

구성 가능한 구문에 대한 설명은 반복하지 않겠지만, 이름 있는 워크플로우가 다시 `take:` 및 `emit:` 블록과 함께 선언되고, 워크플로우 내용이 `main:` 블록 안에 배치된 것을 확인하세요.

### 2.3. 메인 워크플로우 업데이트

두 워크플로우를 모두 사용하도록 `main.nf`를 업데이트합니다:

```groovy title="main.nf" linenums="1"
include { GREETING_WORKFLOW } from './workflows/greeting'
include { TRANSFORM_WORKFLOW } from './workflows/transform'

workflow {
    names = channel.of('Alice', 'Bob', 'Charlie')

    // greeting 워크플로우 실행
    GREETING_WORKFLOW(names)

    // transform 워크플로우 실행
    TRANSFORM_WORKFLOW(GREETING_WORKFLOW.out.timestamped)

    // 결과 확인
    TRANSFORM_WORKFLOW.out.upper.view { "Uppercase: $it" }
    TRANSFORM_WORKFLOW.out.reversed.view { "Reversed: $it" }
}
```

전체 파이프라인을 실행합니다:

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

뒤집힌 파일 중 하나를 확인하면, 인사말의 대문자 버전이 뒤집혀 있는 것을 볼 수 있습니다:

```bash
cat /workspaces/training/side_quests/workflows_of_workflows/work/f0/74ba4a10d9ef5c82f829d1c154d0f6/REVERSED-UPPER-timestamped_Alice-output.txt
```

```console title="Reversed file content"
!ECILA ,OLLEH ]04:50:71 60-30-5202[
```

### 핵심 정리

이제 다음을 수행하는 완전한 파이프라인이 완성되었습니다:

- greeting 워크플로우를 통해 이름 처리
- 타임스탬프가 추가된 인사말을 transform 워크플로우로 전달
- 인사말의 대문자 버전과 뒤집힌 버전 모두 생성

---

## 요약

이 사이드 퀘스트에서는 더 작고 재사용 가능한 구성 요소로 복잡한 파이프라인을 구축할 수 있게 해주는 Nextflow의 강력한 워크플로우 구성 개념을 살펴보았습니다.

이 모듈식 접근 방식은 단일 파이프라인 방식에 비해 여러 가지 장점을 제공합니다:

- 각 워크플로우를 독립적으로 개발, 테스트, 디버그할 수 있습니다
- 워크플로우를 서로 다른 파이프라인에서 재사용할 수 있습니다
- 전체 파이프라인 구조가 더 읽기 쉽고 유지 관리하기 쉬워집니다
- 인터페이스가 일관되게 유지되는 한, 하나의 워크플로우 변경이 다른 워크플로우에 반드시 영향을 미치지는 않습니다
- 필요에 따라 파이프라인의 다른 부분을 실행하도록 진입점을 설정할 수 있습니다

_단, 워크플로우를 호출하는 것이 프로세스를 호출하는 것과 유사해 보이지만, 실제로는 동일하지 않다는 점을 유의해야 합니다. 예를 들어, 크기가 N인 채널로 워크플로우를 호출하여 N번 실행할 수는 없습니다. 크기가 N인 채널을 워크플로우에 전달하고 내부적으로 반복해야 합니다._

이러한 기술을 실제 작업에 적용하면 유지 관리 가능하고 확장 가능한 상태를 유지하면서 복잡한 생물정보학 작업을 처리할 수 있는 더 정교한 Nextflow 파이프라인을 구축할 수 있습니다.

### 핵심 패턴

1.  **워크플로우 구조**: `take:` 및 `emit:` 구문을 사용하여 각 워크플로우의 입력과 출력을 명확하게 정의하고, 구성 요소 간에 잘 정의된 인터페이스를 만들었으며, 워크플로우 로직을 `main:` 블록 안에 배치했습니다.

    ```groovy
    workflow EXAMPLE_WORKFLOW {
        take:
            // 입력 채널은 여기에 선언합니다
            input_ch

        main:
            // 워크플로우 로직은 여기에 작성합니다
            // 프로세스를 실행하고 채널을 조작하는 곳입니다
            result_ch = SOME_PROCESS(input_ch)

        emit:
            // 출력 채널은 여기에 선언합니다
            output_ch = result_ch
    }
    ```

2.  **워크플로우 가져오기**: 두 개의 독립적인 워크플로우 모듈을 만들고 include 구문으로 메인 파이프라인에 가져왔습니다.

    - 단일 워크플로우 가져오기

    ```groovy
    include { WORKFLOW_NAME } from './path/to/workflow'
    ```

    - 여러 워크플로우 가져오기

    ```groovy
    include { WORKFLOW_A; WORKFLOW_B } from './path/to/workflows'
    ```

    - 이름 충돌을 피하기 위한 별칭으로 가져오기

    ```groovy
    include { WORKFLOW_A as WORKFLOW_A_ALIAS } from './path/to/workflow'
    ```

3.  **진입점**: Nextflow는 실행을 시작할 위치를 알기 위해 이름 없는 진입 워크플로우가 필요합니다. 이 진입 워크플로우가 이름 있는 워크플로우를 호출합니다.

    - 이름 없는 워크플로우 (진입점)

    ```groovy
    workflow {
        // 스크립트가 실행될 때의 진입점입니다
        NAMED_WORKFLOW(input_ch)
    }
    ```

    - 이름 있는 워크플로우 (진입 워크플로우에서 호출됨)

    ```groovy
    workflow NAMED_WORKFLOW {
        // 진입 워크플로우에서 호출되어야 합니다
    }
    ```

4.  **데이터 흐름 관리**: 네임스페이스 표기법(`WORKFLOW_NAME.out.channel_name`)을 사용하여 워크플로우 출력에 접근하고 다른 워크플로우로 전달하는 방법을 학습했습니다.

    ```nextflow
    WORKFLOW_A(input_ch)
    WORKFLOW_B(WORKFLOW_A.out.some_channel)
    ```

### 추가 자료

- [Nextflow 워크플로우 문서](https://www.nextflow.io/docs/latest/workflow.html)
- [채널 연산자 참조](https://www.nextflow.io/docs/latest/operator.html)
- [오류 전략 문서](https://www.nextflow.io/docs/latest/process.html#errorstrategy)

---

## 다음 단계

[사이드 퀘스트 메뉴](../)로 돌아가거나 페이지 오른쪽 하단의 버튼을 클릭하여 목록의 다음 주제로 이동하세요.
