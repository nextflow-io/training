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

- [Hello Nextflow](../../hello_nextflow/index.md) 튜토리얼 또는 동급의 입문 과정을 완료해야 합니다.
- 기본적인 Nextflow 개념과 메커니즘(프로세스, 채널, 연산자, 모듈)을 편안하게 사용할 수 있어야 합니다.

---

## 0. 시작하기

#### 교육 코드스페이스 열기

아직 열지 않았다면, [환경 설정](../../envsetup/index.md)에 설명된 대로 교육 환경을 열어야 합니다.

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

편집기가 해당 프로젝트 디렉토리에 집중된 상태로 열립니다.

#### 자료 검토

`modules` 디렉토리에는 프로세스 정의가, `workflows` 디렉토리에는 미리 작성된 두 개의 워크플로우 스크립트가 있으며, 단계적으로 업데이트할 `main.nf` 파일도 확인할 수 있습니다:

```console title="Directory contents"
├── main.nf
├── workflows/
│   ├── greeting.nf              # Standalone greeting workflow (to be made composable)
│   └── transform.nf             # Standalone transform workflow (to be made composable)
└── modules/
    ├── say_hello.nf             # Creates a greeting (from Hello Nextflow)
    ├── say_hello_upper.nf       # Converts to uppercase (from Hello Nextflow)
    ├── timestamp_greeting.nf    # Adds timestamps to greetings
    ├── validate_name.nf         # Validates input names
    └── reverse_text.nf          # Reverses text content
```

`modules/` 디렉토리에는 개별 프로세스 정의가 포함되어 있고, `workflows/` 디렉토리에는 이 사이드 퀘스트에서 사용할 두 개의 미리 작성된 워크플로우 스크립트가 포함되어 있습니다.

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

## 1. 파이프라인에 greeting 워크플로우 추가하기

greeting 워크플로우는 이름을 유효성 검사하고 타임스탬프가 포함된 인사말을 생성합니다.

### 1.1. greeting 워크플로우 검토 및 실행

`workflows/greeting.nf`를 열고 코드를 확인합니다:

```groovy title="workflows/greeting.nf" linenums="1"
include { VALIDATE_NAME } from '../modules/validate_name'
include { SAY_HELLO } from '../modules/say_hello'
include { TIMESTAMP_GREETING } from '../modules/timestamp_greeting'

workflow {
    main:
    names_ch = channel.of('Alice', 'Bob', 'Charlie')

    // 프로세스 연결: 유효성 검사 -> 인사말 생성 -> 타임스탬프 추가
    validated_ch = VALIDATE_NAME(names_ch)
    greetings_ch = SAY_HELLO(validated_ch)
    timestamped_ch = TIMESTAMP_GREETING(greetings_ch)

    publish:
    greetings = greetings_ch
    timestamped = timestamped_ch
}

output {
    greetings {
    }
    timestamped {
    }
}
```

이것은 'Hello Nextflow' 튜토리얼에서 본 것과 동일한 구조를 가진 완전한 단독 실행형 워크플로우입니다.
입력 이름이 하드코딩되어 있고, 세 개의 프로세스를 연결하며, 두 개의 출력을 게시합니다.

모든 것이 정상적으로 작동하는지 실행하여 확인합니다:

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

다른 워크플로우와 구성 가능하게 만들기 위해 몇 가지를 변경해야 합니다.

### 1.2. 워크플로우를 구성 가능하게 만들기

워크플로우를 구성 가능하게 만들려면 네 가지를 변경해야 합니다:
워크플로우에 이름을 부여하고, 입력을 `take:` 블록으로 이동하고, 출력을 `emit:` 블록으로 이동하고,
단독 실행형 `publish:`/`output {}` 블록을 제거합니다(이 블록들은 진입 워크플로우에 속합니다).

각 변경 사항을 하나씩 살펴봅니다.

#### 1.2.1. 워크플로우에 이름 부여하기

상위 워크플로우에서 가져올 수 있도록 워크플로우에 이름을 부여합니다.

=== "후"

    ```groovy title="workflows/greeting.nf" linenums="5" hl_lines="1"
    workflow GREETING_WORKFLOW {
    ```

=== "전"

    ```groovy title="workflows/greeting.nf" linenums="5" hl_lines="1"
    workflow {
    ```

이름이 생기면 워크플로우를 다른 스크립트로 가져올 수 있습니다.

#### 1.2.2. `take:`로 입력 선언하기

하드코딩된 채널 선언을 워크플로우가 기대하는 입력을 선언하는 `take:` 블록으로 교체합니다.
`take:` 블록은 `main:` 앞에 위치하며, `names_ch = channel.of(...)` 줄은 제거합니다.

=== "후"

    ```groovy title="workflows/greeting.nf" linenums="5" hl_lines="2 3 5"
    workflow GREETING_WORKFLOW {
        take:
        names_ch // 이름이 담긴 입력 채널

        main:
        // 프로세스 연결: 유효성 검사 -> 인사말 생성 -> 타임스탬프 추가
        validated_ch = VALIDATE_NAME(names_ch)
        greetings_ch = SAY_HELLO(validated_ch)
        timestamped_ch = TIMESTAMP_GREETING(greetings_ch)
    ```

=== "전"

    ```groovy title="workflows/greeting.nf" linenums="5"
    workflow GREETING_WORKFLOW {
        main:
        names_ch = channel.of('Alice', 'Bob', 'Charlie')

        // 프로세스 연결: 유효성 검사 -> 인사말 생성 -> 타임스탬프 추가
        validated_ch = VALIDATE_NAME(names_ch)
        greetings_ch = SAY_HELLO(validated_ch)
        timestamped_ch = TIMESTAMP_GREETING(greetings_ch)
    ```

`take:` 블록은 채널을 이름으로만 선언합니다. 채널에 무엇이 들어갈지는 상위 워크플로우에서 정의됩니다.

#### 1.2.3. `emit:`으로 출력 선언하기

`publish:` 섹션을 제거하고 `output {}` 블록을 삭제한 후, 출력에 이름을 부여하는 `emit:` 블록으로 교체합니다.

=== "후"

    ```groovy title="workflows/greeting.nf" linenums="14" hl_lines="2 3 4"

        emit:
        greetings = greetings_ch // 원본 인사말
        timestamped = timestamped_ch // 타임스탬프가 추가된 인사말
    }
    ```

=== "전"

    ```groovy title="workflows/greeting.nf" linenums="14"

        publish:
        greetings = greetings_ch
        timestamped = timestamped_ch
    }

    output {
        greetings {
        }
        timestamped {
        }
    }
    ```

`emit:` 블록은 상위 워크플로우가 `GREETING_WORKFLOW.out.greetings` 및 `GREETING_WORKFLOW.out.timestamped`를 통해 접근할 수 있는 이름 있는 출력을 노출합니다.

#### 1.2.4. 결과 확인 및 테스트

세 가지 변경 사항을 모두 적용한 후, 완성된 파일은 다음과 같아야 합니다:

```groovy title="workflows/greeting.nf" linenums="1" hl_lines="5 6 7 9 15 16 17"
include { VALIDATE_NAME } from '../modules/validate_name'
include { SAY_HELLO } from '../modules/say_hello'
include { TIMESTAMP_GREETING } from '../modules/timestamp_greeting'

workflow GREETING_WORKFLOW {
    take:
    names_ch // 이름이 담긴 입력 채널

    main:
    // 프로세스 연결: 유효성 검사 -> 인사말 생성 -> 타임스탬프 추가
    validated_ch = VALIDATE_NAME(names_ch)
    greetings_ch = SAY_HELLO(validated_ch)
    timestamped_ch = TIMESTAMP_GREETING(greetings_ch)

    emit:
    greetings = greetings_ch // 원본 인사말
    timestamped = timestamped_ch // 타임스탬프가 추가된 인사말
}
```

직접 실행해 봅니다:

```bash
nextflow run workflows/greeting.nf
```

??? failure "명령 출력"

    ```console
    N E X T F L O W  ~  version 24.10.0
    Launching `workflows/greeting.nf` [high_brahmagupta] DSL2 - revision: 8f5857af25
    No entry workflow specified
    ```

이것은 **진입 워크플로우(entry workflow)**라는 핵심 개념을 알려줍니다.
Nextflow는 스크립트를 직접 실행할 때 이름 없는 `workflow {}` 블록을 진입점으로 사용합니다.
`GREETING_WORKFLOW`는 이름이 있으므로 Nextflow는 단독으로 실행하는 방법을 알지 못합니다.

이것은 의도된 동작입니다. 구성 가능한 워크플로우는 직접 실행하는 것이 아니라 진입 워크플로우에서 호출하도록 설계되어 있습니다.
해결책은 `main.nf`에 `GREETING_WORKFLOW`를 가져와서 호출하는 진입 워크플로우를 만드는 것입니다.

### 1.3. 메인 워크플로우 업데이트 및 테스트

이제 greeting 워크플로우를 호출하도록 메인 워크플로우를 업데이트합니다.

#### 1.3.1. greeting 워크플로우 가져오기 및 호출

`include` 구문을 추가하고, `GREETING_WORKFLOW`를 호출하도록 워크플로우 본문을 업데이트하며, `publish:`의 `channel.empty()` 플레이스홀더를 교체합니다:

=== "후"

    ```groovy title="main.nf" linenums="1" hl_lines="1 7 8 11"
    include { GREETING_WORKFLOW } from './workflows/greeting'

    workflow {
        main:
        names = channel.of('Alice', 'Bob', 'Charlie')

        // greeting 워크플로우 실행
        GREETING_WORKFLOW(names)

        publish:
        greetings = GREETING_WORKFLOW.out.greetings
    }
    ```

=== "전"

    ```groovy title="main.nf" linenums="1"
    workflow {
        main:
        names = channel.of('Alice', 'Bob', 'Charlie')

        publish:
        greetings = channel.empty()
    }
    ```

진입 워크플로우는 Nextflow가 파이프라인 진입점으로 사용할 수 있도록 이름 없이 유지합니다.

#### 1.3.2. output 블록 업데이트

게시된 인사말을 `greetings/` 하위 디렉토리로 라우팅하는 `path` 지시문을 추가합니다:

=== "후"

    ```groovy title="main.nf" linenums="14" hl_lines="3"
    output {
        greetings {
            path 'greetings'
        }
    }
    ```

=== "전"

    ```groovy title="main.nf" linenums="14" hl_lines="2 3"
    output {
        greetings {
        }
    }
    ```

#### 1.3.3. 워크플로우 실행

워크플로우를 실행하여 정상적으로 작동하는지 테스트합니다:

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
    ```

??? abstract "디렉토리 내용"

    ```console
    results/
    └── greetings
        ├── Alice-output.txt
        ├── Bob-output.txt
        └── Charlie-output.txt
    ```

??? abstract "파일 내용"

    ```console title="results/greetings/Alice-output.txt"
    Hello, Alice!
    ```

인사말 파일이 `results/greetings/`에 게시됩니다.
메인 워크플로우는 `GREETING_WORKFLOW`를 호출하고 그 출력을 `publish:` 섹션에 직접 연결합니다.

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

## 2. 파이프라인에 transformation 워크플로우 추가하기

transform 워크플로우는 타임스탬프가 추가된 인사말에 텍스트 변환을 적용합니다.

### 2.1. 워크플로우 검토 및 실행

`workflows/transform.nf`를 열고 코드를 확인합니다:

```groovy title="workflows/transform.nf" linenums="1"
include { SAY_HELLO_UPPER } from '../modules/say_hello_upper'
include { REVERSE_TEXT } from '../modules/reverse_text'

workflow {
    main:
    input_ch = channel.fromPath('results/timestamped_*.txt')

    // 순서대로 변환 적용
    upper_ch = SAY_HELLO_UPPER(input_ch)
    reversed_ch = REVERSE_TEXT(upper_ch)

    publish:
    upper = upper_ch
    reversed = reversed_ch
}

output {
    upper {
    }
    reversed {
    }
}
```

이 단독 실행형 워크플로우는 `greeting.nf`가 생성한 `results/` 디렉토리의 타임스탬프가 추가된 인사말 파일을 읽어 대문자로 변환한 후 텍스트를 뒤집습니다.

1.1 섹션의 greeting 결과와 함께 정상적으로 작동하는지 실행하여 확인합니다:

```bash
nextflow run workflows/transform.nf
```

??? success "명령 출력"

    ```console
    N E X T F L O W  ~  version 24.10.0
    Launching `workflows/transform.nf` [blissful_curie] DSL2 - revision: 4e7b1c9f02
    executor >  local (6)
    [3e/a14c29] process > SAY_HELLO_UPPER (converting t... [100%] 3 of 3 ✔
    [c8/51b9e3] process > REVERSE_TEXT (reversing UPPER... [100%] 3 of 3 ✔
    ```

`GREETING_WORKFLOW`와 구성 가능하게 만들려면 1.2 섹션과 동일한 세 가지 변경 사항을 적용합니다.

### 2.2. 구성 가능하게 만들기

1.2 섹션과 동일한 세 가지 변경 사항을 적용합니다: 워크플로우에 이름을 부여하고, 하드코딩된 입력을 `take:`로 교체하고, `publish:`/`output {}`을 `emit:`으로 교체합니다.

완성된 파일은 다음과 같아야 합니다:

```groovy title="workflows/transform.nf" linenums="1" hl_lines="4 5 6 8 13 14 15"
include { SAY_HELLO_UPPER } from '../modules/say_hello_upper'
include { REVERSE_TEXT } from '../modules/reverse_text'

workflow TRANSFORM_WORKFLOW {
    take:
    input_ch // 메시지가 담긴 입력 채널

    main:
    // 순서대로 변환 적용
    upper_ch = SAY_HELLO_UPPER(input_ch)
    reversed_ch = REVERSE_TEXT(upper_ch)

    emit:
    upper = upper_ch // 대문자 인사말
    reversed = reversed_ch // 뒤집힌 대문자 인사말
}
```

transform 워크플로우가 이제 구성 가능해졌으며 메인 워크플로우로 가져올 준비가 되었습니다.

### 2.3. 메인 워크플로우 업데이트 및 테스트

이제 transformation 워크플로우를 호출하도록 메인 워크플로우를 업데이트합니다.

#### 2.3.1. transformation 워크플로우 가져오기 및 호출

include 구문을 추가하고, 타임스탬프가 추가된 인사말에 연결된 `TRANSFORM_WORKFLOW` 호출을 추가하며, 두 개의 새로운 `publish:` 항목을 추가합니다:

=== "후"

    ```groovy title="main.nf" linenums="1" hl_lines="2 11 12 16 17"
    include { GREETING_WORKFLOW } from './workflows/greeting'
    include { TRANSFORM_WORKFLOW } from './workflows/transform'

    workflow {
        main:
        names = channel.of('Alice', 'Bob', 'Charlie')

        // greeting 워크플로우 실행
        GREETING_WORKFLOW(names)

        // transform 워크플로우 실행
        TRANSFORM_WORKFLOW(GREETING_WORKFLOW.out.timestamped)

        publish:
        greetings = GREETING_WORKFLOW.out.greetings
        upper = TRANSFORM_WORKFLOW.out.upper
        reversed = TRANSFORM_WORKFLOW.out.reversed
    }
    ```

=== "전"

    ```groovy title="main.nf" linenums="1"
    include { GREETING_WORKFLOW } from './workflows/greeting'

    workflow {
        main:
        names = channel.of('Alice', 'Bob', 'Charlie')

        // greeting 워크플로우 실행
        GREETING_WORKFLOW(names)

        publish:
        greetings = GREETING_WORKFLOW.out.greetings
    }
    ```

이렇게 하면 타임스탬프가 추가된 인사말에 transformation 워크플로우가 실행됩니다.

#### 2.3.2. output 블록 업데이트

`output {}` 블록에 `upper`와 `reversed` 항목을 추가하고, 각각 하위 디렉토리를 위한 `path` 지시문을 설정합니다:

=== "후"

    ```groovy title="main.nf" linenums="20" hl_lines="5 6 7 8 9 10"
    output {
        greetings {
            path 'greetings'
        }
        upper {
            path 'upper'
        }
        reversed {
            path 'reversed'
        }
    }
    ```

=== "전"

    ```groovy title="main.nf" linenums="20" hl_lines="2 3 4 5"
    output {
        greetings {
            path 'greetings'
        }
    }
    ```

이렇게 하면 최종 출력이 적절한 디렉토리에 게시됩니다.

#### 2.3.3. 전체 파이프라인 실행

파이프라인을 실행하여 모든 것이 정상적으로 작동하는지 테스트합니다:

```bash
nextflow run main.nf
```

??? success "명령 출력"

    ```console
    N E X T F L O W  ~  version 24.10.0
    Launching `main.nf` [sick_kimura] DSL2 - revision: 8dc45fc6a8
    executor >  local (15)
    [83/1b51f4] process > GREETING_WORKFLOW:VALIDATE_NAME (validating Alice)  [100%] 3 of 3 ✔
    [68/556150] process > GREETING_WORKFLOW:SAY_HELLO (greeting Alice)        [100%] 3 of 3 ✔
    [de/511abd] process > GREETING_WORKFLOW:TIMESTAMP_GREETING (adding tim... [100%] 3 of 3 ✔
    [cd/e6a7e0] process > TRANSFORM_WORKFLOW:SAY_HELLO_UPPER (converting t... [100%] 3 of 3 ✔
    [f0/74ba4a] process > TRANSFORM_WORKFLOW:REVERSE_TEXT (reversing UPPER... [100%] 3 of 3 ✔
    ```

??? abstract "디렉토리 내용"

    ```console
    results/
    ├── greetings
    │   ├── Alice-output.txt
    │   ├── Bob-output.txt
    │   └── Charlie-output.txt
    ├── reversed
    │   ├── REVERSED-UPPER-timestamped_Alice-output.txt
    │   ├── REVERSED-UPPER-timestamped_Bob-output.txt
    │   └── REVERSED-UPPER-timestamped_Charlie-output.txt
    └── upper
        ├── UPPER-timestamped_Alice-output.txt
        ├── UPPER-timestamped_Bob-output.txt
        └── UPPER-timestamped_Charlie-output.txt
    ```

??? abstract "파일 내용"

    ```console title="results/reversed/REVERSED-UPPER-timestamped_Alice-output.txt"
    !ECILA ,OLLEH ]04:50:71 60-30-5202[
    ```

파이프라인이 처음부터 끝까지 정상적으로 작동합니다: 인사말이 대문자로 변환되고 뒤집혔습니다.

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

[사이드 퀘스트 메뉴](../index.md)로 돌아가거나 페이지 오른쪽 하단의 버튼을 클릭하여 목록의 다음 주제로 이동하세요.
