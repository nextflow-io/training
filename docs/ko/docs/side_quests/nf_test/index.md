# nf-test를 활용한 테스트

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } AI 지원 번역 - [자세히 알아보기 및 개선 제안](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

워크플로우의 모든 부분이 의도한 대로 작동하는지 체계적으로 테스트할 수 있는 능력은 재현성과 장기적인 유지보수에 매우 중요하며, 개발 과정에서도 큰 도움이 됩니다.

테스트가 왜 중요한지 잠시 살펴보겠습니다. 워크플로우를 개발할 때 가장 먼저 하는 일 중 하나는 유효하고 결과를 생성해야 하는 테스트 데이터를 준비하는 것입니다. 파이프라인에 첫 번째 process를 추가하고 입력에 연결하여 작동하게 만듭니다. 그런 다음 모든 것이 제대로 작동하는지 확인하기 위해 테스트 데이터로 실행합니다. 정상적으로 작동하면 다음 process로 넘어가 테스트 데이터를 다시 실행합니다. 만족스러운 파이프라인이 완성될 때까지 이 과정을 반복합니다.

그런 다음 `--skip_process`와 같은 간단한 true 또는 false 매개변수를 추가한다고 가정해 보겠습니다. 이제 파이프라인을 두 번 실행해야 합니다. 각 매개변수로 한 번씩 실행하여 예상대로 작동하는지 확인해야 합니다. 그런데 `--skip_process`가 실제로 process를 건너뛰는지 어떻게 확인할까요? 출력 결과를 직접 확인하거나 로그 파일을 살펴봐야 합니다! 이는 번거롭고 오류가 발생하기 쉽습니다.

파이프라인을 개발하다 보면 매번 수동으로 테스트하기에는 너무 복잡해지고, 속도도 느리며 오류가 발생하기 쉽습니다. 게다가 오류를 발견하더라도 파이프라인의 어느 부분에서 오류가 발생했는지 정확히 파악하기가 매우 어렵습니다. 바로 이때 테스트가 필요합니다.

테스트를 통해 파이프라인의 모든 부분이 예상대로 작동하는지 체계적으로 확인할 수 있습니다. 잘 작성된 테스트가 개발자에게 주는 이점은 매우 큽니다:

- **신뢰성**: 테스트가 전체 파이프라인을 커버하므로, 무언가를 변경해도 다른 부분에 영향을 미치지 않는다는 확신을 가질 수 있습니다.
- **신뢰**: 여러 개발자가 파이프라인을 함께 작업할 때, 다른 개발자가 파이프라인이나 각 컴포넌트를 손상시키지 않았다는 것을 알 수 있습니다.
- **투명성**: 테스트는 파이프라인이 어디서 실패하는지 보여주고 문제를 추적하기 쉽게 만들어 줍니다. 또한 process나 워크플로우를 실행하는 방법을 보여주는 문서 역할도 합니다.
- **속도**: 테스트가 자동화되어 있으므로 매우 빠르게 반복 실행할 수 있습니다. 새로운 버그 도입에 대한 걱정 없이 빠르게 반복 개발할 수 있습니다.

작성할 수 있는 테스트의 종류는 다양합니다:

1. **모듈 수준 테스트**: 개별 process에 대한 테스트
2. **워크플로우 수준 테스트**: 단일 워크플로우에 대한 테스트
3. **파이프라인 수준 테스트**: 파이프라인 전체에 대한 테스트
4. **성능 테스트**: 파이프라인의 속도와 효율성에 대한 테스트
5. **스트레스 테스트**: 극한 조건에서 파이프라인의 성능을 평가하여 한계를 파악하는 테스트

개별 process를 테스트하는 것은 다른 언어의 단위 테스트(unit test)와 유사합니다. 워크플로우나 전체 파이프라인을 테스트하는 것은 다른 언어에서 통합 테스트(integration test)라고 불리는 것과 유사하며, 컴포넌트 간의 상호작용을 테스트합니다.

[**nf-test**](https://www.nf-test.com/)는 모듈, 워크플로우, 파이프라인 수준의 테스트를 작성할 수 있는 도구입니다. 간단히 말해, 파이프라인의 모든 개별 부분이 _독립적으로_ 예상대로 작동하는지 체계적으로 확인할 수 있습니다.

### 학습 목표

이 사이드 퀘스트에서는 nf-test를 사용하여 파이프라인에 대한 워크플로우 수준 테스트와 파이프라인이 호출하는 세 가지 process에 대한 모듈 수준 테스트를 작성하는 방법을 학습합니다.

이 사이드 퀘스트를 마치면 다음 기술들을 효과적으로 활용할 수 있습니다:

- 프로젝트에서 nf-test 초기화하기
- 모듈 수준 및 워크플로우 수준 테스트 생성하기
- 일반적인 유형의 assertion 추가하기
- 스냅샷과 콘텐츠 assertion의 적절한 사용 시점 이해하기
- 전체 프로젝트에 대한 테스트 실행하기

이러한 기술들은 파이프라인 프로젝트에 포괄적인 테스트 전략을 구현하는 데 도움이 되며, 파이프라인을 더욱 견고하고 유지보수하기 쉽게 만들어 줍니다.

### 사전 요구사항

이 사이드 퀘스트를 시작하기 전에 다음을 완료해야 합니다:

- [Hello Nextflow](../hello_nextflow/README.md) 튜토리얼 또는 동급의 입문 과정을 완료해야 합니다.
- 기본적인 Nextflow 개념과 메커니즘(process, 채널, 연산자, 파일 작업, 메타데이터)에 익숙해야 합니다.

---

## 0. 시작하기

#### 교육 코드스페이스 열기

아직 열지 않았다면 [환경 설정](../envsetup/index.md)에 설명된 대로 교육 환경을 열어 주세요.

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

#### 프로젝트 디렉토리로 이동

이 튜토리얼의 파일이 있는 디렉토리로 이동합니다.

```bash
cd side-quests/nf-test
```

VSCode에서 이 디렉토리에 포커스를 설정할 수 있습니다:

```bash
code .
```

#### 자료 검토

메인 워크플로우 파일과 파이프라인의 입력을 담고 있는 `greetings.csv`라는 CSV 파일을 확인할 수 있습니다.

```console title="Directory contents"
.
├── greetings.csv
└── main.nf
```

파일에 대한 자세한 설명은 [Hello Nextflow의 준비 운동](../hello_nextflow/00_orientation.md)을 참조하세요.

테스트할 워크플로우는 [Hello Workflow](../hello_nextflow/03_hello_workflow.md)에서 구축한 Hello 워크플로우의 일부입니다.

??? example "Hello Nextflow 워크플로우는 무엇을 하나요?"

    [Hello Nextflow](../hello_nextflow/index.md) 교육을 수행하지 않은 경우, 이 간단한 워크플로우가 무엇을 하는지 간략히 설명합니다.

    워크플로우는 인사말이 담긴 CSV 파일을 받아 네 가지 연속적인 변환 단계를 실행하고, 재미있는 캐릭터가 인사말을 말하는 ASCII 그림이 담긴 단일 텍스트 파일을 출력합니다.

    네 단계는 별도의 모듈 파일에 저장된 Nextflow process(`sayHello`, `convertToUpper`, `collectGreetings`, `cowpy`)로 구현됩니다.

    1. **`sayHello`:** 각 인사말을 별도의 출력 파일에 씁니다 (예: "Hello-output.txt")
    2. **`convertToUpper`:** 각 인사말을 대문자로 변환합니다 (예: "HELLO")
    3. **`collectGreetings`:** 모든 대문자 인사말을 단일 배치 파일로 수집합니다
    4. **`cowpy`:** `cowpy` 도구를 사용하여 ASCII 아트를 생성합니다

    결과는 `results/` 디렉토리에 게시되며, 파이프라인의 최종 출력(기본 매개변수로 실행 시)은 대문자로 변환된 인사말을 말하는 캐릭터의 ASCII 아트가 담긴 일반 텍스트 파일입니다.

    이 사이드 퀘스트에서는 처음 두 process만 포함하는 Hello 워크플로우의 중간 형태를 사용합니다. <!-- TODO: change this to use the full finished workflow as suggested in https://github.com/nextflow-io/training/issues/735 -->

작업할 부분은 `sayHello`와 `convertToUpper` 두 process로 구성됩니다.
전체 워크플로우 코드는 아래에서 확인할 수 있습니다.

??? example "워크플로우 코드"

    ```groovy title="main.nf"
    /*
    * 파이프라인 매개변수
    */
    params.input_file = "greetings.csv"

    /*
    * echo를 사용하여 'Hello World!'를 표준 출력에 출력
    */
    process sayHello {

        publishDir 'results', mode: 'copy'

        input:
            val greeting

        output:
            path "${greeting}-output.txt"

        script:
        """
        echo '$greeting' > '$greeting-output.txt'
        """
    }

    /*
    * 텍스트 변환 유틸리티를 사용하여 인사말을 대문자로 변환
    */
    process convertToUpper {

        publishDir 'results', mode: 'copy'

        input:
            path input_file

        output:
            path "UPPER-${input_file}"

        script:
        """
        cat '$input_file' | tr '[a-z]' '[A-Z]' > UPPER-${input_file}
        """
    }

    workflow {

        // CSV 파일에서 입력 채널 생성
        greeting_ch = channel.fromPath(params.input_file).splitCsv().flatten()

        // 인사말 출력
        sayHello(greeting_ch)

        // 인사말을 대문자로 변환
        convertToUpper(sayHello.out)
    }
    ```

#### 워크플로우 실행

워크플로우가 예상대로 작동하는지 확인하기 위해 실행합니다.

```bash
nextflow run main.nf
```

```console title="Result of running the workflow"
 N E X T F L O W   ~  version 24.10.2

Launching `main.nf` [soggy_linnaeus] DSL2 - revision: bbf79d5c31

executor >  local (6)
[f7/c3be66] sayHello (3)       | 3 of 3 ✔
[cd/e15303] convertToUpper (3) | 3 of 3 ✔
```

축하합니다! 방금 테스트를 실행했습니다!

"잠깐, 무슨 말인가요? 워크플로우를 실행했고 작동했을 뿐인데, 이게 어떻게 테스트인가요?"

좋은 질문입니다!

방금 일어난 일을 분석해 보겠습니다.

기본 매개변수로 워크플로우를 실행하고, 정상적으로 작동하는 것을 확인했으며, 결과에 만족했습니다. 이것이 바로 테스트의 본질입니다. Hello Nextflow 교육 과정을 진행했다면, 시작점으로 사용하는 워크플로우를 실행하여 모든 것이 올바르게 설정되었는지 확인하는 것으로 매 섹션을 시작했다는 것을 알아챘을 것입니다.

소프트웨어 테스트는 기본적으로 이 과정을 자동화합니다.

#### 과제 검토

여러분의 과제는 nf-test를 사용하여 이 워크플로우에 표준화된 테스트를 추가하는 것입니다. 이를 통해 추가적인 변경이 이루어지더라도 모든 부분이 계속 예상대로 작동하는지 쉽게 확인할 수 있습니다.

<!-- TODO: give a bit more details, similar to how it's done in the Metadata side quest -->

#### 준비 확인 목록

시작할 준비가 되었나요?

- [ ] 이 과정의 목표와 사전 요구사항을 이해했습니다
- [ ] 코드스페이스가 실행 중입니다
- [ ] 작업 디렉토리를 적절히 설정했습니다
- [ ] 워크플로우를 성공적으로 실행했습니다
- [ ] 과제를 이해했습니다

모든 항목을 확인했다면 시작할 준비가 된 것입니다.

---

## 1. `nf-test` 초기화

`nf-test` 패키지는 프로젝트에 대한 테스트 개발을 시작하기 위해 몇 가지 설정을 수행하는 초기화 명령을 제공합니다.

```bash
nf-test init
```

다음과 같은 출력이 생성됩니다:

```bash
🚀 nf-test 0.9.3
https://www.nf-test.com
(c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr

Project configured. Configuration is stored in nf-test.config
```

또한 설정 파일 스텁이 포함된 `tests` 디렉토리가 생성됩니다.

### 1.1. nf-test 생성

`nf-test`에는 nf-test 파일을 구축하는 도구 세트가 포함되어 있어 대부분의 작업을 절약할 수 있습니다. 이는 `generate` 하위 명령어로 제공됩니다. 파이프라인에 대한 테스트를 생성합니다:

```bash
nf-test generate pipeline main.nf
```

```console title="Output"
> nf-test generate pipeline main.nf

Load source file '/workspaces/training/side-quests/nf-test/main.nf'
Wrote pipeline test file '/workspaces/training/side-quests/nf-test/tests/main.nf.test

SUCCESS: Generated 1 test files.
```

이렇게 하면 `tests` 디렉토리 내에 `main.nf.test` 파일이 생성됩니다. 이것이 파이프라인 수준 테스트 파일입니다. `tree tests/`를 실행하면 다음과 같은 내용을 볼 수 있습니다:

```console title="Test directory contents"
tests/
├── main.nf.test
└── nextflow.config
```

`main.nf.test` 파일은 파이프라인 수준 테스트 파일입니다. 파일을 열어 내용을 살펴보겠습니다.

```groovy title="tests/main.nf.test"
nextflow_pipeline {

    name "Test Workflow main.nf"
    script "main.nf"

    test("Should run without failures") {

        when {
            params {
                // define parameters here. Example:
                // outdir = "tests/results"
            }
        }

        then {
            assert workflow.success
        }

    }

}
```

테스트 파일의 구조를 이해하는 데 잠시 시간을 내겠습니다.

`nextflow_pipeline` 블록은 모든 파이프라인 수준 테스트의 진입점입니다. 다음을 포함합니다:

- `name`: 테스트의 이름
- `script`: 파이프라인 스크립트의 경로

`test` 블록은 실제 테스트입니다. 다음을 포함합니다:

- `when`: 테스트가 실행되어야 하는 조건. 파이프라인을 실행하는 데 사용될 매개변수를 포함합니다.
- `then`: 수행해야 할 assertion. 파이프라인의 예상 결과를 포함합니다.

쉽게 말하면, 테스트의 논리는 다음과 같이 읽힙니다:
"이 *파이프라인*에 이 *매개변수*가 제공되면(**When**), 이러한 결과를 기대합니다(**then**)."

이것은 기능적인 테스트가 아니며, 다음 섹션에서 이를 기능적인 테스트로 전환하는 방법을 보여드리겠습니다.

### 테스트 이름에 대한 참고사항

위의 예에서는 파이프라인이 성공적으로 실행되는지만 확인하는 기본 테스트에 적합한 기본 이름 "Should run without failures"를 사용했습니다. 그러나 더 구체적인 테스트 케이스를 추가할 때는 실제로 테스트하는 내용을 나타내는 더 설명적인 이름을 사용해야 합니다. 예를 들어:

- "Should convert input to uppercase" - 특정 기능을 테스트할 때
- "Should handle empty input gracefully" - 엣지 케이스를 테스트할 때
- "Should respect max memory parameter" - 리소스 제약을 테스트할 때
- "Should create expected output files" - 파일 생성을 테스트할 때

좋은 테스트 이름은 다음과 같아야 합니다:

1. "Should"로 시작하여 예상되는 동작을 명확히 합니다
2. 테스트하는 특정 기능이나 시나리오를 설명합니다
3. 테스트가 실패했을 때 어떤 기능이 손상되었는지 알 수 있을 만큼 명확해야 합니다

나중에 더 많은 assertion과 구체적인 테스트 케이스를 추가할 때, 각 테스트가 무엇을 검증하는지 명확히 하기 위해 이러한 설명적인 이름을 사용하겠습니다.

### 1.2. 테스트 실행

테스트를 실행하여 어떤 일이 일어나는지 확인합니다.

```bash
nf-test test tests/main.nf.test
```

```console title="nf-test pipeline fail"
> nf-test test tests/main.nf.test

🚀 nf-test 0.9.3
https://www.nf-test.com
(c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr


Test Workflow main.nf

  Test [693ba951] 'Should run without failures' FAILED (4.652s)

  Assertion failed:

  assert workflow.success
         |        |
         workflow false

  Nextflow stdout:

  ERROR ~ No such file or directory: /workspaces/training/side-quests/nf-test/.nf-test/tests/693ba951a20fec36a5a9292ed1cc8a9f/greetings.csv

   -- Check '/workspaces/training/side-quests/nf-test/.nf-test/tests/693ba951a20fec36a5a9292ed1cc8a9f/meta/nextflow.log' file for details
  Nextflow stderr:

FAILURE: Executed 1 tests in 4.679s (1 failed)
```

테스트가 실패했습니다! 무슨 일이 일어났을까요?

1. nf-test가 `when` 블록의 설정을 사용하여 파이프라인을 그대로 실행하려고 했습니다:

```groovy title="tests/main.nf.test"
when {
    params {
        // define parameters here. Example:
        // outdir = "tests/results"
    }
}
```

2. nf-test가 파이프라인의 상태를 확인하고 `when` 블록과 비교했습니다:

```groovy title="tests/main.nf.test"
then {
    assert workflow.success
}
```

nf-test가 파이프라인 실패를 보고하고 Nextflow의 오류 메시지를 제공한 것을 확인하세요:

```console title="Error"
ERROR ~ No such file or directory: /workspaces/training/side-quests/nf-test/.nf-test/tests/693ba951a20fec36a5a9292ed1cc8a9f/greetings.csv
```

문제가 무엇이었을까요? 파이프라인에는 프로젝트 디렉토리에 greetings.csv 파일이 있습니다. nf-test가 파이프라인을 실행할 때 이 파일을 찾으려 하지만 찾을 수 없습니다. 파일은 있는데 무슨 일이 일어나고 있는 걸까요? 경로를 보면 테스트가 `./nf-test/tests/longHashString/` 경로에서 실행되고 있음을 알 수 있습니다. Nextflow와 마찬가지로 nf-test는 모든 것을 격리된 상태로 유지하기 위해 각 테스트에 대한 새 디렉토리를 생성합니다. 데이터 파일이 그 안에 없으므로 원래 테스트에서 파일 경로를 수정해야 합니다.

테스트 파일로 돌아가서 `when` 블록의 파일 경로를 변경합니다.

테스트에서 파이프라인의 루트를 어떻게 가리킬지 궁금할 것입니다. 이는 일반적인 상황이므로 nf-test에는 작업을 쉽게 만들어 주는 다양한 전역 변수가 있습니다. 전체 목록은 [여기](https://www.nf-test.com/docs/testcases/global_variables/)에서 확인할 수 있으며, 지금은 파이프라인 프로젝트의 루트를 의미하는 `projectDir` 변수를 사용하겠습니다.

=== "후"

    ```groovy title="tests/main.nf.test" linenums="1" hl_lines="3"
    when {
        params {
            input_file = "${projectDir}/greetings.csv"
        }
    }
    ```

=== "전"

    ```groovy title="tests/main.nf.test" linenums="1" hl_lines="3 4"
    when {
        params {
            // define parameters here. Example:
            // outdir = "tests/results"
        }
    }
    ```

테스트를 다시 실행하여 작동하는지 확인합니다.

```bash title="nf-test pipeline pass"
nf-test test tests/main.nf.test
```

```console title="Pipeline passes"
> nf-test test tests/main.nf.test

🚀 nf-test 0.9.3
https://www.nf-test.com
(c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr


Test Workflow main.nf

  Test [1d4aaf12] 'Should run without failures' PASSED (1.619s)


SUCCESS: Executed 1 tests in 1.626s
```

성공입니다! 파이프라인이 성공적으로 실행되고 테스트가 통과됩니다. 몇 번을 실행해도 항상 같은 결과를 얻을 수 있습니다!

기본적으로 Nextflow 출력은 숨겨져 있지만, nf-test가 실제로 워크플로우를 실행하고 있다는 것을 확인하려면 `--verbose` 플래그를 사용할 수 있습니다:

```bash
nf-test test tests/main.nf.test --verbose
```

```console title="Pipeline runs all processes"
> nf-test test tests/main.nf.test

🚀 nf-test 0.9.3
https://www.nf-test.com
(c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr


Test Workflow main.nf

  Test [693ba951] 'Should run without failures'
    > Nextflow 24.10.4 is available - Please consider updating your version to it
    > N E X T F L O W  ~  version 24.10.0
    > Launching `/workspaces/training/side-quests/nf-test/main.nf` [zen_ampere] DSL2 - revision: bbf79d5c31
    > [2b/61e453] Submitted process > sayHello (2)
    > [31/4e1606] Submitted process > sayHello (1)
    > [bb/5209ee] Submitted process > sayHello (3)
    > [83/83db6f] Submitted process > convertToUpper (2)
    > [9b/3428b1] Submitted process > convertToUpper (1)
    > [ca/0ba51b] Submitted process > convertToUpper (3)
    PASSED (5.206s)


SUCCESS: Executed 1 tests in 5.239s
```

### 1.3. Assertion 추가

간단한 확인 방법은 파이프라인이 예상하는 모든 process를 실행하고 있으며 아무것도 자동으로 건너뛰지 않는지 확인하는 것입니다. 파이프라인은 3개의 인사말 각각에 대해 `sayHello`와 `convertToUpper`라는 총 6개의 process를 실행합니다.

파이프라인이 예상한 수의 process를 실행하는지 확인하는 assertion을 테스트에 추가합니다. 또한 테스트 이름을 테스트하는 내용을 더 잘 반영하도록 업데이트합니다.

=== "후"

    ```groovy title="tests/main.nf.test" linenums="1" hl_lines="1 11"
        test("Should run successfully with correct number of processes") {

            when {
                params {
                    input_file = "${projectDir}/greetings.csv"
                }
            }

            then {
                assert workflow.success
                assert workflow.trace.tasks().size() == 6
            }

        }
    ```

=== "전"

    ```groovy title="tests/main.nf.test" linenums="1" hl_lines="1"
        test("Should run without failures") {

            when {
                params {
                    input_file = "${projectDir}/greetings.csv"
                }
            }

            then {
                assert workflow.success
            }

        }
    ```

이제 테스트 이름이 실제로 검증하는 내용을 더 잘 반영합니다. 단순히 파이프라인이 실패 없이 실행되는지뿐만 아니라 예상한 수의 process가 실행되는지도 확인합니다.

테스트를 다시 실행하여 작동하는지 확인합니다.

```bash title="nf-test pipeline pass"
nf-test test tests/main.nf.test
```

```console title="Pipeline passes with assertions"
🚀 nf-test 0.9.3
https://www.nf-test.com
(c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr


Test Workflow main.nf

  Test [1d4aaf12] 'Should run successfully with correct number of processes' PASSED (1.567s)


SUCCESS: Executed 1 tests in 1.588s
```

성공입니다! 파이프라인이 성공적으로 실행되고 테스트가 통과됩니다. 이제 파이프라인의 전체 상태뿐만 아니라 세부 사항도 테스트하기 시작했습니다.

### 1.4. 출력 테스트

출력 파일이 생성되었는지 확인하는 assertion을 테스트에 추가합니다. 결과를 더 쉽게 해석할 수 있도록 설명적인 이름과 함께 별도의 테스트로 추가합니다.

=== "후"

    ```groovy title="tests/main.nf.test" linenums="1" hl_lines="14-33"
        test("Should run successfully with correct number of processes") {

            when {
                params {
                    input_file = "${projectDir}/greetings.csv"
                }
            }

            then {
                assert workflow.success
                assert workflow.trace.tasks().size() == 6
            }

        }

        test("Should produce correct output files") {

            when {
                params {
                    input_file = "${projectDir}/greetings.csv"
                }
            }

            then {
                assert file("$launchDir/results/Bonjour-output.txt").exists()
                assert file("$launchDir/results/Hello-output.txt").exists()
                assert file("$launchDir/results/Holà-output.txt").exists()
                assert file("$launchDir/results/UPPER-Bonjour-output.txt").exists()
                assert file("$launchDir/results/UPPER-Hello-output.txt").exists()
                assert file("$launchDir/results/UPPER-Holà-output.txt").exists()
            }

        }
    ```

=== "전"

    ```groovy title="tests/main.nf.test" linenums="1" hl_lines="14"
        test("Should run successfully with correct number of processes") {

            when {
                params {
                    input_file = "${projectDir}/greetings.csv"
                }
            }

            then {
                assert workflow.success
                assert workflow.trace.tasks().size() == 6
            }

        }
    ```

테스트를 다시 실행하여 작동하는지 확인합니다.

```bash title="nf-test pipeline pass"
nf-test test tests/main.nf.test
```

```console title="Pipeline passes with file assertions"
> nf-test test tests/main.nf.test

🚀 nf-test 0.9.3
https://www.nf-test.com
(c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr


Test Workflow main.nf

  Test [f0e08a68] 'Should run successfully with correct number of processes' PASSED (8.144s)
  Test [d7e32a32] 'Should produce correct output files' PASSED (6.994s)


SUCCESS: Executed 2 tests in 15.165s
```

성공입니다! 파이프라인이 성공적으로 완료되고, 올바른 수의 process가 실행되었으며, 출력 파일이 생성되었기 때문에 테스트가 통과됩니다. 이를 통해 테스트에 설명적인 이름을 제공하는 것이 얼마나 유용한지도 알 수 있습니다.

이것은 시작에 불과합니다. 파이프라인의 세부 사항을 확인하는 assertion을 계속 작성할 수 있지만, 지금은 파이프라인 내부 테스트로 넘어가겠습니다.

### 핵심 정리

파이프라인에 대한 nf-test를 작성하는 방법을 알게 되었습니다.

### 다음 단계

Nextflow process를 테스트하는 방법을 학습합니다.

---

## 2. Nextflow process 테스트

파이프라인의 모든 부분에 대한 테스트를 작성할 필요는 없지만, 테스트가 많을수록 파이프라인에 대해 더 포괄적으로 확인할 수 있고 예상대로 작동한다는 확신을 가질 수 있습니다. 이 섹션에서는 파이프라인의 두 process를 개별 단위로 테스트합니다.

### 2.1. `sayHello` process 테스트

`sayHello` process부터 시작합니다.

`nf-test generate` 명령어를 다시 사용하여 process에 대한 테스트를 생성합니다.

```bash
nf-test generate process main.nf
```

```console title="Output"
> nf-test generate process main.nf

Load source file '/workspaces/training/side-quests/nf-test/main.nf'
Wrote process test file '/workspaces/training/side-quests/nf-test/tests/main.sayhello.nf.test
Wrote process test file '/workspaces/training/side-quests/nf-test/tests/main.converttoupper.nf.test

SUCCESS: Generated 2 test files.
```

지금은 `main.sayhello.nf.test` 파일의 `sayhello` process에 집중합니다.

파일을 열어 내용을 살펴보겠습니다.

```groovy title="tests/main.sayhello.nf.test"
nextflow_process {

    name "Test Process sayHello"
    script "main.nf"
    process "sayHello"

    test("Should run without failures") {

        when {
            params {
                // define parameters here. Example:
                // outdir = "tests/results"
            }
            process {
                """
                // define inputs of the process here. Example:
                // input[0] = file("test-file.txt")
                """
            }
        }

        then {
            assert process.success
            assert snapshot(process.out).match()
        }

    }

}
```

이전과 마찬가지로 테스트 세부 정보로 시작하고 `when` 및 `then` 블록이 이어집니다. 그러나 process에 대한 입력을 정의할 수 있는 추가적인 `process` 블록도 있습니다.

테스트를 실행하여 작동하는지 확인합니다.

```bash title="nf-test pipeline pass"
nf-test test tests/main.sayhello.nf.test
```

```console title="Process test fails"
> nf-test test tests/main.sayhello.nf.test

🚀 nf-test 0.9.3
https://www.nf-test.com
(c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr


Test Process sayHello

  Test [1eaad118] 'Should run without failures' FAILED (4.876s)

  Assertion failed:

  assert process.success
         |       |
         |       false
         sayHello

  Nextflow stdout:

  Process `sayHello` declares 1 input but was called with 0 arguments
  Nextflow stderr:

FAILURE: Executed 1 tests in 4.884s (1 failed)
```

`sayHello` process가 1개의 입력을 선언했지만 0개의 인자로 호출되었기 때문에 테스트가 실패합니다. process에 입력을 추가하여 이를 수정합니다. [Hello Workflow](../hello_nextflow/03_hello_workflow.md)(및 위의 준비 운동 섹션)에서 `sayHello` process가 단일 값 입력을 받는다는 것을 기억하세요. 이를 제공해야 합니다. 또한 테스트 이름을 테스트하는 내용을 더 잘 반영하도록 수정합니다.

=== "후"

    ```groovy title="tests/main.sayhello.nf.test" linenums="1" hl_lines="1 10"
        test("Should run without failures and produce correct output") {

            when {
                params {
                    // define parameters here. Example:
                    // outdir = "tests/results"
                }
                process {
                    """
                    input[0] = "hello"
                    """
                }
            }

            then {
                assert process.success
                assert snapshot(process.out).match()
            }

        }
    ```

=== "전"

    ```groovy title="tests/main.sayhello.nf.test" linenums="1" hl_lines="1 10 11"
        test("Should run without failures") {

            when {
                params {
                    // define parameters here. Example:
                    // outdir = "tests/results"
                }
                process {
                    """
                    // define inputs of the process here. Example:
                    // input[0] = file("test-file.txt")
                    """
                }
            }

            then {
                assert process.success
                assert snapshot(process.out).match()
            }

        }
    ```

테스트를 다시 실행하여 작동하는지 확인합니다.

```console title="nf-test pipeline pass"
> nf-test test tests/main.sayhello.nf.test

🚀 nf-test 0.9.3
https://www.nf-test.com
(c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr


Test Process sayHello

  Test [f91a1bcd] 'Should run without failures and produce correct output' PASSED (1.604s)
  Snapshots:
    1 created [Should run without failures and produce correct output]


Snapshot Summary:
  1 created

SUCCESS: Executed 1 tests in 1.611s
```

성공입니다! `sayHello` process가 성공적으로 실행되고 출력이 생성되었기 때문에 테스트가 통과됩니다.

### 2.2. 테스트로 생성된 스냅샷 확인

`tests/main.sayhello.nf.test` 파일을 보면 assertion 블록에서 `snapshot()` 메서드를 사용하는 것을 확인할 수 있습니다:

```groovy title="tests/main.sayhello.nf.test"
assert snapshot(process.out).match()
```

이는 nf-test에게 `sayHello` process의 출력 스냅샷을 생성하도록 지시합니다. 스냅샷 파일의 내용을 살펴보겠습니다.

```console title="Snapshot file contents"
code tests/main.sayhello.nf.test.snap
```

여기에 출력하지는 않겠지만, process 및 process 출력의 세부 정보가 담긴 JSON 파일을 볼 수 있습니다. 특히 다음과 같은 줄을 확인할 수 있습니다:

```json title="Snapshot file contents"
"0": [
    "hello-output.txt:md5,b1946ac92492d2347c6235b4d2611184"
]
```

이는 명시적으로 테스트하는 `sayHello` process가 생성한 출력을 나타냅니다. 테스트를 다시 실행하면 프로그램은 새 출력이 원래 기록된 출력과 일치하는지 확인합니다. 이는 process 출력이 변경되지 않는지 테스트하는 빠르고 간단한 방법이므로 nf-test가 기본값으로 제공합니다.

!!!warning "경고"

    즉, 원래 실행에서 기록한 출력이 올바른지 확인해야 합니다!

향후 개발 과정에서 코드의 무언가가 변경되어 출력이 달라지면 테스트가 실패하고 변경이 예상된 것인지 아닌지 판단해야 합니다.

- 코드에서 무언가가 손상된 것으로 판명되면 수정해야 하며, 수정된 코드가 테스트를 통과할 것으로 예상됩니다.
- 예상된 변경인 경우(예: 도구가 개선되어 결과가 더 좋아진 경우) 새 출력을 일치시킬 참조로 수락하도록 스냅샷을 업데이트해야 합니다. nf-test에는 이 목적을 위한 `--update-snapshot` 매개변수가 있습니다.

테스트를 다시 실행하면 테스트가 통과되는 것을 확인할 수 있습니다:

```console title="nf-test process pass with snapshot"
> nf-test test tests/main.sayhello.nf.test

🚀 nf-test 0.9.3
https://www.nf-test.com
(c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr


Test Process sayHello

  Test [f91a1bcd] 'Should run without failures and produce correct output' PASSED (1.675s)


SUCCESS: Executed 1 tests in 1.685s
```

성공입니다! `sayHello` process가 성공적으로 실행되고 출력이 스냅샷과 일치하기 때문에 테스트가 통과됩니다.

### 2.3. 스냅샷의 대안: 직접 콘텐츠 Assertion

스냅샷은 출력의 모든 변경 사항을 감지하는 데 유용하지만, 때로는 전체 파일 일치에 너무 엄격하지 않고 특정 콘텐츠를 확인하고 싶을 때가 있습니다. 예를 들어:

- 출력의 일부가 변경될 수 있지만(타임스탬프, 랜덤 ID 등) 특정 핵심 콘텐츠가 반드시 있어야 하는 경우
- 출력에서 특정 패턴이나 값을 확인하고 싶은 경우
- 성공을 구성하는 요소에 대해 테스트를 더 명시적으로 만들고 싶은 경우

다음은 특정 콘텐츠를 확인하도록 테스트를 수정하는 방법입니다:

=== "후"

    ```groovy title="tests/main.sayhello.nf.test" linenums="1" hl_lines="1 5 16 17"
        test("Should run without failures and contain expected greeting") {

            when {
                params {
                    // define parameters here
                }
                process {
                    """
                    input[0] = "hello"
                    """
                }
            }

            then {
                assert process.success
                assert path(process.out[0][0]).readLines().contains('hello')
                assert !path(process.out[0][0]).readLines().contains('HELLO')
            }

        }
    ```

=== "전"

    ```groovy title="tests/main.sayhello.nf.test" linenums="1" hl_lines="1 5 6 17"
        test("Should run without failures and produce correct output") {

            when {
                params {
                    // define parameters here. Example:
                    // outdir = "tests/results"
                }
                process {
                    """
                    input[0] = "hello"
                    """
                }
            }

            then {
                assert process.success
                assert snapshot(process.out).match()
            }

        }
    ```

nf-test는 process 출력을 리스트의 리스트로 인식하므로, `process.out[0][0]`은 이 process의 첫 번째 채널 항목(또는 '방출')의 첫 번째 부분을 가져옵니다.

이 접근 방식은:

- 출력에서 기대하는 내용을 명확하게 합니다
- 출력의 관련 없는 변경에 더 탄력적입니다
- 테스트 실패 시 더 나은 오류 메시지를 제공합니다
- 더 복잡한 검증(regex 패턴, 수치 비교 등)을 허용합니다

테스트를 실행하여 작동하는지 확인합니다.

```bash title="nf-test pipeline pass"
nf-test test tests/main.sayhello.nf.test
```

```console title="Process test fails"
> nf-test test tests/main.sayhello.nf.test

🚀 nf-test 0.9.3
https://www.nf-test.com
(c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr


Test Process sayHello

  Test [58df4e4b] 'Should run without failures and contain expected greeting' PASSED (7.196s)


SUCCESS: Executed 1 tests in 7.208s
```

### 2.4. `convertToUpper` process 테스트

`tests/main.converttoupper.nf.test` 파일을 열어 내용을 살펴보겠습니다:

```groovy title="tests/main.converttoupper.nf.test"
nextflow_process {

    name "Test Process convertToUpper"
    script "main.nf"
    process "convertToUpper"

    test("Should run without failures") {

        when {
            params {
                // define parameters here. Example:
                // outdir = "tests/results"
            }
            process {
                """
                // define inputs of the process here. Example:
                // input[0] = file("test-file.txt")
                """
            }
        }

        then {
            assert process.success
            assert snapshot(process.out).match()
        }

    }

}
```

이것은 `sayHello` process와 유사한 테스트이지만 `convertToUpper` process를 테스트합니다. `sayHello`와 마찬가지로 `convertToUpper` process도 단일 path 입력을 받지만 지정하지 않았기 때문에 이 테스트도 실패할 것입니다.

이제 대문자로 변환하려는 텍스트가 포함된 단일 입력 파일을 convertToUpper process에 제공해야 합니다. 여러 가지 방법이 있습니다:

- 테스트 전용 파일을 생성할 수 있습니다
- 기존 data/greetings.csv 파일을 재사용할 수 있습니다
- 테스트 내에서 즉석으로 생성할 수 있습니다

지금은 파이프라인 수준 테스트에서 사용한 예시를 활용하여 기존 data/greetings.csv 파일을 재사용합니다. 이전과 마찬가지로 테스트하는 내용을 더 잘 반영하도록 테스트 이름을 지정할 수 있지만, 이번에는 (다른 process에서 했던 것처럼) 특정 문자열을 확인하는 대신 콘텐츠를 '스냅샷'으로 남겨두겠습니다.

=== "후"

    ```groovy title="tests/main.converttoupper.nf.test" linenums="1" hl_lines="1 10"
        test("Should run without failures and produce correct output") {

            when {
                params {
                    // define parameters here. Example:
                    // outdir = "tests/results"
                }
                process {
                    """
                    input[0] = "${projectDir}/greetings.csv"
                    """
                }
            }

            then {
                assert process.success
                assert snapshot(process.out).match()
            }

        }
    ```

=== "전"

    ```groovy title="tests/main.converttoupper.nf.test" linenums="1" hl_lines="1 10 11"
        test("Should run without failures") {

            when {
                params {
                    // define parameters here. Example:
                    // outdir = "tests/results"
                }
                process {
                    """
                    // define inputs of the process here. Example:
                    // input[0] = file("test-file.txt")
                    """
                }
            }

            then {
                assert process.success
                assert snapshot(process.out).match()
            }

        }
    ```

테스트를 실행합니다!

```bash title="nf-test pipeline pass"
nf-test test tests/main.converttoupper.nf.test
```

```console title="nf-test process convertToUpper pass"
> nf-test test tests/main.converttoupper.nf.test

🚀 nf-test 0.9.3
https://www.nf-test.com
(c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr


Test Process convertToUpper

  Test [c59b6044] 'Should run without failures and produce correct output' PASSED (1.755s)
  Snapshots:
    1 created [Should run without failures and produce correct output]


Snapshot Summary:
  1 created

SUCCESS: Executed 1 tests in 1.764s
```

`convertToUpper` process에 대한 스냅샷 파일이 `tests/main.converttoupper.nf.test.snap`에 생성되었습니다. 테스트를 다시 실행하면 nf-test가 다시 통과되는 것을 확인할 수 있습니다.

```bash title="nf-test process convertToUpper pass"
nf-test test tests/main.converttoupper.nf.test
```

```console title="nf-test process convertToUpper pass"
> nf-test test tests/main.converttoupper.nf.test

🚀 nf-test 0.9.3
https://www.nf-test.com
(c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr


Test Process convertToUpper

  Test [c59b6044] 'Should run without failures and produce correct output' PASSED (1.798s)


SUCCESS: Executed 1 tests in 1.811s
```

### 핵심 정리

Nextflow process에 대한 테스트를 작성하고 실행하는 방법을 알게 되었습니다.

### 다음 단계

모든 것을 한 번에 테스트하는 방법을 학습합니다!

## 3. 전체 저장소에 대한 테스트 실행

각 컴포넌트에 대해 nf-test를 실행하는 것도 좋지만, 번거롭고 오류가 발생하기 쉽습니다. 모든 것을 한 번에 테스트할 수는 없을까요?

할 수 있습니다!

전체 저장소에 대해 nf-test를 실행합니다.

### 3.1. 전체 저장소에 대해 nf-test 실행

`nf-test test` 명령어를 실행하여 전체 저장소에 대해 nf-test를 실행할 수 있습니다.

```bash
nf-test test .
```

현재 디렉토리에서 모든 것을 실행하기 위해 `.`를 사용합니다. 이렇게 하면 모든 테스트가 포함됩니다!

```console title="nf-test repo pass"
> nf-test test .

🚀 nf-test 0.9.3
https://www.nf-test.com
(c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr


Test Process convertToUpper

  Test [3d26d9af] 'Should run without failures and produce correct output' PASSED (4.155s)

Test Workflow main.nf

  Test [f183df37] 'Should run successfully with correct number of processes' PASSED (3.33s)
  Test [d7e32a32] 'Should produce correct output files' PASSED (3.102s)

Test Process sayHello

  Test [58df4e4b] 'Should run without failures and contain expected greeting' PASSED (2.614s)


SUCCESS: Executed 4 tests in 13.481s
```

확인해 보세요! 단일 명령으로 각 process에 대해 1개, 전체 파이프라인에 대해 2개, 총 4개의 테스트를 실행했습니다. 대규모 코드베이스에서 이것이 얼마나 강력한지 상상해 보세요!

---

## 요약

이 사이드 퀘스트에서는 nf-test의 기능을 활용하여 개별 process에 대한 테스트와 전체 파이프라인에 대한 엔드투엔드 테스트를 생성하고 실행하는 방법을 학습했습니다.
이제 출력 검증의 두 가지 주요 접근 방식인 스냅샷과 직접 콘텐츠 assertion, 그리고 각각을 언제 사용해야 하는지 알게 되었습니다.
또한 테스트를 하나씩 실행하거나 전체 프로젝트에 대해 실행하는 방법도 알게 되었습니다.

이러한 기술을 실제 작업에 적용하면 다음을 보장할 수 있습니다:

- 코드가 예상대로 작동합니다
- 변경 사항이 기존 기능을 손상시키지 않습니다
- 다른 개발자들이 자신 있게 기여할 수 있습니다
- 문제를 빠르게 식별하고 수정할 수 있습니다
- 출력 콘텐츠가 기대치와 일치합니다

### 핵심 패턴

<!-- TODO: Can we add snippets of code below to illustrate? -->

1. 파이프라인 수준 테스트:
   - 기본 성공 테스트
   - process 수 검증
   - 출력 파일 존재 확인
2. process 수준 테스트
3. 출력 검증의 두 가지 접근 방식:
   - 완전한 출력 검증을 위한 스냅샷 사용
   - 특정 콘텐츠 확인을 위한 직접 콘텐츠 assertion 사용
4. 단일 명령으로 저장소의 모든 테스트 실행

### 추가 자료

더 고급 테스트 기능과 모범 사례는 [nf-test 문서](https://www.nf-test.com/)를 참조하세요. 다음을 고려해 볼 수 있습니다:

- 테스트에 더 포괄적인 assertion 추가하기
- 엣지 케이스 및 오류 조건에 대한 테스트 작성하기
- 테스트를 자동으로 실행하도록 지속적 통합(CI) 설정하기
- 워크플로우 및 모듈 테스트와 같은 다른 유형의 테스트 학습하기
- 더 고급 콘텐츠 검증 기술 탐색하기

**기억하세요:** 테스트는 코드가 어떻게 동작해야 하는지에 대한 살아있는 문서입니다. 더 많은 테스트를 작성하고 assertion이 더 구체적일수록 파이프라인의 신뢰성에 대한 확신을 가질 수 있습니다.

---

## 다음 단계

[사이드 퀘스트 메뉴](../)로 돌아가거나 페이지 오른쪽 하단의 버튼을 클릭하여 목록의 다음 주제로 이동하세요.
